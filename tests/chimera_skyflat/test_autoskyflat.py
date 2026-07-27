# SPDX-FileCopyrightText: 2015-present Antonio Kanaan <kanaan@astro.ufsc.br>
# SPDX-License-Identifier: GPL-2.0-or-later

import json
import threading

import numpy as np
import pytest
from astropy.io import fits
from chimera.util.coord import Coord
from chimera.util.position import Position

from chimera_skyflat.controllers.autoskyflat import (
    AutoSkyFlat,
    SkyFlatRunningException,
)
from tests.chimera_skyflat.fakes import (
    COEFFICIENTS,
    ClockEvent,
    FakeCamera,
    FakeClock,
    FakeFilterWheel,
    FakeSite,
    FakeSky,
    FakeTelescope,
)

IDEAL = 25000.0


@pytest.fixture
def coefficients_file(tmp_path):
    path = tmp_path / "skyflat_coefficients.json"
    path.write_text(json.dumps(COEFFICIENTS))
    return path


def build(
    tmp_path,
    coefficients_file,
    alt0=-2.0,
    rate=-0.004,
    factor=1.0,
    real_waits=False,
    **config,
):
    """An AutoSkyFlat wired to the fakes, plus the fakes themselves.

    Waits move the fake clock instead of the wall clock, so a controller
    that waits for the sky to dim makes progress here too.
    """
    clock = FakeClock()
    site = FakeSite(clock, alt0=alt0, rate=rate)
    sky = FakeSky(site, COEFFICIENTS, factor=factor)
    camera = FakeCamera(clock, sky, tmp_path)
    wheel = FakeFilterWheel()
    telescope = FakeTelescope()

    flat = AutoSkyFlat()
    flat["coefficients_file"] = str(coefficients_file)
    flat["ideal_counts"] = IDEAL
    flat["sun_alt_hi"] = -1
    flat["sun_alt_low"] = -20
    flat["exptime_max"] = 300
    flat["tracking"] = True
    for key, value in config.items():
        flat[key] = value

    flat._get_site = lambda: site
    flat._get_cam = lambda: camera
    flat._get_filter_wheel = lambda: wheel
    flat._get_tel = lambda: telescope

    if not real_waits:
        flat._abort = ClockEvent(clock)

    events = []
    flat.expose_complete = lambda *args: events.append(args)

    return flat, dict(
        clock=clock,
        site=site,
        sky=sky,
        camera=camera,
        wheel=wheel,
        telescope=telescope,
        events=events,
    )


#
# sun altitude normalisation
#


def test_sky_level_returns_plain_float(tmp_path):
    """numpy.float64 is not msgspec-serializable: leaking it into the
    expose_complete event payload silently dropped the publication, so
    subscribers (the robobs frame counter) never heard about the frame."""
    path = tmp_path / "flat.fits"
    fits.PrimaryHDU(np.full((8, 8), 1000.0, dtype=np.float32)).writeto(path)

    level = AutoSkyFlat.get_sky_level.func(None, str(path), None)

    assert type(level) is float
    assert level == 1000.0


def test_sun_altitude_from_a_position_is_degrees():
    """Site.sunpos() returns a Position of Coords; np.radians() on a Coord
    raises TypeError, and reading it as radians pinned every flat at
    exptime_max (2026-07-21)."""
    position = Position.from_alt_az(Coord.from_d(-8.5), Coord.from_d(180))

    assert AutoSkyFlat._altitude_in_degrees(position) == pytest.approx(-8.5)


def test_sun_altitude_accepts_the_old_tuple_and_plain_floats():
    assert AutoSkyFlat._altitude_in_degrees((-8.5, 180.0)) == pytest.approx(-8.5)
    assert AutoSkyFlat._altitude_in_degrees(-8.5) == pytest.approx(-8.5)


def test_dusk_and_dawn_come_from_the_sun_not_the_wall_clock(
    tmp_path, coefficients_file
):
    """The fast-forward simulation clock moves the modelled sky without
    moving local time, so localtime().hour cannot tell dusk from dawn."""
    dusk, _ = build(tmp_path, coefficients_file, rate=-0.004)
    dawn, _ = build(tmp_path, coefficients_file, rate=+0.004)

    # the sign of the rate is the dusk/dawn test everywhere in the controller
    assert dusk._sun_track()[1] < 0
    assert dawn._sun_track()[1] > 0


#
# exposure time calculator
#


def test_computed_exposure_reaches_the_ideal_counts(tmp_path, coefficients_file):
    flat, fakes = build(tmp_path, coefficients_file, alt0=-6.0)
    flat._load_coefficients("CLEAR")

    exptime, expected = flat.compute_sky_flat_time()

    # what the sky really delivers in that time, integrated
    delivered = fakes["sky"].counts(fakes["clock"].now, exptime)
    assert expected == pytest.approx(IDEAL, rel=0.02)
    assert delivered == pytest.approx(IDEAL, rel=0.02)


def test_the_calculator_does_not_poll_the_sun_once_per_increment(
    tmp_path, coefficients_file
):
    """A 0.2 s increment over a 300 s frame is 1500 bus round-trips."""
    flat, fakes = build(tmp_path, coefficients_file, alt0=-9.0)
    flat._load_coefficients("CLEAR")

    flat.compute_sky_flat_time()

    assert fakes["site"].sunpos_calls <= 4


def test_dusk_gives_up_when_the_sky_is_past_exptime_max(tmp_path, coefficients_file):
    flat, _ = build(tmp_path, coefficients_file, alt0=-16.0, exptime_max=10)
    flat._load_coefficients("HBETA")

    assert flat.compute_sky_flat_time() is False


def test_dawn_waits_instead_of_giving_up(tmp_path, coefficients_file):
    """At dawn the sky is brightening: waiting gets the exposure under the
    limit, so the calculator must not end the filter."""
    flat, fakes = build(
        tmp_path, coefficients_file, alt0=-16.0, rate=+0.02, exptime_max=30
    )
    flat._load_coefficients("CLEAR")

    exptime, _ = flat.compute_sky_flat_time()

    assert 0 < exptime <= 30


#
# issue 12: the correction loop used to ring
#


@pytest.mark.parametrize("factor", [0.25, 0.5, 2.0, 4.0])
def test_the_correction_converges_without_ringing(tmp_path, coefficients_file, factor):
    """2026-07-22 CLEAR set: 15.7k -> 26.3k -> 40.8k -> 15.5k -> 37.6k ->
    52.1k counts around a 25k target, frame 6 at ~80% of full well. The
    additive unit-gain correction overshot every frame against a sky
    dimming by a factor of a few per frame."""
    flat, fakes = build(tmp_path, coefficients_file, factor=factor)

    taken = flat.get_flats("CLEAR", n_flats=8)

    levels = [level for *_, level in fakes["events"]]
    errors = [level - IDEAL for level in levels]
    assert taken == 8
    # the first frame is blind, the coefficients file is all we know then.
    # From there on: no ringing - the error crosses the target at most once
    # and never grows (the 2026-07-22 set crossed it four times and grew
    # every time)
    crossings = sum(1 for a, b in zip(errors, errors[1:]) if a * b < 0)
    assert crossings <= 1, f"ringing: {levels}"
    assert all(
        abs(b) <= abs(a) + 0.05 * IDEAL for a, b in zip(errors, errors[1:])
    ), f"diverging: {levels}"
    # and it settles on the target
    for level in levels[3:]:
        assert level == pytest.approx(IDEAL, rel=0.25), f"not converged: {levels}"


def test_the_correction_is_damped_not_proportional(tmp_path, coefficients_file):
    flat, _ = build(tmp_path, coefficients_file, correction_damping=0.5)

    # measured four times the expectation: a unit-gain correction would
    # take the whole factor, the damped one takes its square root
    assert flat._correct_model(1.0, 4 * IDEAL, IDEAL)[0] == pytest.approx(2.0)
    assert flat._correct_model(1.0, IDEAL / 4, IDEAL)[0] == pytest.approx(0.5)

    # ... but until a frame lands on target, damping only slows the model
    # down: an uncalibrated model takes the whole ratio
    assert flat._correct_model(1.0, 4 * IDEAL, IDEAL, calibrated=False)[0] == (
        pytest.approx(4.0)
    )


def test_a_single_wild_frame_cannot_run_away_with_the_model(
    tmp_path, coefficients_file
):
    flat, _ = build(tmp_path, coefficients_file)

    # a cloud, a passing satellite, a bad readout: clipped before damping
    assert flat._correct_model(1.0, 1e6 * IDEAL, IDEAL)[0] == pytest.approx(2.0)
    assert flat._correct_model(1.0, 0.0, IDEAL)[0] == 1.0


def test_the_exposure_time_cannot_jump_between_frames(tmp_path, coefficients_file):
    flat, _ = build(tmp_path, coefficients_file, max_exptime_step=1.5)

    exptime, expected = flat._clamp_exptime(90.0, IDEAL, last_exptime=10.0)

    assert exptime == pytest.approx(15.0)
    # the expectation follows the clamp, otherwise the next correction
    # would be fed counts the frame never had a chance to collect
    assert expected == pytest.approx(IDEAL * 15.0 / 90.0)


def test_a_saturated_frame_is_not_published_nor_counted(tmp_path, coefficients_file):
    flat, fakes = build(tmp_path, coefficients_file, factor=4.0, max_counts=1.2 * IDEAL)

    taken = flat.get_flats("CLEAR", n_flats=3)

    levels = [level for *_, level in fakes["events"]]
    assert taken == 3 == len(levels)
    assert all(level <= 1.2 * IDEAL for level in levels)
    # the discarded frames were still exposed, they just did not count
    assert len(fakes["camera"].exposures) > taken


#
# issue 13: filter fallback, both directions
#


def test_dusk_walks_up_to_a_more_sensitive_filter(tmp_path, coefficients_file):
    flat, _ = build(tmp_path, coefficients_file)

    assert flat._next_filter("V", {"V"}, dusk=True) == "R"
    assert flat._next_filter("R", {"R", "V"}, dusk=True) == "B"
    assert flat._next_filter("CLEAR", {"CLEAR"}, dusk=True) is None


def test_dawn_walks_down_to_a_less_sensitive_filter(tmp_path, coefficients_file):
    flat, _ = build(tmp_path, coefficients_file)

    assert flat._next_filter("CLEAR", {"CLEAR"}, dusk=False) == "B"
    assert flat._next_filter("R", {"R", "CLEAR"}, dusk=False) == "V"
    assert flat._next_filter("HBETA", {"HBETA"}, dusk=False) is None


def test_a_saturating_filter_at_dawn_steps_down_instead_of_stopping(
    tmp_path, coefficients_file
):
    """The dawn twin of the dusk fallback: too *short* an exposure used to
    log 'Exposure time too low' and end the sequence with filters to
    spare."""
    flat, fakes = build(
        tmp_path,
        coefficients_file,
        alt0=-2.0,
        rate=+0.004,
        exptime_min=2.0,
    )

    flat.get_flats("CLEAR", n_flats=2)

    assert fakes["wheel"].current != "CLEAR"
    assert [f for f in fakes["wheel"].moves if f != "CLEAR"], fakes["wheel"].moves


def test_the_fallback_can_be_turned_off(tmp_path, coefficients_file):
    flat, _ = build(tmp_path, coefficients_file, filter_fallback=False)

    assert flat._next_filter("V", {"V"}, dusk=True) is None


def test_an_unknown_filter_is_reported_clearly(tmp_path, coefficients_file):
    flat, _ = build(tmp_path, coefficients_file)

    with pytest.raises(Exception, match="NOSUCHFILTER"):
        flat.get_flats("NOSUCHFILTER", n_flats=1)


#
# window, aborts and re-entrancy
#


def test_flats_stop_when_the_sun_leaves_the_window(tmp_path, coefficients_file):
    flat, fakes = build(tmp_path, coefficients_file, alt0=-19.9, rate=-0.05)

    taken = flat.get_flats("CLEAR", n_flats=100)

    assert taken < 100
    assert not fakes["telescope"].tracking


def test_waiting_for_the_window_ends_on_abort(tmp_path, coefficients_file):
    """The wait used to be a plain sleep loop; --stop could not reach it."""
    flat, _ = build(
        tmp_path, coefficients_file, alt0=+10.0, rate=-0.0001, real_waits=True
    )

    done, waiting = threading.Event(), threading.Event()

    class WatchedEvent(threading.Event):
        def wait(self, timeout=None):
            waiting.set()
            return super().wait(timeout)

    flat._abort = WatchedEvent()

    def run():
        flat.get_flats("CLEAR", n_flats=1)
        done.set()

    worker = threading.Thread(target=run, daemon=True)
    worker.start()
    # abort only once the run is parked in the wait, otherwise get_flats
    # clears the abort flag right after it is set
    assert waiting.wait(5)
    flat.abort()

    assert done.wait(5), "abort did not interrupt the sun window wait"


def test_a_second_run_is_refused_while_one_is_in_progress(tmp_path, coefficients_file):
    """The scheduler has been seen forking four concurrent flat programs
    onto one camera."""
    flat, _ = build(tmp_path, coefficients_file)
    started, release = threading.Event(), threading.Event()

    original = flat._take_image

    def blocking_take_image(*args, **kwargs):
        started.set()
        release.wait(5)
        return original(*args, **kwargs)

    flat._take_image = blocking_take_image
    worker = threading.Thread(
        target=lambda: flat.get_flats("CLEAR", n_flats=1), daemon=True
    )
    worker.start()
    assert started.wait(5)

    with pytest.raises(SkyFlatRunningException):
        flat.get_flats("CLEAR", n_flats=1)

    release.set()
    worker.join(timeout=10)


def test_the_scope_is_left_untracked_at_the_end(tmp_path, coefficients_file):
    flat, fakes = build(tmp_path, coefficients_file)

    flat.get_flats("CLEAR", n_flats=2)

    assert not fakes["telescope"].tracking


def test_the_wheel_is_not_moved_when_it_is_already_in_place(
    tmp_path, coefficients_file
):
    """set_filter also applies the filter focus offset, so a redundant call
    costs a focuser move too."""
    flat, fakes = build(tmp_path, coefficients_file)

    flat.get_flats("CLEAR", n_flats=2)

    assert fakes["wheel"].moves == []


def test_a_failed_slew_is_fatal(tmp_path, coefficients_file):
    """Frames taken from wherever the telescope happened to be are not
    flats: the old code logged the failure and exposed anyway."""
    flat, fakes = build(tmp_path, coefficients_file, flat_position_max=0)

    def broken_slew(alt, az):
        raise OSError("mount not responding")

    fakes["telescope"].slew_to_alt_az = broken_slew

    with pytest.raises(Exception, match="flat position"):
        flat.get_flats("CLEAR", n_flats=1)


def test_the_scope_is_not_re_slewed_within_flat_position_max(
    tmp_path, coefficients_file
):
    flat, fakes = build(tmp_path, coefficients_file, flat_position_max=1)
    fakes["telescope"].alt, fakes["telescope"].az = flat["flat_alt"], flat["flat_az"]

    flat.get_flats("CLEAR", n_flats=2)

    assert fakes["telescope"].slews == 0


def test_the_flat_position_check_is_a_great_circle_distance(
    tmp_path, coefficients_file
):
    """Position.angsep() reads its pair as (ra, dec), so on an alt/az
    Position it calls two points 2 degrees apart at the zenith 180 degrees
    apart - and every frame would re-slew."""
    flat, fakes = build(
        tmp_path, coefficients_file, flat_alt=89, flat_az=0, flat_position_max=3
    )
    telescope = fakes["telescope"]

    telescope.alt, telescope.az = 89, 180  # 2 degrees away over the pole
    assert flat._at_flat_position(telescope)

    telescope.alt, telescope.az = 85, 180  # 6 degrees away
    assert not flat._at_flat_position(telescope)
