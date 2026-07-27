# SPDX-FileCopyrightText: 2015-present Antonio Kanaan <kanaan@astro.ufsc.br>
# SPDX-License-Identifier: GPL-2.0-or-later

import datetime as dt
import json
import os
import re
import threading
import time

import numpy as np
from astropy.io import fits
from chimera.controllers.imageserver.imagerequest import ImageRequest
from chimera.core.chimeraobject import ChimeraObject
from chimera.core.event import event
from chimera.core.exceptions import ChimeraException, ProgramExecutionAborted
from chimera.interfaces.camera import Shutter
from chimera.interfaces.telescope import TelescopePierSide
from chimera.util.coord import CoordUtil
from chimera.util.image import Image, ImageUtil

from chimera_skyflat.interfaces.autoskyflat import IAutoSkyFlat

__author__ = "kanaan"

# The sun's altitude is sampled twice, this far apart, and then extrapolated
# linearly over the exposure-time integration below. Sampling it at every
# integration step instead costs one bus round-trip per 0.2 s of modelled
# exposure (~1500 for a 300 s frame); over a 5 min horizon the linear error
# is ~0.01 deg, i.e. ~1% in predicted counts.
SUN_TRACK_BASELINE = 60.0

# Bounds for the multiplicative sky-model correction. A night that needs
# more than this is not a correction, it is the wrong coefficients file.
MODEL_GAIN_MIN = 0.02
MODEL_GAIN_MAX = 50.0

# Largest single-frame correction accepted before damping, so one cloud, one
# satellite trail or one bad readout cannot throw the whole set.
MAX_FRAME_RATIO = 4.0

# How close to ideal_counts a frame has to land for the model to count as
# calibrated - from there on corrections are damped.
CONVERGED_BAND = 0.25


class SkyFlatMaxExptimeReached(ChimeraException):
    """
    Raised when exposure time is longer than exptime_max. See config.
    """


class SkyFlatRunningException(ChimeraException):
    """
    Raised when get_flats is called while another run is in progress.
    """


class AutoSkyFlat(ChimeraObject, IAutoSkyFlat):
    def __init__(self):
        ChimeraObject.__init__(self)
        self._abort = threading.Event()
        self._abort.clear()
        # one run at a time: the same camera cannot serve two flat sets, and
        # the scheduler has been seen forking duplicate programs
        self._run_lock = threading.Lock()
        self.scale = self.slope = self.bias = None
        # the coefficients file, read once per run (see _get_flats)
        self._coefficients = None

    #
    # proxies
    #

    def _get_tel(self):
        return self.get_proxy(self["telescope"])

    def _get_cam(self):
        return self.get_proxy(self["camera"])

    def _get_filter_wheel(self):
        return self.get_proxy(self["filterwheel"])

    def _get_site(self):
        return self.get_proxy(self["site"])

    #
    # sun
    #

    @staticmethod
    def _altitude_in_degrees(sunpos):
        """Sun altitude in DEGREES from whatever Site.sunpos() returned.

        Older cores returned ``(alt, az)`` in degrees; the current one
        returns a ``Position`` whose ``.alt`` is degrees but which also
        unpacks into ``Coord``s. Guessing wrong is what pinned every
        twilight flat at exptime_max on 2026-07-21, and feeding the Coord
        straight to ``np.radians()`` raises TypeError.
        """
        altitude = getattr(sunpos, "alt", None)
        if altitude is None:
            altitude = sunpos[0] if isinstance(sunpos, tuple | list) else sunpos
        return float(altitude)

    def _sun_track(self):
        """(altitude now [deg], rate [deg/s]) from two sun positions.

        The exposure-time calculator needs the rate anyway - it integrates
        the sky forward over the frame - and it doubles as the dusk/dawn
        test: ``rate < 0`` is the sun setting. Local clock hours cannot do
        that job, they are wrong under the fast-forward simulation clock
        (Site.time_speedup), which advances the modelled sky but not the
        wall clock.
        """
        site = self._get_site()
        now = site.ut()
        alt_now = self._altitude_in_degrees(site.sunpos(now))
        alt_later = self._altitude_in_degrees(
            site.sunpos(now + dt.timedelta(seconds=SUN_TRACK_BASELINE))
        )
        return alt_now, (alt_later - alt_now) / SUN_TRACK_BASELINE

    #
    # instrument helpers
    #

    def _set_filter(self, filter_id):
        """Move the wheel only when it is not already on ``filter_id``.

        The controller has to drive the wheel itself: ImageRequest has no
        filter key (it rejects unknown ones), so the camera cannot do it.
        set_filter also applies the configured focus offset, so a redundant
        call costs a focuser move as well as a wheel move.
        """
        if self["filterwheel"] is None:
            return
        fw = self._get_filter_wheel()
        try:
            if fw.get_filter() == filter_id:
                return
        except Exception:
            self.log.debug("Could not read the current filter, setting it anyway.")
        fw.set_filter(filter_id)

    def _take_image(self, exptime, filter_id, request=None):
        cam = self._get_cam()
        self._set_filter(filter_id)

        imrequest = ImageRequest(
            exptime=exptime,
            frames=1,
            shutter=Shutter.OPEN,
            filename=os.path.basename(ImageUtil.make_filename("skyflat-$DATE-$TIME")),
            type="sky-flat",
            compress_format=self["compress_format"],
        )
        if request is not None:
            imrequest.update(request)

        self.log.debug(f"Starting frame. ImageRequest: {imrequest}")
        frames = cam.expose(imrequest)
        self.log.debug("End frame")

        if self._abort.is_set():
            self.log.warning("Aborting exposure!")
            raise ProgramExecutionAborted()

        if not frames:
            raise ChimeraException("Could not take an image")

        image = Image.from_url(frames[0])
        if image is None:
            raise ChimeraException(f"Could not open the image at {frames[0]}")

        if not os.path.exists(image.filename):
            # the image lives on a remote image server: fetch it, we need
            # the pixels to measure the sky level
            t0 = time.time()
            self.log.debug(f"Downloading image from server to {image.filename}")
            if not ImageUtil.download(image):
                raise ChimeraException(
                    f"Error downloading image {image.filename} from {image.http()}"
                )
            self.log.debug(f"Finished download. Took {time.time() - t0:3.2f} seconds")

        return image.filename, image

    def _move_scope(self, tracking=False):
        """Point the scope at the flat field position.

        A failed slew is fatal: flats taken from wherever the telescope
        happened to be are not flats.
        """
        tel = self._get_tel()

        if self["pier_side"] is not None and tel.features("TelescopePier"):
            try:
                side = self["pier_side"]
                if isinstance(side, str):
                    side = TelescopePierSide(side.upper())
                self.log.debug(f"Setting telescope pier side to {side}.")
                tel.set_pier_side(side)
            except Exception:
                self.log.exception("Could not set the pier side, going on without it.")

        if self._at_flat_position(tel):
            self.log.debug(
                f"Telescope is within {self['flat_position_max']} degrees of the "
                f"flat position. Not moving."
            )
        else:
            self.log.debug(
                f"Slewing scope to alt {self['flat_alt']} az {self['flat_az']}."
            )
            try:
                tel.slew_to_alt_az(float(self["flat_alt"]), float(self["flat_az"]))
            except Exception as e:
                raise ChimeraException(
                    f"Could not slew to the flat position "
                    f"(alt {self['flat_alt']} az {self['flat_az']}): {e}"
                ) from e

        if tracking:
            self._start_tracking()
        else:
            self._stop_tracking()

    def _at_flat_position(self, tel):
        """True when the scope is already within flat_position_max."""
        if not self["flat_position_max"]:
            return False
        try:
            alt, az = tel.get_position_alt_az()
            # gcdist, not Position.angsep: angsep reads the pair as
            # (ra, dec), so it comes out wrong for an alt/az Position -
            # 180 deg for two points 2 deg apart at the zenith
            separation = np.degrees(
                CoordUtil.gcdist(
                    (np.radians(float(az)), np.radians(float(alt))),
                    (
                        np.radians(float(self["flat_az"])),
                        np.radians(float(self["flat_alt"])),
                    ),
                )
            )
            return separation < float(self["flat_position_max"])
        except Exception:
            self.log.debug("Could not read the telescope position, slewing anyway.")
            return False

    def _stop_tracking(self):
        try:
            self.log.debug("Skyflat is stopping telescope tracking")
            self._get_tel().stop_tracking()
        except Exception:
            self.log.exception("Error stopping telescope tracking")

    def _start_tracking(self, wait=True, timeout=30.0):
        try:
            self.log.debug("Skyflat is starting telescope tracking")
            tel = self._get_tel()
            tel.start_tracking()
            if not wait:
                return
            # bounded and cancellable: the old unbounded busy-wait hammered
            # the bus and could never be aborted
            deadline = time.time() + timeout
            while time.time() < deadline:
                if tel.is_tracking():
                    return
                if self._abort.wait(0.5):
                    return
            self.log.warning(f"Telescope did not report tracking within {timeout} s.")
        except Exception:
            self.log.exception("Error starting telescope tracking")

    #
    # sky model
    #

    def _read_coefficients(self):
        """The coefficients file, read once per run.

        It is a fit over past nights, it does not change under us, and a
        run re-reads it from scratch (_get_flats drops the cache) so a
        refit lands on the next set without a restart.
        """
        if self._coefficients is None:
            self._coefficients = self.read_coefficients_file(self["coefficients_file"])
        return self._coefficients

    def _load_coefficients(self, filter_id):
        """Point the sky model at ``filter_id``'s scale/slope/bias."""
        coefficients = self._read_coefficients()
        if filter_id not in coefficients:
            raise ChimeraException(
                f"No sky brightness coefficients for filter {filter_id} in "
                f"{self['coefficients_file']} (have: "
                f"{', '.join(sorted(coefficients))})"
            )
        self.scale, self.slope, self.bias = (
            float(value) for value in coefficients[filter_id][:3]
        )
        return self.scale, self.slope, self.bias

    def _sky_rate(self, sun_alt_degrees, model_gain=1.0):
        """Modelled sky count rate (counts/s) at a sun altitude in DEGREES.

        The fit lives in radians (a slope of ~58 only makes sense there),
        while every altitude in this file - the config limits included - is
        in degrees. The conversion used to run the wrong way, which
        underflowed exp() and ran every flat to exptime_max (2026-07-21).

        ``model_gain`` scales the sky term only: the bias is the model's
        floor, not part of the night-to-night sky brightness.
        """
        return self.exp_arg(
            np.radians(sun_alt_degrees), self.scale * model_gain, self.slope, self.bias
        )

    #
    # main entry point
    #

    def get_flats(self, filter_id, n_flats=None, request=None):
        """
        Take flats on filter_id filter.

        * 1 - Wait for the Sun to enter the altitude strip where we can take skyflats
        * 2 - Compute the exposure time from the sky brightness model
        * 3 - Take a flat and measure it
        * 4 - Correct the model with the measurement
        * 5 - Goto 2 until n_flats is reached or the Sun leaves the strip

        :param filter_id: Filter name to take the Flats
        :param n_flats: Number of flats to take. None for maximum on the sun interval.
        :param request: Additional keywords to pass to ImageRequest.
        :return: Number of flats taken.
        """
        if not self._run_lock.acquire(blocking=False):
            raise SkyFlatRunningException(
                "A sky flat run is already in progress on this controller."
            )
        try:
            self._abort.clear()
            return self._get_flats(filter_id, n_flats, request)
        finally:
            self._run_lock.release()

    def _get_flats(self, filter_id, n_flats, request):
        self._coefficients = None  # a run always starts from a fresh read
        self._load_coefficients(filter_id)
        self.log.debug(
            f"Skyflat parameters: n_flats = {n_flats}, filter = {filter_id}, "
            f"scale = {self.scale}, slope = {self.slope}, bias = {self.bias}"
        )

        # move the wheel first: it is the slowest thing between here and the
        # first frame, and the sky is not waiting
        self._set_filter(filter_id)
        self._move_scope(tracking=False)

        if not self._wait_for_sun_window():
            return 0

        i_flat = 0
        model_gain = 1.0
        last_exptime = None
        calibrated = False
        tried_filters = {filter_id}

        while True:
            sun_alt, sun_rate = self._sun_track()
            dusk = sun_rate < 0
            if not self["sun_alt_low"] < sun_alt < self["sun_alt_hi"]:
                self.log.info(
                    f"Sun altitude {sun_alt:.2f} left the flat window "
                    f"[{self['sun_alt_low']}, {self['sun_alt_hi']}]. "
                    f"Done {i_flat} flats on filter {filter_id}."
                )
                break

            if i_flat == n_flats:
                self.log.debug(f"Done {i_flat} flats on filter {filter_id}")
                break

            if self._abort.is_set():
                self.log.warning("Aborting!")
                break

            computed = self.compute_sky_flat_time(model_gain)

            if computed is False:
                # exptime_max reached. At dusk the sky only gets fainter, so
                # this filter is done - but a MORE SENSITIVE one can still
                # reach ideal_counts in time. Walk up the coefficients by
                # scale (V -> R -> CLEAR) instead of ending the whole
                # sequence, which used to abandon twilight with filters to
                # spare. At dawn compute_sky_flat_time waits instead.
                next_filter = self._next_filter(filter_id, tried_filters, dusk=dusk)
                if next_filter is None:
                    break
                self.log.info(
                    f"{filter_id} needs more than {self['exptime_max']} s; "
                    f"switching to {next_filter}."
                )
                filter_id = next_filter
                tried_filters.add(next_filter)
                self._load_coefficients(filter_id)
                self._set_filter(filter_id)
                model_gain, last_exptime, calibrated = 1.0, None, False
                continue

            exptime, expected_counts = computed

            if exptime < float(self["exptime_min"]):
                # sky too bright for this filter
                if dusk:
                    self.log.debug("Exposure time too low. Waiting 5 seconds.")
                    if self._abort.wait(5):
                        break
                    continue
                # at dawn it only gets brighter: step DOWN in sensitivity
                # (CLEAR -> R -> V) instead of ending the sequence.
                next_filter = self._next_filter(filter_id, tried_filters, dusk=False)
                if next_filter is None:
                    self.log.info(
                        "Exposure time too low and no less sensitive filter "
                        "left. Finishing this filter..."
                    )
                    break
                self.log.info(
                    f"{filter_id} saturates in under {self['exptime_min']} s; "
                    f"switching to {next_filter}."
                )
                filter_id = next_filter
                tried_filters.add(next_filter)
                self._load_coefficients(filter_id)
                self._set_filter(filter_id)
                model_gain, last_exptime, calibrated = 1.0, None, False
                continue

            exptime, expected_counts = self._clamp_exptime(
                exptime, expected_counts, last_exptime
            )

            self.log.debug(
                f"Taking sky flat image with exptime = {exptime:.2f} s "
                f"(expecting {expected_counts:.0f} counts, gain {model_gain:.3f})"
            )
            self._move_scope(tracking=self["tracking"])
            filename, image = self._take_image(exptime, filter_id, request=request)
            try:
                sky_level = self.get_sky_level(filename, image)
            finally:
                if image is not None:
                    image.close()
            last_exptime = exptime

            if self["max_counts"] and sky_level > self["max_counts"]:
                # a saturated frame is not a flat: do not publish it, do not
                # count it, just let the correction below shorten the next one
                self.log.warning(
                    f"Sky flat measured {sky_level:.0f} counts, above "
                    f"max_counts = {self['max_counts']}: discarding the frame."
                )
            else:
                i_flat += 1
                self.expose_complete(filter_id, i_flat, exptime, sky_level)

            model_gain, calibrated = self._correct_model(
                model_gain, sky_level, expected_counts, calibrated=calibrated
            )
            self.log.debug(
                f"Done taking image, average counts = {sky_level:.1f} "
                f"(expected {expected_counts:.1f}). New model gain = {model_gain:.3f}"
            )

        self._stop_tracking()
        return i_flat

    def _wait_for_sun_window(self):
        """Hold at the flat position until the sun is inside the strip.

        Returns False when the window is gone for this twilight (or the run
        was aborted), True when it is time to shoot.
        """
        while True:
            sun_alt, sun_rate = self._sun_track()
            dusk = sun_rate < 0
            if self["sun_alt_low"] < sun_alt < self["sun_alt_hi"]:
                return True

            # the window is over, not yet to come: at dusk once the sun is
            # below the strip, at dawn once it is above it
            if dusk and sun_alt < self["sun_alt_low"]:
                self.log.info(
                    f"Finishing flats. Sun altitude {sun_alt:.2f} is below "
                    f"{self['sun_alt_low']}"
                )
                return False
            if not dusk and sun_alt > self["sun_alt_hi"]:
                self.log.info(
                    f"Finishing flats. Sun altitude {sun_alt:.2f} is above "
                    f"{self['sun_alt_hi']}"
                )
                return False

            if self._abort.is_set():
                self.log.warning("Aborting!")
                return False

            self.log.debug(
                f"Sun altitude is {sun_alt:.2f}, waiting for it to be between "
                f"{self['sun_alt_low']} and {self['sun_alt_hi']}"
            )
            if self._abort.wait(5):
                self.log.warning("Aborting!")
                return False

    def _clamp_exptime(self, exptime, expected_counts, last_exptime):
        """Bound the frame-to-frame exposure change.

        Even a well behaved correction should not jump: one bad measurement
        (a cloud, a passing satellite) would otherwise put the next frame
        near full well. The expected counts follow the clamp, so the
        feedback below stays consistent with the frame actually taken.
        """
        step = float(self["max_exptime_step"] or 0)
        if not last_exptime or step <= 1:
            return exptime, expected_counts

        clamped = min(max(exptime, last_exptime / step), last_exptime * step)
        if clamped == exptime:
            return exptime, expected_counts

        self.log.debug(
            f"Clamping exposure time {exptime:.2f} -> {clamped:.2f} s "
            f"(at most x{step} from the last frame)"
        )
        return clamped, expected_counts * clamped / exptime

    def _correct_model(self, model_gain, sky_level, expected_counts, calibrated=True):
        """Fold one measurement into the sky model, gently.

        The old correction was an ADDITIVE count-rate offset applied with
        unit gain. Against a sky dimming by a factor of a few per frame it
        overcorrected every time and rang: 15.7k -> 26.3k -> 40.8k -> 15.5k
        -> 37.6k -> 52.1k counts around a 25k target on 2026-07-22, with
        frame 6 at ~80% of full well.

        The sky's shape in sun altitude is what the coefficients file
        fitted; what varies night to night is its normalisation. So correct
        that instead - multiplicatively - and damp it, which turns the
        oscillation into a geometric convergence.

        Damping is for noise, not for acquisition: while the model has not
        yet put a frame near the target its error IS the normalisation
        error, so take it whole. Damping switches on with the first frame
        that lands inside CONVERGED_BAND, and stays on.

        :return: ``(model_gain, calibrated)``
        """
        if expected_counts <= 0 or sky_level <= 0:
            return model_gain, calibrated

        ratio = float(
            np.clip(sky_level / expected_counts, 1 / MAX_FRAME_RATIO, MAX_FRAME_RATIO)
        )
        on_target = abs(ratio - 1) <= CONVERGED_BAND
        calibrated = calibrated or on_target

        damping = float(self["correction_damping"]) if calibrated else 1.0
        corrected = model_gain * ratio**damping
        return float(np.clip(corrected, MODEL_GAIN_MIN, MODEL_GAIN_MAX)), calibrated

    def compute_sky_flat_time(self, model_gain=1.0):
        """Exposure time that reaches ideal_counts, and the counts expected.

        :param model_gain: Multiplicative correction on the sky model, 1.0
            for the coefficients file as fitted.
        :return: ``(exposure_time, expected_counts)``, or False when the sky
            cannot deliver ideal_counts within exptime_max.

        Requires self.scale, self.slope and self.bias (see
        _load_coefficients).
        """
        increment = float(self["exptime_increment"])
        exptime_max = float(self["exptime_max"])
        ideal_counts = float(self["ideal_counts"])
        n_wait_iter = 0

        while True:
            sun_alt, sun_rate = self._sun_track()
            dusk = sun_rate < 0
            counts = 0.0
            exposure_time = 0.0

            while exposure_time <= exptime_max:
                # the sun keeps moving during the exposure, so integrate the
                # rate forward instead of freezing it at the start
                rate = self._sky_rate(sun_alt + sun_rate * exposure_time, model_gain)
                if counts + rate * increment >= ideal_counts:
                    # finish on a partial step: rounding the exposure down to
                    # a whole increment costs up to exptime_increment of sky,
                    # which is the whole frame when the sky is bright
                    exposure_time += (ideal_counts - counts) / rate
                    self.log.debug(
                        f"Exposure time {exposure_time:.2f} s reaches "
                        f"{ideal_counts:.0f} counts at sun altitude {sun_alt:.2f}"
                    )
                    return exposure_time, ideal_counts
                counts += rate * increment
                exposure_time += increment

            self.log.debug(
                f"Computed exposure: sun altitude {sun_alt:.2f}, counts {counts:.0f}"
            )
            if dusk:
                # dusk: it only gets worse from here
                self.log.warning(
                    f"Computed exposure time {exposure_time:.2f} exceeded the limit "
                    f"of {exptime_max}. Finishing this filter..."
                )
                return False
            if n_wait_iter >= self["max_wait_iter"]:
                self.log.warning(
                    "Maximum number of wait iterations reached. Giving up."
                )
                return False

            # dawn: the sky is brightening, wait for it
            self.log.info(
                f"Computed exposure time {exposure_time:.2f} exceeded the limit of "
                f"{exptime_max}. Waiting 6 sec..."
            )
            n_wait_iter += 1
            if self._abort.wait(6):
                self.log.warning("Aborting!")
                return False

    def get_sky_level(self, filename, image):
        """
        Returns average counts from image
        """
        frame = fits.getdata(filename)
        # plain float: numpy.float64 is not msgspec-serializable, so leaking
        # it into the expose_complete event silently drops the publication
        return float(np.mean(frame))

    def abort(self):
        """Stop the running set as soon as the current frame allows."""
        self._abort.set()
        try:
            self._get_cam().abort_exposure()
        except Exception:
            self.log.exception("Error aborting the current exposure")

    @event
    def expose_complete(self, filter_id, i_flat, exp_time, sky_level):
        """
        Called on exposure completion
        """

    def read_coefficients_file(self, filename):
        with open(os.path.expanduser(filename)) as f:
            coefficients = json.loads(re.sub("#(.*)", "", f.read()))
        return coefficients

    def _next_filter(self, filter_id, tried, dusk=True):
        """The next filter to try when this one cannot reach ideal_counts.

        Sensitivity is the model's own scale term, so the coefficients file
        already ranks the filters (HBETA < V = I < R < B < CLEAR). At dusk
        the sky is fading and we need the next MORE sensitive filter (the
        shortest exposure above the current one); at dawn it is flooding
        and we need the next LESS sensitive one. Returns None once nothing
        is left in that direction.
        """
        if not self["filter_fallback"]:
            return None
        coefficients = self._read_coefficients()
        current = coefficients.get(filter_id)
        if current is None:
            return None

        candidates = [
            (values[0], name)
            for name, values in coefficients.items()
            if name not in tried
            and (values[0] > current[0] if dusk else values[0] < current[0])
        ]
        if not candidates:
            return None
        return min(candidates)[1] if dusk else max(candidates)[1]

    def exp_arg(self, x, scale, slope, bias):
        """Sky brightness model: counts/s at a sun altitude in RADIANS."""
        return scale * np.exp(slope * x) + bias
