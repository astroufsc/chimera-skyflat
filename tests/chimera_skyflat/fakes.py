# SPDX-FileCopyrightText: 2015-present Antonio Kanaan <kanaan@astro.ufsc.br>
# SPDX-License-Identifier: GPL-2.0-or-later

"""Deterministic stand-ins for the site and the instruments.

The sky here is the same exponential the controller models, evaluated on a
clock the camera advances as it exposes, so a whole twilight can be replayed
in a test: the sun sets, the sky dims during each frame, and the measured
counts come out of an actual FITS file.
"""

import datetime as dt
import threading

import msgspec
import numpy as np
from astropy.io import fits

# CLEAR/R/V-like coefficients, same shape as skyflat_coefficients.json
COEFFICIENTS = {
    "CLEAR": [415567.8440142421, 57.8620007247163, 100.0],
    "B": [39974.82501620678, 72.48851256446889, 100.0],
    "R": [17860.060118418703, 42.581608832490545, 100.0],
    "V": [9886.372329954358, 39.29958402578991, 100.0],
    "HBETA": [2441.9884591901236, 52.532639916292844, 100.0],
}


class FakeClock:
    """A clock the fake camera pushes forward as it exposes."""

    def __init__(self, start=None):
        self.now = start or dt.datetime(2026, 7, 22, 20, 55, tzinfo=dt.UTC)

    def advance(self, seconds):
        self.now += dt.timedelta(seconds=seconds)


class ClockEvent(threading.Event):
    """The controller's abort flag, wired to the fake clock.

    Every wait in the controller is `self._abort.wait(seconds)`; here that
    moves the modelled sky forward by those seconds instead of the wall
    clock, so a run waiting for the sky to dim makes progress in a test.
    """

    def __init__(self, clock):
        super().__init__()
        self.clock = clock

    def wait(self, timeout=None):
        if timeout:
            self.clock.advance(timeout)
        return self.is_set()


class FakeSite:
    """Sun altitude linear in time, which is what twilight looks like."""

    def __init__(self, clock, alt0=-2.0, rate=-0.004, azimuth=300.0):
        self.clock = clock
        self.alt0 = alt0
        self.rate = rate
        self.azimuth = azimuth
        self._t0 = clock.now
        self.sunpos_calls = 0

    def ut(self):
        return self.clock.now

    def altitude(self, date=None):
        date = date or self.clock.now
        return self.alt0 + self.rate * (date - self._t0).total_seconds()

    # -- astroufsc/chimera#275/#282 float accessors ---------------------
    # No sunpos() here on purpose: it returns a Position, which the bus
    # cannot encode, so the controller must never call it. A fake without
    # it makes any such call an AttributeError in the tests.

    def sun_altitude(self, date=None):
        self.sunpos_calls += 1
        return self.altitude(date)

    def sun_azimuth(self, date=None):
        return self.azimuth

    def is_dusk(self, date=None):
        return self.rate < 0


class FakeSky:
    """The real sky: the model's shape, scaled by ``factor``.

    ``factor`` is what a coefficients file fitted on another night gets
    wrong, and what the controller's correction has to find.
    """

    def __init__(self, site, coefficients, filter_id="CLEAR", factor=1.0):
        self.site = site
        self.coefficients = coefficients
        self.filter_id = filter_id
        self.factor = factor

    def rate(self, altitude):
        scale, slope, bias = self.coefficients[self.filter_id]
        return self.factor * scale * np.exp(slope * np.radians(altitude)) + bias

    def counts(self, start, exptime, steps=200):
        """Counts collected between ``start`` and ``start + exptime``."""
        step = exptime / steps
        total = 0.0
        for i in range(steps):
            when = start + dt.timedelta(seconds=i * step)
            total += self.rate(self.site.altitude(when)) * step
        return total


class FakeFilterWheel:
    def __init__(self, filters=tuple(COEFFICIENTS), current="CLEAR"):
        self.filters = list(filters)
        self.current = current
        self.moves = []

    def get_filter(self):
        return self.current

    def set_filter(self, filter_id):
        self.moves.append(filter_id)
        self.current = filter_id

    def get_filters(self):
        return tuple(self.filters)


class FakeTelescope:
    def __init__(self, alt=89.0, az=78.0):
        self.alt, self.az = alt, az
        self.tracking = False
        self.slews = 0

    def features(self, interface):
        return False

    def get_position_alt_az(self):
        return self.alt, self.az

    def slew_to_alt_az(self, alt, az):
        self.slews += 1
        self.alt, self.az = alt, az

    def start_tracking(self):
        self.tracking = True

    def stop_tracking(self):
        self.tracking = False

    def is_tracking(self):
        return self.tracking


class FakeCamera:
    """Writes a real FITS frame whose mean is the integrated sky."""

    def __init__(self, clock, sky, tmp_path, readout=0.0, full_well=65535):
        self.clock = clock
        self.sky = sky
        self.tmp_path = tmp_path
        self.readout = readout
        self.full_well = full_well
        self.exposures = []

    def expose(self, request):
        # the bus encodes every request; a numpy scalar in one of them takes
        # out the exposure and every metadata call that follows it
        msgspec.json.encode(dict(request))

        exptime = float(request["exptime"])
        counts = min(
            self.sky.counts(self.clock.now, exptime),
            self.full_well,
        )
        self.exposures.append((exptime, counts))
        self.clock.advance(exptime + self.readout)

        path = self.tmp_path / f"flat-{len(self.exposures):03d}.fits"
        fits.PrimaryHDU(np.full((8, 8), counts, dtype=np.float32)).writeto(path)
        return (f"file://{path}",)

    def abort_exposure(self, readout=True):
        pass
