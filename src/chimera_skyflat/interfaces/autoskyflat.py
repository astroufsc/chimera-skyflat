# SPDX-FileCopyrightText: 2015-present Antonio Kanaan <kanaan@astro.ufsc.br>
# SPDX-License-Identifier: GPL-2.0-or-later

from chimera.core import SYSTEM_CONFIG_DIRECTORY
from chimera.core.event import event
from chimera.core.exceptions import ChimeraException
from chimera.core.interface import Interface
from chimera.util.enum import Enum

__author__ = "kanaan"


class CantPointScopeException(ChimeraException):
    pass


class CanSetScopeButNotThisField(ChimeraException):
    pass


class CantSetScopeException(ChimeraException):
    pass


class Target(Enum):
    CURRENT = "CURRENT"
    AUTO = "AUTO"


class IAutoSkyFlat(Interface):
    __config__ = {
        "telescope": "/Telescope/0",
        "dome": "/Dome/0",
        "camera": "/Camera/0",
        "filterwheel": "/FilterWheel/0",
        "site": "/Site/0",
        # Enable telescope tracking when exposing?
        "tracking": True,
        # Skip the slew when the telescope is already this close to the flat
        # position (degrees). 0 to slew before every frame.
        "flat_position_max": 1,
        # Skyflat altitude. The azimuth is not configurable: flats are shot
        # at the anti-solar point, where the twilight gradient is smallest.
        # Altitude is a site decision - horizon, dome slit, mount limits -
        # so it stays a knob; arXiv:1407.8283 puts the null point at 75.
        "flat_alt": 75,
        # Pier side to take Skyflats on: "EAST", "WEST" or None to leave it
        # to the telescope.
        "pier_side": None,
        # Highest Sun altitude to make Skyflats. (degrees)
        "sun_alt_hi": -5,
        # Lowest Sun altitude to make Skyflats. (degrees)
        "sun_alt_low": -30,
        # Exposure time increment on integration. (seconds)
        "exptime_increment": 0.2,
        # Maximum exposure time. (seconds)
        "exptime_max": 300,
        # Shortest usable exposure, e.g. the camera's own minimum. Below it
        # the sky counts as too bright: wait at dusk, step to a less
        # sensitive filter at dawn. (seconds)
        "exptime_min": 0.2,
        # Maximum number of iterations on the dawn wait loop
        "max_wait_iter": 100,
        # Ideal flat CCD counts.
        "ideal_counts": 25000,
        # Discard (do not publish, do not count) frames above this level,
        # e.g. the detector's saturation. 0 disables the check.
        "max_counts": 0,
        # How hard each measurement corrects the sky model: 1.0 applies the
        # whole measured ratio, 0.5 its square root. Below 1 the correction
        # converges instead of ringing around ideal_counts.
        "correction_damping": 0.5,
        # Largest exposure time ratio between consecutive frames. 0 or 1
        # disables the clamp.
        "max_exptime_step": 2.0,
        # When a filter cannot reach ideal_counts, move to the next more
        # sensitive filter at dusk / less sensitive at dawn instead of
        # ending the sequence.
        "filter_fallback": True,
        "coefficients_file": f"{SYSTEM_CONFIG_DIRECTORY}/skyflat_coefficients.json",
        "compress_format": "NO",
    }

    def get_flats(self, filter_id, n_flats, request):
        """
        Takes sequence of flats, starts taking one frame to determine current level
        Then predicts next exposure time based on exponential decay of sky brightness
        Creates a list of sunZD, intensity.  It should have the right exponential behavior.
        If not exponential raise some flag about sky condition.
        """

    def get_sky_level(self, filename, image):
        """
        Returns average level from image
        """

    def abort(self):
        """
        Aborts the current sky flat sequence, aborting the running exposure.
        """

    @event
    def expose_complete(self, filter_id, i_flat, exp_time, sky_level):
        """
        Called on exposure completion, once per flat frame actually kept.
        """
