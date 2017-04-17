from chimera.core import SYSTEM_CONFIG_DIRECTORY
from chimera.interfaces.telescope import TelescopePierSide

__author__ = 'kanaan'

#! /usr/bin/env python
# -*- coding: iso-8859-1 -*-

# chimera - observatory automation system
# Copyright (C) 2006-2007  P. Henrique Silva <henrique@astro.ufsc.br>

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.


from chimera.core.interface import Interface
from chimera.core.event import event

from chimera.util.enum import Enum
from chimera.core.exceptions import ChimeraException


class CantPointScopeException(ChimeraException):
    pass


class CanSetScopeButNotThisField(ChimeraException):
    pass


class CantSetScopeException(ChimeraException):
    pass

Target = Enum("CURRENT", "AUTO")


class IAutoSkyFlat(Interface):

    __config__ = {"telescope": "/Telescope/0",
                  "dome": "/Dome/0",
                  "camera": "/Camera/0",
                  "filterwheel": "/FilterWheel/0",
                  "site": "/Site/0",
                  "tracking": True,                     # Enable telescope tracking when exposing?
                  "flat_position_max": 1,               # If telescope less than flat_position_max, it does not move prior to expose. (degrees)
                  "flat_alt": 89,                       # Skyflat position - Altitude. (degrees)
                  "flat_az": 78,                        # Skyflat position - Azimuth. (degrees)
                  "pier_side": None,                    # Pier Side to take Skyflat
                  "sun_alt_hi": -5,                     # Lowest Sun Altitude to make Skyflats. (degrees)
                  "sun_alt_low": -30,                   # Highest Sun Altitude to make Skyflats. (degrees)
                  "exptime_increment": 0.2,             # Exposure time increment on integration. (seconds)
                  "exptime_max": 300,                   # Maximum exposure time. (seconds)
                  "max_wait_iter": 100,                 # Maximum number of iterations on wait loop
                  "idealCounts": 25000,                 # Ideal flat CCD counts.
                  "coefficients_file": "%s/skyflats.json" % SYSTEM_CONFIG_DIRECTORY,
                  "compress_format": "NO"
                  }

    def getFlats(self, filter_id, n_flats):
        """
        Takes sequence of flats, starts taking one frame to determine current level
        Then predicts next exposure time based on exponential decay of sky brightness
        Creates a list of sunZD, intensity.  It should have the right exponential behavior.
        If not exponential raise some flag about sky condition.
        """

    def getSkyLevel(self, filename, image):
        """
        Returns average level from image
        """