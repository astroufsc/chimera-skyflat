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
                  "camera": "/Camera/0",
                  "filterwheel": "/FilterWheel/0",
                  "site": "/Site/0",
                  "idealCounts": 25000 # this should be a detector property
                  }

    def getFlats(self,
                 sunInitialZD = -5.0, sunFinalZD = -30.0,
                 filter = 'R',  initialExptime = 1.0,
                 pierSide = 'E',
                 skypar = (2000000., 68., 17.),
                 debug=False):
        """
        Takes sequence of flats. Start taking one frame to determine current level,
        then predicts next exposure time based on exponential decay of sky brightness,
        creates a list of sunZD, intensity.  It should have the right exponential behavior.
        If not exponential raise some flag about sky condition.

        :param sunInitialZD: Initial sun altitude
        :param sunFinalZD: Final sun altitude
        :param filter:
        :param initialExptime:
        :param pierSide: Which side of pier telescope should be positioned?
        :param skypar:
        :param debug:
        :return:
        """

    def getSkyLevel(self, exptime=30, debug=False):
        """
        Takes one sky flat to get current sky counts
        """

    def nextFlatET(self, debug=False):
        """
        Computes exposure time for next exposure
        """

    def setSideOfPier(self, sideOfPier="E", debug=False):
        """
        Puts telescope on right side of pier
        """

    # """
    # @event
    # def pointComplete(self, position, star, frame):
    #     """Raised after every step in the focus sequence with
    #     information about the last step.
    #     """
    # """
