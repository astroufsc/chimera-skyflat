from __future__ import division
import copy

import ntpath
import os
import threading
import time
from datetime import timedelta

import numpy as np
from astropy.io import fits
from chimera.controllers.imageserver.imagerequest import ImageRequest
from chimera.controllers.imageserver.util import getImageServer
from chimera.core.exceptions import ChimeraException, ProgramExecutionAborted
from chimera.interfaces.camera import Shutter
from chimera.util.coord import Coord
from chimera.util.image import ImageUtil
from chimera.util.position import Position

from chimera_skyflat.interfaces.autoskyflat import IAutoSkyFlat

__author__ = 'kanaan'

from chimera.core.chimeraobject import ChimeraObject

class SkyFlatMaxExptimeReached(ChimeraException):
    """
    Raised when exposure time is longer than exptime_max. See config.
    """



class AutoSkyFlat(ChimeraObject, IAutoSkyFlat):
    # normal constructor
    # initialize the relevant variables
    def __init__(self):
        ChimeraObject.__init__(self)
        self._abort = threading.Event()
        self._abort.clear()

    def _getSite(self):
        return self.getManager().getProxy(self["site"])

    def _getTel(self):
        return self.getManager().getProxy(self["telescope"])

    def _getCam(self):
        return self.getManager().getProxy(self["camera"])

    def _getFilterWheel(self):
        return self.getManager().getProxy(self["filterwheel"])

    def _takeImage(self, exptime, filter):

        cam = self._getCam()
        if self["filterwheel"] is not None:
            fw = self._getFilterWheel()
            fw.setFilter(filter)
        self.log.debug("Start frame")
        request = ImageRequest(exptime=exptime, frames=1, shutter=Shutter.OPEN,
                            filename=os.path.basename(ImageUtil.makeFilename("skyflat-$DATE-$TIME")),
                            type='sky-flat')
        self.log.debug('ImageRequest: {}'.format(request))
        frames = cam.expose(request)
        self.log.debug("End frame")

        # checking for aborting signal
        if self._abort.isSet():
            self.log.warning('Aborting exposure!')
            raise ProgramExecutionAborted()

        if frames:
            image = frames[0]
            image_path = image.filename()
            if not os.path.exists(image_path):  # If image is on a remote server, donwload it.

                #  If remote is windows, image_path will be c:\...\image.fits, so use ntpath instead of os.path.
                if ':\\' in image_path:
                    modpath = ntpath
                else:
                    modpath = os.path
                image_path = ImageUtil.makeFilename(os.path.join(getImageServer(self.getManager()).defaultNightDir(),
                                                                 modpath.basename(image_path)))
                t0 = time.time()
                self.log.debug('Downloading image from server to %s' % image_path)
                if not ImageUtil.download(image, image_path):
                    raise ChimeraException('Error downloading image %s from %s' % (image_path, image.http()))
                self.log.debug('Finished download. Took %3.2f seconds' % (time.time() - t0))
            return image_path, image
        else:
            raise Exception("Could not take an image")

    def _moveScope(self, tracking=False):
        """
        Moves the scope, usually to zenith
        """
        tel = self._getTel()
        if tel.getPositionAltAz().angsep(Position.fromAltAz(Coord.fromD(self["flat_alt"]),
                                                            Coord.fromD(self["flat_az"]))).D > self["flat_position_max"]:

            self.log.debug('Telescope is less than {} degrees'.format(self["flat_position_max"]))

        try:
            self.log.debug("Skyflat Slewing scope to alt {} az {}".format(self["flat_alt"], self["flat_az"]))
            tel.slewToAltAz(Position.fromAltAz(self["flat_alt"], self["flat_az"]))
            if tracking:
                self._startTracking()
            else:
                self._stopTracking()
        except:
            self.log.debug("Error moving the telescope")

    def _stopTracking(self):
        """
        stops tracking the telescope
        """
        tel = self._getTel()
        try:
            self.log.debug("Skyflat is stopping telescope tracking")
            tel.stopTracking()
        except:
            self.log.debug("Error stopping the telescope")

    def _startTracking(self, wait=True):
        """
        Starts telescope tracking.
        """
        tel = self._getTel()
        try:
            self.log.debug("Skyflat is restarting telescope tracking")
            tel.startTracking()
            if wait:
                while 1:
                    if tel.isTracking():
                        return
        except:
            self.log.debug("Error starting the telescope")

    def getFlats(self, debug=False):
        """
        DOCME
        :param debug:
        :type self: object
        :param sun_alt_hi:
        :param sun_alt_low:
        :return:
        """
        # 1 - Wait for the Sun to enter the altitude strip where we can take skyflats
        # 2 - Take first image to guess exponential scaling factor
        # 3 - Measure exponential scaling factor
        # 4 - Compute exposure time
        # 5 - Take flat
        # 6 - Goto 2.

        self._abort.clear()

        site = self._getSite()
        pos = site.sunpos()
        self.log.debug(
            'Starting sky flats Sun altitude is {} {} {}'.format(pos.alt.D, self["sun_alt_hi"], self["sun_alt_low"]))

        # Wait with the telescope in the flat position.
        self._moveScope(tracking=False)

        # while the Sun is above or below the flat field strip we just wait
        while pos.alt.D > self["sun_alt_hi"] or pos.alt.D < self["sun_alt_low"]:

            # checking for aborting signal
            if self._abort.isSet():
                self.log.warning('Aborting!')
                return

            # maybe we should test for pos << self.sun_alt_hi K???
            time.sleep(10)
            pos = site.sunpos()
            self.log.debug('Sun altitude is %f waiting to be between %f and %f' % (
            pos.alt.D, self["sun_alt_hi"], self["sun_alt_low"]))

        # now the Sun is in the skyflat altitude strip
        # self._moveScope(tracking=self["tracking"])  # now that the skyflat time arrived move the telescope
        # self.log.debug('Taking image...')
        # sky_level = self.getSkyLevel(exptime=float(self["exptime_default"]))
        # self.log.debug('Mean: %f' % sky_level)
        pos = site.sunpos()
        # self._moveScope(tracking=self["tracking"])  # now that the skyflat time arrived move the telescope


        while self["sun_alt_hi"] > pos.alt.D > self["sun_alt_low"]:  # take flats until the Sun is out of the skyflats regions
            self.log.debug("Initial positions {} {} {}".format(pos.alt.D, self["sun_alt_hi"], self["sun_alt_low"]))
            self._moveScope(tracking=self["tracking"])  # Go to the skyflat pos and shoot!
            expTime = self.computeSkyFlatTime() #sky_level)

            if expTime > 0:
                self.log.debug('Taking sky flat image with exptime = %f' % expTime)
                sky_level = self.getSkyLevel(exptime=expTime)
                self.log.debug('Done taking image, average counts = %f' % sky_level)
            else:
                self.log.debug('Exposure time too low. Waiting 5 seconds.')
                time.sleep(5)

            # checking for aborting signal
            if self._abort.isSet():
                self.log.warning('Aborting!')
                return

            pos = site.sunpos()
            self.log.debug("{} {} {}".format(pos.alt.D,self["sun_alt_hi"],self["sun_alt_low"]))

    def computeSkyFlatTime(self): #, sky_level):
        """
        User specifies:
        :param sky_level - current level of sky counts
        :param pos.alt - initial Sun sun_altitude
        :return exposure time needed to reach sky_level
        Method returns exposureTime
        This computation requires Scale, Slope, Bias
        Scale and Slope are filter dependent, I think they should be part of some enum type variable
        """

        site = self._getSite()
        intCounts = 0.0
        exposure_time = 0
        initialTime = site.ut() #+ timedelta(seconds=20)
        sun_altitude = site.sunpos().alt
        while 1:
            sky_counts = self.expArg(sun_altitude.R, self["Scale"], self["Slope"], self["Bias"])
            initialTime = initialTime + timedelta(seconds=1)
            sun_altitude = site.sunpos(initialTime).alt
            self.log.debug(
                "Computing exposure time {} sky_counts = {} intCounts = {}".format(initialTime, sky_counts, intCounts))
            if intCounts + sky_counts >= self["idealCounts"]:
                self.log.debug("Breaking the Exposure Time Calculation loop. Exposure time {} Computed counts {}"
                               .format(exposure_time,intCounts))
                return float(exposure_time)
            exposure_time += 1
            intCounts += sky_counts
            if exposure_time > self["exptime_max"]:
                self.log.warning("Exposure time exceeded limit of {}".format(self["exptime_max"]))
                return self["exptime_max"]
                # raise SkyFlatMaxExptimeReached("Exposure time is too long")

        return float(exposure_time)

    def expTime(self, sky_level, altitude):
        """
        Computes sky level for current altitude
        Assumes sky_lve
        :param sky_level:
        :param altitude:
        :return:
        """

    def getSkyLevel(self, exptime, debug=False):
        """
        Takes one sky flat and returns average
        Need to get rid of deviant points, mode is best, but takes too long, using mean for now K???
        For now just using median which is almost OK
        """

        self.setSideOfPier("E")  # self.sideOfPier)
        try:
            filename, image = self._takeImage(exptime=exptime, filter=self["filter"])
            frame = fits.getdata(filename)
            img_mean = np.mean(frame)
            image.close()
            return (img_mean)
        except ProgramExecutionAborted:
            return
        except:
            self.log.error("Can't take image")
            raise

    def abort(self):
        self._abort.set()
        cam = copy.copy(self._getCam())
        cam.abortExposure()

    def expArg(self, x, Scale, Slope, Bias):
        return Scale * np.exp(Slope * x) + Bias





