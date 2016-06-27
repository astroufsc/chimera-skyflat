from __future__ import division
import copy
import json
import ntpath
import os
import re
import threading
import time
from datetime import timedelta
import numpy as np
from astropy.io import fits
from chimera.controllers.imageserver.imagerequest import ImageRequest
from chimera.controllers.imageserver.util import getImageServer
from chimera.core.event import event
from chimera.core.exceptions import ChimeraException, ProgramExecutionAborted
from chimera.interfaces.camera import Shutter
from chimera.interfaces.telescope import TelescopePier
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

    def _getDome(self):
        return self.getManager().getProxy(self["dome"])

    def _getCam(self):
        return self.getManager().getProxy(self["camera"])

    def _getFilterWheel(self):
        return self.getManager().getProxy(self["filterwheel"])

    def _takeImage(self, exptime, filter, download=False):

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
            if download and not os.path.exists(image_path):  # If image is on a remote server, donwload it.

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

    def _moveScope(self, tracking=False, pierSide=None):
        """
        Moves the scope, usually to zenith
        """
        tel = self._getTel()
        site = self._getSite()
        self.log.debug('Moving scope to alt %s az %s.' % (self["flat_alt"], self["flat_az"]))
        if tel.getPositionAltAz().angsep(
                Position.fromAltAz(Coord.fromD(self["flat_alt"]), Coord.fromD(self["flat_az"]))).D < self["flat_position_max"]:

            self.log.debug(
                'Telescope is less than {} degrees from flat position. Not moving!'.format(self["flat_position_max"]))
            if tracking and not tel.isTracking():
                tel.startTracking()
            elif not tracking and tel.isTracking():
                tel.stopTracking()
            if pierSide is not None and tel.features(TelescopePier):
                self.log.debug("Setting telescope pier side to %s." % tel.getPierSide().__str__().lower())
                tel.setSideOfPier(self['pier_side'])

            return

        try:
            self.log.debug("Skyflat Slewing scope to alt {} az {}".format(self["flat_alt"], self["flat_az"]))
            tel.slewToRaDec(Position.altAzToRaDec(Position.fromAltAz(Coord.fromD(self["flat_alt"]),
                                                                     Coord.fromD(self["flat_az"])),
                                                  site['latitude'], site.LST()))
            if tracking:
                self._startTracking()
            else:
                self._stopTracking()
        except:
            self.log.debug("Error moving the telescope")

    def _stopTracking(self):
        """
        disables telescope tracking
        """
        tel = self._getTel()
        try:
            self.log.debug("Skyflat is stopping telescope tracking")
            tel.stopTracking()
        except:
            self.log.debug("Error stopping the telescope")

    def _startTracking(self, wait=True):
        """
        enables telescope tracking
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

    def getFlats(self, filter_id, n_flats=None):
        """
        Take flats on filter_id filter.

        * 1 - Wait for the Sun to enter the altitude strip where we can take skyflats
        * 2 - Take first image to guess exponential scaling factor
        * 3 - Measure exponential scaling factor
        * 4 - Compute exposure time
        * 5 - Take flat
        * 6 - Goto 2 until reaches n_flats

        :param filter_id: Filter name to take the Flats
        :param n_flats: Number of flats to take. None for maximum on the sun interval.
        """

        # Read fresh coefficients from file.
        self.scale, self.slope, self.bias = self.readCoefficientsFile(self['coefficients_file'])[filter_id]
        self.log.debug('Skyflat parameters: n_flats = %i, filter = %s, scale = %i, slope = %i, bias = %i' % (
        n_flats, filter_id, self.scale, self.slope, self.bias))

        self._abort.clear()

        site = self._getSite()
        pos = site.sunpos()
        self.log.debug(
            'Starting sky flats Sun altitude is {}. max: {} min: {}'.format(pos.alt.D, self["sun_alt_hi"], self["sun_alt_low"]))

        self.log.debug('Starting dome track.')
        self._getDome().track()

        self.log.debug('Moving to filter {}.'.format(filter_id))
        if self["filterwheel"] is not None:
            fw = self._getFilterWheel()
            fw.setFilter(filter_id)

        # Wait with the telescope in the flat position.
        self._moveScope(tracking=False, pierSide=self['pier_side'])

        # while the Sun is above or below the flat field strip we just wait
        while pos.alt.D > self["sun_alt_hi"] or pos.alt.D < self["sun_alt_low"]:

            # Check if position is outside and return.
            # dusk
            if site.localtime > 12 and pos.alt.D < self["sun_alt_low"]:
                self.log('Finishing flats. Sun position below than {}'.format(self["sun_alt_low"]))
                return
            # dawn
            elif site.localtime < 12 and pos.alt.D > self["sun_alt_hi"]:
                self.log('Finishing flats. Sun position higher than {}'.format(self["sun_alt_low"]))
                return

            # checking for aborting signal
            if self._abort.isSet():
                self.log.warning('Aborting!')
                self._getTel().stopTracking()
                return

            time.sleep(5)
            pos = site.sunpos()
            self.log.debug('Sun altitude is %f waiting to be between %f and %f' % (pos.alt.D, self["sun_alt_hi"], self["sun_alt_low"]))

        pos = site.sunpos()

        # take flats until the Sun is out of the skyflats altitude or n_flats is reached
        i_flat = 0
        correction_factor = 0  #
        while self["sun_alt_hi"] > pos.alt.D > self["sun_alt_low"]:
            if i_flat == n_flats:
                self.log.debug('Done %i flats on filter %s' % (i_flat, filter_id))
                self._stopTracking()
                return

            self.log.debug("Initial positions {} {} {}".format(pos.alt.D, self["sun_alt_hi"], self["sun_alt_low"]))
            self._moveScope(tracking=self["tracking"])  # Go to the skyflat pos and shoot!
            aux = self.computeSkyFlatTime(correction_factor)
            if aux:
                expTime, sky_level_expected = aux
            else:
                self._stopTracking()
                return

            if expTime > 0:
                self.log.debug('Taking sky flat image with exptime = %f' % expTime)
                filename, image = self._takeImage(exptime=expTime, filter=filter_id, download=True)
                i_flat += 1

                sky_level = self.getSkyLevel(filename, image)
                self.exposeComplete(filter_id, i_flat, expTime, sky_level)
                correction_factor += sky_level / expTime - sky_level_expected / expTime
                self.log.debug('Done taking image, average counts = %f. '
                               'New correction factor = %f' % (sky_level, correction_factor))

            else:
                # dusk
                if site.localtime().hour > 12:
                    self.log.debug('Exposure time too low. Waiting 5 seconds.')
                    time.sleep(5)
                # dawn
                else:
                    self.log.debug('Exposure time too low. Finishing this filter...')
                    self._stopTracking()
                    return


            # checking for aborting signal
            if self._abort.isSet():
                self.log.warning('Aborting!')
                self._getTel().stopTracking()
                return

            pos = site.sunpos()
            self.log.debug("{} {} {}".format(pos.alt.D,self["sun_alt_hi"],self["sun_alt_low"]))

    def computeSkyFlatTime(self, correction_factor):
        """
        :param correction_factor: Additive correction factor for the sky counts exponential

        Method returns exposureTime
        This computation requires self.scale, self.slope and self.bias defined.
        """

        site = self._getSite()
        intCounts = 0.0
        exposure_time = 0
        initialTime = site.ut()
        sun_altitude = site.sunpos().alt
        while 1:
            sky_counts = (self.expArg(sun_altitude.R, self.scale, self.slope, self.bias) + correction_factor) * self[
                'exptime_increment']
            initialTime = initialTime + timedelta(seconds=float(self["exptime_increment"]))
            sun_altitude = site.sunpos(initialTime).alt

            if intCounts + sky_counts >= self["idealCounts"]:
                self.log.debug("Breaking the Exposure Time Calculation loop. Sun Altitude {} Exposure time {} "
                               "Computed counts {}".format(sun_altitude.D, exposure_time, intCounts))
                return float(exposure_time), intCounts

            exposure_time += float(self["exptime_increment"])
            intCounts += sky_counts
            if exposure_time > self["exptime_max"]:
                self.log.debug("Computed exposure: Sun Altitude {} time {} sky_counts = {}"
                               " intCounts = {}".format(sun_altitude.D, initialTime, sky_counts, intCounts))
                # dusk
                if site.localtime().hour > 12:
                    self.log.warning(
                        "Exposure time exceeded limit of {}. Finishing this filter...".format(self["exptime_max"]))
                    return False
                else:  # dawn
                    self.log.warning(
                        "Exposure time exceeded limit of {}. Waiting 5 sec...".format(self["exptime_max"]))
                    time.sleep(5)

        return float(exposure_time), intCounts

    def getSkyLevel(self, filename, image):
        """
        Returns average counts from image
        """

        frame = fits.getdata(filename)
        img_mean = np.mean(frame)
        return img_mean

    def abort(self):
        self._abort.set()
        cam = copy.copy(self._getCam())
        cam.abortExposure()

    # @event
    # def exposeBegin(self):
    #     '''
    #     Called on beginning of an expose
    #     '''

    @event
    def exposeComplete(self, filter_id, i_flat, expTime, sky_level):
        '''
        Called on exposuse completion
        '''

    @staticmethod
    def readCoefficientsFile(filename):
        with open(filename) as f:
            coefficients = json.loads(re.sub('#(.*)', '', f.read()))
        return coefficients

    @staticmethod
    def expArg(x, Scale, Slope, Bias):
        return Scale * np.exp(Slope * x) + Bias
