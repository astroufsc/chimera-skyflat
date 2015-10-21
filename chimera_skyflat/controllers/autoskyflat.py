from __future__ import division
import os
from astropy.io import fits
from chimera.controllers.imageserver.util import getImageServer
from chimera.core.exceptions import ChimeraException
from chimera.interfaces.camera import Shutter
from chimera.util.image import ImageUtil
import time
from chimera_skyflat.interfaces.autoskyflat import IAutoSkyFlat

import ntpath

import numpy as np

__author__ = 'kanaan'

from chimera.core.chimeraobject import ChimeraObject

class AutoSkyFlat(ChimeraObject, IAutoSkyFlat):

    # normal constructor
    # initialize the relevant variables
    def __init__(self):
        ChimeraObject.__init__(self)
        self.sideOfPier = "E"
        self.filter="R"
        self.sunInitialZD=3
        self.sunFinalZD=10
        self.defaultExptime=30


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
        frames = cam.expose(exptime=exptime, frames=1, shutter=Shutter.OPEN,
                            filename=os.path.basename(ImageUtil.makeFilename("skyflat-$DATE")))

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


    def _moveScope(self):
        """
        Moves the scope, usually to zenith
        """
        tel = self._getTel()
        # tel.slewToHADec(Position.fromRaDec(ra, dec))


    def getFlats(self, debug=False):
        """
        DOCME
        :param sunInitialZD:
        :param sunFinalZD:
        :return:
        """
        # 1 - Wait for the Sun ZD = sunInitialZD
        # 2 - Take first image to guess exponential scaling factor
        # 3 - Measure exponential scaling factor
        # 4 - Compute exposure time
        # 5 - Take flat
        # 6 - Goto 2.
        site = self._getSite()
        pos = site.sunpos()
        self.log.debug('Sun altitude is %f' % pos.alt)
        self.log.debug('Taking image...')
        sky_level = self.getSkyLevel(exptime=self.defaultExptime)
        self.log.debug('Mean: %f'% sky_level)
        #filename, image = self._takeImage(exptime=initialExptime, filter=filter)
        #frame = fits.getdata(filename)
        #self.log.debug('Mean: %f'% np.mean(frame))


    def getSkyLevel(self, exptime=self.defaultExptime, debug=False):
        """
        Takes one sky flat and returns average
        Need to get rid of deviant points, median is best, but takes too long K???
        """

        self.setSideOfPier("E") # self.sideOfPier)
        # self.move

        try:
            filename, image = self._takeImage(exptime=exptime, filter=self.filter)
            frame = fits.getdata(filename)
            img_mean = np.mean(frame)

            return(img_mean)
        except:
            self.log.error("Can't take image")
            raise





