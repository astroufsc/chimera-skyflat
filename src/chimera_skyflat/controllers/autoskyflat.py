import copy
import json
import os
import re
import threading
import time
from datetime import timedelta
from chimera.util.image import Image
import numpy as np
from astropy.io import fits
from chimera.controllers.imageserver.imagerequest import ImageRequest
from chimera.controllers.imageserver.util import get_image_server
from chimera.core.event import event
from chimera.core.exceptions import ChimeraException, ProgramExecutionAborted
from chimera.interfaces.camera import Shutter
from chimera.interfaces.telescope import TelescopePier
from chimera.util.coord import Coord
from chimera.util.image import ImageUtil
from chimera.util.position import Position
from chimera_skyflat.interfaces.autoskyflat import IAutoSkyFlat

__author__ = "kanaan"

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

    def _get_tel(self):
        return self.get_proxy(self["telescope"])

    def _get_cam(self):
        return self.get_proxy(self["camera"])

    def _get_filter_wheel(self):
        return self.get_proxy(self["filterwheel"])

    def _get_site(self):
        return self.get_proxy("/Site/0")
    
    def _get_dome(self):
        return self.get_proxy(self["dome"])

    def _take_image(self, exptime, filter, download=False, request=None):

        cam = self._get_cam()
        if self["filterwheel"] is not None:
            fw = self._get_filter_wheel()
            fw.set_filter(filter)
        self.log.debug("Start frame")
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
        self.log.debug(f"ImageRequest: {imrequest}")
        frames = cam.expose(imrequest)
        self.log.debug("End frame")

        # checking for aborting signal
        if self._abort.is_set():
            self.log.warning("Aborting exposure!")
            raise ProgramExecutionAborted()

        if frames:
            image = Image.from_url(frames[0])
            if not os.path.exists(image.filename):  # If image is on a remote server, donwload it.

                # #  If remote is windows, image_path will be c:\...\image.fits, so use ntpath instead of os.path.
                # if ':\\' in image_path:
                #     modpath = ntpath
                # else:
                #     modpath = os.path
                # image_path = ImageUtil.make_filename(os.path.join(self["images_dir"], "$LAST_NOON_DATE", modpath.basename(image_path)))
                t0 = time.time()
                self.log.debug(f"Downloading image from server to {image.filename}")
                if not image.download():
                    raise ChimeraException(f"Error downloading image {image.filename} from {image.http()}")
                self.log.debug(f"Finished download. Took {time.time() - t0:3.2f} seconds")
            return image.filename, image
        else:
            raise Exception("Could not take an image")

    def _move_scope(self, tracking=False, pierSide=None):
        """
        Moves the scope, usually to zenith
        """
        tel = self._get_tel()
        site = self._get_site()
        self.log.debug(f"Moving scope to alt {self['flat_alt']} az {self['flat_az']}.")
        # if (
        #     tel.get_position_alt_az()
        #     .angsep(Position.from_alt_az(Coord.from_d(self["flat_alt"]), Coord.from_d(self["flat_az"])))
        #     .degrees
        #     < self["flat_position_max"]
        # ):

        #     self.log.debug(
        #         f"Telescope is less than {self['flat_position_max']} degrees from flat position. Not moving!"
        #     )
        #     if tracking and not tel.isTracking():
        #         tel.startTracking()
        #     elif not tracking and tel.isTracking():
        #         tel.stopTracking()
        #     if pierSide is not None and tel.features(TelescopePier):
        #         self.log.debug(f"Setting telescope pier side to {tel.getPierSide().__str__().lower()}.")
        #         tel.setSideOfPier(self["pier_side"])

        #     return

        try:
            self.log.debug(f"Skyflat Slewing scope to alt {self['flat_alt']} az {self['flat_az']}")
            tel.slew_to_alt_az(self["flat_alt"], self["flat_az"])
            if tracking:
                self._start_tracking()
            else:
                self._stop_tracking()
        except:
            self.log.debug("Error moving the telescope")

    def _stop_tracking(self):
        """
        disables telescope tracking
        """
        tel = self._get_tel()
        try:
            self.log.debug("Skyflat is stopping telescope tracking")
            tel.stop_tracking()
        except:
            self.log.debug("Error stopping the telescope")

    def _start_tracking(self, wait=True):
        """
        enables telescope tracking
        """
        tel = self._get_tel()
        try:
            self.log.debug("Skyflat is restarting telescope tracking")
            tel.start_tracking()
            if wait:
                while True:
                    if tel.is_tracking():
                        return
        except:
            self.log.debug("Error starting the telescope")

    def get_flats(self, filter_id, n_flats=None, request=None):
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
        :param request: Additional keywords to pass to ImageRequest.
        """

        # Read fresh coefficients from file.
        print(self["coefficients_file"] )
        self.scale, self.slope, self.bias = self.read_coefficients_file(self["coefficients_file"])[filter_id]
        self.log.debug(
            f"Skyflat parameters: n_flats = {n_flats}, filter = {filter_id}, scale = {self.scale}, slope = {self.slope}, bias = {self.bias}"
        )

        self._abort.clear()

        site = self._get_site()
        sun_alt, _ = site.sunpos()
        self.log.debug(
            f"Starting sky flats Sun altitude is {sun_alt}. max: {self['sun_alt_hi']} min: {self['sun_alt_low']}"
        )

        # self.log.debug("Starting dome track.")
        # self._get_dome().track()

        self.log.debug(f"Moving to filter {filter_id}.")
        if self["filterwheel"] is not None:
            fw = self._get_filter_wheel()
            fw.set_filter(filter_id)

        # Wait with the telescope in the flat position.
        self._move_scope(tracking=False, pierSide=self["pier_side"])

        # while the Sun is above or below the flat field strip we just wait
        while sun_alt > self["sun_alt_hi"] or sun_alt < self["sun_alt_low"]:

            # Check if position is outside and return.
            # dusk
            if site.localtime().hour > 12 and sun_alt < self["sun_alt_low"]:
                self.log.debug(f"Finishing flats. Sun position below than {self['sun_alt_low']}")
                return
            # dawn
            elif site.localtime().hour < 12 and sun_alt > self["sun_alt_hi"]:
                self.log.debug(f"Finishing flats. Sun position higher than {self['sun_alt_low']}")
                return

            # checking for aborting signal
            if self._abort.is_set():
                self.log.warning("Aborting!")
                self._get_tel().stop_tracking()
                return

            time.sleep(5)
            sun_alt, _ = site.sunpos()
            self.log.debug(
                f"Sun altitude is {sun_alt} waiting to be between {self['sun_alt_hi']} and {self['sun_alt_low']}"
            )

        sun_alt, _ = site.sunpos()

        # take flats until the Sun is out of the skyflats altitude or n_flats is reached
        i_flat = 0
        correction_factor = 0  #
        tried_filters = {filter_id}
        while self["sun_alt_hi"] > sun_alt > self["sun_alt_low"]:
            if i_flat == n_flats:
                self.log.debug(f"Done {i_flat} flats on filter {filter_id}")
                self._stop_tracking()
                return

            self.log.debug(f"Initial positions {sun_alt} {self['sun_alt_hi']} {self['sun_alt_low']}")
            self._move_scope(tracking=self["tracking"])  # Go to the skyflat pos and shoot!
            aux = self.compute_sky_flat_time(correction_factor)
            if aux:
                expTime, sky_level_expected = aux
            else:
                # At dusk the sky only gets fainter, so this filter is done -
                # but a MORE SENSITIVE one can still reach ideal_counts in
                # under exptime_max. Walk up the coefficients by scale
                # (V -> R -> CLEAR) instead of ending the whole sequence,
                # which used to abandon the twilight with filters to spare.
                next_filter = self._next_more_sensitive_filter(filter_id, tried_filters)
                if site.localtime().hour > 12 and next_filter is not None:
                    self.log.info(
                        f"{filter_id} needs more than {self['exptime_max']} s; "
                        f"switching to {next_filter}."
                    )
                    tried_filters.add(next_filter)
                    filter_id = next_filter
                    self.scale, self.slope, self.bias = self.read_coefficients_file(
                        self["coefficients_file"]
                    )[filter_id]
                    if self["filterwheel"] is not None:
                        self._get_filter_wheel().set_filter(filter_id)
                    correction_factor = 0
                    continue
                self._stop_tracking()
                return

            if expTime > 0:
                self.log.debug(f"Taking sky flat image with exptime = {expTime}")
                filename, image = self._take_image(exptime=expTime, filter=filter_id, download=True, request=request)
                i_flat += 1

                sky_level = self.get_sky_level(filename, image)
                self.expose_complete(filter_id, i_flat, expTime, sky_level)
                correction_factor += sky_level / expTime - sky_level_expected / expTime
                self.log.debug(
                    f"Done taking image, average counts = {sky_level}. "
                    f"New correction factor = {correction_factor}"
                )

            else:
                # dusk
                if site.localtime().hour > 12:
                    self.log.debug("Exposure time too low. Waiting 5 seconds.")
                    time.sleep(5)
                # dawn
                else:
                    self.log.debug("Exposure time too low. Finishing this filter...")
                    self._stop_tracking()
                    return

            # checking for aborting signal
            if self._abort.isSet():
                self.log.warning("Aborting!")
                self._getTel().stopTracking()
                return

            sun_alt, _ = site.sunpos()
            self.log.debug(f"{sun_alt} {self['sun_alt_hi']} {self['sun_alt_low']}")

    def compute_sky_flat_time(self, correction_factor):
        """
        :param correction_factor: Additive correction factor for the sky counts exponential

        Method returns exposureTime
        This computation requires self.scale, self.slope and self.bias defined.
        """

        site = self._get_site()
        intCounts = 0.0
        exposure_time = 0
        initial_time = site.ut()
        sun_alt, _ = site.sunpos()
        n_wait_iter = 0
        while 1:
            # The sky model wants the sun altitude in RADIANS (the fitted
            # slope ~58 only makes sense there), but sunpos() returns
            # DEGREES - as sun_alt_hi/sun_alt_low and every other test in
            # this file assume. The old "* 57.30" was the inverse
            # conversion, so the model was evaluated at -471 instead of
            # -0.14 rad: exp() underflowed, sky_counts collapsed to the
            # bias floor and EVERY flat ran to exptime_max regardless of
            # sun altitude (seen live 2026-07-21: 249.8 s twilight flats).
            sky_counts = (
                self.exp_arg(np.radians(sun_alt), self.scale, self.slope, self.bias) + correction_factor
            ) * self["exptime_increment"]
            initial_time = initial_time + timedelta(seconds=float(self["exptime_increment"]))
            sun_alt, _ = site.sunpos(initial_time)
            if intCounts + sky_counts >= self["ideal_counts"]:
                self.log.debug(
                    f"Breaking the Exposure Time Calculation loop. Sun Altitude {sun_alt} Exposure time {exposure_time} "
                    f"Computed counts {intCounts}"
                )
                return float(exposure_time), intCounts

            exposure_time += float(self["exptime_increment"])
            intCounts += sky_counts
            if exposure_time > self["exptime_max"]:
                self.log.debug(
                    f"Computed exposure: Sun Altitude {sun_alt} time {initial_time} sky_counts = {sky_counts}"
                    f" intCounts = {intCounts}"
                )
                # dusk
                if site.localtime().hour > 12:
                    self.log.warning(
                        f"Computed exposure time {exposure_time} exceeded limit of {self['exptime_max']}. "
                        f"Finishing this filter..."
                    )
                    return False
                elif n_wait_iter > self["max_wait_iter"]:
                    self.log.warning("Maximum number of wait iterations reached. Giving up.")
                    return False
                else:  # dawn
                    self.log.info(
                        f"Computed exposure time {exposure_time} exceeded limit of {self['exptime_max']}. Waiting 6 sec..."
                    )
                    time.sleep(6)
                    # Reset initalTime and exposure_time
                    intCounts = 0.0
                    exposure_time = 0
                    initial_time = site.ut()
                    sun_alt, _ = site.sunpos(initial_time)
                    n_wait_iter += 1

            if self._abort.is_set():
                self.log.warning("Aborting!")
                self._get_tel().stop_tracking()
                return False

        return float(exposure_time), intCounts

    def get_sky_level(self, filename, image):
        """
        Returns average counts from image
        """

        frame = fits.getdata(filename)
        img_mean = np.mean(frame)
        return img_mean

    def abort(self):
        self._abort.set()
        cam = copy.copy(self._get_cam())
        cam.abort_exposure()

    # @event
    # def exposeBegin(self):
    #     '''
    #     Called on beginning of an expose
    #     '''

    @event
    def expose_complete(self, filter_id, i_flat, expTime, sky_level):
        """
        Called on exposuse completion
        """

    # @staticmethod
    def read_coefficients_file(self, filename):
        with open(os.path.expanduser(filename)) as f:
            coefficients = json.loads(re.sub("#(.*)", "", f.read()))
        return coefficients

    def _next_more_sensitive_filter(self, filter_id, tried):
        """The filter needing the SHORTEST exposure of those still untried.

        Sensitivity is the model's own scale term, so the coefficients file
        already ranks the filters: HBETA < V = I < R < B < CLEAR. Returns
        None once nothing brighter is left.
        """
        try:
            coefficients = self.read_coefficients_file(self["coefficients_file"])
        except (OSError, ValueError):
            self.log.warning("could not re-read the coefficients file; not switching filter")
            return None

        current = coefficients.get(filter_id)
        if current is None:
            return None

        candidates = [
            (values[0], name)
            for name, values in coefficients.items()
            if name not in tried and values[0] > current[0]
        ]
        return min(candidates)[1] if candidates else None

    # @staticmethod
    def exp_arg(self, x, Scale, Slope, Bias):
        return Scale * np.exp(Slope * x) + Bias
