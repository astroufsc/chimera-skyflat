chimera-skyflat plugin
======================

This is a template plugin for the chimera observatory control system
https://github.com/astroufsc/chimera.

Makes skyflatting easier and less stressful :-)

Usage
-----

Install ``chimera-skyflat``, then configure it on chimera.config and create a json file with the exponential paramenters
on, i.e., ``~/.chimera/skyflats.json``.

Running chimera-skyflats script:
::

    Usage: chimera-skyflat [options]

    Chimera - Observatory Automation System - SkyFlats

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -v, --verbose         Display information while working
      -q, --quiet           Don't display information while working.
                            [default=True]

      skyFlats:
        --auto              Does a sequence of sky flats
        --skyflat=SKYFLAT   Auto Sky Flats
        --sun-high=FILE     Highest Sun altitude
        -n NUMBER, --number=NUMBER
                            Number of skyflats to take on the filter
        -f FILTER, --filter=FILTER
                            Skyflat filter name
        --sun-low=FILE      Lowest Sun altitude

      Client Configuration:
        --config=CONFIG     Chimera configuration file to use.
                            default=/Users/william/.chimera/chimera.config
                            [default=/Users/william/.chimera/chimera.config]
        -P PORT, --port=PORT
                            Port to which the local Chimera instance will listen
                            to. [default=9000]

      Object Paths:
        -C PATH, --controllers-dir=PATH
                            Append PATH to controllers load path. This option
                            could be setted multiple times to add multiple
                            directories. [default=['/Users/william/.virtualenvs/ch
                            imera/lib/python2.7/site-
                            packages/chimera/controllers',
                            '/Users/william/.virtualenvs/chimera/lib/python2.7
                            /site-packages/chimera_pverify/controllers',
                            '/Users/william/.virtualenvs/chimera/lib/python2.7
                            /site-packages/chimera_skyflat/controllers']]

* Example:

::

    $ chimera-skyflat -f R,I -n 3 --auto --sun-hi 0 --sun-lo -12
    Pointing scope to the zenith and waiting for the Sun to reach skyflats altitude range
    Filter: R, Skyflat # 1 - Exposure time: 2.20 seconds - Counts: 1000.00
    Filter: R, Skyflat # 2 - Exposure time: 16.60 seconds - Counts: 1000.00
    Filter: R, Skyflat # 3 - Exposure time: 15.20 seconds - Counts: 1000.00
    Filter: I, Skyflat # 1 - Exposure time: 0.40 seconds - Counts: 1000.00
    Filter: I, Skyflat # 2 - Exposure time: 5.40 seconds - Counts: 1000.00
    Filter: I, Skyflat # 3 - Exposure time: 5.60 seconds - Counts: 1000.00
    Finished.

takes 3 skyflats on filters ``R`` and ``I`` if the sun is between 0 a -12 degrees of altitude.


Installation
------------

::

    pip install -U git+https://github.com/astroufsc/chimera-skyflat.git


Configuration Example
---------------------

Configuration example to be added on ``chimera.config`` file:

::

    controllers:
        - type: AutoSkyFlat
          name: autoskyflat
          site: /Site/0
          telescope: /Telescope/0
          camera: /Camera/0
          filterwheel: /FilterWheel/0
          dome: /Dome/0
          tracking: True                     # Enable telescope tracking when exposing?
          flat_position_max: 1               # If telescope less than flat_position_max it does not move prior to expose. (degrees)
          flat_alt: 89                       # Skyflat position - Altitude. (degrees)
          flat_az: 78                        # Skyflat position - Azimuth. (degrees)
          pier_side: TelescopePierSide.EAST  # Pier Side to take Skyflat
          sun_alt_hi: -5                     # Lowest Sun Altitude to make Skyflats. (degrees)
          sun_alt_low: -30                   # Highest Sun Altitude to make Skyflats. (degrees)
          exptime_increment: 0.2             # Exposure time increment on integration. (seconds)
          exptime_max: 300                   # Maximum exposure time. (seconds)
          idealCounts: 25000                 # Ideal flat CCD counts.
          coefficients_file: ~/.chimera/skyflats.json


``skyflats.json`` file example:

::

    {
        "U": [16500, 70, -14],
        "G": [32002478, 97, 195],
        "R": [355328, 44, 108],
        "I": [41222293, 94, -68],
        "Z": [5985164, 85, 106]
    }

The coefficients on the list are Scale, Slope and Bias from the equation:

``counts_per_sec = scale * exp(slope * sun_altitude) + bias``

Contact
-------

For more information, contact us on chimera's discussion list:
https://groups.google.com/forum/#!forum/chimera-discuss

Bug reports and patches are welcome and can be sent over our GitHub page:
https://github.com/astroufsc/chimera-skyflats/
