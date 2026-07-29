#!/usr/bin/env python
# SPDX-FileCopyrightText: 2015-present Antonio Kanaan <kanaan@astro.ufsc.br>
# SPDX-License-Identifier: GPL-2.0-or-later

import copy
import sys

from chimera.cli.cli import ChimeraCLI, action


class ChimeraAutoSkyFlat(ChimeraCLI):
    def __init__(self):
        ChimeraCLI.__init__(self, "chimera-skyflat", "SkyFlats", 0.1)
        self.add_help_group("SKYFLAT", "skyFlats")
        self.add_controller(
            name="skyflat",
            cls="AutoSkyFlat",
            required=True,
            help_group="SKYFLAT",
            help="Auto Sky Flats",
        )
        self.add_parameters(
            dict(
                name="sun_alt_hi",
                long="sun-high",
                type="float",
                help_group="SKYFLAT",
                help="Highest Sun altitude",
                metavar="DEGREES",
            ),
            dict(
                name="sun_alt_low",
                long="sun-low",
                type="float",
                help_group="SKYFLAT",
                help="Lowest Sun altitude",
                metavar="DEGREES",
            ),
            dict(
                name="filter",
                long="filter",
                short="f",
                required=True,
                type="string",
                help_group="SKYFLAT",
                help="Skyflat filter name, or a comma separated list of them",
            ),
            dict(
                name="number",
                long="number",
                short="n",
                type="int",
                help_group="SKYFLAT",
                help="Number of skyflats to take on each filter",
            ),
            dict(
                name="binning",
                default="1x1",
                help="Apply the selected binning to all frames",
                help_group="SKYFLAT",
            ),
        )

    @action(long="auto", help="Does a sequence of sky flats", help_group="SKYFLAT")
    def do_sequence(self, options):
        """
        Sets variables using command line options
        Take skyflats according to options
        """
        self.out(
            "Pointing scope to the zenith and waiting for the Sun to reach "
            "skyflats altitude range"
        )
        # both limits are optional: only compare them when both were given,
        # comparing None with a float is a TypeError
        if (
            options.sun_alt_low is not None
            and options.sun_alt_hi is not None
            and options.sun_alt_low > options.sun_alt_hi
        ):
            self.exit("sun-low needs to be less than sun-high")
        if options.sun_alt_hi is not None:
            self.skyflat["sun_alt_hi"] = options.sun_alt_hi
        if options.sun_alt_low is not None:
            self.skyflat["sun_alt_low"] = options.sun_alt_low

        num = options.number if options.number else None

        def expose_complete(filter_id, i_flat, exp_time, sky_level):
            self.out(
                f"Filter: {filter_id}, Skyflat # {i_flat} - "
                f"Exposure time: {exp_time:3.2f} seconds - Counts: {sky_level:3.2f}"
            )

        self.skyflat.expose_complete += expose_complete
        try:
            request = {"binning": options.binning}
            for filter_id in options.filter.split(","):
                self.skyflat.get_flats(filter_id.strip(), n_flats=num, request=request)
            self.out("Finished.")
        finally:
            # always unsubscribe: a leaked subscription keeps the dead CLI
            # on the controller's event list
            self.skyflat.expose_complete -= expose_complete

    def __abort__(self):
        self.out("\naborting... ", endl="")

        # copy the Proxy because we are running from a different thread
        skyflat = copy.copy(self.skyflat)
        skyflat.abort()


def main():
    cli = ChimeraAutoSkyFlat()
    cli.run(sys.argv)
    cli.wait()


if __name__ == "__main__":
    main()
