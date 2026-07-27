# SPDX-FileCopyrightText: 2015-present Antonio Kanaan <kanaan@astro.ufsc.br>
# SPDX-License-Identifier: GPL-2.0-or-later

"""The controller driven the way chimera drives it: inside a Manager, over
the bus, against the real Site and the core's fake instruments.

This is what catches the API drift a mocked test cannot see - Site.sunpos()
returning a Position of Coords, ImageUtil.download() replacing
Image.download(), and event payloads that msgspec refuses to encode.
"""

import json
import random
import threading
import time

import pytest
from chimera.controllers.imageserver.imageserver import ImageServer
from chimera.core.bus import Bus
from chimera.core.manager import Manager
from chimera.core.site import Site
from chimera.instruments.fakecamera import FakeCamera
from chimera.instruments.fakefilterwheel import FakeFilterWheel
from chimera.instruments.faketelescope import FakeTelescope

from chimera_skyflat.controllers.autoskyflat import AutoSkyFlat

# a sky bright enough at any sun altitude to keep the test's exposures
# short: the model's constant floor dominates scale * exp(slope * alt)
FLOOD_LIT_SKY = {"R": [1.0, 1.0, 100000.0], "B": [0.5, 1.0, 100000.0]}


@pytest.fixture
def manager():
    bus = Bus(f"tcp://127.0.0.1:{random.randint(20000, 60000)}")
    bus_thread = threading.Thread(
        target=bus.run_forever, name="test-manager-bus", daemon=True
    )
    bus_thread.start()
    assert bus._bus_started.wait(5)

    manager = Manager(bus)
    yield manager

    manager.shutdown()
    bus.shutdown()
    bus_thread.join(timeout=10)


@pytest.fixture
def skyflat(manager, tmp_path):
    coefficients = tmp_path / "coefficients.json"
    coefficients.write_text(json.dumps(FLOOD_LIT_SKY))

    manager.add_class(Site, "test", {"name": "test"})
    manager.add_class(
        ImageServer, "test", {"images_dir": str(tmp_path), "httpd": False}
    )
    manager.add_class(FakeTelescope, "fake")
    manager.add_class(
        FakeCamera, "fake", {"ccd_width": 64, "ccd_height": 64, "use_dss": False}
    )
    manager.add_class(FakeFilterWheel, "fake", {"filters": ["R", "B"]})
    manager.add_class(
        AutoSkyFlat,
        "flats",
        {
            "site": "/Site/test",
            "telescope": "/FakeTelescope/fake",
            "camera": "/FakeCamera/fake",
            "filterwheel": "/FakeFilterWheel/fake",
            "dome": None,
            "coefficients_file": str(coefficients),
            # any sun altitude is inside the window: the test must not
            # depend on the wall clock being at twilight
            "sun_alt_hi": 90,
            "sun_alt_low": -90,
            "ideal_counts": 1000,
            "exptime_min": 0.0,
            "exptime_increment": 0.1,
            "flat_position_max": 0,
        },
    )

    return manager.get_proxy("/AutoSkyFlat/flats")


def test_a_flat_set_runs_over_the_bus(skyflat):
    frames = []

    def expose_complete(filter_id, i_flat, exp_time, sky_level):
        frames.append((filter_id, i_flat, exp_time, sky_level))

    skyflat.expose_complete += expose_complete
    try:
        taken = skyflat.get_flats("R", n_flats=2, request={"binning": "1x1"})

        assert taken == 2
        # the event has to survive msgspec encoding to reach a subscriber:
        # a numpy float in the payload used to drop it silently
        deadline = time.time() + 10
        while len(frames) < 2 and time.time() < deadline:
            time.sleep(0.05)
        assert len(frames) == 2, f"expose_complete never arrived: {frames}"
        assert [f[0] for f in frames] == ["R", "R"]
        assert [f[1] for f in frames] == [1, 2]
        assert all(isinstance(f[3], float) for f in frames)
    finally:
        skyflat.expose_complete -= expose_complete


def test_abort_stops_a_running_set(skyflat):
    skyflat["exptime_max"] = 300
    done = threading.Event()
    taken = []

    def run():
        taken.append(skyflat.get_flats("R", n_flats=1000))
        done.set()

    worker = threading.Thread(target=run, daemon=True)
    worker.start()
    time.sleep(0.5)
    skyflat.abort()

    assert done.wait(30), "the set did not stop after abort()"
