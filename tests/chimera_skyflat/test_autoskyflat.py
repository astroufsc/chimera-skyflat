import numpy as np
from astropy.io import fits

from chimera_skyflat.controllers.autoskyflat import AutoSkyFlat


def test_get_sky_level_returns_plain_float(tmp_path):
    """numpy.float64 is not msgspec-serializable: leaking it into the
    expose_complete event payload silently dropped the publication, so
    subscribers (the robobs frame counter) never heard about the frame."""
    path = tmp_path / "flat.fits"
    fits.PrimaryHDU(np.full((8, 8), 1000.0, dtype=np.float32)).writeto(path)

    # chimera's metaobject wraps methods and rejects a non-instance self:
    # call the wrapped function directly (get_sky_level does not use self)
    level = AutoSkyFlat.get_sky_level.func(None, str(path), None)

    assert type(level) is float
    assert level == 1000.0
