import os

import astropy.io.fits as fits


def test_basic_generation(built_drm_gen):

    os.environ["NUMBA_DISABLE_JIT"] = "1"

    built_drm_gen.to_fits(0, 10, "test", overwrite=True)

    os.environ["NUMBA_DISABLE_JIT"] = "0"


def test_various_options(built_drm_gen):

    os.environ["NUMBA_DISABLE_JIT"] = "0"

    built_drm_gen.set_location(5, 0)

    built_drm_gen.set_location(10, 100)
