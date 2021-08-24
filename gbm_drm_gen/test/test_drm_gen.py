import os
from pathlib import Path

import astropy.io.fits as fits
import numpy as np

from gbm_drm_gen import BgoTTEEdges, DRMGenTTE, NaiTTEEdges
from gbm_drm_gen.utils.package_data import get_path_of_data_file
from gbm_drm_gen.utils.response import OGIPResponse
from gbm_drm_gen import create_rsp2


def test_basic_generation(built_drm_gen):

    os.environ["NUMBA_DISABLE_JIT"] = "1"

    built_drm_gen.to_fits(0, 10, "test", overwrite=True)

    os.environ["NUMBA_DISABLE_JIT"] = "0"


def test_various_options(built_drm_gen):

    os.environ["NUMBA_DISABLE_JIT"] = "0"

    built_drm_gen.set_location(5, 0)

    built_drm_gen.set_location(10, 100)


def test_custom_edges():

    tte_file = get_path_of_data_file(
        "example_data/glg_tte_n6_bn110721200_v00.fit")

    trigdat_file = get_path_of_data_file(
        "example_data/glg_trigdat_all_bn110721200_v01.fit"
    )

    cspec_file = get_path_of_data_file(
        "example_data/glg_cspec_n6_bn110721200_v00.pha")

    edges = NaiTTEEdges.from_log_bins(145)

    n6 = DRMGenTTE(tte_file, trigdat=trigdat_file, mat_type=2,
                   cspecfile=cspec_file, custom_input_edges=edges)

    xx = np.geomspace(5., 50000., 83)

    edges = NaiTTEEdges.from_custom_array(xx)

    n6 = DRMGenTTE(tte_file, trigdat=trigdat_file, mat_type=2,
                   cspecfile=cspec_file, custom_input_edges=edges)



def test_time_change(built_drm_gen_new):

    os.environ["NUMBA_DISABLE_JIT"] = "1"

    built_drm_gen_new.set_time(10)

    os.environ["NUMBA_DISABLE_JIT"] = "0"

    

def test_rsp2_gen(built_drm_gen_new):

    create_rsp2("out.rsp2", built_drm_gen_new, 0, 0, 0, 20, delta_time=1)

    with fits.open("out.rsp2") as f:

        assert f[0].header["DRM_NUM"] == 20

        

    Path("out.rsp2").unlink()




    
