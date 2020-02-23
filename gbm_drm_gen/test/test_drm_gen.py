from gbm_drm_gen import *
from gbm_drm_gen.utils.response import OGIPResponse
from gbm_drm_gen.utils.package_data import get_path_of_data_file
import astropy.io.fits as fits


def test_basic_generation():

    tte_file = get_path_of_data_file("example_data/glg_tte_n6_bn110721200_v00.fit")
    trigdat_file = get_path_of_data_file("example_data/glg_trigdat_all_bn110721200_v01.fit")
    cspec_file = get_path_of_data_file("example_data/glg_cspec_n6_bn110721200_v00.pha")

    n6 = DRMGenTTE(tte_file, trigdat=trigdat_file, mat_type=0, cspecfile=cspec_file)

    n6.to_fits(0, 0, "test", overwrite=True)
