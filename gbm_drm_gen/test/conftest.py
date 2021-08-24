import pytest

from gbm_drm_gen import DRMGenTTE
from gbm_drm_gen.drmgen import DRMGen
from gbm_drm_gen.utils.package_data import get_path_of_data_file
from gbm_drm_gen.utils.response import OGIPResponse



@pytest.fixture(scope="session")
def built_drm_gen() -> DRMGenTTE:

    tte_file = get_path_of_data_file("example_data/glg_tte_n6_bn110721200_v00.fit")
    trigdat_file = get_path_of_data_file(
        "example_data/glg_trigdat_all_bn110721200_v01.fit"
    )
    cspec_file = get_path_of_data_file("example_data/glg_cspec_n6_bn110721200_v00.pha")

    n6 = DRMGenTTE(tte_file, trigdat=trigdat_file, mat_type=2, cspecfile=cspec_file)

    yield n6


@pytest.fixture(scope="session")
def built_drm_gen_new() -> DRMGenTTE:

    tte_file = get_path_of_data_file("example_data/glg_tte_n6_bn110721200_v00.fit")
    trigdat_file = get_path_of_data_file(
        "example_data/glg_trigdat_all_bn110721200_v01.fit"
    )
    cspec_file = get_path_of_data_file("example_data/glg_cspec_n6_bn110721200_v00.pha")

    n6 = DRMGen.from_128_bin_data(det_name="n6" ,trigdat=trigdat_file, mat_type=2, cspecfile=cspec_file)

    yield n6

    
