from gbm_drm_gen import *
from gbm_drm_gen.io.package_utils import get_path_of_data_dir


import os



def test_drm_creation():


    tte_file = os.path.join(get_path_of_data_dir(), 'example_data', 'glg_tte_n6_bn110721200_v00.fit' )
    tdat_file = os.path.join(get_path_of_data_dir(), 'example_data', 'glg_trigdat_all_bn110721200_v01.fit' )
    cspec_file = os.path.join(get_path_of_data_dir(), 'example_data', 'glg_cspec_n6_bn110721200_v00.pha' )
    
    n6 = DRMGenTTE(tte_file,
                   trigdat=tdat_file,
                   mat_type=2,
                   cspecfile=cspec_file)

    n6.set_location(0.,0.)


    n6 = DRMGenTTE(tte_file,
                   trigdat=tdat_file,
                   mat_type=1,
                   cspecfile=cspec_file)

    n6.set_location(0.,0.)


    n6 = DRMGenTTE(tte_file,
                   trigdat=tdat_file,
                   mat_type=0,
                   cspecfile=cspec_file)

    n6.set_location(0.,0.)

