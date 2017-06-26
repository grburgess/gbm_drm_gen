from drmgen import DRMGen
from drmgen_tte import DRMGenTTE
from balrog_drm import BALROG_DRM

import pkg_resources

import h5py


def get_path_of_data_file(data_file):
    file_path = pkg_resources.resource_filename("gbm_drm_gen", 'data/%s' % data_file)

    return file_path





path_to_balrog_db = get_path_of_data_file('balrog_db.h5')


#from gbm_drm_gen.config.gbm_drm_gen_config import gbm_drm_gen_config


__database = h5py.File(path_to_balrog_db,'r')



