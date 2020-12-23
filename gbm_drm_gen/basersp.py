import astropy.io.fits as fits
from glob import glob
import numpy as np

import collections
import re

# from gbm_drm_gen.config.gbm_drm_gen_config import gbm_drm_gen_config

lu = dict(
    n0="NAI_00",
    n1="NAI_01",
    n2="NAI_02",
    n3="NAI_03",
    n4="NAI_04",
    n5="NAI_05",
    n6="NAI_06",
    n7="NAI_07",
    n8="NAI_08",
    n9="NAI_09",
    na="NAI_10",
    nb="NAI_11",
    b0="BGO_00",
    b1="BGO_01",
)

import pkg_resources

import h5py


def get_path_of_data_file(data_file):
    file_path = pkg_resources.resource_filename("gbm_drm_gen", "data/%s" % data_file)

    return file_path


path_to_balrog_db = get_path_of_data_file("balrog_db.h5")


class DetDatabase(object):
    def __init__(self, detector_group):
        """
        :param detector_group:
        """

        self._detector_group = detector_group

        self._load_variables()

    def get_rsp(self, z, az):
        """

        :param z:
        :param az:
        :return:
        """
        # z = str(z)
        # az = str(az)

        return self._rsps["%d_%d" % (az, z)]

        # return self._detector_group["z%s_az%s" %
        #                            (z.zfill(6), az.zfill(6))]  #.value

    def _load_variables(self):

        self._rsps = collections.OrderedDict()

        for key, value in self._detector_group.items():

            if key[0] == "z":

                match = re.match("^z0*(\d+)_az0*(\d+)$", key)
                z, az = map(str, match.groups())

                self._rsps["%s_%s" % (az, z)] = value[()]

        self.at_scat_data = self._detector_group["at_scat_data"][()]
        self.e_in = self._detector_group["e_in"][()]
        self.lat_edge = self._detector_group["lat_edge"][()]
        self.theta_edge = self._detector_group["theta_edge"][()]
        self.phi_edge = self._detector_group["phi_edge"][()]
        self.lat_cent = self._detector_group["lat_cent"][()]
        self.theta_cent = self._detector_group["theta_cent"][()]
        self.phi_cent = self._detector_group["phi_cent"][()]
        self.double_phi_cent = self._detector_group["double_phi_cent"][()]

        self.X = self._detector_group["X"][()]
        self.Y = self._detector_group["Y"][()]
        self.Z = self._detector_group["Z"][()]

        self.Azimuth = self._detector_group["Azimuth"][()]
        self.Zenith = self._detector_group["Zenith"][()]
        self.milliaz = self._detector_group["milliaz"][()]
        self.millizen = self._detector_group["millizen"][()]
        self.LIST = self._detector_group["LIST"][()]
        self.LPTR = self._detector_group["LPTR"][()]
        self.LEND = self._detector_group["LEND"][()]

        self.epx_lo = self._detector_group["epx_lo"][()]
        self.epx_hi = self._detector_group["epx_hi"][()]

        self.ichan = self._detector_group["ichan"][()]
        self.ienerg = self._detector_group["ienerg"][()]
        self.energ_lo = self._detector_group["energ_lo"][()]
        self.energ_hi = self._detector_group["energ_hi"][()]


#_h5_database = h5py.File(path_to_balrog_db, "r")

_all_dets = (
    "n0",
    "n1",
    "n2",
    "n3",
    "n4",
    "n5",
    "n6",
    "n7",
    "n8",
    "n9",
    "na",
    "nb",
    "b0",
    "b1",
)

rsp_database = collections.OrderedDict()


def get_database(det):

    assert det in _all_dets
    
    with h5py.File(path_to_balrog_db, "r") as h5_database:

            db = DetDatabase(h5_database[det])

    return db
