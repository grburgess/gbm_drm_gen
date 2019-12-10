import astropy.io.fits as fits
from glob import glob
import numpy as np

from gbm_drm_gen.config.gbm_drm_gen_config import gbm_drm_gen_config

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


class DetDatabase(object):
    def __init__(self, detector):

        """

        :param detector:
        """
        self._rsp_dict = {}

        self.detector = detector
        self._read_db()

        self._get_tri_info()

        self._read_compressed()

        self._det_atscat_data()

    def get_rsp(self, z, az):

        """

        :param z:
        :param az:
        :return:
        """
        z = str(z)
        az = str(az)

        return self._rsp_dict["z%s_az%s" % (z.zfill(6), az.zfill(6))]

    def _det_atscat_data(self):

        # path = os.environ['BALROG_DB'] + '/test_atscatfile_preeinterp_db002.fits'
        path = (
            gbm_drm_gen_config["gbm drm database location"]
            + "/test_atscatfile_preeinterp_db002.fits"
        )

        if self.detector[0] == "n":
            det_number = 0

        elif self.detector[0] == "b":
            det_number = 1
        else:
            print("Detector name is incorrect!")
            return

        # det_number = int(self.detector[-1])
        at_scat = fits.open(path)

        self.at_scat_data = at_scat["atscat_dbase"].data["AT_SCAT_DATA"][det_number].T

        self.e_in = at_scat["atscat_dbase"].data["E_IN"][det_number]

        self.lat_edge = at_scat["atscat_dbase"].data["LAT_EDGE"][det_number]
        self.theta_edge = at_scat["atscat_dbase"].data["THETA_EDGE"][det_number]
        self.phi_edge = at_scat["atscat_dbase"].data["PHI_EDGE"][det_number]

        self.lat_cent = np.array(
            [
                (self.lat_edge[i] + self.lat_edge[i + 1]) / 2.0
                for i in range(len(self.lat_edge) - 1)
            ]
        )
        self.theta_cent = np.array(
            [
                (self.theta_edge[i] + self.theta_edge[i + 1]) / 2.0
                for i in range(len(self.theta_edge) - 1)
            ]
        )
        self.phi_cent = 180.0 - np.array(
            [
                (self.phi_edge[i] + self.phi_edge[i + 1]) / 2.0
                for i in range(len(self.phi_edge) - 1)
            ]
        )

        self.double_phi_cent = np.array(zip(self.phi_cent, -self.phi_cent))

        at_scat.close()
        del at_scat

    def _get_tri_info(self):

        path = (
            gbm_drm_gen_config["gbm drm database location"]
            + "/InitialGBMDRM_triangleinfo_db001.fits"
        )

        trifile = fits.open(path)

        self.X = trifile["TRIANGULATION"].data["X"]
        self.Y = trifile["TRIANGULATION"].data["Y"]
        self.Z = trifile["TRIANGULATION"].data["Z"]

        self.Azimuth = trifile["TRIANGULATION"].data["Azimuth"][0]
        self.Zenith = trifile["TRIANGULATION"].data["Zenith"][0]
        self.milliaz = trifile["TRIANGULATION"].data["milliaz"][0]
        self.millizen = trifile["TRIANGULATION"].data["millizen"][0]
        self.LIST = trifile["TRIANGULATION"].data["LIST"][0]
        self.LPTR = trifile["TRIANGULATION"].data["LPTR"][0]
        self.LEND = trifile["TRIANGULATION"].data["LEND"][0]

        trifile.close()
        del trifile

    def _read_db(self):

        path = (
            gbm_drm_gen_config["gbm drm database location"]
            + "/GBMDRMdb002/"
            + lu[self.detector]
            + "/"
        )

        self.all_leafs = glob(path + "glg_leaf_%s_z*" % self.detector)

        self._chop_leaf()

        self._construct_dict()

    def _read_compressed(self):

        #        path = os.environ['BALROG_DB']+'/GBMDRMdb002/'+lu[self.detector]+'/'

        fits_file = fits.open(self.all_leafs[0])

        self.epx_lo = fits_file["ECOMPRESS"].data["E_MIN"]
        self.epx_hi = fits_file["ECOMPRESS"].data["E_MAX"]

        fits_file.close()
        del fits_file

    def _construct_dict(self):

        for key, filename in zip(self.dict_names, self.all_leafs):
            self._rsp_dict[key] = self._get_spec_rsp(filename)

    def _chop_leaf(self):

        dict_names = []

        for leaf in self.all_leafs:
            leaf_name = leaf.split("/")[-1]
            _, _, _, z, az, _ = leaf_name.split("_")
            dict_name = "%s_%s" % (z, az)
            dict_names.append(dict_name)

        self.dict_names = dict_names

    def _get_spec_rsp(self, leaf):

        rsp = fits.open(leaf)

        self.ichan = rsp["SPECRESP MATRIX"].header["DETCHANS"]
        self.ienerg = rsp["SPECRESP MATRIX"].header["NUMEBINS"]
        self.energ_lo = rsp["SPECRESP MATRIX"].data["ENERG_LO"]
        self.energ_hi = rsp["SPECRESP MATRIX"].data["ENERG_HI"]
        n_grp = rsp["SPECRESP MATRIX"].data["N_GRP"]
        fchan = rsp["SPECRESP MATRIX"].data["F_CHAN"]
        nchan = rsp["SPECRESP MATRIX"].data["N_CHAN"]

        tmp1 = fchan
        tmp2 = nchan

        matrix = np.zeros((self.ienerg, self.ichan))

        for fcs, ncs, i in zip(tmp1, tmp2, range(self.ienerg)):
            colIndx = 0

            for fc, nc in zip(fcs, ncs):
                matrix[i, fc - 1 : fc + nc - 1] = rsp["SPECRESP MATRIX"].data["MATRIX"][
                    i
                ][colIndx : colIndx + nc]
                colIndx += nc

        rsp.close()
        del rsp
        return matrix
