import at_scat
import ftran

import numpy as np

try:
    from threeML.utils.response import InstrumentResponse

except(ImportError):
    from gbm_drm_gen.utils.response import InstrumentResponse

from gbm_drm_gen.basersp import rsp_database
from gbm_drm_gen.utils.geometry import ang2cart, is_occulted

lu = [
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
]


class DRMGen(object):

    def __init__(self,
                 quaternions,
                 sc_pos,
                 det_number,
                 ebin_edge_in,
                 mat_type=0,
                 ebin_edge_out=None,
                 occult=True):

        """
        A generic GBM DRM generator. This can be inherited for specific purposes.
        It takes as input various spacecraft files to figure out geometry. The user can
        either supply custom output energy edges or they can be read from a cspec file.

        The great benefit here is the ability to make custom input edges if desired. The
        user needs to have the GBM database pointed at by the environment variable:
        BALROG_DB

        Additionally, matrix folding routines are supplied for spectral fitting. One can also
        simulate Poisson rates from supplied spectral models.

        :param trigdat: the path to a trigdat file
        :param det_number: the number (0-13) of the detector to be used
        :param ebin_edge_in: an array of energy **edges** lenght 1 longer than the num of bins
        :param mat_type: 0=direct 1=scattered 2=direct+scattered
        :param time: time relative to trigger to pull spacecraft position
        :param cspecfile: the cspecfile to pull energy output edges from
        :param ebin_edge_out: an array of output edges

        """

        # Attach inputs to object

        self._matrix_type = mat_type

        self._nobins_in = len(ebin_edge_in) - 1

        # If you want to use occulting
        self._occult = occult

        ##### Initiate database loading:

        self._det_number = det_number

        self._database = rsp_database[lu[det_number]]
        self._ein = np.zeros(self._nobins_in, dtype=np.float32)
        energ_lo = self._database.energ_lo
        energ_hi = self._database.energ_hi
        self._ein[: self._database.ienerg] = energ_lo
        self._ein[self._database.ienerg] = energ_hi[-1]

        self._out_edge = ebin_edge_out
        self._in_edge = ebin_edge_in

        self._nobins_out = len(self._out_edge) - 1

        self._quaternions = quaternions
        self._sc_pos = sc_pos

        self._compute_spacecraft_coordinates()

    @property
    def ebounds(self):

        return self._out_edge

    @property
    def monte_carlo_energies(self):

        return self._in_edge

    @property
    def matrix(self):

        return self._drm.T

    def to_3ML_response(self, ra, dec):

        self.set_location(ra, dec)

        response = InstrumentResponse(
            self.matrix, self.ebounds, self.monte_carlo_energies
        )

        return response

    def to_fits(self, ra, dec, filename, overwrite):

        response = self.to_3ML_response(ra, dec)

        split_filename = filename.split(".")

        if len(split_filename) > 1:

            filename = "".join(split_filename)

        response.to_fits(
            "%s_%s.rsp" % (filename, lu[self._det_number]), "GLAST", "GBM", overwrite
        )

    def set_location(self, ra, dec):
        """
        Set the Ra and Dec of the DRM to be built. This invokes DRM generation as well.

        :param ra: ra in degrees
        :param dec: dec in degrees
        """

        self.ra = ra
        self.dec = dec

        if self._occult:
            if is_occulted(ra, dec, self._sc_pos):
                self._drm = self._occulted_DRM

            else:

                # get the spacecraft coordinates
                az, el = self._get_coords(ra, dec)

                # build the DRM
                self._drm = self._make_drm(az, el, self._geo_az, self._geo_el)
        else:
            # get the spacecraft coordinates
            az, el = self._get_coords(ra, dec)

            # build the DRM
            self._drm = self._make_drm(az, el, self._geo_az, self._geo_el)

        # go ahead and transpose it for spectal fitting, etc.
        # self._drm_transpose = self._drm.T

    def set_time(self, time):

        self._time = time

        self._sc_quaternions_updater()

    def _sc_quaternions_updater(self):

        raise NotImplementedError("implemented in subclass")

    def _compute_spacecraft_coordinates(self):
        """
        GBM geometry calculations
        """

        self._scx = np.zeros(3)
        self._scy = np.zeros(3)
        self._scz = np.zeros(3)

        geodir = np.zeros(3)

        self._scx[0] = (
            self._quaternions[0] ** 2
            - self._quaternions[1] ** 2
            - self._quaternions[2] ** 2
            + self._quaternions[3] ** 2
        )
        self._scx[1] = 2.0 * (
            self._quaternions[0] * self._quaternions[1]
            + self._quaternions[3] * self._quaternions[2]
        )
        self._scx[2] = 2.0 * (
            self._quaternions[0] * self._quaternions[2]
            - self._quaternions[3] * self._quaternions[1]
        )
        self._scy[0] = 2.0 * (
            self._quaternions[0] * self._quaternions[1]
            - self._quaternions[3] * self._quaternions[2]
        )
        self._scy[1] = (
            -self._quaternions[0] ** 2
            + self._quaternions[1] ** 2
            - self._quaternions[2] ** 2
            + self._quaternions[3] ** 2
        )
        self._scy[2] = 2.0 * (
            self._quaternions[1] * self._quaternions[2]
            + self._quaternions[3] * self._quaternions[0]
        )
        self._scz[0] = 2.0 * (
            self._quaternions[0] * self._quaternions[2]
            + self._quaternions[3] * self._quaternions[1]
        )
        self._scz[1] = 2.0 * (
            self._quaternions[1] * self._quaternions[2]
            - self._quaternions[3] * self._quaternions[0]
        )
        self._scz[2] = (
            -self._quaternions[0] ** 2
            - self._quaternions[1] ** 2
            + self._quaternions[2] ** 2
            + self._quaternions[3] ** 2
        )

        geodir[0] = -self._scx.dot(self._sc_pos)
        geodir[1] = -self._scy.dot(self._sc_pos)
        geodir[2] = -self._scz.dot(self._sc_pos)

        denom = np.sqrt(geodir.dot(geodir))

        geodir /= denom

        geo_az = np.arctan2(geodir[1], geodir[0])

        if geo_az < 0.0:
            geo_az += 2 * np.pi
        while geo_az > 2 * np.pi:
            geo_az -= 2 * np.pi

        geo_el = np.arctan2(np.sqrt(geodir[0] ** 2 + geodir[1] ** 2), geodir[2])

        self._geo_el = 90 - np.rad2deg(geo_el)

        self._geo_az = np.rad2deg(geo_az)

        # Also setup the occulted matrix
        self._occulted_DRM = np.zeros((self._nobins_in, self._nobins_out))

    def _get_coords(self, ra, dec):
        source_pos_sc = np.zeros(3)
        source_pos = ang2cart(ra, dec)

        source_pos_sc[0] = self._scx.dot(source_pos)
        source_pos_sc[1] = self._scy.dot(source_pos)
        source_pos_sc[2] = self._scz.dot(source_pos)

        el = np.arccos(source_pos_sc[2])
        az = np.arctan2(source_pos_sc[1], source_pos_sc[0])

        if az < 0.0:
            az += 2 * np.pi
        el = 90 - np.rad2deg(el)
        az = np.rad2deg(0.0 + az)

        return [az, el]

    def _make_drm(self, src_az, src_el, geo_az, geo_el):
        """
        The DRM generator. Should not be invoked by the user except via the set location member funtion.
        :param src_az:
        :param src_el:
        :param geo_az:
        :param geo_el:
        :return:
        """

        final_drm = np.zeros((self._nobins_in, self._nobins_out), dtype=np.float32)

        ## SKY Interpolation

        rlon = src_az
        rlat = src_el

        sf = np.arctan(1.0) / 45.0
        plat = src_el * sf
        plon = src_az * sf
        P = np.array(
            [np.cos(plat) * np.cos(plon), np.cos(plat) * np.sin(plon), np.sin(plat)],
            np.float32,
        )

        N = len(self._database.LEND)
        b1, b2, b3, i1, i2, i3 = ftran.trfind(
            0,
            P,
            N,
            self._database.X,
            self._database.Y,
            self._database.Z,
            self._database.LIST,
            self._database.LPTR,
            self._database.LEND,
        )

        mat1 = self._database.get_rsp(
            self._database.millizen[i1 - 1], self._database.milliaz[i1 - 1]
        )
        mat2 = self._database.get_rsp(
            self._database.millizen[i2 - 1], self._database.milliaz[i2 - 1]
        )
        mat3 = self._database.get_rsp(
            self._database.millizen[i3 - 1], self._database.milliaz[i3 - 1]
        )
        ## Interpolator on triangle

        sum = b1 + b2 + b3
        b1n = b1 / sum
        b2n = b2 / sum
        b3n = b3 / sum

        out_matrix = b1n * mat1 + b2n * mat2 + b3n * mat3

        n_tmp_phot_bin = 2 * self._nobins_in + self._nobins_in % 2
        tmp_phot_bin = np.zeros(n_tmp_phot_bin, dtype=np.float32)
        tmp_phot_bin[::2] = self._in_edge[:-1]
        tmp_phot_bin[1::2] = 10 ** (
            (np.log10(self._in_edge[:-1]) + np.log10(self._in_edge[1:])) / 2.0
        )

        #### Atmospheric scattering
        if self._matrix_type == 1 or self._matrix_type == 2:

            # these change each time the source position is changed
            theta_geo = 90.0 - geo_el
            phi_geo = geo_az
            theta_source = 90 - rlat
            phi_source = rlon

            ## Get new coordinates in the proper space
            gx, gy, gz, sl = ftran.geocords(
                theta_geo, phi_geo, theta_source, phi_source
            )
            lat = 180.0 - np.rad2deg(
                np.arccos(sl[0] * gz[0] + sl[1] * gz[1] + sl[2] * gz[2])
            )
            atscat_diff_matrix = np.zeros((n_tmp_phot_bin - 1, self._nobins_out))
            if lat <= self._database.lat_edge[-1] and (
                lat < self._database.lat_cent[-1]
            ):
                coslat_corr = np.abs(np.cos(np.deg2rad(lat)))

                if lat <= self._database.lat_cent[0]:
                    il_low = 0
                    il_high = 1
                    l_frac = 0.0
                else:
                    idx = np.where(self._database.lat_cent < lat)[0][-1]
                    il_low = idx - 1
                    il_high = idx
                    l_frac = 1.0 - (lat - self._database.lat_cent[idx]) / (
                        self._database.lat_cent[idx + 1] - self._database.lat_cent[idx]
                    )

                # We now have to loop over all the at scat data
                # This could be sped up by moving this to FORTRAN
                tmp_out = np.zeros(
                    (
                        len(self._database.theta_cent)
                        * len(self._database.double_phi_cent)
                        * 2,
                        self._database.ienerg,
                        self._nobins_out,
                    )
                )

                tmp_out = at_scat.get_at_scat(
                    gx,
                    gy,
                    gz,
                    il_low,
                    il_high,
                    l_frac,
                    self._nobins_out,
                    self._out_edge,
                    self._database,
                )

                tmp_out *= coslat_corr
                atscat_diff_matrix = ftran.atscat_highres_ephoton_interpolator(
                    tmp_phot_bin, self._ein, tmp_out
                )

        ###################################

        new_epx_lo, new_epx_hi, diff_matrix = ftran.highres_ephoton_interpolator(
            tmp_phot_bin,
            self._ein,
            out_matrix,
            self._database.epx_lo,
            self._database.epx_hi,
            self._database.ichan,
            n_tmp_phot_bin,
        )

        binned_matrix = ftran.echan_integrator(
            diff_matrix, new_epx_lo, new_epx_hi, self._database.ichan, self._out_edge
        )

        if self._matrix_type == 1:
            binned_matrix = atscat_diff_matrix

        if self._matrix_type == 2:
            binned_matrix[:-1, :] += atscat_diff_matrix

        # Integrate photon edge with trapazoid

        final_drm[:-1, :] = (
            binned_matrix[::2, :][:-1, :] / 2.0
            + binned_matrix[1::2, :][:-1, :]
            + binned_matrix[2::2, :] / 2.0
        ) / 2.0

        return final_drm

    def _at_scat(
        self, itheta, theta_u, iphi, phi_u, gx, gy, gz, il_low, il_high, l_frac
    ):
        """

        :param itheta:
        :param theta_u:
        :param iphi:
        :param phi_u:
        :param gx:
        :param gy:
        :param gz:
        :param il_low:
        :param il_high:
        :param l_frac:
        :return:
        """
        dirx, diry, dirz, az, el = ftran.geo_to_space(theta_u, phi_u, gx, gy, gz)

        sf = np.arctan(1.0) / 45.0
        plat = el * sf
        plon = az * sf
        P = np.array(
            [np.cos(plat) * np.cos(plon), np.cos(plat) * np.sin(plon), np.sin(plat)],
            np.float32,
        )

        # Find a new interpolated matrix
        ist = 0

        N = len(self._database.LEND)
        b1, b2, b3, i1, i2, i3 = ftran.trfind(
            ist,
            P,
            N,
            self._database.X,
            self._database.Y,
            self._database.Z,
            self._database.LIST,
            self._database.LPTR,
            self._database.LEND,
        )
        i_array = np.array([i1, i2, i3])

        dist1 = ftran.calc_sphere_dist(
            az, 90.0 - el, self._database.Azimuth[i1 - 1], self._database.Zenith[i1 - 1]
        )
        dist2 = ftran.calc_sphere_dist(
            az, 90.0 - el, self._database.Azimuth[i2 - 1], self._database.Zenith[i2 - 1]
        )
        dist3 = ftran.calc_sphere_dist(
            az, 90.0 - el, self._database.Azimuth[i3 - 1], self._database.Zenith[i3 - 1]
        )

        dist_array = np.array([dist1, dist2, dist3])
        i1 = i_array[np.argmin(dist_array)]

        tmpdrm = self._database.get_rsp(
            self._database.millizen[i1 - 1], self._database.milliaz[i1 - 1]
        )

        # intergrate the new drm
        direct_diff_matrix = ftran.echan_integrator(
            tmpdrm,
            self._database.epx_lo,
            self._database.epx_hi,
            self._database.ichan,
            self._out_edge,
        )

        # Now let FORTRAN add the at scat to the direct
        tmp_out = ftran.sum_at_scat(
            direct_diff_matrix,
            self._database.at_scat_data[:, :, il_low, itheta, iphi],
            self._database.at_scat_data[:, :, il_high, itheta, iphi],
            l_frac,
        )

        return tmp_out

        ###################################################  ################
