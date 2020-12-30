# import at_scat
# import ftran

import astropy.io.fits as fits
import astropy.units as u
import gbmgeometry
import numba as nb
import numpy as np
from numba import prange

try:
    from threeML.utils.OGIP.response import InstrumentResponse

except:
    from responsum import InstrumentResponse

from gbm_drm_gen.basersp_numba import get_database
from gbm_drm_gen.matrix_functions import (atscat_highres_ephoton_interpolator,
                                          calc_sphere_dist, echan_integrator,
                                          geo_to_space, geocoords,
                                          highres_ephoton_interpolator, trfind)
from gbm_drm_gen.utils.geometry import ang2cart, is_occulted

det_name_lookup = {
    "NAI_00": 0,
    "NAI_01": 1,
    "NAI_02": 2,
    "NAI_03": 3,
    "NAI_04": 4,
    "NAI_05": 5,
    "NAI_06": 6,
    "NAI_07": 7,
    "NAI_08": 8,
    "NAI_09": 9,
    "NAI_10": 10,
    "NAI_11": 11,
    "BGO_00": 12,
    "BGO_01": 13,
}


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
    def __init__(
        self,
        quaternions,
        sc_pos,
        det_number,
        ebin_edge_in,
        mat_type=0,
        ebin_edge_out=None,
        occult=True,
    ):
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

        # Initiate database loading:

        self._det_number = det_number

        self._database_nb = get_database(lu[det_number])
        self._ein = np.zeros(self._nobins_in, dtype=np.float32)
        energ_lo = self._database_nb.energ_lo
        energ_hi = self._database_nb.energ_hi
        self._ein[: self._database_nb.ienerg] = energ_lo
        self._ein[self._database_nb.ienerg] = energ_hi[-1]

        self._out_edge = ebin_edge_out
        self._in_edge = ebin_edge_in

        self._nobins_out = len(self._out_edge) - 1

        self._quaternions = quaternions
        self._sc_pos = sc_pos

        self._compute_spacecraft_coordinates()

    # @classmethod
    # def from_tte(cls)

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

        self.set_location(ra, dec, use_numba=True)

        response = InstrumentResponse(
            self.matrix, self.ebounds, self.monte_carlo_energies
        )

        return response

    def to_3ML_response_direct_sat_coord(self, az, el):

        self.set_location_direct_sat_coord(az, el)

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
            "%s_%s.rsp" % (filename, lu[self._det_number]
                           ), "GLAST", "GBM", overwrite
        )

    def set_location(self, ra, dec, use_numba=True):
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

            self._drm = self._make_drm_numba(
                az, el, self._geo_az, self._geo_el)

        # go ahead and transpose it for spectal fitting, etc.
        # self._drm_transpose = self._drm.T
        # go ahead and transpose it for spectal fitting, etc.
        # self._drm_transpose = self._drm.T

    def set_location_direct_sat_coord(self, az, el):
        """
        Set the AZ and EL in satellite coordinates of the DRM to be built. This invokes DRM generation as well.

        :param az: az in degrees
        :param el: el in degrees
        """

        self.ra = az
        self.dec = el

        if self._occult:
            if is_occulted(az, el, self._sc_pos):
                self._drm = self._occulted_DRM

            else:
                # build the DRM
                self._drm = self._make_drm_numba(
                    az, el, self._geo_az, self._geo_el)
        else:
            # build the DRM
            self._drm = self._make_drm_numba(
                az, el, self._geo_az, self._geo_el)

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

        geo_el = np.arctan2(
            np.sqrt(geodir[0] ** 2 + geodir[1] ** 2), geodir[2])

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

    def _make_drm_numba(self, src_az, src_el, geo_az, geo_el):

        # move outside loop

        n_tmp_phot_bin = 2 * self._nobins_in + self._nobins_in % 2

        tmp_phot_bin = np.zeros(n_tmp_phot_bin)
        tmp_phot_bin[::2] = self._in_edge[:-1]
        tmp_phot_bin[1::2] = 10 ** (
            (np.log10(self._in_edge[:-1]) + np.log10(self._in_edge[1:])) / 2.0
        )

        return _build_drm(
            src_az,
            src_el,
            geo_az,
            geo_el,
            nobins_in=self._nobins_in,
            nobins_out=self._nobins_out,
            Azimuth=self._database_nb.Azimuth,
            Zenith=self._database_nb.Zenith,
            grid_points_list=self._database_nb.grid_points_list,
            milliaz=self._database_nb.milliaz,
            millizen=self._database_nb.millizen,
            in_edge=self._in_edge,
            lat_edge=self._database_nb.lat_edge,
            lat_cent=self._database_nb.lat_cent,
            theta_cent=self._database_nb.theta_cent,
            phi_cent=self._database_nb.phi_cent,
            double_phi_cent=self._database_nb.double_phi_cent,
            ienerg=self._database_nb.ienerg,
            out_edge=self._out_edge,
            ein=self._ein,
            epx_lo=self._database_nb.epx_lo,
            epx_hi=self._database_nb.epx_hi,
            ichan=self._database_nb.ichan,
            matrix_type=self._matrix_type,
            rsps=self._database_nb.rsps,
            n_tmp_phot_bin=n_tmp_phot_bin,
            tmp_phot_bin=tmp_phot_bin,
            at_scat_data=self._database_nb.at_scat_data,
        )


@nb.njit(fastmath=True, parallel=False)
def _build_drm(
    src_az,
    src_el,
    geo_az,
    geo_el,
    nobins_in,
    nobins_out,
    Azimuth,
    Zenith,
    grid_points_list,
    milliaz,
    millizen,
    in_edge,
    lat_edge,
    lat_cent,
    theta_cent,
    phi_cent,
    double_phi_cent,
    ienerg,
    out_edge,
    ein,
    epx_lo,
    epx_hi,
    ichan,
    matrix_type,
    rsps,
    n_tmp_phot_bin,
    tmp_phot_bin,
    at_scat_data,
):

    final_drm = np.zeros((nobins_in, nobins_out))

    # SKY Interpolation

    dtr = np.pi / 180.0

    rlon = src_az
    rlat = src_el

    sf = np.arctan(1.0) / 45.0
    plat = src_el * sf
    plon = src_az * sf
    P = np.array(
        [np.cos(plat) * np.cos(plon), np.cos(plat)
         * np.sin(plon), np.sin(plat)]
    )

    b1, b2, b3, i1, i2, i3 = trfind(P, grid_points_list)

    i1 += 1
    i2 += 1
    i3 += 1

    # mod here

    mat1 = rsps[milliaz[i1 - 1] + "_" + millizen[i1 - 1]]
    mat2 = rsps[milliaz[i2 - 1] + "_" + millizen[i2 - 1]]
    mat3 = rsps[milliaz[i3 - 1] + "_" + millizen[i3 - 1]]

    # Interpolator on triangle

    sum = b1 + b2 + b3
    b1n = b1 / sum
    b2n = b2 / sum
    b3n = b3 / sum

    # need to do a loop here

    out_matrix = np.empty((mat1.shape[0], mat1.shape[1]))

    for i in range(mat1.shape[0]):
        for j in range(mat1.shape[1]):
            out_matrix[i, j] = b1n * mat1[i, j] + \
                b2n * mat2[i, j] + b3n * mat3[i, j]

    # move outside loop
    # n_tmp_phot_bin = 2 * self._nobins_in + self._nobins_in % 2
    # tmp_phot_bin = np.zeros(n_tmp_phot_bin, dtype=np.float32)
    # tmp_phot_bin[::2] = self._in_edge[:-1]
    # tmp_phot_bin[1::2] = 10 ** (
    #     (np.log10(self._in_edge[:-1]) + np.log10(self._in_edge[1:])) / 2.0
    # )

    # Atmospheric scattering
    if matrix_type == 1 or matrix_type == 2:

        # these change each time the source position is changed
        theta_geo = 90.0 - geo_el
        phi_geo = geo_az
        theta_source = 90 - rlat
        phi_source = rlon

        # Get new coordinates in the proper space
        gx, gy, gz, sl = geocoords(
            theta_geo, phi_geo, theta_source, phi_source)

        lat = 180.0 - np.rad2deg(
            np.arccos(sl[0] * gz[0] + sl[1] * gz[1] + sl[2] * gz[2])
        )
        atscat_diff_matrix = np.zeros((n_tmp_phot_bin - 1, nobins_out))
        if lat <= lat_edge[-1] and (lat < lat_cent[-1]):
            coslat_corr = np.abs(np.cos(np.deg2rad(lat)))

            if lat <= lat_cent[0]:
                il_low = 0
                il_high = 1
                l_frac = 0.0
            else:
                idx = np.where(lat_cent < lat)[0][-1]
                il_low = idx - 1
                il_high = idx
                l_frac = 1.0 - (lat - lat_cent[idx]) / (
                    lat_cent[idx + 1] - lat_cent[idx]
                )

            num_theta = len(theta_cent)
            num_phi = len(double_phi_cent)

            theta_u = theta_cent
            phi_u = double_phi_cent

            # loop over all the fucking atm matrices

            tmp_out = np.zeros((ienerg, nobins_out))

            # itr = 0
            for i in range(num_theta):
                for j in range(num_phi):
                    for k in range(1):

                        dirx, diry, dirz, az, el = geo_to_space(
                            theta_u[i], phi_u[j, k], gx, gy, gz)
                        sf = np.arctan(1.0) / 45.0
                        plat = el * sf
                        plon = az * sf
                        P = np.array(
                            [
                                np.cos(plat) * np.cos(plon),
                                np.cos(plat) * np.sin(plon),
                                np.sin(plat),
                            ]
                        )

                        # Find a new interpolated matrix
                        ist = 0

                        b1, b2, b3, i1, i2, i3 = trfind(P, grid_points_list)

                        i1 += 1
                        i2 += 1
                        i3 += 1

                        #                        N = len(LEND)

                        i_array = [i1, i2, i3]

                        dist1 = calc_sphere_dist(
                            az, 90.0 - el, Azimuth[i1 - 1], Zenith[i1 - 1], dtr
                        )
                        dist2 = calc_sphere_dist(
                            az, 90.0 - el, Azimuth[i2 - 1], Zenith[i2 - 1], dtr
                        )
                        dist3 = calc_sphere_dist(
                            az, 90.0 - el, Azimuth[i3 - 1], Zenith[i3 - 1], dtr
                        )

                        dist_array = np.array([dist1, dist2, dist3])

                        i1 = i_array[np.argmin(dist_array)]
                        #                        print(tmpdrm.dtype)

                        # intergrate the new drm
                        # direct_diff_matrix =

                        msum(np.ascontiguousarray(at_scat_data[..., il_low, i, j]),
                             np.ascontiguousarray(
                                 at_scat_data[..., il_high, i, j]),
                             np.ascontiguousarray(echan_integrator(rsps[milliaz[i1 - 1] + "_" + millizen[i1 - 1]],
                                                                   epx_lo,
                                                                   epx_hi,
                                                                   ichan,
                                                                   out_edge
                                                                   )),
                             tmp_out,
                             l_frac,
                             ienerg,
                             nobins_out
                             )

            tmp_out *= coslat_corr
            atscat_diff_matrix = atscat_highres_ephoton_interpolator(
                tmp_phot_bin,  ein, tmp_out)

#    return atscat_diff_matrix
    ###################################

    new_epx_lo, new_epx_hi, diff_matrix = highres_ephoton_interpolator(
        tmp_phot_bin, ein, out_matrix, epx_lo, epx_hi, 64,
    )

    binned_matrix = echan_integrator(
        diff_matrix, new_epx_lo, new_epx_hi, ichan, out_edge
    )

    if matrix_type == 1:
        binned_matrix = atscat_diff_matrix

    if matrix_type == 2:
        binned_matrix[:-1, :] += atscat_diff_matrix

    # Integrate photon edge with trapazoid

    final_drm[:-1, :] = (
        binned_matrix[::2, :][:-1, :] / 2.0
        + binned_matrix[1::2, :][:-1, :]
        + binned_matrix[2::2, :] / 2.0
    ) / 2.0

    return final_drm


@nb.njit(fastmath=True, parallel=True)
def msum(this_data_lo, this_data_hi, direct_diff_matrix, tmp_out, l_frac, ienerg, nobins_out):
    for ii in prange(ienerg):
        for jj in prange(nobins_out):
            for kk in prange(ienerg):

                tmp_out[ii, jj] += (
                    this_data_lo[ii, kk] * l_frac
                    + this_data_hi[ii, kk]
                    * (1 - l_frac)
                ) * direct_diff_matrix[kk, jj]
