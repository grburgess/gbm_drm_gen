from typing import Optional

import astropy.io.fits as fits
import astropy.units as u
import gbmgeometry
import numba as nb
import numpy as np
from numba import prange, types
# TODO Remove this
from numba.typed import Dict
from responsum import InstrumentResponse

from gbm_drm_gen.basersp_numba import (get_database,
                                       get_trigdat_precalc_database)
from gbm_drm_gen.input_edges import (BgoTTEEdges, InputEdges, NaiTTEEdges,
                                     trigdat_edges, trigdat_out_edge,
                                     tte_edges)
from gbm_drm_gen.matrix_functions import (atscat_highres_ephoton_interpolator,
                                          calc_sphere_dist, closest,
                                          echan_integrator, geo_to_space,
                                          geocoords,
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

det_name_lookup2 = {
    "n0": "NAI_00",
    "n1": "NAI_01",
    "n2": "NAI_02",
    "n3": "NAI_03",
    "n4": "NAI_04",
    "n5": "NAI_05",
    "n6": "NAI_06",
    "n7": "NAI_07",
    "n8": "NAI_08",
    "n9": "NAI_09",
    "na": "NAI_10",
    "nb": "NAI_11",
    "b0": "BGO_00",
    "b1": "BGO_01"}



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
        position_interpolator: gbmgeometry.PositionInterpolator,
        det_number: str,
        ebin_edge_in,
        mat_type: int=0,
        ebin_edge_out=None,
        occult:bool=True,
        time:float=0,
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
        #######################################################################################

        # this checks if the output edges do not change
        
        self._trigdat = True
        self._trigdat_mask = np.array([], dtype=int)

        if det_number < 12:
            trigdat_edges = np.array(
                [3.4, 10., 22., 44., 95., 300., 500., 800., 2000.])
        else:
            trigdat_edges = np.array(
                [150., 400., 850., 1500., 3000., 5500., 10000., 20000., 50000.])

        for i, (e_l, e_h) in enumerate(zip(ebin_edge_out[:-1], ebin_edge_out[1:])):
            elow_index = np.argwhere(np.isclose(e_l,
                                                trigdat_edges,
                                                rtol=0.01))
            if len(elow_index) == 0:
                self._trigdat = False
            elif np.isclose(e_h, trigdat_edges[elow_index[0, 0]+1], rtol=0.01):
                self._trigdat_mask = np.append(self._trigdat_mask,
                                               elow_index[0, 0])
            else:
                self._trigdat = False

        if self._trigdat:
            self._database_precalc_trigdat = get_trigdat_precalc_database(
                lu[det_number], self._trigdat_mask)

        else:
            # dummy
            self._database_precalc_trigdat = get_trigdat_precalc_database(lu[det_number], [
                                                                          0])
            self._trigdat_mask = None

        #######################################################################################
        self._ein = np.zeros(self._nobins_in, dtype=np.float32)
        energ_lo = self._database_nb.energ_lo
        energ_hi = self._database_nb.energ_hi
        self._ein[: self._database_nb.ienerg] = energ_lo
        self._ein[self._database_nb.ienerg] = energ_hi[-1]

        self._out_edge = ebin_edge_out
        self._in_edge = ebin_edge_in

        self._nobins_out = len(self._out_edge) - 1


        self._position_interpolator = position_interpolator

        self.set_time(time)
        

    @classmethod
    def from_128_bin_data(
            cls,
            det_name,
            time: float = 0.0,
            cspecfile: Optional[str] = None,
            trigdat: Optional[str] = None,
            poshist: Optional[str] = None,
            T0: Optional[float] = None,
            mat_type: int = 0,
            custom_input_edges: Optional[InputEdges] = None,
            occult: bool = False,
    ):
        """
        A TTE/CSPEC specific drmgen already incorporating the standard input edges. Output edges are obtained
        from lib import funs the input cspec file. Spacecraft position is read from the TTE file. For further details see the 
        generic reader (DRMGen).

        :param det_name: either NAI_{**} or BGO_{**} or n* b*
        :param trigdat: the path to a trigdat file
        :param mat_type: 0=direct 1=scattered 2=direct+scattered
        :param time: time relative to trigger to pull spacecraft position or MET if using a poshist file
        :param cspecfile: the cspecfile to pull energy output edges from
        :param poshist: read a poshist file
        """


        if det_name not in det_name_lookup:

            if det_name not in det_name_lookup2:

                raise RuntimeError(f"{det_name} is not valid")

            else:

                det_name = det_name_lookup2[det_name]

        det_number = det_name_lookup[det_name]

        if det_number > 11:
            # BGO

            if custom_input_edges is None:

                in_edge = tte_edges["bgo"]

            else:

                assert isinstance(
                    custom_input_edges, BgoTTEEdges), f"custom edges are not an instance of BgoTTEEdges!"

                in_edge = custom_input_edges.edges

        else:

            if custom_input_edges is None:

                in_edge = tte_edges["nai"]

            else:

                assert isinstance(
                    custom_input_edges, NaiTTEEdges), f"custom edges are not an instance of NaiTTEEdges!"

                in_edge = custom_input_edges.edges

        # Create the out edge energies
        with fits.open(cspecfile) as f:
            out_edge = np.zeros(129, dtype=np.float32)
            out_edge[:-1] = f["EBOUNDS"].data["E_MIN"]
            out_edge[-1] = f["EBOUNDS"].data["E_MAX"][-1]

        out_edge = out_edge

        if trigdat is not None:

            try:
                position_interpolator = gbmgeometry.PositionInterpolator.from_trigdat(
                    trigdat_file=trigdat
                )

            except:

                position_interpolator = gbmgeometry.PositionInterpolator.from_trigdat_hdf5(
                    trigdat_file=trigdat
                )

            gbm = gbmgeometry.GBM(
                position_interpolator.quaternion(time),
                position_interpolator.sc_pos(time) * u.km,
            )

        elif poshist is not None:

            try:

                position_interpolator = gbmgeometry.PositionInterpolator.from_poshist(
                    poshist_file=poshist, T0=T0
                )

            except:

                position_interpolator = gbmgeometry.PositionInterpolator.from_poshist_hdf5(
                    poshist_file=poshist, T0=T0
                )

            gbm = gbmgeometry.GBM(
                position_interpolator.quaternion(time),
                position_interpolator.sc_pos(time) * u.m,
            )

        else:

            raise RuntimeError("No trigdat or posthist file used!")

        out = cls(
            position_interpolator=position_interpolator,
            det_number=det_number,
            ebin_edge_in=in_edge,
            mat_type=mat_type,
            ebin_edge_out=out_edge,
            occult=occult,
            time=time
        )

        out._gbm = gbm

        return out

    @classmethod
    def from_trigdat(cls,
                     trigdat_file: str,
                     det,
                     mat_type=0,
                     tstart=0,
                     tstop=0.0,
                     time=0.0,
                     occult=False):

        """TODO describe function

        :param cls: 
        :type cls: 
        :param trigdat_file: 
        :type trigdat_file: str
        :param det: 
        :type det: 
        :param mat_type: 
        :type mat_type: 
        :param tstart: 
        :type tstart: 
        :param tstop: 
        :type tstop: 
        :param time: 
        :type time: 
        :param occult: 
        :type occult: 
        :returns: 

        """

        matrix_type = mat_type

        
        det_number = det

        time = time

        maxen = 140

        
        # Setup the input side energy edges
        if det > 11:

            in_edge = trigdat_edges["bgo"]
            out_edge = trigdat_out_edge["bgo"]

        else:

            in_edge = trigdat_edges["nai"]
            out_edge = trigdat_out_edge["nai"]


        try:
            position_interpolator = gbmgeometry.PositionInterpolator.from_trigdat(
                trigdat_file=trigdat_file
            )

        except:

            position_interpolator = gbmgeometry.PositionInterpolator.from_trigdat_hdf5(
                trigdat_file=trigdat_file
                )

        out = cls(
            position_interpolator=position_interpolator,
            det_number=det_number,
            ebin_edge_in=in_edge,
            mat_type=matrix_type,
            ebin_edge_out=out_edge,
            occult=occult,
            time = time
        )

        out._tstart = tstart
        out._tstop = tstop

        return out



    @property
    def postion_interpolator(self) -> gbmgeometry.PositionInterpolator:

        return self._position_interpolator

    @property
    def current_met(self) -> float:

        return self._position_interpolator.met(self._time)

    def met_at(self, t: float) -> float:

        """
        returns the MET at a given relative time

        :param t: 
        :type t: float
        :returns: 

        """
        return self._position_interpolator.met(t)
    

        
    
    @property
    def ebounds(self):

        return self._out_edge

    @property
    def monte_carlo_energies(self):

        return self._in_edge

    @property
    def matrix(self):

        return self._drm.T

    def to_3ML_response(self, ra, dec, coverage_interval=None):
        """
        create an instrument reponse object for 3ML
        
        :param ra: 
        :type ra: 
        :param dec: 
        :type dec: 
        :returns: 

        """
        self.set_location(ra, dec)

        out = self.matrix

        if not np.all(np.isfinite(out)):

            for i, j in zip(np.where(np.isnan(out))[0], np.where(np.isnan(out))[1]):

                out[i, j] = 0.

        response = InstrumentResponse(
            out, self.ebounds, self.monte_carlo_energies, coverage_interval=coverage_interval
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

    def set_location(self, ra, dec) -> None:
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
                az, el = self._get_coords(ra, dec)
                self._drm = self._make_drm_numba(
                    az, el, self._geo_az, self._geo_el)
        else:
            # get the spacecraft coordinates
            az, el = self._get_coords(ra, dec)

            self._drm = self._make_drm_numba(
                az, el, self._geo_az, self._geo_el)

        # go ahead and transpose it for spectal fitting, etc.
        # self._drm_transpose = self._drm.T
        # go ahead and transpose it for spectal fitting, etc.
        # self._drm_transpose = self._drm.T

    def set_location_direct_sat_coord(self, az, el) -> None:
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

        self._compute_spacecraft_coordinates()
        

    def _sc_quaternions_updater(self):

        self._quaternions = self._position_interpolator.quaternion(self._time)

        self._sc_pos = self._position_interpolator.sc_pos(self._time)

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
            trigdat_precalc_rsps=self._database_precalc_trigdat.rsps,
            trigdat=self._trigdat,
            trigdat_mask=self._trigdat_mask
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
    trigdat_precalc_rsps,
    trigdat,
    trigdat_mask,
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
        ################
        theta_geo = 90.0 - geo_el
        phi_geo = geo_az
        theta_source = 90 - rlat
        phi_source = rlon

        # Get new coordinates in the proper space
        # the np.rad2deg were missing!
        gx, gy, gz, sl = geocoords(
            np.deg2rad(theta_geo), np.deg2rad(phi_geo), np.deg2rad(theta_source), np.deg2rad(phi_source))

        lat = 180.0 - np.rad2deg(
            np.arccos(sl[0] * gz[0] + sl[1] * gz[1] + sl[2] * gz[2])
        )
        atscat_diff_matrix = np.zeros((n_tmp_phot_bin - 1, nobins_out))
        if lat <= lat_edge[-1] and (lat < lat_cent[-1]):
            coslat_corr = np.abs(np.cos(np.deg2rad(lat)))

            if lat <= lat_cent[0]:
                il_low = 0
                il_high = 1
                # TODO Changes this. Was 0.
                l_frac = 1.0
            else:
                idx = np.where(lat_cent < lat)[0][-1]
                # TODO Changed this. This was il_low = idx -1 and il_high = idx
                il_low = idx
                il_high = idx+1
                l_frac = 1.0 - (lat - lat_cent[idx]) / (
                    lat_cent[idx + 1] - lat_cent[idx]
                )

            num_theta = len(theta_cent)
            num_phi = len(double_phi_cent)

            theta_u = theta_cent
            phi_u = double_phi_cent

            # loop over all the fucking atm matrices

            tmp_out = np.zeros((num_theta, num_phi, 2, ienerg, nobins_out))
            sf = np.arctan(1.0) / 45.0

            at_scat(tmp_out, num_theta, num_phi, theta_u, phi_u, gx, gy, gz,
                    sf, grid_points_list, trigdat_precalc_rsps, rsps, milliaz, millizen, epx_lo, epx_hi, ichan, out_edge,
                    at_scat_data, il_low, il_high, l_frac, trigdat)
            tmp_out2 = np.zeros((ienerg, nobins_out))
            for i in range(num_theta):
                for j in range(num_phi):
                    for k in range(2):
                        tmp_out2 += tmp_out[i, j, k]

            tmp_out2 *= coslat_corr

            atscat_diff_matrix = atscat_highres_ephoton_interpolator(
                tmp_phot_bin,  ein, tmp_out2)

        ##############
        """
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
                l_frac = 0.0 #TODO Not 1?
            else:
                idx = np.where(lat_cent < lat)[0][-1]
                #TODO Check this. This was il_low = idx -1 and il_high = idx
                il_low = idx
                il_high = idx + 1
                l_frac = 1.0 - (lat - lat_cent[idx]) / (
                    lat_cent[idx + 1] - lat_cent[idx]
                )

            num_theta = len(theta_cent)
            num_phi = len(double_phi_cent)

            theta_u = theta_cent
            phi_u = double_phi_cent

            # loop over all the fucking atm matrices

            tmp_out = np.zeros((ienerg, nobins_out))
            sf = np.arctan(1.0) / 45.0

            for i in prange(num_theta):
                for j in prange(num_phi):
                    #TODO: Check this, this was range(1), but I think 2 is needed to also cover
                    # the negative phi values
                    for k in prange(2):

                        dirx, diry, dirz, az, el = geo_to_space(
                            theta_u[i], phi_u[j, k], gx, gy, gz)
                        plat = el * sf
                        plon = az * sf
                        P = np.array(
                            [
                                np.cos(plat) * np.cos(plon),
                                np.cos(plat) * np.sin(plon),
                                np.sin(plat),
                            ]
                        )

                        # Find the closest direct matrix grid point and use it. No interpolation... should be accurate enough
                        i1 = closest(P, grid_points_list)

                        if new:
                            if trigdat:
                                #TODO why i1-1? makes no sense to me...
                                rsps_echan_integrator = trigdat_precalc_rsps[milliaz[i1] + "_" + millizen[i1]]

                            else:
                                rsps_echan_integrator = echan_integrator(rsps[milliaz[i1] + "_" + millizen[i1]],
                                                                         epx_lo,
                                                                         epx_hi,
                                                                         ichan,
                                                                         out_edge
                                )

                            msum(at_scat_data[..., il_low, i, j],
                                 at_scat_data[..., il_high, i, j],
                                 np.ascontiguousarray(rsps_echan_integrator),
                                 tmp_out,
                                 l_frac,
                            )
                            #msum(np.ascontiguousarray(at_scat_data[..., il_low, i, j]),
                            #     np.ascontiguousarray(
                            #         at_scat_data[..., il_high, i, j]),
                            #     np.ascontiguousarray(rsps_echan_integrator),
                            #     tmp_out,
                            #     l_frac,
                            #     ienerg,
                            #     nobins_out
                            #)

                        else:
                            msum(np.ascontiguousarray(at_scat_data[..., il_low, i, j]),
                                 np.ascontiguousarray(
                                     at_scat_data[..., il_high, i, j]),
                                 np.ascontiguousarray(echan_integrator(rsps[milliaz[i1] + "_" + millizen[i1]],
                                                                       epx_lo,
                                                                       epx_hi,
                                                                       ichan,
                                                                       out_edge
                                 )),
                                 tmp_out,
                                 l_frac,
                            )

            tmp_out *= coslat_corr
            atscat_diff_matrix = atscat_highres_ephoton_interpolator(
                tmp_phot_bin,  ein, tmp_out)
        """

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
        # TODO: Why :-1 ?
        binned_matrix[:-1, :] += atscat_diff_matrix
    # Integrate photon edge with trapazoid
    final_drm[:-1, :] = (
        binned_matrix[::2, :][:-1, :] / 2.0
        + binned_matrix[1::2, :][:-1, :]
        + binned_matrix[2::2, :] / 2.0
    ) / 2.0
    return final_drm


@nb.njit(parallel=True, fastmath=True)
def at_scat(tmp_out, num_theta, num_phi, theta_u, phi_u, gx, gy, gz,
            sf, grid_points_list, trigdat_precalc_rsps, rsps, milliaz, millizen, epx_lo, epx_hi, ichan, out_edge,
            at_scat_data, il_low, il_high, l_frac, trigdat):
    for i in prange(num_theta):
        for j in prange(num_phi):
            # TODO: Changed this. This was range(1), but I think 2 is needed to also cover
            # the negative phi values
            for k in prange(2):
                dirx, diry, dirz, az, el = geo_to_space(
                    theta_u[i], phi_u[j, k], gx, gy, gz)
                plat = el * sf
                plon = az * sf
                P = np.array(
                    [
                        np.cos(plat) * np.cos(plon),
                        np.cos(plat) * np.sin(plon),
                        np.sin(plat),
                    ]
                )

                # Find the closest direct matrix grid point and use it. No interpolation... should be accurate enough

                # TODO: Difference, this was wrong in the old way. calc_sphere_dist(az, 90.-el, database.Azimuth[...], database.Zenith[...])
                # must be calc_sphere_dist(az, el, database.Azimuth[...], 90-database.Zenith[...]. Gave wrong grid point
                # as closest grid point (gave angles>20 degree)
                i1 = closest(P, grid_points_list)

                if trigdat:
                    # i1
                    rsps_echan_integrator = trigdat_precalc_rsps[milliaz[i1] +
                                                                 "_" + millizen[i1]]

                else:
                    rsps_echan_integrator = echan_integrator(rsps[milliaz[i1] + "_" + millizen[i1]],  # i1
                                                             epx_lo,
                                                             epx_hi,
                                                             ichan,
                                                             out_edge
                                                             )

                msum(at_scat_data[..., il_low, i, j],
                     at_scat_data[..., il_high, i, j],
                     np.ascontiguousarray(rsps_echan_integrator),
                     tmp_out[i, j, k],
                     l_frac,
                     )


@nb.njit(fastmath=True, parallel=False)
def msum(this_data_lo, this_data_hi, direct_diff_matrix, tmp_out, l_frac):
    tmp = this_data_lo*l_frac+this_data_hi*(1-l_frac)
    tmp_out += np.dot(tmp, direct_diff_matrix)
    # for ii in prange(ienerg):
    #    for jj in prange(nobins_out):
    #        for kk in prange(ienerg):
    #            tmp_out[ii, jj] += (
    #                this_data_lo[ii, kk] * l_frac
    #                + this_data_hi[ii, kk]
    #                * (1 - l_frac)
    #            ) * direct_diff_matrix[kk, jj]
