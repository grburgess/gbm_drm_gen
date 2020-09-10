__author__ = "fkunzweiler"
import numpy as np
import astropy.io.fits as fits

import gbmgeometry

import astropy.units as u

from gbm_drm_gen.drmgen import DRMGen

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


class DRMGenCTIME(DRMGen):
    """
    TODO: Fix docstring
    A TTE/CSPEC specific drmgen already incorporating the standard input edges. Output edges are obtained
    from the input cspec file. Spacecraft position is read from the TTE file. For further details see the 
    generic reader (DRMGen).

    :param trigdat: the path to a trigdat file
    :param det: the number (0-13) of the detector to be used
    :param mat_type: 0=direct 1=scattered 2=direct+scattered
    :param time: time relative to trigger to pull spacecraft position or MET if using a poshist file
    :param cspecfile: the cspecfile to pull energy output edges from
    :param poshist: read a poshist file
    """

    def __init__(
        self,
        ctime_file=None,
        det_name=None,
        time=0.0,
        trigdat=None,
        poshist=None,
        T0=None,
        mat_type=0,
        occult=False,
    ):

        self._occult = occult

        self._time = time

        self._matrix_type = mat_type

        if det_name is None:

            assert ctime_file is not None

            with fits.open(ctime_file) as f:

                det_name = f["PRIMARY"].header["DETNAM"]

        else:
            assert det_name in det_name_lookup, "must use a valid detector name"

        self._det_number = det_name_lookup[det_name]

        if self._det_number > 11:
            # BGO
            # TODO: FIX BGO
            self._in_edge = np.array(
                [
                ],
                dtype=np.float32,
            )

        else:
            self._in_edge = np.array(
                np.logspace(0.5, 3.7, 301),
                dtype=np.float32,
            )

        # Create the out edge energies
        with fits.open(ctime_file) as f:
            out_edge_start = f["EBOUNDS"].data["E_MIN"]
            out_edge_stop = f["EBOUNDS"].data["E_MAX"][-1]

        out_edge = np.append(out_edge_start, out_edge_stop)

        self._out_edge = out_edge

        if trigdat is not None:

            self._position_interpolator = gbmgeometry.PositionInterpolator.from_trigdat(
                trigdat_file=trigdat
            )

            self._gbm = gbmgeometry.GBM(
                self._position_interpolator.quaternion(time),
                self._position_interpolator.sc_pos(time) * u.km,
            )

        elif poshist is not None:

            self._position_interpolator = gbmgeometry.PositionInterpolator.from_poshist(
                poshist_file=poshist, T0=T0
            )

            self._gbm = gbmgeometry.GBM(
                self._position_interpolator.quaternion(time),
                self._position_interpolator.sc_pos(time) * u.m,
            )

        else:

            raise RuntimeError("No trigdat or posthist file used!")

        self._sc_quaternions_updater()

    def _sc_quaternions_updater(self):

        quaternions = self._position_interpolator.quaternion(self._time)

        sc_pos = self._position_interpolator.sc_pos(self._time)

        super(DRMGenCTIME, self).__init__(
            quaternions=quaternions,
            sc_pos=sc_pos,
            det_number=self._det_number,
            ebin_edge_in=self._in_edge,
            mat_type=self._matrix_type,
            ebin_edge_out=self._out_edge,
            occult=self._occult,
        )
