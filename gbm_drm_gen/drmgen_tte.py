__author__ = "grburgess"
from typing import Any, Dict, List, Optional

import astropy.io.fits as fits
import astropy.units as u
import gbmgeometry
import numpy as np

from gbm_drm_gen.drmgen import DRMGen
from gbm_drm_gen.input_edges import (BgoTTEEdges, InputEdges, NaiTTEEdges,
                                     tte_edges)

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


class DRMGenTTE(DRMGen):
    def __init__(
        self,
        tte_file: Optional[str] = None,
        det_name: Optional[str] = None,
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

        :param tte_file: the TTE file the data are associated to (used if det name NOT given)
        :param det_name: either NAI_{**} or BGO_{**} or n* b*
        :param trigdat: the path to a trigdat file
        :param mat_type: 0=direct 1=scattered 2=direct+scattered
        :param time: time relative to trigger to pull spacecraft position or MET if using a poshist file
        :param cspecfile: the cspecfile to pull energy output edges from
        :param poshist: read a poshist file
        """

        
        if det_name is None:

            assert tte_file is not None

            with fits.open(tte_file) as f:

                det_name = f["PRIMARY"].header["DETNAM"]

        else:

            # check the dector name lookups

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

                self._in_edge = custom_input_edges.edges

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

            self._gbm = gbmgeometry.GBM(
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

            self._gbm = gbmgeometry.GBM(
                position_interpolator.quaternion(time),
                position_interpolator.sc_pos(time) * u.m,
            )

        else:

            raise RuntimeError("No trigdat or posthist file used!")

        super(DRMGenTTE, self).__init__(
            position_interpolator=position_interpolator,
            det_number=det_number,
            ebin_edge_in=in_edge,
            mat_type=mat_type,
            ebin_edge_out=out_edge,
            occult=occult,
            time=time
        )

