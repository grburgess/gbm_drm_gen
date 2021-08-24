import gbmgeometry
import numpy as np

from .drmgen import DRMGen
from .input_edges import trigdat_edges, trigdat_out_edge


class DRMGenTrig(DRMGen):
    def __init__(
            self,
            trigdat_file: str,
            det,
            mat_type=0,
            tstart=0,
            tstop=0.0,
            time=0.0,
            occult=False,
    ):
        """
        Inherited drmgen from the TTE version. Builds 8-channel RSPs for trigdat data using 140 input edges
        :param trigdat: a TrigReader object with the background and source selection already exists
        :param det: a number corresponding to the GBM detector to be used
        :param mat_type: the type of matrix to produce: 0=direct 1=scattered 2=direct+scattered
        :param time: the time of the spacecraft position to use
        :param occult: (bool) occult points blocked by the Earth
        """

        matrix_type = mat_type

        self._tstart = tstart
        self._tstop = tstop

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


            super(DRMGenTrig, self).__init__(
                position_interpolator=position_interpolator,
                det_number=det_number,
                ebin_edge_in=in_edge,
                mat_type=matrix_type,
                ebin_edge_out=out_edge,
                occult=occult,
                time = time
        )



################
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
