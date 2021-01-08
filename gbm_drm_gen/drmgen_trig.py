import numpy as np

from .drmgen import DRMGen
from .input_edges import trigdat_edges, trigdat_out_edge


class DRMGenTrig(DRMGen):
    def __init__(
        self,
        quaternions,
        sc_pos,
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

        self._matrix_type = mat_type

        self._tstart = tstart
        self._tstop = tstop

        self._det_number = det

        self._time = time

        maxen = 140

        self._occult = occult
        # Setup the input side energy edges
        if det > 11:

            self._in_edge = trigdat_edges["bgo"]
            self._out_edge = trigdat_out_edge["bgo"]

        else:

            self._in_edge = trigdat_edges["nai"]
            self._out_edge = trigdat_out_edge["nai"]


        self._all_quats = quaternions
        self._all_sc_pos = sc_pos

        self._sc_quaternions_updater()

        #
        # super(DRMGenTrig, self).__init__(quaternions=quaternions,
        #                                 sc_pos=sc_pos,
        #                                 det_number=det,
        #                                 ebin_edge_in=self._in_edge,
        #                                 mat_type=mat_type,
        #                                 ebin_edge_out=self._out_edge)

    def _sc_quaternions_updater(self):

        condition = np.logical_and(
            self._tstart <= self._time, self._time <= self._tstop
        )

        quaternions = self._all_quats[condition][0]
        sc_pos = self._all_sc_pos[condition][0]

        super(DRMGenTrig, self).__init__(
            quaternions=quaternions,
            sc_pos=sc_pos,
            det_number=self._det_number,
            ebin_edge_in=self._in_edge,
            mat_type=self._matrix_type,
            ebin_edge_out=self._out_edge,
            occult=self._occult,
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
