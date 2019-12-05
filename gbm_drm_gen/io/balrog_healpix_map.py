import healpy as hp
import numpy as np
import copy


class BALROGHealpixMap(object):
    def __init__(self, analysis_results, nside=32):

        # assuming there is only one source!
        self._point_source_name = analysis_results.optimized_model.get_point_source_name(
            0
        )

        self._nside = nside

        self._bras = analysis_results.get_variates(
            "%s.position.ra" % self._point_source_name
        ).samples
        self._bdecs = analysis_results.get_variates(
            "%s.position.dec" % self._point_source_name
        ).samples

        self._generate_map()

    def _generate_map(self):

        hpix = hp.pixelfunc.ang2pix(self._nside, self.ra_healpix, self.dec, lonlat=True)

        newmap = np.zeros(hp.nside2npix(self._nside))

        for px in hpix:
            newmap[px] += 1.0 / float(self._bras.shape[0])

        self._map = newmap

    def write_map(self, filename):

        hp.write_map(filename, self._map, coord="C")

    @property
    def map(self):

        return self._map

    @property
    def ra(self):

        return self._bras

    @property
    def ra_healpix(self):

        # healpix has RA that goes from -180 to 180
        idx = self._bras >= 180.0

        bras_new = copy.copy(self._bras)

        bras_new[idx] = self._bras[idx] - 360

        return bras_new

    @property
    def dec(self):

        return self._bdecs
