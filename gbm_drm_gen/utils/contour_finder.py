import numpy as np
from .general_utils import check_power_of_two, pix_to_sky
import healpy as hp


class ContourFinder(object):
    def __init__(self, healpix_map, nside=512):
        """
        modeled from giacomov


        """

        self._nside = nside

        if not check_power_of_two(self._nside):
            raise RuntimeError("nside must be a power of 2.")

        # hpx_orig, header = hp.read_map(healpix_map, h=True, verbose=False)

        # Use power=-2 so the sum is still 1

        self._map = healpix_map  # hp.read_map(healpix_map)

        # self._map = hp.pixelfunc.ud_grade(healpix_map, nside_out=self._nside, power=-2)

        # Check that the total probability is still equal to the input map

        # assert abs(np.sum(healpix_map) - np.sum(self._map)) < 1e-3, "Total probability after resize has not been kept!"

    @property
    def map(self):
        """Return the resized map"""

        return self._map

    @property
    def nside(self):
        """Return the nside of the re-sized map"""
        return self._nside

    @property
    def pixel_size(self):
        """Return average side of pixel in degrees for the re-sized map"""

        # Formula from the manual of the function HEALPIXWINDOW

        arcmin_to_degree = 1.0 / 60.0

        return np.sqrt(3.0 / np.pi) * 3600.0 / self._nside * arcmin_to_degree

    def find_contour(self, containment_level=0.9):
        """Return the *indexes* of the pixels in the map within the given containment level"""

        # Get indexes sorting the array in decreasing order

        index_revsort = np.argsort(self._map)[::-1]

        # Compute the cumulative sum of the probabilities in the ordered "space"

        cumsum = np.cumsum(self._map[index_revsort])

        # Find out which pixels are within the containment level requested in the ordered "space"

        idx_prime = cumsum <= containment_level

        # Convert back to the "unordered space"

        idx = index_revsort[idx_prime]

        assert (
            abs(np.sum(self._map[idx]) - containment_level) < 1e-2
        ), "Total prob. within containment is too far from requested value"

        return idx

    def get_sky_coordinates(self, indexes):
        ra, dec = pix_to_sky(indexes, self.nside)

        return ra, dec
