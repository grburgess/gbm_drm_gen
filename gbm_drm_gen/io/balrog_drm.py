try:
    from threeML.utils.OGIP.response import InstrumentResponse

except:
    from responsum import InstrumentResponse

import numpy as np
from gbm_drm_gen.utils.geometry import get_ang, ang2cart

    
class BALROG_DRM(InstrumentResponse):
    def __init__(self, drm_generator, ra, dec):
        """
        
        :param drm_generator: BALROG DRM generator
        :param ra: RA of the source
        :param dec: DEC of the source
        """

        self._drm_generator = drm_generator

        self._drm_generator.set_location(ra, dec)
        self._min_dist = np.deg2rad(.5) 
        

        super(BALROG_DRM, self).__init__(
            self._drm_generator.matrix,
            self._drm_generator.ebounds,
            self._drm_generator.monte_carlo_energies,
        )

        self._cache = []


    def _cache_key_to_value(self, key):
        pass

    def _check_cache(self, ra, dec):

        max_sep = 1e99

        found_key = None
        

        this_cart = ang2cart(ra, dec)
        
        for v in self._cache[::-1]:

            if get_ang(this_cart, v["cart"]) < self._min_dist:

                return v["matrix"]

        self._drm_generator.set_location(ra, dec)

        self._cache.append(dict(cart=this_cart, matrix=self._drm_generator.matrix))

        return self._drm_generator.matrix
    
        
    def set_location(self, ra, dec, cache=False):
        """
        Set the source location
        :param ra: 
        :param dec: 
        :return: 
        """
        if not cache:

            self._drm_generator.set_location(ra, dec)

            self._matrix = self._drm_generator.matrix
            self._matrix_transpose = self._matrix.T

        else:
            self._matrix  = self._check_cache(ra, dec)

            
    def set_time(self, time):
        """
        set the time of the source
        :param time: 
        :return: 
        """

        self._drm_generator.set_time(time)

        self._matrix = self._drm_generator.matrix
