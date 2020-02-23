try:
    from threeML.utils.response import InstrumentResponse

except:
    from gbm_drm_gen.utils.response import InstrumentResponse


class BALROG_DRM(InstrumentResponse):
    def __init__(self, drm_generator, ra, dec):
        """
        
        :param drm_generator: BALROG DRM generator
        :param ra: RA of the source
        :param dec: DEC of the source
        """

        self._drm_generator = drm_generator

        self._drm_generator.set_location(ra, dec)

        super(BALROG_DRM, self).__init__(
            self._drm_generator.matrix,
            self._drm_generator.ebounds,
            self._drm_generator.monte_carlo_energies,
        )

    def set_location(self, ra, dec):
        """
        Set the source location
        :param ra: 
        :param dec: 
        :return: 
        """

        self._drm_generator.set_location(ra, dec)

        self._matrix = self._drm_generator.matrix

    def set_time(self, time):
        """
        set the time of the source
        :param time: 
        :return: 
        """

        self._drm_generator.set_time(time)

        self._matrix = self._drm_generator.matrix
