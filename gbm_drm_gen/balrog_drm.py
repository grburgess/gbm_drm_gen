from threeML.plugins.OGIP.response import InstrumentResponse




class BALROG_DRM(InstrumentResponse):

    def __init__(self, drm_generator, ra, dec):


        self._drm_generator = drm_generator

        self._drm_generator.set_location(ra,dec)

        super(BALROG_DRM, self).__init__(self._drm_generator.matrix,
                                         self._drm_generator.ebounds,
                                         self._drm_generator.monte_carlo_energies)

    def set_location(self, ra, dec):

        self._drm_generator.set_location(ra, dec)

        self._matrix = self._drm_generator.matrix



    def set_time(self, time):

        self._drm_generator.set_time(time)

        self._matrix = self._drm_generator.matrix