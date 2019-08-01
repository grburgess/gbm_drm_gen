from threeML.plugins.DispersionSpectrumLike import DispersionSpectrumLike
from astromodels.functions.priors import Uniform_prior, Cosine_Prior
from balrog_drm import BALROG_DRM


import numpy as np


def ra_calc(x,y):
    return np.rad2deg(np.arctan2(y,x))


def dec_calc(z):
    return -(np.rad2deg(np.arccos(z)) - 90.  )






class BALROGLike(DispersionSpectrumLike):
    def __init__(self,
                 name,
                 observation,
                 drm_generator=None,
                 background=None,
                 time=0,
                 free_position=True,
                 verbose=True):
        """
        BALROGLike is a general plugin for fitting GBM spectra and locations at the same time


        :param name: plugin name
        :param observation: observed spectrum
        :param drm_generator: the drm generator for this 
        :param background: background spectrum
        :param time: time of the observation
        :param free_position: keep the position free
        :param verbose: the verbosity level of the plugin
        """

        self._free_position = free_position



        if drm_generator is None:

            # If a generator is not supplied
            # then make sure that there is a
            # balrog response

            assert isinstance(observation.response, BALROG_DRM), 'The response associated with the observation is not a BALROG'


        else:

            # here we will reset the response
            # this is violating the fact that
            # the response is provate

            balrog_drm =  BALROG_DRM(drm_generator,0.,0.)

            observation._rsp = balrog_drm



        super(BALROGLike, self).__init__(name, observation, background,
                                         verbose)

        # only on the start up

        self._rsp.set_time(time)

    def set_model(self, likelihoodModel):
        """
        Set the model and free the location parameters


        :param likelihoodModel:
        :return: None
        """

        # set the standard likelihood model

        super(BALROGLike, self).set_model(likelihoodModel)

        # now free the position
        # if it is needed

        if self._free_position:

            if self._verbose:
                print('Freeing the position of %s and setting priors' % self.name)

            for key in self._like_model.point_sources.keys():
                pass
#                self._like_model.point_sources[key].position.ra.free = True
#                self._like_model.point_sources[key].position.dec.free = True

        # self._like_model.point_sources[
        #     key].position.ra.prior = Uniform_prior(
        #     lower_bound=0., upper_bound=360)
        # self._like_model.point_sources[
        #     key].position.dec.prior = Cosine_Prior(
        #     lower_bound=-90., upper_bound=90)

        ra = self._like_model.point_sources[key].position.ra.value
        dec = self._like_model.point_sources[key].position.dec.value

        self._rsp.set_location(ra, dec)

    def get_model(self):

        # Here we update the GBM drm parameters which creates and new DRM for that location
        # we should only be dealing with one source for GBM

        # update the location

        if self._free_position:

            # assumes that the is only one point source which is how it should be!
            # ra, dec = self._like_model.get_point_source_position(0)  

            xyz = np.array([self._like_model.x.value,
                            self._like_model.y.value,
                            self._like_model.z.value])

            x, y, z = xyz/np.linalg.norm(xyz)

            ra = ra_calc(x,y)

            dec = dec_calc(z)

            if ra < 0:
                ra +=360.

            
            

            
            self._rsp.set_location(ra, dec)

        return super(BALROGLike, self).get_model()

    @classmethod
    def from_spectrumlike(cls, spectrum_like, time, drm_generator=None,free_position=True):
        """
        Generate a BALROGlike from an existing SpectrumLike child
        
        
        :param spectrum_like: the existing spectrumlike
        :param time: the time to generate the RSPs at
        :param drm_generator: optional BALROG DRM generator
        :param free_position: if the position should be free
        :return: 
        """

        return cls(spectrum_like.name, spectrum_like._observed_spectrum,drm_generator,
                   spectrum_like._background_spectrum, time, free_position,
                   spectrum_like._verbose)
