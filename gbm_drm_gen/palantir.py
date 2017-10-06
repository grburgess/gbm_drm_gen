import numpy as np
import gbmgeometry
import matplotlib.pyplot as plt
from matplotlib.projections.geo import GeoAxes
import healpy as hp

import astropy.units as u

import collections


from balrog_healpix_map import BALROGHealpixMap


class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi
    Shifts labelling from -180,180 to 0-360"""
    extrapi=False

    def __init__(self,*a,**aa):
        if 'extrapi' in aa:
            self.extrapi=aa['extrapi']
        del aa['extrapi']
        super(ThetaFormatterShiftPi,self).__init__(*a,**aa)

    def __call__(self, x, pos=None):
        if self.extrapi:
            x+=np.pi
        if x != 0:
            x *= -1
        if x < 0:
            x += 2*np.pi
        return GeoAxes.ThetaFormatter.__call__(self, x, pos)


def check_power_of_two(num):
    """Check that the given number is a power of 2"""

    return num != 0 and ((num & (num - 1)) == 0)


def pix_to_sky(idx, nside):
    """Convert the pixels corresponding to the input indexes to sky coordinates (RA, Dec)"""

    theta, phi = hp.pix2ang(nside, idx)

    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    return ra, dec


def sky_to_pix(ra, dec, nside):
    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)
    ipix = hp.ang2pix(nside, theta, phi)
    return ipix


_det_kwargs = {'colorbar': True,
               'cmap':'viridis',
               'show_earth':True,
               'width':18,
               'fov': 60
               }


class Palantir(object):
    def __init__(self, result, nside=64, trigdat=None, time=0.):


        """

        View teh 3ML BALROG results

        :param result: 3ML BALROG result
        :param nside: power of 2
        :param trigdat: optional trigdat files
        :param time: time of observation
        """

        self._balrog_map = BALROGHealpixMap(result, nside=nside)

        self._nside = nside

        self._map = self._balrog_map.map

        self._trigdat = trigdat

        if trigdat is not None:

            self._position_interpolator = gbmgeometry.PositionInterpolator(trigdat=trigdat)

            self._gbm = gbmgeometry.GBM(self._position_interpolator.quaternion(time),
                                        self._position_interpolator.sc_pos(time) * u.km)

        else:

            self._position_interpolator = None

        self._extra_maps = collections.OrderedDict()


    def add_healpix_map(self,healpix_map,name):


        self._extra_maps[name] = healpix_map


    def skymap(self, *dets, **kwargs):

        """

        :param dets: optional detectors
        :param cmap: colormap for location
        :param show_earth: show the earth points
        :param width: no idea
        :param fov: FOV of GBM detectors
        :return: figure
        """
        vmin = 0.
        vmax = self._map.max()

        xsize = 2000
        ysize = xsize / 2.

        for kw, val in kwargs.iteritems():

            _det_kwargs[kw] = val


        theta = np.linspace(np.pi, 0, ysize)


        phi = np.linspace(-np.pi, np.pi, xsize)
        longitude = np.radians(np.linspace(-180, 180, xsize))
        latitude = np.radians(np.linspace(-90, 90, ysize))

        # project the map to a rectangular matrix xsize x ysize
        PHI, THETA = np.meshgrid(phi, theta)
        grid_pix = hp.pixelfunc.ang2pix(self._nside, THETA, PHI)

        grid_map = self._map[grid_pix]

        # for width in [ 8.8]:
        # for width in [18., 12., 8.8]:

        fig, ax = plt.subplots(subplot_kw=dict(projection='mollweide'))

        # rasterized makes the map bitmap while the labels remain vectorial
        # flip longitude to the astro convention
        image = ax.pcolormesh(longitude[::-1],
                              latitude,
                              grid_map,
                              vmin=vmin,
                              vmax=vmax,
                              rasterized=True,
                              cmap=_det_kwargs['cmap'], zorder=-32)





        for det in dets:


            ra, dec = self._get_detector_map(det, _det_kwargs['fov'])

            # pix = np.where(detector_map == 1)[0]
            # ra, dec = pix_to_sky(pix,self._nside)

            
            # need to shift the ra because.... dumb

            idx =ra >180.

            ra[idx] -=360.



            x = np.radians(ra)
            y = np.radians(dec)

            idx = y.argsort()

            ax.plot(x[idx], y[idx], '.', markersize=4)

        if _det_kwargs['show_earth']:
            earth_points = self._gbm.get_earth_points(False)


            ra = earth_points.ra.value

            # need to shift the ra because.... dumb

            idx = ra > 180.

            ra[idx] -= 360.



            ax.plot(np.radians(ra),
                    earth_points.dec.radian,
                    '.',
                    color='k',
                    alpha=.5,
                    zorder=-31)

        unit = r'$p$'


        for name, extra_map in self._extra_maps.iteritems():

            contour_worker = ContourFinder(extra_map,self._nside)




            for level in [.68,.95]:

                contour_idx = contour_worker.find_contour(containment_level=level)

                ra, dec = contour_worker.get_sky_coordinates(contour_idx)

                x = np.radians(ra - 180)
                y = np.radians(dec)

                idx = y.argsort()

                ax.plot(x[idx], y[idx], '.', markersize=3)



            #plt.contour(longitude[::-1], latitude, region, colors='jet', levels=levels)



        # colorbar
        if _det_kwargs['colorbar']:

            cb = fig.colorbar(image,
                              orientation='horizontal',
                              shrink=.6, pad=0.05,
                              ticks=[vmin, vmax])

            cb.ax.xaxis.set_label_text(unit)
            cb.ax.xaxis.labelpad = -8
            # workaround for issue with viewers, see colorbar docstring
            cb.solids.set_edgecolor("face")



            # graticule

        shift = 0.

        lon_g = ax.set_longitude_grid(60.)

        ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60., extrapi=False))

        width = _det_kwargs['width']

        if width < 10:
            lat_g = ax.set_latitude_grid(45)
            lon_g = ax.set_longitude_grid_ends(90)

        for g in ax.xaxis.get_gridlines() + ax.yaxis.get_gridlines():
            g.set_linestyle("dotted")
            g.set_color("black")
            g.set_alpha(0.5)

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10, color='k')
        ax.grid(True)

        return fig

    def _get_detector_map(self, det,fov):



        cir = self._gbm.detectors[det].get_fov(fov)
        ra = cir[0]
        dec = cir[1]

        return ra, dec


    def get_flux_upper_bound(self,ene_min,ene_max, percentile=95.):

        flux = self._result.get_point_source_flux(ene_min,
                                                  ene_max,
                                                  get_distribution=True)

        # now get the flux argument

        ub = np.percentile(flux['distribution'][0],percentile)

        return ub * u.Unit('erg/(cm2 s)')



    @property
    def healpix_map(self):

        return self._map

    @property
    def balrog_map(self):

        return self._balrog_map



class ContourFinder(object):
    def __init__(self, healpix_map, nside=128):
        """
        modeled from giacomov


        """

        self._nside = nside

        if not check_power_of_two(self._nside):
            raise RuntimeError("nside must be a power of 2.")

        #hpx_orig, header = hp.read_map(healpix_map, h=True, verbose=False)

        # Use power=-2 so the sum is still 1

        self._map = hp.pixelfunc.ud_grade(healpix_map, nside_out=self._nside, power=-2)

        # Check that the total probability is still equal to the input map

        assert abs(np.sum(healpix_map) - np.sum(self._map)) < 1e-3, "Total probability after resize has not been kept!"

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

        idx_prime = (cumsum <= containment_level)

        # Convert back to the "unordered space"

        idx = index_revsort[idx_prime]

        assert abs(np.sum(self._map[
                              idx]) - containment_level) < 1e-2, "Total prob. within containment is too far from requested value"

        return idx

    def get_sky_coordinates(self, indexes):
        ra, dec = pix_to_sky(indexes, self.nside)

        return ra, dec














