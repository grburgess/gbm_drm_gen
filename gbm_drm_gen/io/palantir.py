import numpy as np
import gbmgeometry
import matplotlib.pyplot as plt

import healpy as hp

import astropy.units as u

import collections


from .balrog_healpix_map import BALROGHealpixMap
from gbm_drm_gen.utils.general_utils import ThetaFormatterShiftPi, theta_shifter
from gbm_drm_gen.utils.contour_finder import ContourFinder


_det_kwargs = {
    "colorbar": True,
    "cmap": "viridis",
    "show_earth": True,
    "width": 18,
    "fov": 60,
}


class Palantir(object):
    def __init__(self, result, nside=64, trigdat=None, poshist=None, time=0.0):

        """

        View the 3ML BALROG results

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

            self._position_interpolator = gbmgeometry.PositionInterpolator(
                trigdat=trigdat
            )

            self._gbm = gbmgeometry.GBM(
                self._position_interpolator.quaternion(time),
                self._position_interpolator.sc_pos(time) * u.km,
            )

        elif poshist is not None:

            self._position_interpolator = gbmgeometry.PositionInterpolator(
                poshist=poshist
            )

            self._gbm = gbmgeometry.GBM(
                self._position_interpolator.quaternion(time),
                self._position_interpolator.sc_pos(time) * u.km,
            )

        else:

            self._position_interpolator = None

        self._extra_maps = collections.OrderedDict()
        self._extra_cmap = collections.OrderedDict()

    def add_healpix_map(self, healpix_map, name, cmap):

        self._extra_maps[name] = healpix_map
        self._extra_cmap[name] = cmap

    def skymap(self, *dets, **kwargs):

        """

        :param dets: optional detectors
        :param cmap: colormap for location
        :param show_earth: show the earth points
        :param width: no idea
        :param fov: FOV of GBM detectors
        :return: figure
        """
        vmin = 0.0
        vmax = self._map.max()

        xsize = 2000
        ysize = xsize / 2.0

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

        fig, ax = plt.subplots(subplot_kw=dict(projection="mollweide"))

        # rasterized makes the map bitmap while the labels remain vectorial
        # flip longitude to the astro convention
        image = ax.pcolormesh(
            longitude[::-1],
            latitude,
            grid_map,
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
            cmap=_det_kwargs["cmap"],
            zorder=-32,
        )

        for det in dets:

            ra, dec = self._get_detector_map(det, _det_kwargs["fov"])

            # need to shift the ra because.... dumb

            idx = ra > 180.0

            ra[idx] -= 360.0

            x = np.radians(ra)
            y = np.radians(dec)

            theta_shifter(x)

            idx = y.argsort()

            ax.plot(x[idx], y[idx], ".", markersize=4)

        if _det_kwargs["show_earth"]:

            earth_points = self._gbm.get_earth_points(False)

            ra = earth_points.ra.value

            # need to shift the ra because.... dumb

            idx = ra > 180.0

            ra[idx] -= 360.0

            x = np.radians(ra)
            y = np.radians(earth_points.dec)

            theta_shifter(x)

            ax.plot(x, y, ".", color="k", alpha=0.5, zorder=-31)

        unit = r"$p$"

        for name, extra_map in self._extra_maps.iteritems():

            # nside = hp.npix2nside(len(extra_map))
            #
            #
            #
            # grid_pix = hp.pixelfunc.ang2pix(nside, THETA, PHI)
            #
            # grid_map = extra_map[grid_pix]
            #
            # idx = grid_map>0.
            #
            # vmin = min(grid_map[idx])
            # vmax = extra_map.max()
            #
            #
            #
            # # rasterized makes the map bitmap while the labels remain vectorial
            # # flip longitude to the astro convention
            # image = ax.pcolormesh(longitude[idx][::-1],
            #                       latitude[idx],
            #                       grid_map[idx],
            #                       vmin=vmin,
            #                       vmax=vmax,
            #                       rasterized=True,
            #                       cmap=self._extra_cmap[name],alpha=.6)
            #

            contour_worker = ContourFinder(extra_map, hp.npix2nside(len(extra_map)))

            for level in [0.68]:

                contour_idx = contour_worker.find_contour(containment_level=level)

                ra, dec = contour_worker.get_sky_coordinates(contour_idx)

                idx = ra > 180.0

                ra[idx] -= 360.0

                x = np.radians(ra)
                y = np.radians(dec)

                idx = y.argsort()

                theta_shifter(x)

                ax.plot(x[idx], y[idx], ".", markersize=3, alpha=0.5)

            # plt.contour(longitude[::-1], latitude, region, colors='jet', levels=levels)

        # colorbar
        if _det_kwargs["colorbar"]:

            cb = fig.colorbar(
                image,
                orientation="horizontal",
                shrink=0.6,
                pad=0.05,
                ticks=[vmin, vmax],
            )

            cb.ax.xaxis.set_label_text(unit)
            cb.ax.xaxis.labelpad = -8
            # workaround for issue with viewers, see colorbar docstring
            cb.solids.set_edgecolor("face")

            # graticule

        shift = 0.0

        lon_g = ax.set_longitude_grid(60.0)

        ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60.0, extrapi=False))

        width = _det_kwargs["width"]

        if width < 10:
            lat_g = ax.set_latitude_grid(45)
            lon_g = ax.set_longitude_grid_ends(90)

        for g in ax.xaxis.get_gridlines() + ax.yaxis.get_gridlines():
            g.set_linestyle("dotted")
            g.set_color("black")
            g.set_alpha(0.5)

        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=10, color="k")
        ax.grid(True)

        return fig

    def _get_detector_map(self, det, fov):

        cir = [[ra, dec] for ra, dec in self._gbm.detectors[det].get_fov(fov)][0]
        ra = cir[0]
        dec = cir[1]

        return ra, dec

    def get_flux_upper_bound(self, ene_min, ene_max, percentile=95.0):

        flux = self._result.get_point_source_flux(
            ene_min, ene_max, get_distribution=True
        )

        # now get the flux argument

        ub = np.percentile(flux["distribution"][0], percentile)

        return ub * u.Unit("erg/(cm2 s)")

    @property
    def healpix_map(self):

        return self._map

    @property
    def balrog_map(self):

        return self._balrog_map
