import numpy as np
import healpy as hp
from matplotlib.projections.geo import GeoAxes


class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi
    Shifts labelling from -180,180 to 0-360"""

    extrapi = False

    def __init__(self, *a, **aa):
        if "extrapi" in aa:
            self.extrapi = aa["extrapi"]
        del aa["extrapi"]
        super(ThetaFormatterShiftPi, self).__init__(*a, **aa)

    def __call__(self, x, pos=None):
        if self.extrapi:
            x += np.pi

        # theta_shifter(x)
        if x != 0:
            x *= -1
        if x < 0:
            x += 2 * np.pi

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


def theta_shifter(x):

    x[x != 0] *= -1
    # x[x<0] += 2*np.pi
