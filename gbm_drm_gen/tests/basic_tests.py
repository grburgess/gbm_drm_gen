from gbm_drm_gen import *
from threeML.plugins.OGIP.response import OGIPResponse

import astropy.io.fits as fits

n6 = DRMGenTTE('glg_tte_n6_bn110721200_v00.fit',
               trigdat='glg_trigdat_all_bn110721200_v01.fit',
               mat_type=2,
               cspecfile='glg_cspec_n6_bn110721200_v00.pha')

n6.to_fits(ra=0,dec=0,'test',overwrite=True)

ogip_rsp = OGIPResponse('test_n6.rsp')
