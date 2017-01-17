from numpy import array, matrix, mean, zeros, where, arange, sqrt, arctan2, cos, sin, arccos
import numpy as np
import astropy.io.fits as fits

import ftran
import at_scat
from .detdatabase import DetDataBase
from ._geometry import ang2cart, is_occulted


class DRMGen(object):
    def __init__(self, trigdat, det, ebin_edge_in, mat_type=0, time=0, cspecfile=None, ebin_edge_out=None):
        """
        A generic GBM DRM generator. This can be inherited for specific purposes.
        It takes as input various spacecraft files to figure out geometry. The user can
        either supply custom output energy edges or they can be read from a cspec file.

        The great benefit here is the ability to make custom input edges if desired. The
        user needs to have the GBM database pointed at by the environment variable:
        BALROG_DB

        Additionally, matrix folding routines are supplied for spectral fitting. One can also
        simulate Poisson rates from supplied spectral models.

        :param trigdat: the path to a trigdat file
        :param det: the number (0-13) of the detector to be used
        :param ebin_edge_in: an array of energy **edges** lenght 1 longer than the num of bins
        :param mat_type: 0=direct 1=scattered 2=direct+scattered
        :param time: time relative to trigger to pull spacecraft position
        :param cspecfile: the cspecfile to pull energy output edges from
        :param ebin_edge_out: an array of output edges

        """

        # Attach inputs to object

        self.mat_type = mat_type
        self.in_edge = ebin_edge_in

        self.nobins_in = len(ebin_edge_in) - 1

        # If you want to use occulting
        self._occult = False

        ##### Initiate database loading:

        self.database = DetDataBase(lu[det])
        self.ein = np.zeros(self.nobins_in, dtype=np.float32)
        energ_lo = self.database.energ_lo
        energ_hi = self.database.energ_hi
        self.ein[:self.database.ienerg] = energ_lo
        self.ein[self.database.ienerg] = energ_hi[-1]

        edgeFlag = False
        if cspecfile != None:
            # Create the out edge energies
            tmp = fits.open(cspecfile)
            self.out_edge = np.zeros(129, dtype=np.float32)
            self.out_edge[:-1] = tmp['EBOUNDS'].data['E_MIN']
            self.out_edge[-1] = tmp['EBOUNDS'].data['E_MAX'][-1]

            tmp.close()
            edgeFlag = True

        elif ebin_edge_out != None:

            self.out_edge = ebin_edge_out

            edgeFlag = True

        if not edgeFlag:
            print "No output edges set! Use a cspecfile or your own!"
            return
        self.nobins_out = len(self.out_edge) - 1
        self.numEnergyBins = self.nobins_in
        self.numDetChans = self.nobins_out
        self.det = det

        # Some setup for matrix folding
        self.photonE = np.array(zip(self.in_edge[:-1], self.in_edge[1:]))
        self.chanWidth = self.photonE[:, 1] - self.photonE[:, 0]
        self.channelE = np.array(zip(self.out_edge[:-1], self.out_edge[1:]))
        self.chanMin = self.channelE[:, 0]
        self.chanMax = self.channelE[:, 1]
        self.meanChan = np.array(map(np.mean, self.channelE))
        self.meanPhtE = np.array(map(np.mean, self.photonE))
        # Initialize the energy selection to all MC channels
        self._energySelection = np.where(np.array([True] * len(self.meanChan)))[0].tolist()

        # Space craft stuff from TRIGDAT
        trigdat = fits.open(trigdat)
        trigtime = trigdat['EVNTRATE'].header['TRIGTIME']
        tstart = trigdat['EVNTRATE'].data['TIME'] - trigtime
        tstop = trigdat['EVNTRATE'].data['ENDTIME'] - trigtime
        condition = np.logical_and(tstart <= time, time <= tstop)

        self.qauts = trigdat['EVNTRATE'].data['SCATTITD'][condition][0]
        self.sc_pos = trigdat['EVNTRATE'].data['EIC'][condition][0]
        trigdat.close()

        self._scCoord(self.qauts, self.sc_pos)

        self._CreatePhotonEvalEnergies()

    def set_location(self, ra, dec):
        """
        Set the Ra and Dec of the DRM to be built. This invokes DRM generation as well.

        :param ra: ra in degrees
        :param dec: dec in degrees
        """

        self.ra = ra
        self.dec = dec

        if self._occult:
            if is_occulted(ra, dec, self.sc_pos):
                self.drm = self._occulted_DRM

            else:

                # get the spacecraft coordinates
                az, el = self._get_coords(ra, dec)

                # build the DRM
                self.drm = self._make_drm(az,
                                          el,
                                          self.geo_az,
                                          self.geo_el)
        else:
            # get the spacecraft coordinates
            az, el = self._get_coords(ra, dec)

            # build the DRM
            self.drm = self._make_drm(az,
                                      el,
                                      self.geo_az,
                                      self.geo_el)

        # go ahead and transpose it for spectal fitting, etc.
        self.drmTranspose = self.drm.T

    def get_energy_selection(self):
        """
        Obtain the user specified energy selection indices

        :returns energy selection: indices of output edges selected
        """
        return self._energySelection

    def get_model_cnts(self, ignore=False):
        """
        Get the model counts from forward folding. If the ignore flag is set then
        the user energy selection is ignored.

        :rtype: object
        :param ignore: bool to ignore user energy selection
        :returns: An array of model counts
        """

        self._ConvolveMatrix()
        if not ignore:
            modelCnts = self.counts[self._energySelection]

        else:
            modelCnts = self.counts
        return modelCnts

    def set_model(self, model):
        """
        Set the function to be used for DRM convolution
        :param model: a spectral model
        """
        self.model = model

    def set_params(self, params):
        """
        Set the parameters of the spectral model and invoke folding

        :param params: a [list] of params
        """
        self.params = params
        self._EvalModel()

    def select_energies(self, selection):
        """
        An array or nested array of energies (keV) is passed
        that specfies the energy selection(s) that will be used
        in the fitting procedure. This selection is memorized for
        plotting purposes later.

        :param selection: [10-900] or [[10-30],[40-900]] in KeV
        """

        # Make sure you've got an array
        selection = array(selection)
        self.selMins = []
        self.selMaxs = []

        # If the selection is a simple one
        if len(selection.shape) == 1:

            # Find out what the corresponding channels are
            tmp = map(self._GetChannel, selection)

            # Create a boolean array the length of all the channels
            tt = [False] * len(self.meanChan)

            # For all the good channels, flip the bits
            tt[tmp[0]:tmp[1] + 1] = [True] * (tmp[1] - tmp[0] + 1)

            # Record the max and min energies for plotting later
            self.emin = min(selection)
            self.emax = max(selection)

            # Record all the mins and maxes for plotting later
            self.selMins.append(self.chanMin[tmp[0]])
            self.selMaxs.append(self.chanMax[tmp[1]])

        # If instead we have a more complex selection...
        elif len(selection.shape) == 2:

            # Find the nested channel selections
            tmp = array(map(lambda x: [self._GetChannel(x[0]), self._GetChannel(x[1])], selection))

            # Create a boolean array the length of all the channels
            tt = [False] * len(self.meanChan)

            # For each selection, flip the bits
            for x in tmp:
                tt[x[0]:x[1] + 1] = [True] * (x[1] + 1 - x[0])

            # Record all the good selection min and maxes
            self.selMins = self.chanMin[tmp[:, 0]]
            self.selMaxs = self.chanMax[tmp[:, 1]]

            self.emin = selection.min()
            self.emax = selection.max()

        tt = array(tt)

        tt = where(tt)[0].tolist()

        # Save the energy selections
        self._energySelection = tt

    def sim_spectrum(self, params):
        """
        Simulates a spectrum with Poisson noise using the specified parameters.
        :param params: [list] of spectral parameters
        :returns: and array of simulated model counts corresponding to the params

        """

        self.set_params(params)
        meanCnts = self.get_model_cnts(ignore=True)

        simRates = array(map(self._PoissonRate, meanCnts))

        return simRates.T[0]

    def effective_area(self,type='total'):
        """




        :param type:
        :return:
        """


        if type == 'photon':

            return array(self.drm.sum(axis=1))

        if type=='chan':

            return array(self.drmTranspose.sum(axis=1))

        if type == 'total':

            return  self.drm.sum()


    def _scCoord(self, sc_quat, sc_pos):
        """
        GBM geometry calculations
        """

        self.scx = np.zeros(3)
        self.scy = np.zeros(3)
        self.scz = np.zeros(3)

        geodir = np.zeros(3)
        #        source_pos=zeros(3)
        #        source_pos_sc=zeros(3)

        # xx = sc_quat[0]**2
        # yy = sc_quat[1]**2
        # zz = sc_quat[2]**2
        # ww = sc_quat[3]**2
        # xy = sc_quat[0]*sc_quat[1]
        # yz = sc_quat[1]*sc_quat[2]
        # xz = sc_quat[0]*sc_quat[2]
        # wx = sc_quat[0]*sc_quat[3]
        # wy = sc_quat[1]*sc_quat[3]
        # wz = sc_quat[2]*sc_quat[3]



        # self.scx[0] = (xx - yy - zz + ww)
        # self.scx[1] = 2.0 * (xy + wz)
        # self.scx[2] = 2.0 * (xz - wy)
        # self.scy[0] = 2.0 * (xy - wz)
        # self.scy[1] = (-xx + yy - zz + ww)
        # self.scy[2] = 2.0 * (yz + wx)
        # self.scz[0] = 2.0 * (xz + wy)
        # self.scz[1] = 2.0 * (yz - wx)
        # self.scz[2] = (-xx - yy + zz + ww)


        self.scx[0] = (sc_quat[0] ** 2 - sc_quat[1] ** 2 - sc_quat[2] ** 2 + sc_quat[3] ** 2)
        self.scx[1] = 2.0 * (sc_quat[0] * sc_quat[1] + sc_quat[3] * sc_quat[2])
        self.scx[2] = 2.0 * (sc_quat[0] * sc_quat[2] - sc_quat[3] * sc_quat[1])
        self.scy[0] = 2.0 * (sc_quat[0] * sc_quat[1] - sc_quat[3] * sc_quat[2])
        self.scy[1] = (-sc_quat[0] ** 2 + sc_quat[1] ** 2 - sc_quat[2] ** 2 + sc_quat[3] ** 2)
        self.scy[2] = 2.0 * (sc_quat[1] * sc_quat[2] + sc_quat[3] * sc_quat[0])
        self.scz[0] = 2.0 * (sc_quat[0] * sc_quat[2] + sc_quat[3] * sc_quat[1])
        self.scz[1] = 2.0 * (sc_quat[1] * sc_quat[2] - sc_quat[3] * sc_quat[0])
        self.scz[2] = (-sc_quat[0] ** 2 - sc_quat[1] ** 2 + sc_quat[2] ** 2 + sc_quat[3] ** 2)

        geodir[0] = -self.scx.dot(sc_pos)
        geodir[1] = -self.scy.dot(sc_pos)
        geodir[2] = -self.scz.dot(sc_pos)

        denom = sqrt(geodir.dot(geodir))

        geodir /= denom

        geo_az = np.arctan2(geodir[1], geodir[0])

        if (geo_az < 0.0):
            geo_az += 2 * np.pi
        while (geo_az > 2 * np.pi):
            geo_az -= 2 * np.pi

        geo_el = np.arctan2(sqrt(geodir[0] ** 2 + geodir[1] ** 2), geodir[2])

        self.geo_el = 90 - np.rad2deg(geo_el)

        self.geo_az = np.rad2deg(geo_az)

        # Also setup the occulted matrix
        self._occulted_DRM = np.zeros((self.nobins_in, self.nobins_out))

    def _get_coords(self, ra, dec):
        source_pos_sc = zeros(3)
        source_pos = ang2cart(ra, dec)

        source_pos_sc[0] = self.scx.dot(source_pos)
        source_pos_sc[1] = self.scy.dot(source_pos)
        source_pos_sc[2] = self.scz.dot(source_pos)

        el = arccos(source_pos_sc[2])
        az = arctan2(source_pos_sc[1], source_pos_sc[0])

        if (az < 0.0):
            az += 2 * np.pi
        el = 90 - np.rad2deg(el)
        az = np.rad2deg(0. + az)

        return [az, el]

    def _make_drm(self, src_az, src_el, geo_az, geo_el):
        """
        The DRM generator. Should not be invoked by the user except via the set location member funtion.
        :param src_az:
        :param src_el:
        :param geo_az:
        :param geo_el:
        :return:
        """

        final_drm = np.zeros((self.nobins_in, self.nobins_out), dtype=np.float32)

        ## SKY Interpolation

        rlon = src_az
        rlat = src_el

        sf = np.arctan(1.) / 45.
        plat = src_el * sf
        plon = src_az * sf
        P = np.array([np.cos(plat) * np.cos(plon),
                      np.cos(plat) * np.sin(plon),
                      np.sin(plat)], np.float32)

        N = len(self.database.LEND)
        b1, b2, b3, i1, i2, i3 = ftran.trfind(0,
                                              P,
                                              N,
                                              self.database.X,
                                              self.database.Y,
                                              self.database.Z,
                                              self.database.LIST,
                                              self.database.LPTR,
                                              self.database.LEND)

        mat1 = self.database.get_rsp(self.database.millizen[i1 - 1],
                                     self.database.milliaz[i1 - 1])
        mat2 = self.database.get_rsp(self.database.millizen[i2 - 1],
                                     self.database.milliaz[i2 - 1])
        mat3 = self.database.get_rsp(self.database.millizen[i3 - 1],
                                     self.database.milliaz[i3 - 1])
        ## Interpolator on triangle




        sum = b1 + b2 + b3
        b1n = b1 / sum
        b2n = b2 / sum
        b3n = b3 / sum

        out_matrix = b1n * mat1 + b2n * mat2 + b3n * mat3

        n_tmp_phot_bin = 2 * self.nobins_in + self.nobins_in % 2
        tmp_phot_bin = np.zeros(n_tmp_phot_bin, dtype=np.float32)
        tmp_phot_bin[::2] = self.in_edge[:-1]
        tmp_phot_bin[1::2] = 10 ** ((np.log10(self.in_edge[:-1]) + np.log10(self.in_edge[1:])) / 2.)

        #### Atmospheric scattering
        if self.mat_type == 1 or self.mat_type == 2:
            theta_geo = 90. - geo_el
            phi_geo = geo_az
            theta_source = 90 - rlat
            phi_source = rlon

            ## Get new coordinates in the proper space
            gx, gy, gz, sl = ftran.geocords(theta_geo, phi_geo, theta_source, phi_source)
            lat = 180. - np.rad2deg(np.arccos(sl[0] * gz[0] + sl[1] * gz[1] + sl[2] * gz[2]))
            atscat_diff_matrix = np.zeros((n_tmp_phot_bin - 1, self.nobins_out))
            if lat <= self.database.lat_edge[-1] and (lat < self.database.lat_cent[-1]):
                coslat_corr = np.abs(np.cos(np.deg2rad(lat)))

                if lat <= self.database.lat_cent[0]:
                    il_low = 0
                    il_high = 1
                    l_frac = 0.0
                else:
                    idx = np.where(self.database.lat_cent < lat)[0][-1]
                    il_low = idx - 1
                    il_high = idx
                    l_frac = 1.0 - (lat - self.database.lat_cent[idx]) / (
                        self.database.lat_cent[idx + 1] - self.database.lat_cent[idx])

                # We now have to loop over all the at scat data
                # This could be sped up by moving this to FORTRAN
                tmp_out = np.zeros((len(self.database.theta_cent) * len(self.database.double_phi_cent) * 2,
                                    self.database.ienerg, self.nobins_out))

                NUMBER_OF_PROCESSES = 2

                #                pool = Pool(nodes=NUMBER_OF_PROCESSES)

                # its=[]
                # theatas=[]
                # ips=[]
                # phis=[]
                # gxs=[]
                # gys=[]
                # gzs=[]
                # ils=[]
                # ihs=[]
                # lfs=[]
                # for itheta, theta_u in enumerate(self.database.theta_cent):
                #     for iphi, dp in enumerate(self.database.double_phi_cent):
                #         for phi_u in dp: # get the postive and negative

                #             its.append(itheta)
                #             theatas.append(theta_u)
                #             ips.append(iphi)
                #             phis.append(phi_u)
                #             gxs.append(gx)
                #             gys.append(gy)
                #             gzs.append(gz)
                #             ils.append(il_low)
                #             ihs.append(il_high)
                #             lfs.append(l_frac)


                # tmp_out = self.pool.map(self.at_scat,its,theatas,ips,phis,gxs,gys,gzs,ils,ihs,lfs)
                # tmp_out=np.array(tmp_out)

                # TASKS = []
                # itr = 0
                # for itheta, theta_u in enumerate(self.database.theta_cent):
                #     for iphi, dp in enumerate(self.database.double_phi_cent):
                #         for phi_u in dp: # get the postive and negative

                #             TASKS.append((self.at_scat,(itheta,theta_u,iphi,phi_u,gx,gy,gz,il_low,il_high,l_frac)))

                # task_queue = mp.Queue()
                # done_queue = mp.Queue()

                # for task in TASKS:
                #     task_queue.put(task)

                # for i in range(NUMBER_OF_PROCESSES):
                #     mp.Process(target=self.worker, args=(task_queue, done_queue)).start()

                # for i in range(len(TASKS)):
                #     tmp_out[i]=done_queue.get()

                # for i in range(NUMBER_OF_PROCESSES):
                #     task_queue.put('STOP')


                # itr = 0
                # for itheta, theta_u in enumerate(self.database.theta_cent):
                #     for iphi, dp in enumerate(self.database.double_phi_cent):
                #         for phi_u in dp: # get the postive and negative

                #             tmp_out[itr]= self._at_scat(itheta,theta_u,iphi,phi_u,gx,gy,gz,il_low,il_high,l_frac)
                #             itr+=1

                # tmp_out = tmp_out.sum(axis=0)
                tmp_out = at_scat.get_at_scat(gx, gy, gz, il_low, il_high, l_frac, self.nobins_out, self.out_edge,
                                              self.database)

                tmp_out *= coslat_corr
                atscat_diff_matrix = ftran.atscat_highres_ephoton_interpolator(tmp_phot_bin,
                                                                               self.ein,
                                                                               tmp_out)

        ###################################


        new_epx_lo, new_epx_hi, diff_matrix = ftran.highres_ephoton_interpolator(tmp_phot_bin,
                                                                                 self.ein,
                                                                                 out_matrix,
                                                                                 self.database.epx_lo,
                                                                                 self.database.epx_hi,
                                                                                 self.database.ichan,
                                                                                 n_tmp_phot_bin)

        binned_matrix = ftran.echan_integrator(diff_matrix,
                                               new_epx_lo,
                                               new_epx_hi,
                                               self.database.ichan,
                                               self.out_edge)

        if self.mat_type == 1:
            binned_matrix = atscat_diff_matrix

        if self.mat_type == 2:
            binned_matrix[:-1, :] += atscat_diff_matrix

        # Integrate photon edge with trapazoid




        final_drm[:-1, :] = (binned_matrix[::2, :][:-1, :] / 2. + binned_matrix[1::2, :][:-1, :] + binned_matrix[2::2,
                                                                                                   :] / 2.) / 2.

        return final_drm


##### MATRIX FOLDING AND SPECTRAL FITTING



    def _at_scat(self, itheta, theta_u, iphi, phi_u, gx, gy, gz, il_low, il_high, l_frac):
        """

        :param itheta:
        :param theta_u:
        :param iphi:
        :param phi_u:
        :param gx:
        :param gy:
        :param gz:
        :param il_low:
        :param il_high:
        :param l_frac:
        :return:
        """
        dirx, diry, dirz, az, el = ftran.geo_to_space(theta_u, phi_u, gx, gy, gz)

        sf = np.arctan(1.) / 45.
        plat = el * sf
        plon = az * sf
        P = np.array([np.cos(plat) * np.cos(plon),
                      np.cos(plat) * np.sin(plon),
                      np.sin(plat)], np.float32)

        # Find a new interpolated matrix
        ist = 0

        N = len(self.database.LEND)
        b1, b2, b3, i1, i2, i3 = ftran.trfind(ist,
                                              P,
                                              N,
                                              self.database.X,
                                              self.database.Y,
                                              self.database.Z,
                                              self.database.LIST,
                                              self.database.LPTR,
                                              self.database.LEND)
        i_array = np.array([i1, i2, i3])

        dist1 = ftran.calc_sphere_dist(az, 90. - el,
                                       self.database.Azimuth[i1 - 1],
                                       self.database.Zenith[i1 - 1])
        dist2 = ftran.calc_sphere_dist(az, 90. - el,
                                       self.database.Azimuth[i2 - 1],
                                       self.database.Zenith[i2 - 1])
        dist3 = ftran.calc_sphere_dist(az, 90. - el,
                                       self.database.Azimuth[i3 - 1],
                                       self.database.Zenith[i3 - 1])

        dist_array = np.array([dist1, dist2, dist3])
        i1 = i_array[np.argmin(dist_array)]

        tmpdrm = self.database.get_rsp(self.database.millizen[i1 - 1],
                                       self.database.milliaz[i1 - 1])

        # intergrate the new drm
        direct_diff_matrix = ftran.echan_integrator(tmpdrm,
                                                    self.database.epx_lo,
                                                    self.database.epx_hi,
                                                    self.database.ichan,
                                                    self.out_edge)

        # Now let FORTRAN add the at scat to the direct
        tmp_out = ftran.sum_at_scat(direct_diff_matrix,
                                    self.database.at_scat_data[:, :, il_low, itheta, iphi],
                                    self.database.at_scat_data[:, :, il_high, itheta, iphi],
                                    l_frac)

        return tmp_out

    def _CreatePhotonEvalEnergies(self):
        """
        Creates fast evaluation energies. Adopted from the MFIT routine.
        """

        resFrac511 = 0.2  # lifted from RMFIT!
        resExp = -0.15  # lifted from RMFIT!

        binCenter = np.array(map(np.mean, self.photonE))

        resFrac = resFrac511 * (binCenter / 511.) ** resExp
        resFWHM = binCenter * resFrac
        numEchans = np.ones(len(binCenter))

        self.lowEval = self.chanWidth < resFWHM / 2.
        self.lowEvalWhere = where(self.lowEval)[0].tolist()

        self.medEval = self.chanWidth >= resFWHM / 2.
        self.medEvalWhere = where(self.medEval)[0].tolist()

        numEchans[self.medEval] = 3.
        self.highEval = self.chanWidth / 2. >= resFWHM / 3.
        self.highEvalWhere = where(self.highEval)[0].tolist()
        numEchans[self.highEval] = 7.

        self.lowEne = binCenter[self.lowEval]
        self.medEne = np.array(map(lambda x, y: [x - 0.333333 * y, x, x + 0.333333 * y]
                                   , binCenter[self.medEval],
                                   self.chanWidth[self.medEval]))
        self.highEne = np.array(map(lambda x, y: [x - 0.5 * y,
                                                  x - 0.333333 * y,
                                                  x - 0.16667 * y,
                                                  x,
                                                  x + 0.16667 * y,
                                                  x + 0.333333 * y, x - 0.5 * y]
                                    , binCenter[self.highEval]
                                    , self.chanWidth[self.highEval]))

    def _EvalModel(self):
        tmpCounts = zeros(len(self.photonE))

        lowRes = self.model(self.lowEne, *self.params)
        medRes = array(map(lambda x: sum(self.model(x, *self.params)) / 3., self.medEne))

        hiRes = array(map(lambda x: sum(self.model(x, *self.params)) / 7., self.highEne))

        tmpCounts[self.lowEval] = lowRes
        tmpCounts[self.medEval] = medRes
        tmpCounts[self.highEval] = hiRes

        self.vec = tmpCounts * self.chanWidth

    def _GetChannel(self, energy):
        """
        Private function that finds the channel for a given energy.
        """

        if energy < self.chanMin[0]:
            return 0
        elif energy > self.chanMax[-1]:
            return len(self.chanMax) - 1

        ch = 0
        for lo, hi in zip(self.chanMin, self.chanMax):

            if energy >= lo and energy <= hi:
                return ch
            else:
                ch += 1

    def _ConvolveMatrix(self):
        self.counts = np.dot(self.drmTranspose, self.vec)

    ######################################################
    def _PoissonRate(self, meanCnts):
        numPhotons = np.random.poisson(meanCnts, 1)
        return numPhotons

    @staticmethod
    def worker(input, output):
        for func, args in iter(input.get, 'STOP'):
            result = func(*args)
            output.put(result)







        ###################################################  ################


lu = ['n0', "n1", 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']
