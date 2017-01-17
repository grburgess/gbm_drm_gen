__author__ = "J. Michael Burgess"

import numpy as np
import warnings
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Grid

import gbmgeometry as gg

import astropy.coordinates as coord

from _background_fitter import BackgroundFitter


import pandas as pd

from IPython.display import HTML, display

from ._likelihood import li_and_ma


class TrigReader(object):
    """
    This class reads a GBM trigdat file and performs background fitting, source selection, and plotting.
    It is also fed to the Balrog when performing localization with trigdat data.

    :param triddat_file: string that is the path to the trigdat file you wish ot read
    :param fine: optional argument to use trigdat fine resolution data. Defaults to False
    """

    def __init__(self, trigdat_file, fine=False, time_resolved=False):

        self._backgroundexists = False
        self._sourceexists = False
        self._time_resolved = time_resolved
        # Read the trig data file and get the appropriate info

        trigdat = fits.open(trigdat_file)
        self._filename = trigdat_file
        self.out_edge_bgo = np.array([150., 400.0, 850.0, 1500.0, 3000.0, 5500.0, 10000.0, 20000.0, 50000.0],
                                     dtype=np.float32)
        self.out_edge_nai = np.array([3.4, 10.0, 22.0, 44.0, 95.0, 300.0, 500.0, 800.0, 2000.], dtype=np.float32)
        self.binwidth_bgo = self.out_edge_bgo[1:] - self.out_edge_bgo[:-1]
        self.binwidth_nai = self.out_edge_nai[1:] - self.out_edge_nai[:-1]

        # Get the times
        evntrate = "EVNTRATE"

        self.trigtime = trigdat[evntrate].header['TRIGTIME']
        self.tstart = trigdat[evntrate].data['TIME'] - self.trigtime
        self.tstop = trigdat[evntrate].data['ENDTIME'] - self.trigtime

        self.rates = trigdat[evntrate].data['RATE']

        # Now fix the shape because GBM stored this shit
        # so only fucking IDL can read it. Jesus Christ!
        num_times = len(self.tstart)
        self.rates = self.rates.reshape(num_times, 14, 8)

        # Obtain the positional information
        self.qauts = trigdat[evntrate].data['SCATTITD']  # [condition][0]
        self.sc_pos = trigdat[evntrate].data['EIC']  # [condition][0]

        # Get the flight software location
        self._fsw_ra = trigdat["PRIMARY"].header["RA_OBJ"]
        self._fsw_dec = trigdat["PRIMARY"].header["DEC_OBJ"]
        self._fsw_err = trigdat['PRIMARY'].header['ERR_RAD']

        # Clean up
        trigdat.close()

        # Sort out the high res times because they are dispersed with the normal
        # times.


        # The delta time in the file.
        # This routine is modeled off the procedure in RMFIT.
        myDelta = self.tstop - self.tstart
        self.tstart[myDelta < .1] = np.round(self.tstart[myDelta < .1], 4)
        self.tstop[myDelta < .1] = np.round(self.tstop[myDelta < .1], 4)

        self.tstart[~(myDelta < .1)] = np.round(self.tstart[~(myDelta < .1)], 3)
        self.tstop[~(myDelta < .1)] = np.round(self.tstop[~(myDelta < .1)], 3)

        if fine:

            # Create a starting list of array indices.
            # We will dumb then ones that are not needed

            all_index = range(len(self.tstart))

            # masks for all the different delta times and
            # the mid points for the different binnings
            temp1 = myDelta < .1
            temp2 = np.logical_and(myDelta > .1, myDelta < 1.)
            temp3 = np.logical_and(myDelta > 1., myDelta < 2.)
            temp4 = myDelta > 2.
            midT1 = (self.tstart[temp1] + self.tstop[temp1]) / 2.
            midT2 = (self.tstart[temp2] + self.tstop[temp2]) / 2.
            midT3 = (self.tstart[temp3] + self.tstop[temp3]) / 2.

            # Dump any index that occurs in a lower resolution
            # binning when a finer resolution covers the interval
            for indx in np.where(temp2)[0]:
                for x in midT1:
                    if self.tstart[indx] < x < self.tstop[indx]:
                        try:

                            all_index.remove(indx)
                        except:
                            pass

            for indx in np.where(temp3)[0]:
                for x in midT2:
                    if self.tstart[indx] < x < self.tstop[indx]:
                        try:

                            all_index.remove(indx)

                        except:
                            pass
            for indx in np.where(temp4)[0]:
                for x in midT3:
                    if self.tstart[indx] < x < self.tstop[indx]:
                        try:

                            all_index.remove(indx)
                        except:
                            pass

            all_index = np.array(all_index)
        else:

            # Just deal with the first level of fine data
            all_index = np.where(myDelta > 1.)[0].tolist()

            temp1 = np.logical_and(myDelta > 1., myDelta < 2.)
            temp2 = myDelta > 2.
            midT1 = (self.tstart[temp1] + self.tstop[temp1]) / 2.

            for indx in np.where(temp2)[0]:
                for x in midT1:
                    if self.tstart[indx] < x < self.tstop[indx]:

                        try:

                            all_index.remove(indx)

                        except:
                            pass

            all_index = np.array(all_index)

        # Now dump the indices we do not need
        self.tstart = self.tstart[all_index]
        self.tstop = self.tstop[all_index]
        self.qauts = self.qauts[all_index]
        self.sc_pos = self.sc_pos[all_index]
        self.rates = self.rates[all_index, :, :]

        # Now we need to sort because GBM may not have done this!

        sort_mask = np.argsort(self.tstart)
        self.tstart = self.tstart[sort_mask]
        self.tstop = self.tstop[sort_mask]
        self.qauts = self.qauts[sort_mask]
        self.sc_pos = self.sc_pos[sort_mask]
        self.rates = self.rates[sort_mask, :, :]

        self._pos_interp = gg.PositionInterpolator(trigdat=trigdat_file)

    def view_detector_movement(self):
        """

        :return:
        """

        time_gird = np.linspace(-50, 200, 25)

        gbm = gg.GBM(self._pos_interp.quaternion(-100))

        for t in time_gird:
            gbm.set_quaternion(self._pos_interp.quaternion(t))

            gbm.plot_pointing()

        gbm.plot_pointing(point=self.get_fsw_location())

    def get_fsw_location(self):
        """
        Returns the ra and dec determined by the FSW
        :return: ra,dec of FSW
        """

        return coord.SkyCoord(ra=self._fsw_ra, dec=self._fsw_dec, frame='icrs', unit='deg')

    def select_source(self, *intervals):
        """
        Select the source intervals that we want to use. An arbitrary number of intervals can be
        used. 

        :param intervals: 2-D lists of times to use
        :param time_resolved: whether or not to use separate time selections
        """

        # make a mask containg all the intervals we wish to use
        self._sourceexists = True

        mean_times = []
        all_intervals = []
        for interval in intervals:
            mask = np.logical_and(self.tstart >= interval[0], self.tstop <= interval[1])
            all_intervals.append(mask)

            mean_times.append(np.mean(interval))

        self.mean_times = mean_times

        # If there are multiple masks:
        if not self._time_resolved:
            source_mask = all_intervals[0]
            if len(all_intervals) > 1:
                for mask in all_intervals[1:]:
                    source_mask = np.logical_or(source_mask, mask)

            self.source_mask = source_mask

            # Compute the counts per channel in each detector (rate * exposure)
            # Each interval is summed over. This produces a matrix [num_dets,num_channels]
            exposure = (self.tstop[source_mask] - self.tstart[source_mask])  # need deadtime correction!!!
            self.exposure = exposure

            # print "Total Exposure:"
            # print exposure.sum()

            tmp = []
            for i, rate in enumerate(self.rates[source_mask, :, :]):
                tmp2 = rate * exposure[i]
                tmp.append(tmp2)
            tmp = np.array(tmp)

            self.counts_per_channel = tmp.sum(
                axis=0)  # (self.rates[source_mask,:,:] * exposure.reshape( (1,len(exposure))) ).sum(axis=0)
            self.rates_per_channel = self.rates[source_mask, :, :].sum(axis=0)

            # Now compute the background counts and errors
            self.bkg_counts = np.zeros((14, 8))
            self.bkg_err = np.zeros((14, 8))

            for i, det_poly in enumerate(self._all_polys):

                # for each channel, sum over the intervals
                for j, poly in enumerate(det_poly):

                    bkgcnt = 0.
                    bkgerr = 0.

                    for interval in zip(self.tstart[source_mask], self.tstop[source_mask]):
                        bkgcnt += poly.integral(interval[0], interval[1])
                        bkgerr += (poly.integralError(interval[0], interval[1])) ** 2

                    self.bkg_counts[i, j] = bkgcnt
                    self.bkg_err[i, j] = np.sqrt(bkgerr)

        elif self._time_resolved:

            # If we want to instead use multiple intervals
            # we will scroll through them and add them separately


            self.num_intervals = len(all_intervals)
            self.source_mask = all_intervals
            self.counts_per_channel = []
            self.rates_per_channel = []
            self.bkg_counts = np.zeros((self.num_intervals, 14, 8))
            self.bkg_err = np.zeros((self.num_intervals, 14, 8))
            self.exposure = []
            for k, source_mask in enumerate(all_intervals):

                exposure = (self.tstop[source_mask] - self.tstart[source_mask])  # need deadtime correction!!!
                self.exposure.append(exposure)

                # print "Total Exposure:"
                # print exposure.sum()

                tmp = []
                for i, rate in enumerate(self.rates[source_mask, :, :]):
                    tmp2 = rate * exposure[i]
                    tmp.append(tmp2)
                tmp = np.array(tmp)

                self.counts_per_channel.append(tmp.sum(axis=0))
                self.rates_per_channel.append(self.rates[source_mask, :, :].sum(axis=0))

                for i, det_poly in enumerate(self._all_polys):

                    # for each channel, sum over the intervals
                    for j, poly in enumerate(det_poly):

                        bkgcnt = 0.
                        bkgerr = 0.

                        for interval in zip(self.tstart[source_mask], self.tstop[source_mask]):
                            bkgcnt += poly.integral(interval[0], interval[1])
                            bkgerr += (poly.integralError(interval[0], interval[1])) ** 2

                        self.bkg_counts[k, i, j] = bkgcnt
                        self.bkg_err[k, i, j] = np.sqrt(bkgerr)

    def view_spectrum(self, detector):
        """
        View the spectrum of a detector for a specified selection
        :param detector: the detector number to view
        """
        if not self._backgroundexists and self._sourceexists:
            print "Error: Must select a source!"
            return

        if detector < 12:

            edges = np.array(zip(self.out_edge_nai[:-1], self.out_edge_nai[1:]))

            width = self.out_edge_nai[1:] - self.out_edge_nai[:-1]
        else:
            edges = np.array(zip(self.out_edge_bgo[:-1], self.out_edge_bgo[1:]))

            width = self.out_edge_bgo[1:] - self.out_edge_bgo[:-1]

        fig = plt.figure(666)
        ax = fig.add_subplot(111)

        # detectors = np.array(detectors)

        # total_rates = self.counts_per_channel[detectors,:].sum(axis=0)
        total_counts = self.counts_per_channel[detector, :]  # .sum(axis=0)
        total_bkg = self.bkg_counts[detector, :]
        self.Step(ax, edges, total_counts / width)
        self.Step(ax, edges, total_bkg / width, col='r')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("MC Energy (keV)")
        ax.set_ylabel("Counts")

    def view_detector(self, detector, tmin=-100., tmax=100., chan=-1):
        """
        Plots the lightcurves of the selected detector.
        If a background fit has been made, then it will also plot this

        :param chan:
        :param detector: the detector number to view
        :param tmin: The minimum time of the x-axis
        :param tmax: The maximum time of the x-axis

        """
        if detector < 12:
            binwidth = self.binwidth_nai
        else:
            binwidth = self.binwidth_bgo

        fig = plt.figure(777)
        ax = fig.add_subplot(111)

        tbins = np.array(zip(self.tstart, self.tstop))

        if self._sourceexists:

            if self._time_resolved:

                for source_mask, color in zip(self.source_mask, ["blue", "green", "purple", "orange", "black",
                                                                 "cyan"]):  # Add more colors via a color map

                    tbins_selected = np.array(zip(self.tstart[source_mask], self.tstop[source_mask]))

                    if chan < 0:
                        self.Step(ax, tbins_selected, self.rates[source_mask, detector, :].sum(axis=1), col=color,
                                  fill=True)
                    else:
                        self.Step(ax, tbins_selected, self.rates[source_mask, detector, chan], col=color, fill=True)
                ax.set_xlim(tmin, tmax)

            else:

                tbins_selected = np.array(zip(self.tstart[self.source_mask], self.tstop[self.source_mask]))

                if chan < 0:
                    self.Step(ax, tbins_selected, self.rates[self.source_mask, detector, :].sum(axis=1), col="green",
                              fill=True)
                else:
                    self.Step(ax, tbins_selected, self.rates[self.source_mask, detector, chan], col="green", fill=True)
                ax.set_xlim(tmin, tmax)

        if chan < 0:
            self.Step(ax, tbins, self.rates[:, detector, :].sum(axis=1))
        else:
            self.Step(ax, tbins, self.rates[:, detector, chan])
        ax.set_xlim(tmin, tmax)

        if self._backgroundexists:

            if chan < 0:
                det_poly = self._all_polys[detector]
                bkg_rates = []
                bkg = []
                for tb in tbins:
                    tmpbkg = 0.
                    for poly in det_poly:
                        tmpbkg += poly.integral(tb[0], tb[1]) / (tb[1] - tb[0])
                    bkg.append(tmpbkg)
                ax.plot(self.meantime, bkg, 'r')
            else:
                det_poly = self._all_polys[detector]
                bkg_rates = []
                bkg = []
                for tb in tbins:
                    tmpbkg = 0.
                    tmpbkg += det_poly[chan].integral(tb[0], tb[1]) / (tb[1] - tb[0])
                    bkg.append(tmpbkg)
                ax.plot(self.meantime, bkg, 'r')
        ax.set_xlim(tmin, tmax)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Rate (cnts/s)")

    def view_lightcurves(self, tmin=-100, tmax=100, chan_min=3, chan_max=4):
        """
        Plots the lightcurves in each detector. If a background fit has been made, then it will also plot this

        :param chan_min: min chan to sum over
        :param chan_max: max chan to sum over
        :param tmin: The minimum time of the x-axis
        :param tmax: The maximum time of the x-axis

        """

        fig = plt.figure(123, (10, 10))
        grid = Grid(fig,
                    111,
                    nrows_ncols=(7, 2),
                    axes_pad=0.6,
                    share_y=False,
                    label_mode='all')

        tbins = np.array(zip(self.tstart, self.tstop))

        condition = np.logical_and(self.tstart >= tmin,
                                    self.tstop<=tmax)



        # If there is a background we wnat to grab the max above background
        # else we just use the total rate. This way the user has a better idea
        # of which detectors are brighter

        if self._backgroundexists:

            li_ma_significance = []
            maxes = []

            for i in range(12):
                bkg = []
                for tb in tbins[condition]:
                    tmpbkg = 0.
                    for poly in self._all_polys[i][chan_min:chan_max]:
                        tmpbkg += poly.integral(tb[0], tb[1]) / (tb[1] - tb[0])
                    bkg.append(tmpbkg)

                bkg = np.array(bkg)

                tot = (self.rates[condition, i, chan_min:chan_max].sum(axis=1)).max()

                S = li_and_ma(tot, bkg)



                li_ma_significance.append(S.max())

                maxes.append(S.max())

            maxes = np.array(maxes)

            maxes /= maxes.max()


        else:

            maxes = [(self.rates[condition, det, chan_min:chan_max]).sum(axis=1).max() for det in range(12)]

        max_index = np.argsort(maxes)

        color_itr = np.linspace(0, 1, 12)

        for detector in range(14):
            if detector < 12:

                if self._backgroundexists:
                    color = plt.cm.winter_r(maxes[detector])




                else:
                    color = plt.cm.winter_r(color_itr[max_index[detector]])

                self.Step(grid[detector], tbins, self.rates[:, detector, chan_min:chan_max].sum(axis=1), col=color)
            else:
                self.Step(grid[detector], tbins, self.rates[:, detector, chan_min:chan_max].sum(axis=1))

            grid[detector].set_xlim(tmin, tmax)

            if self._backgroundexists and detector < 12:
                grid[detector].set_title("Detector:{} $\sigma$:{:.2f}".format(detector, li_ma_significance[detector]))

            else:

                grid[detector].set_title("Detector:{}".format(detector))

        if self._backgroundexists:

            bkg_rates = []

            for i, det_poly in enumerate(self._all_polys):

                bkg = []
                for tb in tbins:
                    tmpbkg = 0.
                    for poly in det_poly[chan_min:chan_max]:
                        tmpbkg += poly.integral(tb[0], tb[1]) / (tb[1] - tb[0])
                    bkg.append(tmpbkg)
                grid[i].plot(self.meantime, bkg, 'r')
                grid[i].set_xlim(tmin, tmax)

        grid[-1].set_xlabel("Time (s)")
        grid[-2].set_xlabel("Time (s)")
        for i in range(0, 13, 2):
            grid[i].set_ylabel("Rate (cnts/s)")

    def get_quats(self, time):

        condition = np.logical_and(self.tstart <= time, time <= self.tstop)

        return self.qauts[condition][0]

    def get_sc_pos(self, time):

        condition = np.logical_and(self.tstart <= time, time <= self.tstop)

        return self.sc_pos[condition][0]

    def fit_background(self, bkg_intervals=[]):

        self._backgroundexists = True
        allbkgmasks = []

        for bkgsel in bkg_intervals:
            mask = np.logical_and(self.tstart >= bkgsel[0], self.tstop <= bkgsel[1])
            allbkgmasks.append(mask)

        backgroundmask = allbkgmasks[0]

        # If there are multiple masks:
        if len(allbkgmasks) > 1:
            for mask in allbkgmasks[1:]:
                backgroundmask = np.logical_or(backgroundmask, mask)

        self.backgroundmask = backgroundmask
        # Now we will find the the best poly order unless the use specified one
        # The total cnts (over channels) is binned to 1 sec intervals



        self.meantime = np.array(map(np.mean, np.array(zip(self.tstart, self.tstop))))

        self._all_polys = []
        # Fit all the bakcgrounds
        for i in range(14):
            self._all_polys.append(self._fit_bkg_per_detector(i))

    def determine_trigger_characteristics(self, time=0.,tmin=-10,tmax=100):
        """
        Determines brightest detectors, legal pairs, etc.
        :return:
        """

        if not self._backgroundexists:
            print "No BKG set!"
            return

        chan_min, chan_max = 3, 4

        tbins = np.array(zip(self.tstart, self.tstop))

        # filter times

        condition = np.logical_and(self.tstart >= tmin,
                                    self.tstop<=tmax)

        # Assuming a background is fit, we want to get the detectors with the most significance


        li_ma_significance = []
        maxes = []

        for i in range(12):
            bkg = []
            for tb in tbins[condition]:
                tmpbkg = 0.
                for poly in self._all_polys[i][chan_min:chan_max]:
                    tmpbkg += poly.integral(tb[0], tb[1]) / (tb[1] - tb[0])
                bkg.append(tmpbkg)

            bkg = np.array(bkg)

            tot = (self.rates[condition, i, chan_min:chan_max].sum(axis=1)).max()

            S = np.nan_to_num(li_and_ma(tot, bkg))

            li_ma_significance.append(S.max())

            maxes.append(S.max())

        maxes = np.array(maxes)

        li_ma_significance = np.array(li_ma_significance)

        # Which detectors have at least 5 sigma?

        five_sigma_detectors = li_ma_significance > 5.

        sort_id = np.argsort(li_ma_significance[five_sigma_detectors])

        # legal sets from the table

        legal_sets = np.array([[1, 3],
                               [0, 2, 5],
                               [1, 5, 10],
                               [0, 4, 5],
                               [3, 5, 8],
                               [1, 2, 3, 4],
                               [7, 9],
                               [6, 8, 11],
                               [4, 7, 11],
                               [6, 10, 11],
                               [2, 9, 11],
                               [7, 8, 9, 10]])

        # get the sets that have detectors with a lot of counts

        five_sigma_sets = legal_sets[five_sigma_detectors]

        self.test = five_sigma_sets
        # Get the intersection of the detectors

        # intersection = functools.reduce(np.intersect1d, five_sigma_sets)

        # Now we want to see which detectors see the flight software position
        # fov = 60 for now

        gbm = gg.GBM(quaternion=self._pos_interp.quaternion(time))

        good_dets = gbm.get_good_detectors(self.get_fsw_location(), fov=60.)

        good_dets = np.array([detLU[d] for d in good_dets])

        seen = []

        for d in np.arange(12)[five_sigma_detectors]:
            if d in good_dets:
                seen.append("Yes")
            else:
                seen.append("No")

        seen = np.array(seen)

        # now get the angular separation from FSW to each
        # detector

        keys = [detLU2[k] for k in np.arange(12)[five_sigma_detectors]]

        centers = gbm.get_centers(keys=keys)

        seps = np.array([c.separation(self.get_fsw_location()).value for c in centers])

        def color_FSW_det_yellow(val):

            color = []

            for v in val:

                if v == "Yes":

                    color.append('background-color: #FFF26C')
                else:

                    color.append('background-color: #FF6F58')

            return color

        # Make a DF of 5 sigma detects

        five_sigma_df = pd.DataFrame({"Detectors": np.arange(12)[five_sigma_detectors][sort_id],
                                      "Sigma": li_ma_significance[five_sigma_detectors][sort_id],
                                      "Sees FSW Loc": seen[sort_id],
                                      "FSW Separation": seps[sort_id]})

        # intersect_dets_df = pd.DataFrame({"Detectors": intersection})

        display(
            five_sigma_df.style.set_table_styles(styles).background_gradient(cmap='winter_r',subset="FSW Separation").apply(color_FSW_det_yellow, subset='Sees FSW Loc').set_caption(
                "5 Sigma Detectors"))

        # display(intersect_dets_df.style.apply(color_FSW_det_red).set_caption("5 Sigma Legal Set"))

        print "Yellow indicates near FSW location"

    @staticmethod
    def Step(ax, tBins, y, col='k', lw=1, ls='-', fill=False):

        x = []
        newY = []
        lastX = -9999
        for t, v in zip(tBins, y):

            if lastX != -9999:

                if lastX != t[0]:
                    newY.append(0)
                    newY.append(0)
                    x.append(lastX)
                    x.append(t[0])

            x.append(t[0])
            newY.append(v)
            x.append(t[1])
            newY.append(v)
            lastX = t[1]
        if fill:
            ax.fill_between(x, 0, newY, color=col, alpha=.6, linestyle=ls, linewidth=lw)
        else:
            ax.plot(x, newY, color=col, linewidth=lw, linestyle=ls)

        newY = np.array(newY)
        mask = newY > 0.
        minY = min(newY[mask]) - .05 * min(newY[mask])
        maxY = max(newY) + .05 * max(newY)
        ax.set_ylim(bottom=minY, top=maxY)

    def _fit_bkg_per_detector(self, detector):

        background_fitter = BackgroundFitter()
        rates = self.rates[:, detector, :]
        optimalPolGrade = background_fitter._fitGlobalAndDetermineOptimumGrade(rates.sum(axis=1)[self.backgroundmask],
                                                                               self.meantime[self.backgroundmask])

        polys = []

        for chan in rates.T:
            thispolynomial, cstat = background_fitter._fitChannel(chan[self.backgroundmask],
                                                                  self.meantime[self.backgroundmask], optimalPolGrade)
            polys.append(thispolynomial)

        return polys


detLU = dict(n0=0, n1=1, n2=2, n3=3, n4=4, n5=5, n6=6, n7=7, n8=8, n9=9, na=10, nb=11, b0=12, b1=13)
detLU2 = ('n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1')


## Table Styles


def hover(hover_color="#7AFC93"):
    return dict(selector="tr:hover",
                props=[("background-color", "%s" % hover_color)])


styles = [
    hover(),
    dict(selector="th", props=[("font-size", "110%"),
                               ("text-align", "center"),
                               ('color', '#EE226D'),
                               ('background-color', '#1C0E4B')]),
    dict(selector="td", props=[("font-size", "100%"),
                               ("text-align", "center"),
                               ('color', 'k')
                               ]),
    dict(selector="caption", props=[("caption-side", "top")])]
