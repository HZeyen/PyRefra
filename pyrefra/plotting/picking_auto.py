# -*- coding: utf-8 -*-
"""
Semi-automatic and automatic picking and line fitting.

Contains the following Class:
    AutoPickMixin
"""

import numpy as np
from PyQt5 import QtCore
from obspy.signal import trigger
from .regression import best_lines, l1_regres


class AutoPickMixin:
    """Semi-automatic and automatic picking and line fitting."""

    def findNearest2ndDerivative(self, v, n0, nwin, nav):
        """
        Find position of maximum second derivqtive nearest to a reference
        sample

        Parameters
        ----------
        v : 1D Numpy array float
            Data vector in which to search
        n0 : int
            number of sample around which point has to be searched
        nwin : int
            Length of window to each side of n0 where to search given in number
            of samples
        nav : int
            number of samples to be summed on each side of sample point to
            calculate an approximation of derivative

        For all samples n0-nwin <= i <= n0+nwin, the difference of the average
        slopes of values v[i:i+nav] and v[i-nav:i] is calculated, giving the
        first derivatives. Then the position is searched where the difference
        of successive first derivatives is maximum.

        Returns
        -------
        di : int
            number of samples before (negative) or after (positive)sample n0
            where second derivative is maximum

        """
        dif = []
        im = []
        n = nav+1
        x = np.arange(n)
        sx = np.sum(x)
        sx2 = sx*sx
        sxx = n*np.dot(x, x)
        for i in range(-nwin, nwin+1):
            ii = n0+i
            im.append(i)
            y1 = v[ii-nav:ii+1]
            y2 = v[ii:ii+nav+1]
            a1 = max((n*np.dot(x, y1)-sx*np.sum(y1))/(sxx-sx2), 0.)
            a2 = max((n*np.dot(x, y2)-sx*np.sum(y2))/(sxx-sx2), 0.)
            dif.append(a2-a1)
        im = np.array(im, dtype=int)
        dif = np.array(dif)
        di = im[np.argmax(dif)]
        return di

    def corrPick(self):
        """
        Function does semi-automatic picking using cross-correlation
        First, the user must make a manual pick (without uncertainty) on a
        good, representative trace. The "source function" which has to be
        correlated with the other traces is taken rather arbitrarily as the
        signal of the chosen trace 10ms on both sides of the manual pick.
        These 10ms are generally ok for near-surface data. if another value
        should be chosen, search for the line:
            t_half = 0.01
        and replace 0.01 by another time (in seconds).
        Then the function searches for all other traces the maximum of the
        cross-correlation function. The data of each trace to be compared with
        the reference trace are not the full trace but only a length of +/-20
        ms (as for the reference trace) around a guessed time. There are two
        possible ways of guessing this position:

        If there are already picks from other shots nearby, search for the
        nearest pick having the same offset (supposing that the arrival times
        change only slowly from one shot point to the next). If there are no
        picks available, take the same time as found for the last pick
        (initially, the pick of the reference trace, then the one from the
        last trace treated). For small offsets (actually up to <5m, very low
        apparent velocities are normal, which is why there, the times are
        change by +/-5ms depending on whether the offset increases or
        decreases. If this should be changed, search the following lines:

        off_lim = 5
        dt_pick = 0.005

        and replace off_lim by the distance limit (in meters) up to which very
        low velocities are supposed and dt_pick by the time shift per trace
        (in seconds) for thes short offsets. The picks are only searched for
        the traces visible in the actual zoom. This may be convenient if data
        are noisy such that the user may first pick near-offset traces, then
        filter the data and pick far-offsets. Uncertainty limits are fixed to 2
        samples for raw data and to 4 samples for filtered data. If you want to
        change this do it in the if structure:

        if self.filtered:
            umin = 4*self.dt
        else:
            umin = 2*self.dt

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "corr_pick"
        self.setHelp(self.cpick_text)
        figure = self.figs[self.fig_plotted]

        def onPress(event):
            """
            Come here if a mouse button has been pressed
            """
            self.d_time = 0
# Set some default values explained further up in the explanation of the
# function
            if self.main.utilities.filtered:
                self.umin = 2*self.unc_base*self.data.dt
            else:
                self.umin = self.unc_base*self.data.dt
# If the left mouse button was pressed (manual pick set), come here
# Search the trace where the pick has been defined
            if event.button == 1:
                ret = self.searchTrace(event.xdata)
                if not ret:
                    return
                trace = self.itrace
                sht = self.traces.shot[trace]
                stn = self.traces.receiver[trace]
# Set the reference time the the pick position
                time_ref = event.ydata
                self.time_ref = time_ref
                print(f"\nshot {sht+1}, receiver {stn+1}, "
                      + f"picked time {time_ref*1000:0.2f}ms")
                self.traces.npick[trace] = 1
                try:
                    self.traces.pick_times[trace][0] = time_ref
                    self.traces.pick_times_min[trace][0] = \
                        self.traces.pick_times[trace][0]-self.umin
                    self.traces.pick_times_max[trace][0] = \
                        self.traces.pick_times[trace][0]+self.umin
                except (KeyError, IndexError, NameError):
                    self.traces.pick_times[trace].append(time_ref)
                    self.traces.pick_times_min[trace].append(
                        self.traces.pick_times[trace][0]-self.umin)
                    self.traces.pick_times_max[trace].append(
                        self.traces.pick_times[trace][0]+self.umin)
                self.n_pick_ref = int(np.round((time_ref-self.data.t0)
                                               / self.data.dt, 0))
                trace_ref = self.n_tr
# Calculate the width of the zone to find maxima in the correlation function
# For this, search the first maximum behind the manual pick. The distance
# between the pick and this first maximum is taken as halfwidth of the
# maximum to be searched (the full width corresponds thus approximately to
# half a wavelenght of the signal)
                mapo, _, _, _ = self.main.utilities.min_max(
                    self.v[trace_ref, self.n_pick_ref:], half_width=10)
# half_width is distance of next maximum from pickes position in number of
# samples
                self.half = max(mapo[0], 10)*2
                print(f"Half-width: {self.half}")
                n_ref1 = 0
                # nd = int(np.round(self.half/self.data.dt, 0))
                n_ref1 = max(0, self.n_pick_ref-self.half)
                n_max = self.v.shape[1]
                n_ref2 = min(self.n_pick_ref+self.half, n_max)
# tr_ref will contain the "source signal"
                tr_ref = self.v[trace_ref, n_ref1:n_ref2]
                self.tr_ref = (tr_ref-np.mean(tr_ref))/np.std(tr_ref)
                n_ref = len(tr_ref)
                self.add_samp = n_ref1+n_ref
                off0 = self.traces.offset[trace]
                off_ref = off0
                self.off_ref = abs(off0)
                dtrac = abs(self.traces.offset[1]-self.traces.offset[0])
                if self.off_ref > 5.:
                    vmin = 2000.
                else:
                    vmin = 1000.
                self.shift_samp = int(dtrac/(vmin*self.data.dt))
                dm = 0
                dm_last = 0
# Loop over the traces left of the reference trace
                if trace_ref > 0:
                    for i in range(trace_ref-1, -1, -1):
                        if i < np.array(self.tr, dtype=int).min():
                            break
                        ddm = abs(dm-dm_last)
                        dm_last = dm
                        off_last = dm
                        if ddm > 0:
                            self.shift_samp = ddm
                        else:
                            if abs(off0 > 5.):
                                vmin = 2000.
                            else:
                                vmin = 1000.
                            self.shift_samp = int(dtrac/(vmin*self.data.dt))
                        dm, off0 = self.search_corr_pick(i, dm, off0)
                        if dm is None:
                            dm = dm_last
                            off0 = off_last
# Loop over the traces right of the reference trace
                if trace_ref < self.actual_number_traces:
                    off0 = off_ref
                    dm = 0
                    for i in range(trace_ref+1, self.actual_number_traces):
                        if i < np.array(self.tr, dtype=int).min():
                            break
                        dm_last = dm
                        off_last = dm
                        dm, off0 = self.search_corr_pick(i, dm, off0)
                        if dm is None:
                            dm = dm_last
                            off0 = off_last

# When all traces inside the zoom have been treated, leave the function
                self.end = True
# if any other than left button was pressed, cancel correlation picking
            else:
                print("\nCorrelation picking cancelled")
                self.setHelp(self.main_text)
                self.main.function = "main"
                return

# Search picks done earlier in nearby gathers an plot a line connecting them
        self.x_pk_help, self.t_pk_help, self.i_tr_help = self.searchNearPicks()
# Plot found picks as dotted green line
        x1 = self.x_pk_help[self.x_pk_help >= self.x_zoom_min]
        t1 = self.t_pk_help[self.x_pk_help >= self.x_zoom_min]
        x_pk = x1[x1 <= self.x_zoom_max]
        t_pk = t1[x1 <= self.x_zoom_max]
        if len(x_pk):
            self.axes[self.fig_plotted].plot(x_pk, t_pk, c="g", ls="--")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)

        self.end = False
        self.coor_x = []
        self.coor_y = []
        figure = self.figs[self.fig_plotted]
        self.line, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y,
                                                      animated=True)
        self.cidpress = self.line.figure.canvas.mpl_connect(
            'button_press_event', onPress)

# Wait for mouse event
        while self.end is False:
            QtCore.QCoreApplication.processEvents()

# Redraw record section with new picks
        self.PlotPicks.setEnabled(True)
        self.Tomography.setEnabled(True)
        self.MovePicks.setEnabled(True)
        self.drawNew(True)
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
        if len(x_pk):
            self.axes[self.fig_plotted].plot(x_pk, t_pk, c="g", ls="--")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)
        self.setHelp(self.main_text)
        self.traces.storePicks()
        self.main.function = "main"

    def search_corr_pick(self, ntrac, dm, off0):
        """
        Search picks of one trace nearest to the pick of the trace beside or of
        a pick done on another shot at the same offset and the nearest cmp.

        Parameters
        ----------
        ntrac : int
            Number of trace for which the pick is to be found.
        dm : int
            Shift of last calculated trace with respect to reference trace in
            number of samples.
        off0 : float
            Offset of the last calculated trace.

        Returns
        -------
        dm : int
            Shift of the pick with respect to the reference trace in number of
            samples. If the signal arrives earlier than in the reference trace,
            dm is negative. If trace does not have data, dm is None
        off : float
            offset of the actual trace.

        """
        trace = self.actual_traces[ntrac]
        off = self.traces.offset[trace]
        off_a = abs(off)
        # off_r = self.off_ref
        n0 = abs(int(self.data.t0/self.data.dt))
        n_act = self.n_pick_ref+dm
        n_max = self.v.shape[1]
# If offset < 1cm, set pick time to zero
        if off_a < 0.01:
            self.traces.npick[trace] = 1
            try:
                self.traces.pick_times[trace][0] = 0.
                self.traces.pick_times_min[trace][0] =\
                    -2*self.data.dt
                self.traces.pick_times_max[trace][0] =\
                    2*self.data.dt
            except (KeyError, IndexError, NameError):
                self.traces.pick_times[trace].append(0.)
                self.traces.pick_times_min[trace].append(
                    -2*self.data.dt)
                self.traces.pick_times_max[trace].append(
                    2*self.data.dt)
                dm = n0-self.n_pick_ref
            return dm, off
        nt_help = np.where(self.i_tr_help == ntrac)[0]
        help_flag = False
# Set window of the signal to be correlated with the reference signal (tr2)
# The window is defined as the position of the last pick (trace beside the
# actual one) +/- (self.shift_sam+self.half). self.half is defined in function
# corrPick as the half-period of the signal
        if nt_help > 0:
            n_act = int((
                self.t_pk_help[nt_help[0]]-self.data.t0)
                / self.data.dt)
            help_flag = True
            dm = n_act-self.n_pick_ref
        doff = off_a-abs(off0)
        if doff < 0:
            n_sam1 = int(max(0, n_act-self.shift_samp-self.half))
            n_sam2 = int(min(n_act+self.shift_samp+self.half, n_max))
        else:
            n_sam1 = int(max(0, n_act-self.shift_samp))
            n_sam2 = int(min(n_act+self.shift_samp+self.half, n_max))
        # t1 = (n_sam1*self.data.dt+self.data.t0)*1000.
        # t2 = (n_sam2*self.data.dt+self.data.t0)*1000.
        tr2 = self.v[ntrac, n_sam1:n_sam2]
        s = np.std(tr2)
# If trace has no data (std == 0) skip picking
        if np.isclose(s, 0.):
            return None, None
        tr2 = (tr2-np.mean(tr2))/s
# Do correlation
        cor = np.correlate(tr2, self.tr_ref, 'full')
# Find maximum of correlation function. If nearby picks exist,
# search relative correlation maximum nearest to that time
        max_pos, max_val, _, _ = \
            self.main.utilities.min_max(cor, half_width=int(self.half/4))
        if len(max_pos > 0):
            max_pos += 1+n_sam1-self.add_samp
            index = np.where(max_pos+self.n_pick_ref >= n0)[0]
            max_pos = max_pos[index]
        if len(max_pos) == 0:
            i_max_pos = np.argmax(cor)+1+n_sam1-self.add_samp
        else:
            if help_flag:
                i_max_pos = max_pos[np.argmin(np.abs(max_pos-dm))]
# If not, use absolute maximum
            else:
                if doff < 0:
                    index = np.where(max_pos-dm <= int(self.half/2))[0]
                else:
                    index = np.where(max_pos-dm >= -int(self.half/2))[0]
                if len(index) > 0:
                    mpos = max_pos[index]
                else:
                    mpos = max_pos-dm
                i_max_pos = mpos[np.argmin(np.abs(mpos-dm))]
        isamp = int((self.time_ref-self.data.t0)/self.data.dt)+i_max_pos
        # t3 = (isamp*self.data.dt+self.data.t0)*1000.
        a = self.v[ntrac, isamp]
# Find next local maximum going from the found point to larger times until
# the end of the window. The maximum is defind as the first point where the
# amplitude is larger than the 10 samples before and after.
        imx = 0
        for i in range(isamp, n_sam2):
            if self.v[ntrac, i] == self.v[ntrac, i-10:i+11].max():
                imx = i
                break
# If no local maximum is found, take the position of the maximum amplitude
# within the window
        if imx == 0:
            imx = np.argmax(self.v[ntrac, n_sam1:n_sam2])+n_sam1
        # t4 = (imx*self.data.dt+self.data.t0)*1000.
        # print(f"     a:{a}, amax:{amx}")
# From the found local maximum go backward until the amplitude is smaller than
# 10% of the amplitude of the found maximum
        amx = self.v[ntrac, imx]*0.1
        if abs(a) > amx:
            for i in range(imx, n_sam1-20, -1):
                if self.v[ntrac, i] < amx:
                    # t5 = (i*self.data.dt+self.data.t0)*1000.
                    # print(f"     tpk={t5:0.2f}")
                    break
            isamp = i
        # t6 = (isamp*self.data.dt+self.data.t0)*1000.
# Find point nearby where the second derivative is maximum
        nwin = int(self.half/4)
        nwin2 = nwin*2
        der2 = self.v[ntrac, isamp+5:isamp+nwin2+6]+self.v[ntrac, isamp-nwin2-5:isamp-4]-\
            2*self.v[ntrac, isamp-nwin:isamp+nwin+1]
        isamp += np.argmax(der2)-nwin
        dm = isamp - self.n_pick_ref
        # t7 = (isamp*self.data.dt+self.data.t0)*1000.
        # print(f"trace {ntrac}: offset: {off:0.2f}; search:("
        #       f"{t1:0.2f}-{t2:0.2f}), t[isamp]={t3:0.2f}, tmax={t4:0.2f}, "
        #       f"tpk1={t6:0.2f}, tpk2={t7:0.2f}")
        self.traces.pick_times[trace].\
            append(isamp*self.data.dt+self.data.t0)
        self.traces.npick[trace] += 1
        self.traces.pick_times_min[trace].\
            append(self.traces.pick_times[trace][-1]-self.umin)
        self.traces.pick_times_max[trace].\
            append(self.traces.pick_times[trace][-1]+self.umin)
        return dm, off

    def Sta_Lta(self):
        """
        Calculate first arrival times using Sta-Lta method

        Returns
        -------
        None.

        """

        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "sta_lta_pick"
        visible_flag = True
# set figure (simply shorter variable name)
        figure = self.figs[self.fig_plotted]
        results, okButton = self.main.dialog(
            ["Sta length [ms]", "Lta length [ms]",
             "Only visible traces (y/n)"],
            ["e", "e", "e"], [2, 100, "y"], "Sta-Lta parameter input")
        trig = []

        if okButton is False:
            print("Picking cancelled")
            return
        S_time = float(results[0])/1000
        L_time = float(results[1])/1000
        if results[2].lower() == "n":
            visible_flag = False
        nsta = int(S_time/self.data.dt)
        nlta = int(L_time/self.data.dt)
        if visible_flag:
            n1 = self.tr[0]
            n2 = self.tr[-1]+1
        else:
            n1 = 0
            n2 = np.size(self.v, 0)
        for i in range(n1, n2):
            trace = self.actual_traces[i]
            npi = self.traces.npick[trace]
# If a pick exists for this trace, skip search (the idea is that for some
# traces, picks may have been placed before manually or calling Sta_Lta with
# other settings
            if npi > 0:
                continue
# Check whether the trace has data
            if np.isclose(np.max(np.abs(self.v[i, :])), 0.):
                continue
            t_trig_min = abs(self.traces.offset[trace])/self.v_max_trigger
            n_trig_min = int((t_trig_min - self.data.t0)/self.data.dt)
            trg = trigger.recursive_sta_lta(self.v[i, :], nsta, nlta)
            trg[0:n_trig_min] = 0
            indices = np.arange(len(trg))
            j50 = indices[trg > max(trg)*0.5][0]
            j25 = indices[trg > max(trg)*0.45][0]
            j75 = indices[trg > max(trg)*0.55][0]
            dminus = j25-j50
            dplus = j75-j50
            dav = int((dminus+dplus)/2)
# Search next maximum and go back to near zero crossing
            alim = self.v[i, j50]*0.1
            for k in range(j50, 0, -1):
                if self.v[i, k] < alim:
                    j50 = k
                    # j25 = k+dminus
                    # j75 = k+dplus
                    j25 = k-dav
                    j75 = k+dav
                    break
            t50 = self.data.time[j50]
            t25 = min(self.data.time[j25], t50-0.0005)
            t75 = max(self.data.time[j75], t50+0.0005)
            self.traces.npick[trace] += 1
            self.traces.pick_times[trace].append(t50)
            self.traces.pick_times_max[trace].append(t75)
            self.traces.pick_times_min[trace].append(t25)
            trig.append(t50)
        self.background = figure.canvas.copy_from_bbox(figure.bbox)
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
        self.Tomography.setEnabled(True)
        self.PlotPicks.setEnabled(True)
        self.MovePicks.setEnabled(True)
        self.main.function = "main"
        self.traces.storePicks()

    def ampPick(self, data, istart=0, iend=0, half_width=7):
        """
        Function locates first arrivals by analysing full amplitudes.
        First all maxima and minima are determined.
        For each maximum (minimum), the function calculates first an amplitude
        parameter corresponding to amp_pos = max-(min1+min2) for maxima and
        amp_neg = max1+max2-min for minima,
        where 1 means the next extreme value before the actual one and
        2 the next extreme value after the actual one.
        Then a vector is calculated with the relative amplitudes, i.e.
        dmax[i] = amp_pos[i]/amp[i-1] and dmin = amp_neg[i]/amp_neg[i-1]
        The pick is located in a first try near the highest value of dmax,
        located later than istart, precisely at the position of this maximum
        minus half the distance between the corrresponding maximum position and
        the following minimum position. If the maximum value of dmin is earlier
        than that of dmax and its value is larger than the one of dmax, then
        the pick is switched to the position of maximum dmin, shifted again by
        half the distance between this minimum and the following maximum.
        Then as a last test, the function searches for a (near) zero-crossing
        between the chosen pick position and the corresponding maximum
        (minimum) of dmax (dmin) If there is one, it will be the final pick, if
        not the chosen pick is maintained.

        Parameters
        ----------
        data : numpy 1D array FLOAT
            Contains the data vector
        istart : Int, optional
            Number of sample of the earliest accepted pick position.
            If istart is 0, the minimum accepted pick time is 0 (i.e. -t0)
            The default is 0.
        iend : Int, optional
            Number of the last analysed sample.
            The default is 0.
            If iend == 0, all data are used from istart to the end
        half_width : Int, optional
            Used for searching the extreme positions.
            It corresponds to the number of samples analyzed to both sides of
            every data sample. See function min_max in class Utilities
            The default is 7. This values has been found to give the best
            results.

        Returns
        tpk : Float
            Time of the found pick [s]
        amp : Float
            Amplitude at pick position

        """

        N = len(data)
    # Check starting and end positions
        ie = int(iend)
        if ie == 0:
            ie = N
        zero_pos = self.nt_0
        if istart > 0:
            zero_pos = istart
# Search relative maxima and minima from beginning of data set to iend
        max_pos, max_val, min_pos, min_val =\
            self.main.utilities.min_max(data[:ie], half_width)
# Prepare vectors for amplitude parameters
# sums_pos will contain the sums for relative maxima. i.e.
# maximum - (minimum earler+minimum later)
# sums_neg will contain the sums for relative minima. i.e.
# (maximum earler+maximum later)-minimum
        lx = len(max_pos)
        ln = len(min_pos)
        sums_pos = np.zeros(lx)
        sums_neg = np.zeros(ln)
# Calculate the sums
# First relative minimum comes earlier than first relative maximum
# In this case, calculate sums_pos from the first maximum on. For later
# search of the pick for minima, the following maximum will have
# an index higher by 1 than the value of the positionof the mminimum (iadd_pos)
        if min_pos[0] < max_pos[0]:
            ipos = 0
            iadd_pos = 1
            while ipos < lx-1 and ipos < ln-1:
                sums_pos[ipos] = np.abs(max_val[ipos]-(min_val[ipos]
                                        + min_val[ipos+1]))
                if sums_pos[ipos] == 0:
                    sums_pos[ipos] = np.nan
                ipos += 1
            ipos = 1
            iadd_neg = 0
            while ipos < lx-1 and ipos < ln-1:
                sums_neg[ipos] = np.abs(-min_val[ipos]+(max_val[ipos-1]
                                        + max_val[ipos]))
                if sums_neg[ipos] == 0:
                    sums_neg[ipos] = np.nan
                ipos += 1
# First relative minimum comes later than first relative maximum. In this case,
# calculate sums_pos from the second maximum on. For later search of the pick
# for maxima, the following minimum will have an index higher by 1 than the
# value of the position of the maximum (iadd_neg)
        else:
            ipos = 1
            iadd_pos = 0
            while ipos < lx and ipos < ln:
                sums_pos[ipos] = np.abs(max_val[ipos]-(min_val[ipos-1]
                                        + min_val[ipos]))
                if sums_pos[ipos] == 0:
                    sums_pos[ipos] = np.nan
                ipos += 1
            ipos = 0
            iadd_neg = 1
            while ipos < lx-1 and ipos < ln-1:
                sums_neg[ipos] = np.abs(-min_val[ipos]+(max_val[ipos]
                                        + max_val[ipos+1]))
                if sums_neg[ipos] == 0:
                    sums_neg[ipos] = np.nan
                ipos += 1
# Calculate amplitude relations for maxima (dsums_pos) and minima (dsums_neg)
        dsums_pos = np.zeros_like(sums_pos)
        dsums_neg = np.zeros_like(sums_neg)
# For possible plotting purposes, the amplitude coefficients are normalized
# to the data maximum. This may be eliminated without changing the result
        max_where = np.where(max_pos > zero_pos)
        min_where = np.where(min_pos > zero_pos)
        if len(sums_pos) < 4 or len(sums_neg) < 4 or np.size(max_where) == 0\
                or np.size(min_where) == 0:
            return np.nan, np.nan
        dsums_pos[2:-2] = sums_pos[3:-1]/sums_pos[1:-3]
        dsums_neg[2:-2] = sums_neg[3:-1]/sums_neg[1:-3]
# Also the following normalization may be eliminated, it is only for
#   plotting purposes
#    mx_sums = max(np.max(dsums_pos),np.max(dsums_neg))
#    dsums_pos = dsums_pos*maxi/mx_sums
#    dsums_neg = dsums_neg*maxi/mx_sums
# Search the first relative maximum (minimum) after startimg position
        ini_pos = np.min(max_where)
        ini_neg = np.min(min_where)
# Search the pick position corresponding to the maximum value of the
#   amplitude relations for the maxima.
        nmx_pos = ini_pos
        if iadd_pos > 0:
            nmx_pos = np.argmax(dsums_pos[ini_pos:])+ini_pos
        else:
            nmx_pos = np.argmax(dsums_pos[ini_pos:])+ini_pos
# Search the pick position corresponding to the maximum value of the
#   amplitude relations for the minima.
        nmx_neg = ini_neg
        if iadd_neg > 0:
            nmx_neg = np.argmax(dsums_neg[ini_neg:])+ini_neg
        else:
            nmx_neg = np.argmax(dsums_neg[ini_neg:])+ini_neg
# The position of the relations for maxima is shifted by half
#   the distance between this maximum relation and the next relative minimum
#   (i.e., in principle, a quarter of a wavelength)
        pos = int(max_pos[nmx_pos]-(min_pos[nmx_pos+iadd_pos] -
                                    max_pos[nmx_pos])/2)
# The position of the relations for minima is shifted by half the
#   distance between this maximum relation and the next relative maximum
#   (i.e., in principle, a quarter of a wavelength)
        neg = int(min_pos[nmx_neg]-(max_pos[nmx_neg+iadd_neg] -
                                    min_pos[nmx_neg])/2)
        tpos = pos*self.data.dt+self.data.t0
        tneg = neg*self.data.dt+self.data.t0
# Now check which of the two pick positions to choose:
#   normally, the one of the maxima is used, supposing that most traces
#   have a positive first arrival (it will halp the program if the input
#   data are signed in this sense).
#   If, however, the pick calculated from the minima is earlier and has
#   in addition a higher relationship indicator take that position.
        if tneg < tpos and dsums_neg[nmx_neg] > dsums_pos[nmx_pos]:
            tpk = tneg
            npk1 = neg
            npk2 = int(min_pos[nmx_neg])
# Check wether there is a zero crossing or near zero crossing between
#   the shifted position of the minima pick and the corresponding minimum.
#   If one is found, this will be the final pick position.
            a_comp = min_val[nmx_neg]/100
            for i in range(npk2, npk1-1, -1):
                if np.abs(data[i]) < np.abs(a_comp) or a_comp*data[i] <= 0:
                    break
# Similar procededure as above if the pick is obtained from relations of
#   maxima
        else:
            tpk = tpos
            npk1 = pos
            npk2 = max_pos[nmx_pos]
            a_comp = max_val[nmx_pos]/100
            for i in range(npk2, npk1-1, -1):
                if np.abs(data[i]) < np.abs(a_comp) or a_comp*data[i] <= 0:
                    break
    # Calculate and return final pick time
        tpk = (i-1)*self.data.dt+self.data.t0
        amp = data[i]
        return tpk, amp

    def ampPicks(self):
        """
        Function searches first arrivals for a full record section. Fist,
        function Amp_pick is used to make a first guess. Then nlines regression
        lines are fitted through the obtained picks and picks being 1.5*sigma
        apart from the regression lines are relocated to the nearest maximum
        to the corresponding regression line. Traces with zero amplitude
        will have numpy.nan value for arrival time.

        Parameters
        ----------
        data : Numpy.float 2D array [ndata,ntraces]
            Array contains teh data of one record section, each trace in a row
        x : Numpy float 1D array [ntraces]
            Signed offset of trace position [m]
        dt : Float
            Time step of samples [s]
        t0 : Float
            Time of first sample [s].
        v_max : Float
            Maximum expected velocity [m/s].
            No picks are accepted before time = abs(x)/v_max
            If v_max==0, the limit for all traces is set to t=0
        sig_len : Float
            Length of signal to be treated after first accepted time [s]
        half_width : Int, optional
            Used for searching the extreme positions.
            It corresponds to the number of samples analyzed to both sides of
            every data sample. See function min_max in class Utilities
            The default is 5. This values has been found to give the best
            results.
        nlines : Int, optional
            Number of lines to be fitted to first-guess picks. The first line
            is a direct wave, i.e., it passes through the origin. the other
            lines are refracted waves. The slope has to increase for every
            line.
            See function Best_lines in Seg2_utility_routines
            The default is 3, minimum 2.
        dist1 : Float, optional
            For the calculation of the direct wave two lines are fitted, the
            one for the direct wave and the first refracted wave. For this,
            not all picks are used, but only those with a distance smaller
            than dist1 [m]
            The default is 15.

        Returns
        -------
        t_pick : Numpy.float 1D array [ntraces]
            Array contains the calculated arrival times for each trace [s].
            If no arrival time has been calculated, the value is nan.
        a_picks : Numpy.float 1D array [ntraces]
            Amplitude at the position of the calculated pick.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "amp_pick"
# set figure (simply shorter variable name)
        figure = self.figs[self.fig_plotted]
        self.keys_held = set()
        results, okButton = self.main.dialog(
            ["Half_width for max", "Number of lines fitted",
             "Initial distance [m]", "signal length [ms]",
             "Maximum velocity [m/s]"], ["e", "e", "e", "e", "e"],
            [5, 3, 15, 50, self.v_max_trigger], "Amp_Pick parameter input")

        if okButton is False:
            self.main.function = "main"
            return
        half_width = int(results[0])
        nlines = int(results[1])
        dist1 = float(results[2])
        sig_len = float(results[3])
        sig_len /= 1000
        self.v_max_trigger = float(results[4])
        visible_flag = False
        ntrac = self.actual_number_traces
        t_pick = np.empty(ntrac)
        t_pick.fill(np.nan)
        a_pick = np.zeros(ntrac)
        n_end = np.zeros(ntrac, dtype=int)
        n_len = int(sig_len/self.data.dt)
        nt0 = int(-self.data.t0/self.data.dt)
        if visible_flag:
            n1 = self.tr[0]
            n2 = self.tr[-1]+1
        else:
            n1 = 0
            n2 = ntrac
        for j in range(n1, n2):
            i = j
            itry = 0
            data = self.v[i, :]
            if np.std(self.v[i, :]) > 0:
                if self.v_max_trigger > 0:
                    n_start = int((np.abs(self.x[self.indices[i]]) /
                                   self.v_max_trigger-self.data.t0) /
                                  self.data.dt)
                else:
                    n_start = nt0
                n_end[j] = n_start+n_len
                ne = n_end[j]
                dne = n_end[j]-n_start
                while np.isnan(t_pick[j]):
                    itry += 1
                    t_pick[j], a_pick[j] = self.ampPick(data, n_start, ne,
                                                        half_width=half_width)
                    ne = min(ne+dne, len(data))
        for _ in range(2):
            x = np.array(self.x)[self.indices]
            _, _, _, _, _, y_points_regres, _, _ =\
                self.bestLines(x, t_pick, nlines, dist1, 1)
            dist = 1.5*np.nanstd(y_points_regres-t_pick)
            for i in range(ntrac):
                if i < n1 or i > n2:
                    continue
                if np.abs(t_pick[i]-y_points_regres[i]) > dist:
                    nt_reg = int((y_points_regres[i]-self.data.t0) /
                                 self.data.dt)
                    max_pos, max_val, min_pos, _ =\
                        self.main.utilities.min_max(data[:n_end[i]], 5)
                    nmx_pos = len(max_pos)
                    nmn_pos = len(min_pos)
                    dpos = np.zeros(nmx_pos)
                    if max_pos[0] > min_pos[0]:
                        for k in range(nmx_pos):
                            if k > nmn_pos-2:
                                break
                            dpos = int((min_pos[k+1]-max_pos[k])/2)
                    else:
                        for k in range(1, nmx_pos):
                            if k > nmn_pos-1:
                                break
                            dpos = int((min_pos[k]-max_pos[k])/2)
                    mxp = max_pos-dpos
                    k = np.argmin(np.abs(mxp-nt_reg))
                    nt_pik = mxp[k]
                    for kk in range(max_pos[k], mxp[k], -1):
                        if data[kk] < max_val[k]/100 or\
                                data[kk]*max_val[k] < 0:
                            nt_pik = kk
                            break
                    t_pick[i] = nt_pik*self.data.dt+self.data.t0
                    a_pick[i] = data[nt_pik]
        if self.main.utilities.filtered:
            unc = 2*self.unc_base*self.data.dt
        else:
            unc = self.unc_base*self.data.dt
        for j in range(ntrac):
            if j < n1 or j > n2:
                continue
            if np.isnan(t_pick[j]):
                continue
            i = j
            trace = self.actual_traces[i]
            self.traces.npick[trace] += 1
            self.traces.pick_times[trace].append(t_pick[j])
            self.traces.pick_times_max[trace].append(t_pick[j]+unc)
            self.traces.pick_times_min[trace].append(t_pick[j]-unc)
        self.background = figure.canvas.copy_from_bbox(figure.bbox)
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
        self.Tomography.setEnabled(True)
        self.PlotPicks.setEnabled(True)
        self.MovePicks.setEnabled(True)
        self.traces.storePicks()
        self.main.function = "main"

    def bestLines(self, x_points, y_points, n_lines=2, xoff1=0, Lnorm=2,
                  origin=True, refra=True):
        """
        Function to calculate best fitting series of n_lines lines through
        points (x_points, y_points).
        First line is supposed to pass through the origin.
        Slopes have to decrease monotonously
        The function calculates first the unique absolute positions of
        x_points, which means that the lines are supposed to be symmetric
        around the origin. The function then searches the optimum points where
        to break slopes If more than 2 lines have to be fitted, it is possible
        to give a maximum distance for the points to be fitted for the two
        first lines.

        Parameters
        ----------
        x_points : numpy.float 1D array
            Contains X coordinates of points to be fitted
        y_points : numpy.float 1D array
            Contains Y coordinates of points to be fitted
        n_lines : int (optional, default: 2)
            Number of lines to be fitted
        xoff1 : float (optional, default: 0 (equivalent to all data))
            If given, for the fit of the first line, only points with offsets
            up to xoff1 are used.
        Lnorm : int (optional, default: 2)
            norm used for fitting (1 or 2)
        origin : bool(optional, default True)
            if True, first line must go through origin.
        refra : bool (optional, default True)
            if True, slopes are calculated for seismic refactions, i.e. slopes
                must be positive and decrease from one line to the next

        Returns
        -------
        best_off : numpy.float 1D array [n_lines+1]
            offsets of line breaks. First line goes from best_off[0] until
            best_off[1] inclusive second line goes from best_off[1] exclusive
            to best_off[2] inclusive etc.
        best_slope : numpy.float 1D array [n_lines]
            optimum slope of each line
        best_intercept : numpy.float 1D array [n_lines]
            optimum intercept of each line (first line has best_inter[0]=0)
        x_regres : numpy.float 1D array [number of unique offsets]
            unique offsets in increasing order
        y_regres : numpy.float 1D array [number of unique offsets]
            regression line y values for all unique offsets, ordered by offsets
        y_points_regres : numpy.float 1D array [number of data points]
            regression y coordinates at all positions x_points
        medians : Don't know anymore what this is...

        """

        return best_lines(x_points, y_points, n_lines, xoff1, Lnorm, origin,
                          refra)

    def L1_regres(self, x, y, par_ini=[0, 0]):
        """
        Function calculates linear regression line over points[x,y] using an L1
        norm

        Parameters
        ----------
        x : numpy.float 1D array
            Contains X coordinates of points to be fitted
        y : numpy.float 1D array
            Contains X coordinates of points to be fitted
        par_ini : float 1D list with one or two values optional, default: [0,0]
            Initial guess (usually obtained from a linear regression fit)
            if len(par_ini)==1, it is supposed that the given value is an
            initial guess of the slope and that the regression line passes
            through the origin if len(par_ini)==2, first parameter is slope,
            second intercept and both parameters are optimized

        Returns
        -------
        param : numpy.float 1D array , same size as par_ini
            optimized parameters
            if len(par_ini)==1, param contains only one value, the slope
            if len(par_ini)==2, param contains [slope,intercept]

        """

        return l1_regres(x, y, par_ini)
