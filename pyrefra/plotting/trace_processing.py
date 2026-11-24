# -*- coding: utf-8 -*-
"""
Gain, muting, polarity and display of original data.

Contains the following Class:
    ProcessingMixin
"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from .canvas_manager import newWindow


class ProcessingMixin:
    """Gain, muting, polarity and display of original data."""

    def original(self):
        """
        Plot original data of all shots (i.e. eliminate filtering and/or
        muting)

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.data.st = copy.deepcopy(self.data.st_ori)
        self.v_set = False
        self.drawNew(True)
        self.v_set = True
# Change help text to self module
        self.setHelp(self.main_text)
        self.main.utilities.filtered = False
        self.main.utilities.filtered_all = False

    def originalScreen(self):
        """
        Plot original data of the shot shown on screen
        (i.e. eliminate filtering and/or muting)

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        for t in self.actual_traces:
            nf = self.traces.file[t]
            nt = self.traces.trace[t]
            self.data.st[nf][nt].data = np.copy(self.data.st_ori[nf][nt].data)
        self.v_set = False
        self.drawNew(True)
        self.v_set = True
        self.main.utilities.filtered = False
# Change help text to self module
        self.setHelp(self.main_text)

    def changeAmp(self, sens):
        """
        Change plotting amplitude of data set
        The value is selftained for the whole run until called
        again, also when changing data set.

        Parameters
        ----------
        sens : int or float
            If > 0 increase amplitude by a factor 2,
            if < 0 decrease amplitude by a factor 1/sqrt(2)
        """
        self.main.function = "change_amp"
        if sens > 0:
            self.amp_plt *= 2
        else:
            self.amp_plt /= np.sqrt(2)
        print("Amp:", self.amp_plt)
        self.drawNew(True)
# Change help text to self module
        self.setHelp(self.main_text)
        self.main.function = "main"

    def tNorm(self):
        """
        Function set gain as trace normalized
        """
        answer = self.main.test_function()
        if not answer:
            return
        self.gain = "tnorm"
        self.trace_norm = True
        self.t_Gain.setEnabled(True)
        self.t_Norm.setEnabled(False)
        self.d_Gain.setEnabled(True)
        self.Agc.setEnabled(True)
        self.t_Gain.setChecked(False)
        self.t_Norm.setChecked(True)
        self.d_Gain.setChecked(False)
        self.Agc.setChecked(False)
        self.drawNew(True)
# Change help text to self module
        self.setHelp(self.main_text)

    def tGain(self):
        """
        Function set gain as time**time_gain
        """
        answer = self.main.test_function()
        if not answer:
            return
        results, okButton = self.main.dialog(
            ["Exponent for time gain"], ["e"], ["2"], "Time gain")
        if okButton is False:
            print("Time gain cancelled")
        else:
            self.time_gain = float(results[0])
            self.gain = "time"
            self.trace_norm = True
            self.t_Gain.setEnabled(True)
            self.t_Norm.setEnabled(True)
            self.d_Gain.setEnabled(True)
            self.Agc.setEnabled(True)
            self.t_Gain.setChecked(True)
            self.t_Norm.setChecked(False)
            self.d_Gain.setChecked(False)
            self.Agc.setChecked(False)
            self.drawNew(True)
# Change help text to self module
            self.setHelp(self.main_text)

    def dGain(self):
        """
        Function set gain as distance**dist_gain
        """
        answer = self.main.test_function()
        if not answer:
            return
        results, okButton = self.main.dialog(
            ["Exponent for distance gain"], ["e"], ["2"], "Distance gain")
        if okButton is False:
            print("Distance gain cancelled")
        else:
            self.dist_gain = float(results[0])
            self.gain = "dist"
            self.trace_norm = False
            self.t_Gain.setEnabled(True)
            self.t_Norm.setEnabled(True)
            self.d_Gain.setEnabled(True)
            self.Agc.setEnabled(True)
            self.t_Gain.setChecked(False)
            self.t_Norm.setChecked(False)
            self.d_Gain.setChecked(True)
            self.Agc.setChecked(False)
            self.drawNew(True)
# Change help text to self module
            self.setHelp(self.main_text)

    def AGC(self):
        """
        Function set gain as AGC
        """
        answer = self.main.test_function()
        if not answer:
            return
        results, okButton = self.main.dialog(
            ["AGC window length [ms]"], ["e"], ["100"], "AGC gain")
        if okButton is False:
            print("Distance gain cancelled")
        else:
            self.agc_window = float(results[0])*0.001
            self.n_agc_win = int(self.agc_window/self.data.dt)
            self.gain = "agc"
            self.trace_norm = True
            self.t_Gain.setEnabled(True)
            self.t_Norm.setEnabled(True)
            self.d_Gain.setEnabled(True)
            self.Agc.setEnabled(True)
            self.t_Gain.setChecked(False)
            self.t_Norm.setChecked(False)
            self.d_Gain.setChecked(False)
            self.Agc.setChecked(True)
            self.drawNew(True)
# Change help text to self module
            self.setHelp(self.main_text)

    def phasePlot(self):
        """
        Toggles flag for plotting phase angles for multi-component data

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        if self.phase_plot:
            self.phase_plot = False
        else:
            self.phase_plot = True

    def phaseCalc(self):
        """
        Plot of phase angles for multi-component data

        Returns
        -------
        None.

        """
        traces = np.array(self.actual_traces)[self.tr]
        shots = self.traces.shot[traces]
        receivers = self.traces.receiver[traces]
        xx = self.x_ori[self.tr]
        xr, xr_ind = np.unique(xx, return_index=True)
        sr = np.vstack((shots, receivers)).T[xr_ind]
        trac = np.array(traces)[xr_ind]
        xx = []
        v = []
        for i in range(len(trac)):
            sht = sr[i, 0]
            rec = sr[i, 1]
            pos = self.geom.rec_dict[rec]["unique"]
            traces = []
            components = []
            for j, t in enumerate(self.geom.pos_dict[pos]["rec"]):
                if t in receivers:
                    traces.append(t)
                    components.append(self.geom.pos_dict[pos]["comp"][j])
            traces = np.array(traces, dtype=int)
            ntr = len(traces)
            if ntr < 2:
                continue
            if ntr == 2:
                t1 = self.traces.sht_rec_dict[(sht, traces[0])]
                t2 = self.traces.sht_rec_dict[(sht, traces[1])]
                f1 = self.traces.file[t1]
                f2 = self.traces.file[t2]
                tf1 = self.traces.trace[t1]
                tf2 = self.traces.trace[t2]
                v1 = self.data.st[f1][tf1].data
                v2 = self.data.st[f2][tf2].data
                v1 = np.abs(scipy.signal.hilbert(v1))
                v2 = np.abs(scipy.signal.hilbert(v2))
                if components[0] == "Z" or components[0] == "V":
                    ang = np.arctan(v1/v2)
#                    ang = v2/v1
                else:
                    ang = np.arctan(v2/v1)
#                    ang = v1/v2
                # amin = np.quantile(ang,0.01)
                # amax = np.quantile(ang,0.99)
                # ang = np.clip(ang,amin,amax)
                xx.append(xr[i])
                v.append(list(ang))
            else:
                t1 = self.traces.sht_rec_dict[(sht, traces[0])]
                t2 = self.traces.sht_rec_dict[(sht, traces[1])]
                t3 = self.traces.sht_rec_dict[(sht, traces[2])]
                f1 = self.traces.file[t1]
                f2 = self.traces.file[t2]
                f3 = self.traces.file[t3]
                tf1 = self.traces.trace[t1]
                tf2 = self.traces.trace[t2]
                tf3 = self.traces.trace[t3]
                if components[0] == "Z" or components[0] == "V":
                    v1 = self.data.st[f1][tf1].data
                    v2 = self.data.st[f2][tf2].data
                    v3 = self.data.st[f3][tf3].data
                elif components[1] == "Z" or components[1] == "V":
                    v1 = self.data.st[f2][tf2].data
                    v2 = self.data.st[f1][tf1].data
                    v3 = self.data.st[f3][tf3].data
                else:
                    v1 = self.data.st[f3][tf3].data
                    v2 = self.data.st[f1][tf1].data
                    v3 = self.data.st[f2][tf3].data
                v1 = np.abs(scipy.signal.hilbert(v1))
                v2 = np.abs(scipy.signal.hilbert(v2))
                v3 = np.abs(scipy.signal.hilbert(v3))
                a = v2+v3
                ang = np.arctan(v1/a)
                xx.append(xr[i])
                v.append(list(ang))
                ang = np.arctan(v3/v2)
                xx.append(xr[i])
                v.append(list(ang))
        xx = np.array(xx)
        v = np.array(v)
        return xx, v

    def traceMute(self, QselfWindow):
        """
        Eliminate one trace which is chosen by a mouse click, i.e. set
        amplitudes to zero. The choice is selftained for the whole run
        until called again, also when changing data set.
        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "trace_mute"
        self.finish = False
        self.picked = False
        self.setHelp(self.trace_mute_text)

        def onPress(event):
            if event.button == 1:
                ret = self.searchTrace(event.xdata)
                if not ret:
                    return
                trace = self.itrace
                nf = self.traces.file[trace]
                nt = self.traces.trace[trace]
                if np.isclose(self.traces.amplitudes[trace], 0.):
                    self.traces.amplitudes[trace] = self.general_sign
                    try:
                        fac = self.data.receiver_corr_dict[self.p_r]["amp"]
                        if fac != 0:
                            self.traces.amplitudes[trace] *= fac
                    except (KeyError, IndexError, NameError):
                        pass
                    self.v[self.n_tr, :] = self.data.st[nf][nt].data *\
                        self.general_sign
                    self.axes[self.fig_plotted].plot(
                        self.x[self.n_tr], self.time[self.nt_mn+5], "*g")
                    print("trace ", self.n_tr, " recovered")
                else:
                    self.traces.amplitudes[trace] = 0
                    self.v[self.n_tr, :] *= 0
                    self.axes[self.fig_plotted].plot(
                        self.x[self.n_tr], self.time[self.nt_mn+5], "*r")
                    print("trace ", self.n_tr, " muted")
                self.canvas.draw()
            elif event.button == 3:
                print("finish trace_mute")
                self.lin.set_animated(False)
                self.finish = True
                self.lin.figure.canvas.mpl_disconnect(self.cidpress)
                self.drawNew(True)
# Change help text to self module
                self.setHelp(self.main_text)
                self.main.function = "main"
        self.x_coor = []
        self.y_coor = []
        self.lin, = self.axes[self.fig_plotted].plot(self.x_coor, self.y_coor,
                                                     animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect(
            'button_press_event', onPress)

    def muteAir(self, QselfWindow):
        """
        Routine zeroes data around the air wave arrival for all traces

        Parameters
        ----------
        QselfWindow

        Returns
        -------
        modified traces

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "Mute_air"
        results, okButton = self.main.dialog(
            ["Air velocity [m/s]", "Mute-time on each side [ms]",
             "Width sin**2 taper [ms]"], ["e", "e", "e"], [340, 10, 3],
            "Mute air wave")
        if okButton is False:
            print("Air mute cancelled")
            return
        air_speed = float(results[0])
        mute_half = float(results[1])/1000.
        mute = 2*mute_half
        sin_width = float(results[2])/1000.
        nmute = int(mute/self.data.dt)
        nsin = int(sin_width/self.data.dt)
        nfac = int(nmute+2*nsin)
        factor = np.zeros(nfac)
        factor[0:nsin] = (np.cos(0.5*np.pi/nsin*np.arange(nsin)))**2
        factor[nfac-1:nfac-nsin-1:-1] = factor[0:nsin]
        for i in range(self.traces.number_of_traces):
            nf = self.traces.file[i]
            nt = self.traces.trace[i]
            tair = np.abs(self.traces.offset[i])/air_speed
            tini = tair - mute_half
            nini = int((tini-self.data.t0)/self.data.dt-nsin)
            nend = nini+nfac
            nini_fac = 0
            nend_fac = nfac
            if nini < 0:
                nini_fac = -nini
                nini = 0
            elif nend > self.data.nsamp:
                nend_fac = nfac-(nend-self.data.nsamp)
                nend = self.data.nsamp
            self.data.st[nf][nt].data[nini:nend] *= factor[nini_fac:nend_fac]
        self.main.function = "main"
        self.v_set = False
        self.drawNew(True)
        self.setHelp(self.main_text)

    def muteBefore(self):
        """
        Mute data at times smaller than a line drawn interactively

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.muteTime(-1)

    def muteAfter(self):
        """
        Mute data at times larger than a line drawn interactively

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.muteTime(1)

    def muteTime(self, sign):
        """
        Function zeroes all samples arriving earlier than a user-traced line.
        The line is defined by a series of points given by the user through
        left clicks. The last point of the line is marked with a right click.

        Only the traces crossed by the line are treated, which makes it
        possible to treat e.g. just the (noisy) beginning of a single trace.

        Before zeroing, the function seraches the sample nearest to the traced
        line which has a relative minimum value or is nearest to zero crossing.
        Zeroing is done from this sample on in order to avoid a large
        amplitude step at the cutting point.

        Parameters
        ----------
        sign : int
            If >0: data are muted after mouse position, else before
        """
        self.main.function = "Mute_before"
        self.setHelp(self.mute_before_text)
# Define the line before which the samples are zeroed
        _, _, = self.followLine(False, 1, 1)
# (xline,yline) are the distance and time coordinates of the points defining
# the line
        xline = np.array(self.coor_x)
        yline = np.array(self.coor_y)
# Sort the points in increasing X order
        line_index = np.argsort(xline)
        ntraces = np.arange(len(self.x))
# Treat one segment after the other
        for i in range(len(line_index)-1):
            i0 = line_index[i]
            i1 = line_index[i+1]
# If too a point is double clicked, skip the second occurrence
            if i0 > 0:
                if xline[i0] == xline[line_index[i-1]]:
                    continue
# Search traces lying inside a line segment
            npos = ntraces[self.x >= xline[i0]]
# Not sure why the folloing line does not work, replaced by the following for
# loop
            nt = []
            for j in npos:
                if self.x[j] > xline[i1]:
                    break
                nt.append(j)
# ntreat contains the numbers of all traces found within the segment
            ntreat = np.array(nt, dtype=int)
# Treat the corresponding traces
# First find time where the line crosses the trace and the sample number
            for j in ntreat:
                nt_cut = yline[i0]+(self.x[j]-xline[i0]) /\
                    (xline[i1]-xline[i0])*(yline[i1]-yline[i0])
                ncut = int((nt_cut-self.data.t0)/self.data.dt)
# Find the sample after and including the line cut position which has a
# (relative) minimum amplitude
                test = False
                k = ncut
                while test is False:
                    if np.abs(self.v[j, k]) <= np.abs(self.v[j, k-1]):
                        if np.abs(self.v[j, k]) <= np.abs(self.v[j, k+1]):
                            test = True
                            break
                    k += 1
                kpos = k-1
# Find the sample before the line cut position which has a (relative)
# minimum amplitude
                test = False
                k = ncut-1
                while test is False:
                    if np.abs(self.v[j, k]) <= np.abs(self.v[j, k-1]):
                        if np.abs(self.v[j, k]) <= np.abs(self.v[j, k+1]):
                            test = True
                            break
                    k -= 1
                kneg = k+1
# Out of the two minimum positions chose the one nearer to the cutting point
                if ncut-kneg > kpos-ncut:
                    ncut = kpos
                else:
                    ncut = kneg
# Set all values before the found position to zero
                if sign < 0:
                    self.v[j, :ncut] = 0
                else:
                    self.v[j, ncut:] = 0
                trace = self.actual_traces[j]
                nf = self.traces.file[trace]
                nt = self.traces.trace[trace]
                self.data.st[nf][nt].data = self.v[j, :]
# Plot the treated record section
        self.main.function = "main"
        self.drawNew(True)
        self.setHelp(self.main_text)

    def traceSign(self, QselfWindow):
        """
        Change polarity of one trace which is chosen by a mouse click
        (the value is selftained for the whole run until called again,
          also when changing data set)
        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "trace_sign"
        self.finish = False
        self.picked = False
        self.setHelp(self.trace_sign_text)

        def onPress(event):
            if event.button == 1:
                ret = self.searchTrace(event.xdata)
                if not ret:
                    return
                self.traces.amplitudes[self.itrace] *= -1
                self.v[self.n_tr, :] *= -1
                print("trace ", self.n_tr, " changed sign")
# Rewrite help text
                self.setHelp(self.trace_sign_text)
                self.picked = True
                self.axes[self.fig_plotted].plot(
                    self.x[self.n_tr], self.data.time[self.nt_mn+5], "*r")
                self.canvas.draw()
            elif event.button == 3:
                print("finish trace_sign")
                self.lin.set_animated(False)
                self.finish = True
                self.lin.figure.canvas.mpl_disconnect(self.cidpress)
                self.drawNew(True)
# Change help text to self module
                self.setHelp(self.main_text)
                self.main.function = "main"
        self.x_coor = []
        self.y_coor = []
        self.lin, = self.axes[self.fig_plotted].plot(self.x_coor, self.y_coor,
                                                     animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect(
            'button_press_event', onPress)

    def changeSign(self):
        """
        Change polarity of all traces
        (the sign is selftained for the whole run until called again,
          also when changing data set)
        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "change_sign"
        self.traces.amplitudes *= -1
        self.data.general_sign *= -1
        self.v *= -1
        self.drawNew(True)
# Change help text to self module
        self.setHelp(self.main_text)
        self.main.function = "main"

    def animateLine(self):
        """
        Plots animation of wave evolution in time along the geophone line.

        Amplitudes are clipped at 99% quantile and smoothed with 3-point
        moving average. To change clipping, change value of clip_level.

        Animation goes from first to last sample of actual time zoom.

        Input parameters:
            None

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
# chose data to be plotted and set clipping to 99% of extreme amplitudes
        nd_start = self.zooms[self.i_zooms][2]
        nd_end = self.zooms[self.i_zooms][3]
        clip_level = 0.99
        vt = self.v[:, nd_start:nd_end].transpose()
        for i in range(vt.shape[0]):
            vt[i, 1:-1] = (vt[i, :-2]+vt[i, 1:-1]+vt[i, 2:])/3
        t = self.data.t0+self.data.dt*np.arange(nd_start, nd_end)
        vmax = np.quantile(vt, clip_level)
        vt[np.abs(vt) > vmax] = vmax*np.sign(vt[np.abs(vt) > vmax])
# Define new figure
#        self.drawNew(False)
#        fig = self.figs[self.fig_plotted]
#        ax = self.axes[self.fig_plotted]
        self.w_anim = newWindow("Animation")
        fig = self.w_anim.fig
        ax = fig.subplots()
        ax.set_xlim(self.x.min(), self.x.max())
        ax.set_ylim(vt.min(), vt.max())
        ax.set_xlabel("Offset [m]", fontsize=18)
        ax.set_ylabel("Amplitude [n.u.]", fontsize=18)
        ax.tick_params(axis='both', labelsize=18)
        self.w_anim.show()
        self.setHelp(self.animate_text)
# Set matplotlib animated mode
        plt.ion()
        line = ax.plot(self.x, vt[0, :])[0]
        fig.canvas.draw()
        fig.canvas.flush_events()
        for i, tt in enumerate(t):
            line.set_ydata(vt[i, :])
            ax.set_title(f"Time {tt*1000:0.2f}ms", fontsize=20)
            fig.canvas.draw()
            fig.canvas.flush_events()
        plt.ioff()
