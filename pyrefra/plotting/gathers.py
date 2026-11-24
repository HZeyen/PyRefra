# -*- coding: utf-8 -*-
"""
Construction and plotting of shot, file, receiver and distance gathers.

Contains the following Class:
    GatherMixin
"""

import numpy as np
from PyQt5 import QtWidgets
from matplotlib.figure import Figure
from .view_state import GatherType


class GatherMixin:
    """Construction and plotting of shot, file, receiver and distance gathers."""

    def changeFig(self, item):
        """
        After having clicked on a plottable data set on the right side of the
        screen this function is called to change the plot.

        Parameters
        ----------
        item : int (?)
            index of the data set to be added, automatically assigned in line
            "self.mplfigs.itemClicked.connect"
        """
        self.v_set = False
        text = item.text()
# keep name of plotted data set
        self.fig_plotted = self.plot_names.index(text)
# Draw only a frame
        self.drawNew(False)
# Plot record section, depending on which type of plot is wanted
# (shot section, file section, record section, distance section)
# Keep zoom in time, but plot all traces
        if self.sg_flag:
            key = self.sh_plot[self.fig_plotted]
            self.zooms[self.i_zooms][:2] = self.gatherXLimits(key)
            self.fig_plotted_S = self.fig_plotted
            self.plotShot(self.axes[self.fig_plotted], self.fig_plotted)
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
# Keep zoom in time, but plot all traces
        elif self.fg_flag:
            key = self.fig_plotted
            self.zooms[self.i_zooms][:2] = self.gatherXLimits(key)
            self.fig_plotted_F = self.fig_plotted
            self.plotFile(self.axes[self.fig_plotted], self.fig_plotted)
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
# Keep zoom in time, but plot all traces. Offset of shots with respect to
# receivers has reverse sign than receiver vs. offset, therefore the negative
# sign and inversion of min and max
        elif self.rg_flag:
            key = int(text.split()[1])-1
            self.zooms[self.i_zooms][:2] = self.gatherXLimits(key)
            self.fig_plotted_R = self.fig_plotted
            nrec = np.where(self.receivers == key)[0][0]
            self.plotReceiver(self.axes[self.fig_plotted], self.receivers[nrec])
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
# Keep zoom in time, but plot all traces. Distances are plotted with respect to
# common depth point coordinate
        else:
            key = None
            self.fig_plotted_D = self.fig_plotted
            self.zooms[self.i_zooms][:2] = self.gatherXLimits(key)
            self.plotDistance(self.axes[self.fig_plotted],
                              self.distances[self.fig_plotted])
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
        if self.i_zooms == 0:
            self.i_zooms = 1
            try:
                self.zooms[1] = self.zooms[0]
            except IndexError:
                self.zooms.append(self.zooms[0])
        self.zooms[self.i_zooms][:2] = self.gatherXLimits(key)
# Plot picks if available
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
# plot help line for self module
        self.setHelp(self.main_text)
        self.v_set = True
        self.main.function = "main"
        if self.main.utilities.filtered_all:
            self.main.utilities.filtered = True
        else:
            self.main.utilities.filtered = False
        try:
            self.main.utilities.mgr
        except NameError:
            pass
        else:
            self.Change_colors.setEnabled(True)

    def gatherXLimits(self, key=None):
        """
        x-axis limits that include every trace of one gather in the actual
        ordering (offset or station number).

        Parameters
        ----------
        key : int, optional
            Shot key (shot gather), file index (file gather) or receiver key
            (receiver gather). Not used for distance gathers.

        Returns
        -------
        tuple of two floats
        """
        traces = self.traces
        if self.dg_flag:
            return traces.xcdp.min(), traces.xcdp.max()
        if self.order_number:
            if self.sg_flag:
                numbers = traces.sht_pt_dict[key]["receiver"]
            elif self.fg_flag:
                numbers = traces.receiver[self.files.file_dict[key]["traces"]]
            else:
                numbers = traces.rec_pt_dict[key]["shot"]
            return min(numbers), max(numbers)+1
# Offset of shots with respect to receivers has reverse sign than receiver
# vs. offset, therefore the negative sign and inversion of min and max
        if self.rg_flag:
            return -traces.off_max, -traces.off_min
        return traces.off_min, traces.off_max

    def toggleNumber(self):
        self.order_number = not self.order_number
        self.NumberOrder.toggle()
        if not self.order_number:
            self.x_text = "Offset [m]"
        elif self.rg_flag:
            self.x_text = "Shot point number"
        else:
            self.x_text = "Receiver number"
        if self.sg_flag:
            key = self.actual_shot
        elif self.fg_flag:
            key = self.fig_plotted
        else:
            key = self.traces.receiver[self.actual_traces[0]]
        self.zooms[self.i_zooms][:2] = self.gatherXLimits(key)
        self.v_set = False
        self.drawNew(True)

    def chooseComponent(self):
        """
        If several components have been recorded, choose the one to be plotted
        or to plot all components

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        labels = list(self.geom.types)+["All"]
        n = len(labels)
        results, okButton = self.main.dialog([labels], ["r"], [n],
                                             "Coose geophone component")
        if okButton is False:
            print("No component chosen")
            return
        self.plotComponent = labels[int(results[0])]
        print(f"Plot {self.plotComponent}-component")
        self.v_set = False
        self.drawNew(True)

    def plotRG(self):
        """
        Prepare all plots for the different existing receiver gathers

        Returns
        -------
        None.

        """
# Set flags correctly
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "plot_RG"
        self.mplfigs.clear()
        self.v_set = False
        self.gather_type = GatherType.RECEIVER
        self.figs = []
        self.axes = []
        self.plot_names = []
        self.keys_held = set()
# Find all unique receiver numbers. For every single one, a plot gather is
# prepared
        self.receivers = np.sort(np.array(list(self.traces.rec_pt_dict.keys()),
                                          dtype=int))
# Loop over all found receivers
        for i, r in enumerate(self.receivers):
            if r in self.data.receiver_corr_dict:
                if self.data.receiver_corr_dict[r]["amp"] == 0.:
                    continue
            self.figs.append(Figure())
            self.axes.append(self.figs[-1].subplots())
# plot data of the first receiver to the screen
            if i == 0:
                self.zooms[self.i_zooms][:2] = self.gatherXLimits(r)
                self.plotReceiver(self.axes[-1], r)
# Add the figure to the list of figures in the right-side tool bar
            self.addFig(f"Receiver {r+1}", self.figs[-1])
            self.plot_names.append(f"Receiver {r+1}")
# Erase actual plot in central widget and plot the new one
        self.rmMPL()
        self.addMPL(self.figs[0])
        self.fig_plotted = 0
        self.fig_plotted_R = 0
# Plot picks
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
        if self.PlotCalculatedTimes.isChecked():
            self.plotCalcPicks()
        self.main.function = "main"

    def plotReceiver(self, ax, irec):
        """
        Plot receiver gather.

        Parameters
        ----------
        ax : Matplotlib Axes object
            Axis where to plot the receiver gather
        irec : int
            Number of receiver to be plotted
        """
        self.time = self.data.t0+np.arange(int(self.data.nsamp))*self.data.dt
        self.x = []
        self.x_off = []
        self.x_ori = []
        self.tr = []
        self.x_plot = []
        self.stdev = []
# Search traces having been recorded by receiver irec
        ntraces = len(np.where(self.traces.receiver == irec)[0])
        self.traces.plotted[:] = False
        if self.v_set is False:
            self.v = np.zeros((ntraces, self.data.nsamp))
        self.v_norm = np.zeros((ntraces, self.data.nsamp))
        self.actual_traces = []
        self.actual_number_traces = -1
        if self.order_number:
            x_text = "Shot point number"
        else:
            x_text = "Shot offset [m]"
# Loop over all traces of the receiver gather
        j = -1
        for i, isht in enumerate(self.traces.rec_pt_dict[irec]["shot"]):
            j += 1
            if self.plotComponent not in (self.traces.component[irec], "All"):
                j -= 1
                continue
            ifile = self.traces.rec_pt_dict[irec]["file"][i]
            itrace = self.traces.rec_pt_dict[irec]["trace"][i]
            ntr = self.traces.sht_rec_dict[(isht, irec)]
            self.actual_traces.append(ntr)
            self.traces.plotted[ntr] = True
            if self.order_number:
                xx = float(isht+1)
                self.x.append(xx)
                self.x_ori.append(xx)
                self.x_off.append(-self.traces.offset[ntr])
            else:
                xx = -self.traces.offset[ntr]
                self.x.append(xx)
                self.x_off.append(xx)
                self.x_ori.append(xx)
# If a trace exists already at the actual position, shift the first one by
# 0.2m, the second one by +0.2m so they may be distinguished on the screen
            if j > 0:
                if self.x[-1] == self.x[-2]:
                    self.x[-1] += self.shift
                    self.x[-2] -= self.shift
            self.actual_number_traces += 1
# Copy data into array self.v
# v_set = True means that this work has already been done earlier
            if self.v_set is False:
                self.v[j, :] = self.data.st[ifile][itrace].data *\
                    self.traces.amplitudes[ntr]
                if self.nt_0 > 0:
                    self.v[j, :] -= \
                        np.mean(self.v[self.actual_number_traces, :self.nt_0])
# Normalize data by the standard deviation of each trace
            s = np.std(self.v[j, :])
            self.stdev.append(s)
            if s > 0:
                self.v_norm[j, :] = self.v[j, :]/s
            else:
                self.v_norm[j, :] = 0.
# If trace position is within the zoom limits, include it in list of traces to
#    be plotted (self.tr)
            if (self.x[-1] >= self.zooms[self.i_zooms][0]-self.shift and
                self.x[-1] <= self.zooms[self.i_zooms][1]+self.shift) or\
                    self.i_zooms == 0:
                self.tr.append(i)
                self.x_plot.append(self.x[-1])
        self.x = np.array(self.x)
        self.x_off = np.array(self.x_off)
        self.x_ori = np.array(self.x_ori)
        self.x_plot = np.array(self.x_plot)
        self.tr_plot = np.array(self.tr_plot)
        n_traces = len(self.x)
        self.v = self.v[:n_traces, :]
        self.v_norm = self.v_norm[:n_traces, :]
        self.actual_number_traces += 1
        self.indices = np.argsort(self.x)
        self.ind_plot = np.argsort(self.x_plot)
        self.tr_plot = np.array(self.tr)[self.ind_plot]
        self.x_plot = self.x_plot[self.ind_plot]
        self.stdev = np.array(self.stdev)
# Seismogram does the plotting
        if self.main.config:
            text_t = self.main.title+": "
        else:
            text_t = ""
        text_t += f"Receiver {irec+1}"
        if self.plotComponent != "All":
            text_t += f"; {self.plotComponent}-component"
        self.seismogram(ax, self.time, self.x, self.v, fill=True,
                        amp=self.amp_plt, traces=self.tr,
                        nt_min=self.nt_mn, nt_max=self.nt_mx,
                        text_x=x_text, text_t=text_t)
# Change_colors is only meant for plot of tomography result, so deactivate it
        self.Change_colors.setEnabled(False)
        self.main.function = "main"

    def plotDG(self):
        """
        Prepare all plots for the different existing offset gathers

        Returns
        -------
        None.

        """
# Set flags correctly
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "plot_DG"
        self.v_set = False
        self.mplfigs.clear()
        self.gather_type = GatherType.DISTANCE
        self.figs = []
        self.axes = []
        self.plot_names = []
        self.keys_held = set()

# Find all unique distances. For every single one, a plot gather is prepared
        self.distances = np.sort(np.unique(np.abs(self.traces.off_round)))
# Since typîcally, the preparation of distance gathers is rather time
# consuming, activate a progress bar
        self.progressBar = QtWidgets.QProgressBar(self)
        self.mplvl.addWidget(self.progressBar)
        self.progressBar.show()
        self.progressBar.setValue(0)
# Loop over all found offsets
# prepare a figure for each offset
        for i, d in enumerate(self.distances):
            self.figs.append(Figure())
            self.axes.append(self.figs[-1].subplots())
# plot the smallest found offset to the screen
            if i == 0:
                self.zooms[self.i_zooms][0] = np.min(self.traces.xcdp)
                self.zooms[self.i_zooms][1] = np.max(self.traces.xcdp)
                self.plotDistance(self.axes[-1], d)
# Add the figure to the list of figures in the right-side tool bar
            self.addFig(f"Distance {d}", self.figs[-1])
            self.plot_names.append(f"Distance {d}")
# Advance the progress bar
            completed = int((i+1)/len(self.distances)*100)
            self.progressBar.setValue(completed)
# Erase actual plot in central widget and plot the new one
        self.rmMPL()
        self.addMPL(self.figs[0])
# Remove the progress bar
        self.progressBar.setValue(0)
        self.mplvl.removeWidget(self.progressBar)
        self.progressBar.close()
        self.fig_plotted = 0
        self.fig_plotted_D = 0
# Plot picks
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
        if self.PlotCalculatedTimes.isChecked():
            self.plotCalcPicks()
        self.main.function = "main"

    def plotDistance(self, ax, dist):
        """
        Plot a distance gather

        Parameters
        ----------
        ax : matplotlib Axes
            Axes into which to plot the distance gather
        dist : float
            Offset to be plotted

        Returns
        -------
        None.

        """
        self.actual_dist = dist
# Search traces having an absolute offset equal to dist
        itraces = np.where(np.abs(self.traces.off_round) == dist)[0]
# CDP positions are used as X-coordinate
        indices = np.argsort(self.traces.ncdp[itraces])
        traces = itraces[indices]
        xx = self.traces.xcdp[traces]
        self.actual_traces = []
        self.traces.plotted[:] = False
        self.actual_number_traces = 0
        nrecplt = len(xx)
        self.x = []
        self.x_ori = []
        self.tr = []
        self.x_plot = []
        self.stdev = []
# Loop over all traces of the distance gather
        j = -1
        for i in range(nrecplt):
            j += 1
            if self.plotComponent not in (self.traces.component[traces[i]],
                                          "All"):
                j -= 1
                continue
            self.x.append(xx[i])
            self.x_ori.append(xx[i])
            self.actual_traces.append(traces[i])
            self.traces.plotted[traces[i]] = True
            self.actual_number_traces += 1
# If a trace exists already at the actual position (reciprocal positions of
#    shot and receiver), shift the first one by -0.2m, the second one by +0.2m
#    so they may be distinguished on the screen
            if j > 0:
                if self.x[-1] == self.x[-2]:
                    self.x[-1] += self.shift
                    self.x[-2] -= self.shift
# If trace position is within the zoom limits, include it in list of traces to
#    be plotted (self.tr)
            if (self.x[-1] >= self.zooms[self.i_zooms][0]-self.shift and
                self.x[-1] <= self.zooms[self.i_zooms][1]+self.shift) or\
                    self.i_zooms == 0:
                self.tr.append(j)
                self.x_plot.append(self.x[-1])
        self.x = np.array(self.x)
        self.x_off = np.copy(self.x)
        self.x_ori = np.array(self.x_ori)
        self.x_plot = np.array(self.x_plot)
        self.tr = np.array(self.tr, dtype=int)
        n_traces = len(self.x)
        self.v = self.v[:n_traces, :]
        self.v_norm = self.v_norm[:n_traces, :]
        ns = self.data.nsamp
# Copy data into array self.v
# v_set = True means that this work has already been done earlier
        if self.v_set is False:
            self.v = np.zeros((nrecplt, ns))
        self.v_norm = np.zeros((nrecplt, ns))
        for i in range(nrecplt):
            ifile = self.traces.file[traces[i]]
            itrace = self.traces.trace[traces[i]]
            if self.v_set is False:
                self.v[i, :] = self.data.st[ifile][itrace].data *\
                    self.traces.amplitudes[traces[i]]
                if self.nt_0 > 0:
                    self.v[i, :] -= np.mean(self.v[i, :self.nt_0])
# Normalize data by the standard deviation of each trace
            s = np.std(self.v[i, :])
            self.stdev.append(s)
            if s > 0:
                self.v_norm[i, :] = self.v[i, :]/s
            else:
                self.v_norm[i, :] = 0.
        self.indices = np.argsort(self.x)
        self.ind_plot = np.argsort(self.x_plot)
        self.tr_plot = np.copy(self.tr[self.ind_plot])
        self.x_plot = self.x_plot[self.ind_plot]
        self.stdev = np.array(self.stdev)
# Seismogram does the plotting
        if self.main.config:
            text_t = self.main.title+": "
        else:
            text_t = ""
        text_t += f"Distance {dist} m"
        if self.plotComponent != "All":
            text_t += f"; {self.plotComponent}-component"
        self.seismogram(ax, self.time, self.x, self.v, fill=True,
                        amp=self.amp_plt, traces=self.tr,
                        nt_min=self.nt_mn, nt_max=self.nt_mx,
                        text_x="Midpoint position [m]", text_t=text_t)
# Change_colors is only meant for plot of tomography result, so deactivate it
        self.Change_colors.setEnabled(False)
        self.main.function = "main"

    def plotSG(self, isht=0):
        """
        Prepare all plots for the different existing shot gathers

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
# Set flags correctly
        self.main.function = "shot_gather"
        self.v_set = False
        self.gather_type = GatherType.SHOT
        self.mplfigs.clear()
        self.figs = []
        self.axes = []
        self.plot_names = []
        self.sh_plot = []
        self.stdev = []
        self.keys_held = set()
# Loop over all shot points
        for i, sh in enumerate(self.traces.sht_pt_dict):
            self.sh_plot.append(sh)
            self.figs.append(Figure())
            self.axes.append(self.figs[-1].subplots())
            self.traces.sht_pt_dict[sh]["axes"].append(i)
# Set zoom limits in x such that all traces are plotted to the screen
            if i == 0:
                if self.order_number:
                    self.zooms[self.i_zooms][0] = min(self.traces.receiver)
                    self.zooms[self.i_zooms][1] = max(self.traces.receiver)+1
                else:
                    self.zooms[self.i_zooms][0] = self.traces.off_min
                    self.zooms[self.i_zooms][1] = self.traces.off_max
# Plot data of first shot point to the screen
                self.plotShot(self.axes[-1], i)
# Add the figure to the list of figures in the right-side tool bar
            self.addFig(f"Shot {sh+1}", self.figs[-1])
            self.plot_names.append(f"Shot {sh+1}")
# Erase actual plot in central widget and plot the new one
        self.rmMPL()
        self.addMPL(self.figs[isht])
        self.fig_plotted = 0
        self.fig_plotted_S = 0
# Plot picks
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
        if self.PlotCalculatedTimes.isChecked():
            self.plotCalcPicks()
        self.main.function = "main"

    def plotShot(self, ax, sh0):
        """
        Plot a shot point gather

        Parameters
        ----------
        ax : matplotlib Axes
            Axes into which to plot the shot gather
        sh0 : int
            Shot point number to be plotted

        Returns
        -------
        None.

        """
        sh = self.sh_plot[sh0]
        ns = self.data.nsamp
# Search traces that recorded shot point sh0
        ntraces = len(np.where(self.traces.shot == sh)[0])
        if ntraces == 0:
            return
        self.x = []
        self.x_off = []
        self.x_ori = []
        self.x_plot = []
        self.tr = []
        self.actual_traces = []
        self.traces.plotted[:] = False
        self.actual_number_traces = 0
        self.stdev = []
        if self.v_set is False:
            self.v = np.zeros((ntraces, ns))
        self.v_norm = np.zeros((ntraces, ns))
        self.actual_shot = sh
        j = -1
# Loop over all traces of the shot point gather
        if self.order_number:
            x_text = "Receiver number"
        else:
            x_text = "Offset [m]"
        for i, nt in enumerate(self.traces.sht_pt_dict[sh]["trace"]):
            j += 1
            ifile = self.traces.sht_pt_dict[sh]["file"][i]
            irec = self.traces.sht_pt_dict[sh]["receiver"][i]
            ntr = self.traces.sht_rec_dict[(sh, irec)]
            if self.plotComponent not in (self.traces.component[ntr], "All"):
                j -= 1
                continue
            self.actual_traces.append(ntr)
            self.traces.plotted[ntr] = True
            self.actual_number_traces += 1
            if self.order_number:
                # xx = float(nt)
                xx = float(irec+1)
                self.x.append(xx)
                self.x_ori.append(xx)
                self.x_off.append(self.traces.offset[ntr])
            else:
                xx = self.traces.offset[ntr]
                self.x.append(xx)
                self.x_ori.append(xx)
                self.x_off.append(xx)
# If a trace exists already at the actual position, shift the first one by
# -0.2m, the second one by +0.2m so they may be distinguished on the screen
            if j > 0:
                if self.x[-1] == self.x[-2]:
                    self.x[-1] += self.shift
                    self.x[-2] -= self.shift
# Copy data into array self.v
# v_set = True means that this work has already been done earlier
            if self.v_set is False:
                try:
                    self.v[j, :] = self.data.st[ifile][nt].data *\
                                  self.traces.amplitudes[ntr]
# If there is a problem filling self.v, it is usually due to variable trace
# lengths between shot points. In this case, since v.shape[1] is defined as the
# length of the first found trace, one of two possibilities exist: If the
# actual trace is too long, cut it at the end; if the actual trace is shorter,
# fill v with zeros at the end.
                except (KeyError, IndexError, NameError, ValueError):
                    nv = self.v.shape[1]
                    nd = len(self.data.st[ifile][nt].data)
                    if nv > nd:
                        self.v[j, :nd] = self.data.st[ifile][nt].data *\
                            self.traces.amplitudes[ntr]
                    else:
                        self.v[j, :] = self.data.st[ifile][nt].data[:nv] *\
                            self.traces.amplitudes[ntr]
                if self.nt_0 > 0:
                    self.v[j, :] -= np.mean(self.v[j, :self.nt_0])
# Normalize data by the standard deviation of each trace
            s = np.std(self.v[j, :])
            self.stdev.append(s)
            if s > 0:
                self.v_norm[j, :] = self.v[j, :]/s
            else:
                self.v_norm[j, :] = 0.
# If trace position is within the zoom limits, include it in list of traces to
#    be plotted (self.tr)
            if (self.x[-1] >= self.zooms[self.i_zooms][0]-self.shift and
                self.x[-1] <= self.zooms[self.i_zooms][1]+self.shift) or\
                    self.i_zooms == 0:
                self.tr.append(j)
                self.x_plot.append(self.x[-1])
        self.x = np.array(self.x)
        self.x_off = np.array(self.x_off)
        self.x_ori = np.array(self.x_ori)
        self.x_plot = np.array(self.x_plot)
        n_traces = len(self.x)
        self.v = self.v[:n_traces, :]
        self.v_norm = self.v_norm[:n_traces, :]
        self.indices = np.argsort(self.x)
        self.ind_plot = np.argsort(self.x_plot)
        self.tr_plot = np.copy(np.array(self.tr))
        self.tr_plot = self.tr_plot[self.ind_plot]
        self.x_plot = self.x_plot[self.ind_plot]
        self.stdev = np.array(self.stdev)
# Seismogram does the plotting
        if self.main.config:
            text_t = self.main.title+": "
        else:
            text_t = ""
        text_t += f"Shot point {sh+1}"
        if self.phase_plot:
            xx, vv = self.phaseCalc()
            text_t += ", Inclination/declination"
            tt = np.arange(len(xx))
            self.seismogram(ax, self.data.time, xx, vv, fill=True,
                            amp=self.amp_plt, traces=tt, nt_min=self.nt_mn,
                            nt_max=self.nt_mx, text_x=x_text,
                            text_t=text_t)
        else:
            if self.plotComponent != "All":
                text_t += f"; {self.plotComponent}-component"
            self.seismogram(ax, self.data.time, self.x, self.v, fill=True,
                            amp=self.amp_plt, traces=self.tr,
                            nt_min=self.nt_mn, nt_max=self.nt_mx,
                            text_x=x_text, text_t=text_t)
# The following lines are only for testing purpose, to show potential picks
# from false colour plots
        # try:
        #     self.axes[self.fig_plotted].plot(self.x,self.main.utilities.pk,
        #                                      "m")
        #     self.axes[self.fig_plotted].set_ylim([self.time_plt_min,
        #                                           self.time_plt_max])
        #     renderer = self.figs[self.fig_plotted].canvas.renderer
        #     self.axes[self.fig_plotted].draw(renderer)
        #     self.canvas.blit(self.axes[self.fig_plotted].bbox)
        # except:
        #     pass
# Change_colors is only meant for plot of tomography result, so deactivate it
        self.Change_colors.setEnabled(False)
        self.main.function = "main"

    def plotFG(self, isht=0):
        """
        Prepare all plots for the different files

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
# Set flags correctly
        self.main.function = "shot_gather"
        self.v_set = False
        self.gather_type = GatherType.FILE
        self.mplfigs.clear()
        self.figs = []
        self.axes = []
        self.plot_names = []
        self.stdev = []
        self.keys_held = set()
# Loop over all recorded files
        for i in range(self.files.file_count):
            self.figs.append(Figure())
            self.axes.append(self.figs[-1].subplots())
# Set zoom limits in x such that all traces are plotted to the screen
            if i == 0:
                if self.order_number:
                    self.zooms[self.i_zooms][0] = min(self.traces.receiver)
                    self.zooms[self.i_zooms][1] = max(self.traces.receiver)+1
                else:
                    self.zooms[self.i_zooms][0] = self.zooms[0][0]
                    self.zooms[self.i_zooms][1] = self.zooms[0][1]
# Plot data of first file to the screen
                self.plotFile(self.axes[-1], i)
# Add the figure to the list of figures in the right-side tool bar
            self.addFig(f"File {self.files.numbers[i]}", self.figs[-1])
            self.plot_names.append(f"File {self.files.numbers[i]}")
# Erase actual plot in central widget and plot the new one
        self.rmMPL()
        self.addMPL(self.figs[isht])
        self.fig_plotted = 0
        self.fig_plotted_F = 0
# Plot picks
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
        if self.PlotCalculatedTimes.isChecked():
            self.plotCalcPicks()
        self.main.function = "main"

    def plotFile(self, ax, isht):
        """
        Plot data of file isht into axis ax
        """
        self.x = []
        self.x_off = []
        self.x_ori = []
        self.tr = []
        self.x_plot = []
        self.actual_traces = []
        self.traces.plotted[:] = False
        self.actual_number_traces = 0
        self.stdev = []
# Set a few plotting parameters
        ntr = len(self.data.st[isht])
        nsamp = len(self.data.st[isht][0].data)
        sht = int(self.data.st[isht][0].stats.seg2['SOURCE_STATION_NUMBER'])

        if self.v_set is False:
            self.v = np.zeros((ntr, nsamp))
        self.v_norm = np.zeros((ntr, nsamp))
        self.actual_shot = self.files.numbers[isht]
        if self.order_number:
            x_text = "Receiver number"
        else:
            x_text = "Offset [m]"
        j = -1
# Loop over all traces of the file
        for i, t in enumerate(self.files.file_dict[isht]["traces"]):
            j += 1
            if self.plotComponent not in (self.traces.component[t], "All"):
                j -= 1
                continue
            self.actual_traces.append(t)
            self.traces.plotted[t] = True
            self.actual_number_traces += 1
            if self.order_number:
                xx = float(self.traces.receiver[t]+1)
                self.x.append(xx)
                self.x_ori.append(xx)
                self.x_off.append(self.traces.offset[t])
            else:
                xx = self.traces.offset[t]
                self.x.append(xx)
                self.x_off.append(xx)
                self.x_ori.append(xx)
# If a trace exists already at the actual position, shift the first one by
# -0.2m, the second one by +0.2m so they may be distinguished on the screen
            if j > 0:
                if self.x[-1] == self.x[-2]:
                    self.x[-1] += self.shift
                    self.x[-2] -= self.shift
# Copy data into array self.v
# v_set = True means that this work has already been done earlier
            if self.v_set is False:
                self.v[j, :] = self.data.st[isht][i].data *\
                    self.traces.amplitudes[t]
                if self.nt_0 > 0:
                    self.v[j, :] -= np.mean(self.v[j, :self.nt_0])
# Normalize data by the standard deviation of each trace
            s = np.std(self.v[j, :])
            self.stdev.append(s)
            if s > 0:
                self.v_norm[j, :] = self.v[j, :]/s
            else:
                self.v_norm[j, :] = 0.
# If trace position is within the zoom limits, include it in list of traces to
#    be plotted (self.tr)
            if (self.x[-1] >= self.zooms[self.i_zooms][0]-self.shift and
                self.x[-1] <= self.zooms[self.i_zooms][1]+self.shift) or\
                    self.i_zooms == 0:
                self.tr.append(j)
                self.x_plot.append(self.x[-1])
        self.x = np.array(self.x)
        self.x_off = np.array(self.x_off)
        self.x_ori = np.array(self.x_ori)
        self.x_plot = np.array(self.x_plot)
        self.tr = np.array(self.tr, dtype=int)
        n_traces = len(self.x)
        self.v = self.v[:n_traces, :]
        self.v_norm = self.v_norm[:n_traces, :]
        self.indices = np.argsort(self.x)
        self.ind_plot = np.argsort(self.x_plot)
        self.tr_plot = np.copy(self.tr[self.ind_plot])
        self.x_plot = self.x_plot[self.ind_plot]
        self.stdev = np.array(self.stdev)
# Seismogram does the plotting
        if self.main.config:
            text_t = self.main.title+": "
        else:
            text_t = ""
        text_t += f"File {self.files.numbers[isht]}, shot point {sht}"
        if self.plotComponent != "All":
            text_t += f"; {self.plotComponent}-component"
        self.seismogram(ax, self.data.time, self.x, self.v, fill=True,
                        amp=self.amp_plt, traces=self.tr,
                        nt_min=self.nt_mn, nt_max=self.nt_mx,
                        text_x=x_text, text_t=text_t)
# Change_colors is only meant for plot of tomography result, so deactivate it
        self.Change_colors.setEnabled(False)
        self.main.function = "main"
