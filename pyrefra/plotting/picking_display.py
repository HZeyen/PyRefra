# -*- coding: utf-8 -*-
"""
Display and lookup of picks in the data window.

Contains the following Class:
    PickDisplayMixin
"""

import sys
import numpy as np
from PyQt5 import QtWidgets
from .canvas_manager import newWindow


class PickDisplayMixin:
    """Display and lookup of picks in the data window."""

    def picksPlot(self, fig, ax, col="r"):
        """
        Plot all measured picks into the actual gather
        Uses function pickPlot to plot each pick found

        Parameters:
        ----------
        ax:  axes where to plot picks
        col: str
            color to be used for pcik plot. Default : red ("r")

        Returns
        -------
        None.

        """
# Check which trace of the actual gather is the first one plotted on the
# screen. Traces actually plotted on the screen ares stored in self.tr, the
# ones of the actual gather are stored in self.actual_traces.
        if np.sum(np.array(
                self.traces.npick, dtype=int)[self.actual_traces]) > 0:
            fil = []
            for itr in self.tr:
                ntr = self.actual_traces[itr]
                fil.append(self.traces.file[ntr])
            fil = np.array(fil, dtype=int)
            f_nr = np.unique(fil)
            cols = ["red", "lime", "blue", "cyan", "yellow", "magenta", "orange"]
            nchar = len(cols[0])
            for c in cols[1:]:
                nchar = max(nchar, len(c))
            form = f"<U{nchar}"
            c = np.zeros_like(fil, dtype=form)
            if self.sg_flag:
                for i, f in enumerate(f_nr):
                    c[fil == f] = cols[i]
            else:
                c[:] = cols[0]
            for i, itr in enumerate(self.tr):
                ntr = self.actual_traces[itr]
                if self.traces.npick[ntr] > 0:
                    for j in range(self.traces.npick[ntr]):
                        self.pickPlot(itr, ntr, j, c[i])
            renderer = fig.canvas.renderer
            ax.draw(renderer)
            self.PlotPicks.setEnabled(True)
            self.MovePicks.setEnabled(True)
            self.Tomography.setEnabled(True)
            self.Pseudo_velocity.setEnabled(True)
        else:
            pass

    def pickPlot(self, i, ntr, j, col, v=True, h=True):
        """
        Plot one single pick

        Parameters
        ----------
        i : int
            number of trace in gather
        ntr : int
            number of trace in list of all traces
        j : int
            number of pick for that trace
        col : string
            Color to be used for plot
        v : bool
            If True, vertical uncertainty bar is plotted
        h : bool
            If True, horizontal bar is plotted at the time of the pick

        Returns
        -------
        None.

        """
        t = self.traces.pick_times[ntr][j]
# Picks outside the actual time zoom are not plotted
        if self.time_plt_min <= t <= self.time_plt_max:
            if h:
                x_c = []
                x_c.append(self.x[i]-0.3)
                x_c.append(self.x[i]+0.3)
                y_c = [t, t]
                self.line, = self.axes[self.fig_plotted].plot(x_c, y_c, col)
# If plot pick is called before the uncertainty of a travel time has been
# picked, i.e. before the picking has been finished, add standard uncertainty
            if v:
                try:
                    tmn = self.traces.pick_times_min[ntr][j]
                    tmx = self.traces.pick_times_max[ntr][j]
                except (KeyError, IndexError, NameError):
                    try:
                        self.traces.pick_times_min[ntr][j] = t-2*self.data.dt
                        self.traces.pick_times_max[ntr][j] = t+2*self.data.dt
                        tmn = self.traces.pick_times_min[ntr][j]
                        tmn = self.traces.pick_times_max[ntr][j]
                    except (KeyError, IndexError, NameError):
                        return
                x_c = [self.x[i], self.x[i]]
                y_c = [max(tmn, self.time_plt_min), min(tmx, self.time_plt_max)]
                self.line, = self.axes[self.fig_plotted].plot(x_c, y_c, col)

    def plotPickSection(self):
        """
        Opens new window and plots all picks vs receiver poition. Picks of
        different shots get different colors (7 base colors, except white).

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        cols = ["blue", "lime", "red", "cyan", "magenta", "yellow", "orange",
                "black", "gray", "lightgray", "brown", "skyblue", "gold",
                "teal", "limegreen", "pink"]
        marks = ["+", "x", ".", "d", "v", "^"]
        labels = ["Plot picks by",
                  ["Trace position","Offset"],
                  "Define colors by",
                  ["shot nr", "ray direction"],
                  "Connect points with line"]
        types = ["l", "r", "l", "r", "c"]
        values = [None, 1, None, 1, None]
        results, okButton = self.main.dialog(
            labels, types, values, "Plot picks")
        if okButton is False:
            print("Plot all picks method cancelled")
            return
        ordering = int(results[1])
        coloring = int(results[3])
        connect = int(results[4]) > -1
        dangle = 20.
        self.w_picks = newWindow("All picks")
        self.ax = self.w_picks.fig.subplots()
        npk = np.array(self.traces.npick)
        if coloring:
            shts = self.traces.shot
            recs = self.traces.receiver
            values = np.zeros_like(shts, dtype=int)
            colors = np.arange(int(180./dangle))
            for i, s in enumerate(shts):
                dx = self.geom.sht_dict[s]["x"]-self.geom.rec_dict[recs[i]]["x"]
                dy = self.geom.sht_dict[s]["y"]-self.geom.rec_dict[recs[i]]["y"]
                values[i] = np.arctan(dy/(dx+sys.float_info.min))*180./np.pi
                if values[i] < 0.:
                    values[i] += 180.
                values[i] = int(values[i]/dangle)
                # values[i] = min(values[i], 17)
            values = values[npk > 0]
        else:
            values = self.traces.shot[npk > 0]
            colors = np.unique(values)
        if ordering:
            x = np.abs(np.array(self.traces.offset)[npk > 0])
            xtext = "Offset [m]"
            name = "offset"
        else:
            x = np.array(self.traces.receiver_pos)[npk > 0]
            xtext = "Trace position [m]"
            name = "receiver"
        t = np.array([tt[0] for tt in
                      np.array(self.traces.pick_times[:],
                               dtype=object)[npk > 0]])
        icol = -1
        imark = 0
        nleg_lines = 0
        for i, s in enumerate(colors):
            xx = x[values == s]
            tt = t[values == s]
            if len(xx) == 0:
                continue
            icol += 1
            nleg_lines += 1
            if icol == len(cols):
                icol = 0
                imark += 1
            mark = marks[imark]
            col = cols[icol]
            if coloring:
                label = f"{i*dangle}-{(i+1)*dangle}°"
            else:
                label = f"Shot {s+1}"
            if connect:
                self.ax.plot(xx, tt, color=col, marker=mark, label=label)
            else:
                self.ax.scatter(xx, tt, marker=mark, c=col, label=label)
        self.ax.set_xlabel(xtext, fontsize=18)
        self.ax.set_ylabel("Time [s]", fontsize=18)
        npks = np.sum(npk)
        self.ax.set_title(f"Picked arrival times ({npks} picks)", fontsize=20)
        self.ax.tick_params(axis='both', labelsize=18)
        ax_xmin, ax_xmax = self.ax.get_xlim()
        ax_ymin, ax_ymax = self.ax.get_ylim()
        self.ax.text(ax_xmin, ax_ymax+(ax_ymax-ax_ymin)*0.01,
                     self.main.dir_start,
                     horizontalalignment="left",
                     verticalalignment="bottom", fontsize=18)
        self.ax.text(ax_xmax, ax_ymax+(ax_ymax-ax_ymin)*0.01,
                     self.main.dir_end,
                     horizontalalignment="right",
                     verticalalignment="bottom", fontsize=18)
        self.ax.grid(axis="both", which="major", linestyle="-")
        self.ax.grid(axis="y", which="minor", linestyle="--")
        self.ax.minorticks_on()
        ncolumns = int((nleg_lines-1)/20)+1
        self.ax.legend(ncols=ncolumns)
        self.w_picks.show()
        self.w_picks.fig.savefig(f"all_picks_{name}.png")

    def plotCalcPicks(self):
        """
        Plot calculated travel times, see function Traces.readCalcPicks for
        further explanations

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        if not self.traces.calc_picks:
            print("\nNo calculated travel times available\n")
            return
        xc = []
        tc = []
        if self.dg_flag:
            for t in self.actual_traces:
                if np.isnan(self.traces.calc_t[t]):
                    continue
                xc.append(self.traces.xcdp[t])
                tc.append(self.traces.calc_t[t])
        else:
            for t in self.actual_traces:
                if np.isnan(self.traces.calc_t[t]):
                    continue
                if self.rg_flag:
                    xc.append(-self.traces.offset[t])
                else:
                    xc.append(self.traces.offset[t])
                tc.append(self.traces.calc_t[t])
# Next lines are important for file gather to avoid picks plotted outside of
# the plotted traces
        xc = np.array(xc)
        tc = np.array(tc)
        tc = tc[(xc <= max(self.x)) & (xc >= min(self.x))]
        xc = xc[(xc <= max(self.x)) & (xc >= min(self.x))]
        nsort = np.argsort(xc)
        xx = xc[nsort]
        yy = tc[nsort]
# The following try-except construction avoids problems in case no calculated
# travel times are found for a gather. This is mainly the case for distance
# gathers at zero distance, since there, no times are calculated.
        k = len(xx)
        try:
            if k < 1:
                print("\nNo calculated travel times available")
            else:
                xmn = self.zooms[self.i_zooms][0]
                xmx = self.zooms[self.i_zooms][1]
                xxx = xx[xx >= xmn]
                yyy = yy[xx >= xmn]
                xx = xxx[xxx <= xmx]
                yy = yyy[xxx <= xmx]
                self.line, = self.axes[self.fig_plotted].plot(xx, yy, "g")
                self.axes[self.fig_plotted].set_ylim([self.time_plt_min,
                                                      self.time_plt_max])
                renderer = self.figs[self.fig_plotted].canvas.renderer
                self.axes[self.fig_plotted].draw(renderer)
                self.canvas.blit(self.axes[self.fig_plotted].bbox)
        except Exception as error:
            if k < 1:
                print("\nNo calculated travel times available")
            else:
                print(f"\nProblems plotting {k} calculated travel times: "
                      + f"{error}.")

    def searchTrace(self, x):
        """
        Search the trace nearest to a mouse click
        Input:
                x (float): position of mouse click
        Output:
                self.p_s: shot number of selected trace
                self.p_r: receiver number of selected trace
                self.n_tr: number of trace on actual screen
        """
        if x is None:
            return False
        dmin = 1e10
        self.n_tr = 0
        for i in range(self.actual_number_traces):
            dist = np.abs(self.x[i]-x)
            if dist < dmin:
                dmin = dist
                self.n_tr = i
        self.itrace = self.actual_traces[self.n_tr]
        self.p_s = self.traces.shot[self.itrace]
        self.p_r = self.traces.receiver[self.itrace]
        return True

    def searchPick(self, x, t=0):
        """
        Search the pick nearest to a mouse click
        Input:
                x (float): x-position of mouse click
                t (float): time position of mouse click
        """
# Search trace containing a pick nearest to x position of mouse click
# Save actual trace number in case click is on a trace without pick
        try:
            itr = self.itrace
        except NameError:
            itr = -1
            for i, it in enumerate(self.actual_traces):
                if self.traces.npick[it] > 0:
                    itr = it
                    break
            if itr < 0:
                _ = QtWidgets.QMessageBox.warning(
                    None, "Warning", "No pick exists for this gather.\n",
                    QtWidgets.QMessageBox.Close)
                return False
        ret = self.searchTrace(x)
        if not ret:
            return
# If there are several picks for the trace, search pick nearest to time of
# mouse position
        if self.traces.npick[self.itrace] > 0:
            dmin = np.abs(self.traces.pick_times[self.itrace][0]-t)
            self.ipk = 0
            if self.traces.npick[self.itrace] > 1:
                for i in range(1, self.traces.npick[self.itrace]):
                    dist = np.abs(self.traces.pick_times[self.itrace][i]-t)
                    if dist < dmin:
                        dmin = dist
                        self.ipk = i
# If clicked trace has no pick, stay on last trace
        else:
            self.itrace = itr
            self.ipk = 0
        return True

    def searchNearPicks(self):
        """
        At the beginning of manual picking, the function is called to search
        for the gather plotted on the screen existing picks done at nearby
        gathers with the same offsets for each trace. They are used as hint
        for the position of picks on the actual gather.

        Returns
        -------
        x_pk_near : numpy float
            1D array with maximum length = number of traces on screen. Contains
            the signed offset of corresponding trace
        t_pk_near : numpy float
            1D array with maximum length = number of traces on screen. Contains
            the picked traveltime measured for the offset abs(x_pk_near) in
            nearby gather
        i_tr_near : numpy int
            1D array with maximum length = number of traces on screen. Contains
            the picked traveltime measured for the offset abs(x_pk_near) in
            nearby gather

        If in other gathers no pick has been measured for the offset
        corresponding to a trace on the screen, the arrays do not contain any
        value for this trace.
        """
        x_pk_near = []
        t_pk_near = []
        i_tr_near = []
        for i in range(self.actual_number_traces):
            n_tr = self.indices[i]
            trace = self.actual_traces[n_tr]
            p_r = self.traces.receiver[trace]
            if self.dg_flag:
                x_pk = abs(self.traces.offset[trace])
            else:
                x_pk = self.x[n_tr]
            similar = np.where(np.abs(self.traces.offset-abs(x_pk)) < 0.5)
            sim_t = []
            sim_s = []
            sim_r = []
            sim_cdp = []
            d_trac = []
            for s in similar[0]:
                if self.traces.npick[s] < 1:
                    continue
                if s == trace:
                    continue
                sim_t.append(s)
                sim_s.append(self.traces.shot[s])
                sim_r.append(self.traces.receiver[s])
                sim_cdp.append(self.traces.xcdp[s])
                d_trac.append(abs(trace-s))
            if len(sim_t) < 1:
                continue
            sim_t = np.array(sim_t, dtype=int)
            sim_s = np.array(sim_s, dtype=int)
            sim_r = np.array(sim_r, dtype=int)
            sim_cdp = np.array(sim_cdp, dtype=float)
            d_trac = np.array(d_trac, dtype=float)
            if self.sg_flag or self.fg_flag:
                itr = np.argsort(abs(sim_cdp-self.traces.xcdp[trace]))
                sim_t = sim_t[itr]
                sim_cdp = sim_cdp[itr]
                if len(sim_t) > 1:
                    for k, cdp in enumerate(sim_cdp[1:]):
                        if cdp != sim_cdp[0]:
                            break
                        if d_trac[k] < d_trac[0]:
                            d_trac[0] = d_trac[k]
                            sim_t[0] = sim_t[k]
                t_pk = None
                for i in range(min(len(itr), 6)):
                    try:
                        t_pk = self.traces.pick_times[sim_t[0]][0]
                        d = abs(self.traces.xcdp[sim_t[0]] -
                                self.traces.xcdp[trace])
                        break
                    except (KeyError, IndexError, NameError):
                        continue
                if not t_pk or d > 10:
                    continue
                t_pk_near.append(t_pk)
                x_pk_near.append(self.x[n_tr])
                i_tr_near.append(n_tr)
            elif self.rg_flag:
                itr = np.argsort(np.abs(sim_r-p_r))
                t_pk = None
                for i in range(min(len(itr), 6)):
                    try:
                        t_pk = self.traces.pick_times[sim_t[i]][0]
                        break
                    except (KeyError, IndexError, NameError):
                        continue
                if not t_pk:
                    continue
                t_pk_near.append(t_pk)
                x_pk_near.append(self.x[n_tr])
                i_tr_near.append(n_tr)
        return np.array(x_pk_near), np.array(t_pk_near), \
            np.array(i_tr_near, dtype=int)
