# -*- coding: utf-8 -*-
"""
Rendering helpers for the seismogram window that do not depend on window state.

Contains the following functions:
    agc_calc
    set_ticks
"""

import numpy as np
from PyQt5 import QtWidgets


def agc_calc(data, n_agc_win):
    """
    Calculate AGC gain for one trace

    Parameters
    ----------
    data : Numpy array one dimension; float
        Data to be transformed
    n_agc_win : int
        length of AGC window in number of samples

    Returns
    -------
    data : Numpy array with same shape and type as data
        AGC gained data

    """
    x = np.abs(data)
    nsamp = len(data)
    for j in range(nsamp):
        i1 = max(j-n_agc_win, 0)
        i2 = min(j+n_agc_win, nsamp)
        m = np.max(x[i1:i2])
        if m > 0:
            data[j] = data[j]/m
    return data


def set_ticks(minval, maxval, ntick=5, dtick=5.):
    """
    Function

    Parameters
    ---------s-
    minval : float
        Minimum value of data to be represented (first tick will be
                at a position larger or equal to this value)
    maxval : float
        Maximum value of data to be represented (last tick will be
                at a position smaller or equal to this value)
    ntick : int, optional
        Approximate number of ticks. The default is 5.
    dtick : float, optional
        Ticks will be placed at multiples of dtick such that the number of
        ticks is approximately nticks. The default is 5.

    Returns
    -------
    ticks : numpy array
        contains the positions of the ticks

    """
    mxv = max(minval, maxval)
    mnv = min(minval, maxval)
    d_ticks = (mxv-mnv)/ntick
    mult = round(d_ticks/dtick)
    mult = max(mult, 1)
    d_ticks = dtick*mult
    min_tick = np.ceil(mnv/d_ticks)*d_ticks
    max_tick = int(mxv/d_ticks)*d_ticks
    ticks = np.arange(min_tick, max_tick+d_ticks/2, d_ticks)
    return ticks


class RenderMixin:
    """Seismogram drawing of the data window."""

    def agcCalc(self, data, n_agc_win):
        """
        Calculate AGC gain for one trace

        Parameters
        ----------
        data : Numpy array one dimension; float
            Data to be transformed
        n_agc_win : int
            length of AGC window in number of samples

        Returns
        -------
        data : Numpy array with same shape and type as data
            AGC gained data

        """
        return agc_calc(data, n_agc_win)

    def seismogram(self, ax, time, x_pos, data, nt_min=0, nt_max=0,
                   traces=None, text_x="Distance [m]", text_y="Time [s]",
                   text_t="", amp=0.667, fill=False):
        """
        plots a seismogram section

        Parameters
        ----------
        ax : Matplotlib Axes object
            axis where to plot the seismogram
        time : numpy 1D float array with length nsamples
            time of each sample
        x_pos : numpy 1D float array
            position of geophones
        data : numpy 1D float array with length nsamples
            measured amplitudes
        nt_min : int, optional; default: 0
            number of first sample to be plotted (counting starts at 0!)
        nt_max : int, optional; default: 0
            number of last sample to be plotted. If 0: plot all samples
        traces list of int, optional; default: None
            list of traces to be plotted. If None: plot all traces
        text_x : str, optional; default: "Distance [m]"
            Text to be written on x-axis
        text_y : str, optional; default: "Time [s]"
            Text to be written on y-axis
        amp : float, optional; default: 0.667
            Amplitude of maximum value of a trace in units of x-axis
        fill : bool, optional; default: False
            Fill positive part of traces
        """
        try:
            tab_col = -1
            self.actual_axis = ax
            self.line_colors=[]
            if nt_max == 0:
                nt_max = len(time)
            t = time[nt_min:nt_max]
            ntrace = len(x_pos)
            if ntrace == 0:
                _ = QtWidgets.QMessageBox.warning(
                    None, "Warning",
                    f"No data in {text_t}\nChoose other gather",
                    QtWidgets.QMessageBox.Close)
                return
            nsamp = len(t)
            if nt_max == 0:
                nt_max = nsamp
            if np.array(traces).any() is None or len(traces) == 0:
                traces = np.arange(ntrace, dtype=int)
            else:
                traces = np.array(traces, dtype=int)
            xmin = x_pos[traces].min()
            xmax = x_pos[traces].max()
            nt = len(x_pos[traces])-1
            dx = (xmax-xmin)/(nt-1)
            self.x_zoom_min = xmin
            self.x_zoom_max = xmax
            self.t_zoom_min = nt_min*self.data.dt+self.data.t0
            self.t_zoom_max = nt_max*self.data.dt+self.data.t0
            xmin -= 1
            xmax += 1
            d = np.copy(data[:, nt_min:nt_max])
            dmax = np.zeros(ntrace)
            if np.std(d) > 0.:
                if self.gain == "tnorm" or (self.gain == "dist" and
                                            self.dg_flag):
                    dmax = np.abs(d).max(axis=1)
                elif self.gain == "time":
                    fac = np.zeros(len(t))
                    fac[t > 0] = t[t > 0]**self.time_gain
                    for i in traces:
                        d[i, :] = d[i, :]*fac
                    dmax = np.abs(d).max(axis=1)
                elif self.gain == "dist":
                    fac = np.abs(x_pos)**self.dist_gain
                    for i in traces:
                        d[i, :] = d[i, :]*fac[i]
                    dmax = np.ones(ntrace)*np.max(np.abs(d))
                elif self.gain == "agc":
                    for i in traces:
                        d[i, :] = self.agcCalc(d[i, :], self.n_agc_win)
                    dmax = np.abs(d).max(axis=0)
            if np.isclose(dmax.max(), 0.):
                _ = QtWidgets.QMessageBox.warning(
                    None, "Warning",
                    "All traces are empty\nChose other gather",
                    QtWidgets.QMessageBox.Close)
                return
            for i in traces:
                if np.isclose(dmax[i], 0.):
                    continue
                dd = d[i, :]/dmax[i]*amp
                dd[dd > dx] = dx
                dd[dd < -dx] = -dx
                x = np.array(x_pos[i]+dd)
                x = np.clip(x, x-2*dx, x+2*dx)
                x = np.clip(x, xmin, xmax)
                tab_col += 1
                if tab_col > 9:
                    tab_col = 0
                ax.plot(x, time[nt_min:nt_max], c=self.tab_colors[tab_col])
                self.line_colors.append(self.tab_colors[tab_col])
                if fill:
                    nulls = np.ones_like(t)*x_pos[i]
                    ax.fill_betweenx(t, nulls, x, where=nulls <= x)
            ax.grid(which="major", axis="y")
            ax.tick_params(axis='both', labelsize=18)
            ax.set_xlabel(text_x, fontsize=18)
            ax.set_ylabel(text_y, fontsize=18)
            ax.set_title(text_t, fontsize=20)
            ax_xmin, ax_xmax = ax.get_xlim()
            ax_ymin, ax_ymax = ax.get_ylim()
            ax.text(ax_xmin, ax_ymax+(ax_ymax-ax_ymin)*0.01,
                    self.main.dir_start, horizontalalignment="left",
                    verticalalignment="bottom", fontsize=18)
            ax.text(ax_xmax, ax_ymax+(ax_ymax-ax_ymin)*0.01, self.main.dir_end,
                    horizontalalignment="right", verticalalignment="bottom",
                    fontsize=18)
        except Exception as error:
            print(f"Error in seismogram: {error}.")

    def set_ticks(self, minval, maxval, ntick=5, dtick=5.):
        """
        Function

        Parameters
        ---------s-
        minval : float
            Minimum value of data to be represented (first tick will be
                    at a position larger or equal to this value)
        maxval : float
            Maximum value of data to be represented (last tick will be
                    at a position smaller or equal to this value)
        ntick : int, optional
            Approximate number of ticks. The default is 5.
        dtick : float, optional
            Ticks will be placed at multiples of dtick such that the number of
            ticks is approximately nticks. The default is 5.

        Returns
        -------
        ticks : numpy array
            contains the positions of the ticks

        """
        return set_ticks(minval, maxval, ntick, dtick)
