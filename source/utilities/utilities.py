import os
from os.path import exists
from copy import deepcopy
from datetime import datetime, date
import numpy as np
from PyQt5 import QtWidgets, QtCore
import statsmodels.api as sm
from obspy.signal import trigger
import obspy.signal.filter as obs_filter
import scipy.stats
import scipy.signal
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path
from matplotlib import tri
from mpl_toolkits.axes_grid1 import make_axes_locatable
import colorcet as cc
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
from scipy.signal import hilbert, medfilt
from sklearn.linear_model import LinearRegression

from ..plotting.refraPlot import newWindow


class Utilities:
    """
    Contains utility methods for PyRefra:
        min_max
        tauP
        pModel
        envel
        falseColour
        secondDerivative
        sta_lta
        akaike
        max_min_amplitudes
        spectrum
        filterAll
        filterTrace
            onPress
        ffilter
        frequencyFilter
        FK_filt
        airWaveFilter
        velocityFilter
        v_nmo
        inversion
            vel_scale
        invCol
        prepareSOFI2D
        atten_amp
    """

    def __init__(self, main, files, data, traces, geom, window):
        self.main = main
        self.files = files
        self.data = data
        self.traces = traces
        self.geom = geom
        self.window = window
        self.pmodel = []
        self.thick_l = []
        self.thick_r = []
        self.vels_l = []
        self.vels_r = []
        self.depths_l = []
        self.depths_r = []
        self.tints_l = []
        self.tints_r = []
        self.pk = []
        self.coor_x = []
        self.coor_y = []
        self.x_coor = []
        self.y_coor = []
        self.v_red = 400.
        self.fmax_plot = 400
        self.kmax_plot = 1000
        self.fmin = 0.
        self.fmax = 0.
        self.high_cut_flag = False
        self.low_cut_flag = False
        self.filtered = False
        self.filtered_all = False
        self.all_freq_filter = False
        self.finish = False
        self.picked = False
        self.neg_flag = False
        self.mod_v1 = 800.
        self.mod_v2 = 4000.
        self.mod_h = 3.
# ax_tt will contain the plot of the attenuation of maximum amplitudes
        self.ax_att = None
        self.w_tau = None
        self.w_env = None
        self.w_fcol = None
# w_amp is the window into which ax_amp is plotted
        self.w_amp = None
# ax_amp will contain the animated plot of the wave amplitudes
        self.ax_amp = None
        self.w_pseudo = None
        self.nd_start = 0
        self.nd_end = 0
        self.max_lag = 0
        self.lin = None
        self.cidpress = None
        self.nd_start_slta = 0
        self.tick_size_mod = 16
        self.tick_size_sec = 12
        self.scheme = None
        self.px = None
        self.gx = []
        self.sx = []
        self.sens = []
        self.offsets = None
# The following variables are only used in function inversion, however, their
# value must be kept after leaving the function in case "C" is pressed later in
# order to recreate the image of the tomography result. Most of the variables
# are initialized here only with dummy values.
#
# w_tomo is the window into which the figure is plotted
        self.w_tomo = None
# figinv is the figure into which the tomography result is plotted
        self.figinv = None
# gs will contain the subplots
        self.gs = None
# ax_mod will contain the final model plot
        self.ax_mod = None
# ax_start will contain the plot of the starting model
        self.ax_start = None
# ax_rays will contain the plot of the coverage and rays of the final model
        self.ax_rays = None
# ax_tt will contain the plot of measured traval times
        self.ax_tt = None
# ax_diff will contain the plot of the differences between calculated and
# measured times
        self.ax_diff = None
# ax_av_diff will contain the plot of the average misfits of shots and
# receivers
        self.ax_av_diff = None
# ax_chi will contain the plot of the evolution of chi**2
        self.ax_chi = None
# zmax will contain the maximum depth of the inversion model
        self.zmax = 0.
# xax_min and xax_max will contain the limits of the x_axis for model plots
        self.xax_min = 0.
        self.xax_max = 0.
# model_colors contains the color maps that may be used for model plots
        self.model_colors = ["Special P (cyan=1500)", "Special S (cyan=500)",
                             "rainbow", "viridis", "seismic"]
# smooth, s_fact and zSmooth contain default smoothing parameters for the
#    inversion. Values may be modified interactively
        self.smooth = 200
        self.s_fact = 0.8
        self.zSmooth = 0.2
# maxiter contins the maximum number of iterations allowed (0 = no limit)
        self.maxiter = 0
# vmin and vmax contain the velocities at the top and the bottomof the initial
# model
        self.vmin = 200
        self.vmax = 4000
# vmin_limit and vmax_limit contain the limits of the velocity search space
        self.vmin_limit = 150
        self.vmax_limit = 6000
# v_scale_min and v_scale_max contain the limits of the colorscale for model
# plot
        self.v_scale_min = 200.
        self.v_scale_max = 6000.
# zmax_plt is the maximum depth to which the models are plotted
# (may be different from zmax)
        self.zmax_plt = 20.
# If ray_flag=True plot rays of last iteration onto plot of final model
        self.rays_flag = False
# mgr contains the pygimli manager
        self.mgr = None
# startmodel contains the starting model
        self.startModel = None
# endModel contains the final inverted model
        self.endModel = None
# mesh_coor_all contains the center coordinates of all model triangles
        self.mesh_coor_all = []
# mesh_coor only those where cover is not zero
        self.mesh_coor = []
# cell_rays contains the length of rays for every cell and every iteration
        self.cell_rays = []
# ncover is the sum of cell_rays over iterations
        self.ncover = []
# cover contains the coverage of all triangles (summed ray lengths) divided by
# cell size
        self.cover = []
# path is the path where pygimli writes the results
        self.path = ""
# p_aim is the path where results are copied at the end
# (path is too complicated)
        self.p_aim = ""
# cov_txt is the title of the coverage plot
        self.cov_txt = ""
# dat contains the measured arrival times used for the inversion
        self.dat = []
# calc contins the calculated arrival times from the final model
        self.calc = []
# min_end and max_end are the minimum and maximum velocities of the final model
        self.min_end = 0
        self.max_end = 0
# q1_start and q2_start are 1% and 99% quantiles of starting model velocities
        self.q1_start = 0.01
        self.q2_start = 0.01
# q1_end and q2_end are 1% and 99% quantiles of final model
        self.q1_end = 0.99
        self.q2_end = 0.99
# levels are color scale levels for model plot
        self.levels = []
# levs are annotated levels for color scale of model plot
        self.levs = []
# tick_x and tick_y are annotated ticks of model axes
        self.tick_x = []
        self.tick_y = []
# triang contains triangles for interpolation of model plot
        self.triang = []
# mask contains all triangles to be masked (coverage is zero)
        self.mask = []

    def min_max(self, data, half_width=3):
        """
        Find all relative minima and maxima in a data vector. A maximum is
        found if a value at position i of the vector is larger than or equal to
        all other values in a range [i-half_width:i+half_with] and at the same
        time strictly larger than all values of one side. Sometimes, e.g. if
        seismograms are saturated, a slope exists on one side, but the values
        are constant on the other side. The one-sided test is to avoid that a
        flat curve with 2*half_width constant values is also considered as
        local maximum. Equivalent scheme for definition of a minimum.
        In addition, the function reduces possible local maxima and minima such
        that a maximum is always followed by a minimum and vice versa.
        If several local maxima follow each other (i.e. no wide enough local
        minimum exists between them), the program searches the strongest one of
        the subsequent maxima or, if several equal maximum values exists, takes
        as position of the maximum the center point between those multiple
        maxima (again for saturated curves).

        Parameters
        ----------
        data : numpy float vectore
            data vector to be analysed
        half_width: integer, optional. Defult: 3
                value giving the number of samples analyzed to
                both sides of every data sample. The first and last
                "half_width" samples of the vector are not analyzed.

        Returns
        -------
        max_pos : numpy integer vector
            All position numbers where there is a relative maximum in vector
            "data"
        max_val : numpy float vector
            Values at maximum positions
        min_pos : numpy integer vector
            All position numbers where there is a relative minimum in vector
            "data"
        min_val: numpy float vector
            Values at minimum positions
        """
        N = len(data)
        NN = np.arange(N, dtype='int')
# extreme_pos (extreme_neg) will contain the sum of all values <= (>=) the
# central value
# half_extreme_pos (half_extreme_neg) will contain the maximum of the
# number of values < (>) the central value on the left and the right side
# A maximu (minimum) is found if extreme_xxx[i]==(2*half_width+1) and if
# half_extreme_xxx[i]==half_width.
        extreme_pos = np.zeros(N)
        half_extreme_pos = np.zeros(N)
        extreme_neg = np.zeros(N)
        half_extreme_neg = np.zeros(N)
# Define vectors with ones have lengths of the full or half sliding window
        test = np.ones(2*half_width+1)
        half_test = np.ones(half_width)
        width = 2*half_width+1
# Start loop over data points
# Sum of neigbouring points for which value[i] <= (>=) value[test_point]
        for k in range(half_width, N-half_width):
            extreme_pos[k] =\
                sum(test[data[k]-data[k-half_width:k+half_width+1] >= 0])
            extreme_neg[k] =\
                sum(test[data[k]-data[k-half_width:k+half_width+1] <= 0])
# Sum of neighbouring values to the left (half1) and right (half2) < value
# [test_point]
            half1 = sum(half_test[data[k]-data[k-half_width:k] > 0])
            half2 = sum(half_test[data[k]-data[k+1:k+half_width+1] > 0])
            half_extreme_pos[k] = max(half1, half2)
# Sum of neighbouring values to the left (half1) and right (half2) > value
# [test_point]
            half1 = sum(half_test[data[k]-data[k-half_width:k] < 0])
            half2 = sum(half_test[data[k]-data[k+1:k+half_width+1] < 0])
            half_extreme_neg[k] = max(half1, half2)
# Search all points that fulfill the criteria for local maximum and minimum
        max_pos = NN[(extreme_pos == width) & (half_extreme_pos == half_width)]
        max_val = data[max_pos]
        min_pos = NN[(extreme_neg == width) & (half_extreme_neg == half_width)]
        min_val = data[min_pos]
        del extreme_pos, extreme_neg, half_extreme_pos, half_extreme_neg
# mx_sig is a vector with length equal to number of found maxima with +1
# mn_sig is a vector with length equal to number of found maxima with -1
#   These vectors will be used to know which position is a maximum, which one
#   a minimum, once all extrema are concatenated in a single vector in order to
#   intercalate maxima and minima and to find places where multiple maxima or
#   minima follow each other
        mx_sig = np.ones(len(max_pos))
        mn_sig = -np.ones(len(min_pos))
# Concatenate positions, values and signs of maxima and minia into a single
#   vector for each of them
        signs = np.concatenate((mx_sig, mn_sig))
        positions = np.concatenate((max_pos, min_pos))
        values = np.concatenate((max_val, min_val))
    # Order the concatenated vectors by positions
        iord = np.argsort(positions)
        pord = positions[iord]
        vord = values[iord]
        sord = signs[iord]
        ls = len(sord)
# Prepare lists that will contain positions, values and signs of alternating
#   extreme values (avoiding having several maxima (minima) following each
#   other without a minumum (maximum) between them).
        pos = []
        val = []
        sig = []
        i = 1
# Start loop over concatenated extreme positions
# If sign of position [i] is different from position [i-1] accept position
#   [i-1] into a new list
        while i < ls:
            if sord[i] != sord[i-1]:
                pos.append(pord[i-1])
                val.append(vord[i-1])
                sig.append(sord[i-1])
            if i == ls-1:
                if sord[i] != sord[i-1]:
                    pos.append(pord[i])
                    val.append(vord[i])
                    sig.append(sord[i])
                    i += 1
                break
            else:
                i += 1
                continue
# if sign of position i is the same as the one of position i-1 search for next
#   position where sign changes
            i1 = i+1
            for i1 in range(i+1, ls):
                if sord[i] != sord[i1]:
                    break
            if i1 < i:
                break

# Search maximum values of the positions having the same sign
#   the chosen position is the average position of all equal maximum (minimum)
#   values. If one of the relative maxima (minima) has the strongest value, its
#   position and value will be copied into the new list.
            if sord[i] > 0:
                mx = np.where(vord[i:i1] == max(vord[i:i1]))
                mpos = int(np.mean(pord[i:i1][mx]))
                pos.append(mpos)
                val.append(max(vord[i:i1]))
                sig.append(sord[i])
            else:
                mx = np.where(vord[i:i1] == min(vord[i:i1]))
                mpos = int(np.mean(pord[i:i1][mx]))
                pos.append(mpos)
                val.append(min(vord[i:i1]))
                sig.append(sord[i])
            i = i1+1
        del max_pos, max_val, min_pos, min_val, iord, pord, vord, sord
    # Transform lists to numpy arrays
        pos = np.array(pos)
        val = np.array(val)
        sig = np.array(sig)
    # Separate again relative maxima from relative minima
        max_val = val[sig > 0]
        max_pos = pos[sig > 0]
        min_val = val[sig < 0]
        min_pos = pos[sig < 0]
        del pos, val, sig, positions, signs, values
        return max_pos, max_val, min_pos, min_val

    def tauP(self):
        """
        Function calculates tau-pi transform and plots result into actual
        window For the moment, this is just for plotting purpose, no use of
        this data is made.

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "tau_p"
        results, okButton = self.main.dialog(
            ["min velocity [m/s]",
             "max velocity [m/s]",
             "velocity step [m/s]"],
            ["e", "e", "e"],
            ["100", "3100", "50"], "tau_p parameters")
        if okButton is False:
            print("Tau-p calculation cancelled")
            self.main.function = "main"
            return
        vmin = float(results[0])
        vmax = float(results[1])
        dv = float(results[2])
        velocities = np.arange(vmin, vmax, dv)
        tau = np.arange(0, self.data.tmax/2, self.data.dt)
        v_tau = np.zeros((len(velocities), len(tau)))
        if self.window.dg_flag:
            print("Tau-p transform can only be done for shot, file or "
                  + "receiver gather")
            return
        offset = np.abs(self.traces.offset[self.traces.plotted])
# Since calculation takes time, install a progress bar that reports progress on
#    calculated velocities
        progressBar = QtWidgets.QProgressBar(None)
        self.window.mplvl.addWidget(progressBar)
        progressBar.show()
        progressBar.setValue(0)
        nmax = np.size(self.window.v, 1)
        for iv, v in enumerate(velocities):
            sample = np.array((offset[offset >= 0]/v -
                              self.data.t0)/self.data.dt, dtype='int')
            for it in range(len(tau)):
                for i, s in enumerate(sample):
                    if s < nmax:
                        v_tau[iv, it] += self.window.v_norm[i, s]
                sample += 1
            completed = int((iv+1)/len(velocities)*100)
            progressBar.setValue(completed)
        progressBar.setValue(0)
        self.window.mplvl.removeWidget(progressBar)
        progressBar.close()
        if self.w_tau:
            self.w_tau.close()
            self.w_tau = None
        self.w_tau = newWindow("TauP")
        ax = self.w_tau.fig.subplots()
        ax.pcolormesh(velocities, tau, v_tau.T, cmap=plt.cm.jet,
                      shading='gouraud')
        ax.set_xlabel("Velocity[m/s]", fontsize=18)
        ax.set_ylabel("Intercept time [s]", fontsize=18)
        if self.window.fg_flag:
            ax.set_title(
                "Tau_P, file "
                + f"{self.files.file_numbers[self.window.fig_plotted]}",
                fontsize=20)
        elif self.window.sg_flag:
            ax.set_title(
                "Tau_P, shot "
                + f"{self.traces.shot[self.window.actual_traces[0]]+1}",
                fontsize=20)
        elif self.window.rg_flag:
            ax.set_title(
                "Tau_P, receiver "
                + f"{self.traces.receiver[self.window.actual_traces[0]]+1}",
                fontsize=20)
        ax.tick_params(axis='both', labelsize=18)
        self.w_tau.show()
        self.main.function = "main"

    def pModel(self, model="P"):
        """

        Function determines line slopes and intercept times to calculate
        velocities and interface depths

        Input:
            model: str, optional
            may be "P" for P-wave model (default) or "S" for S-wave model
            Used only for contructing the name of the output file.

        """
        answer = self.main.test_function()
        if not answer:
            return
        if not model:
            model = "P"
        if self.window.dg_flag:
            print("Measuring velocities makes no sense in a distance gather")
            print(f"{model}_Model canceled")
            return
        self.main.function = "P_model"
        self.window.setHelp(self.window.p_model_text)
        n_lay_l = 0
        n_lay_r = 0
        finish = False
        self.thick_l = []
        self.thick_r = []
        self.vels_l = []
        self.vels_r = []
        self.depths_l = []
        self.depths_r = []
        self.tints_l = []
        self.tints_r = []
        x = np.zeros(2)
        t = np.zeros(2)
        figure = self.window.figs[self.window.fig_plotted]
        renderer = figure.canvas.renderer
# copy background picture
        while not finish:
            background = figure.canvas.copy_from_bbox(figure.bbox)
            self.window.followLine(True, n_lay_l, n_lay_r)
            if len(self.window.coor_x) <= 0:
                print(f"{model}-wave model cancelled")
                self.main.function = "main"
                return
            x[0] = self.window.coor_x[0]
            x[1] = self.window.coor_x[1]
            t[0] = self.window.coor_y[0]
            t[1] = self.window.coor_y[1]
# Make sure plotted lines stay within zoomed window
            xmn = self.window.x.min()
            xmx = self.window.x.max()
            nmn = self.window.zooms[self.window.i_zooms][2]
            nmx = self.window.zooms[self.window.i_zooms][3]
            tmn = nmn*self.main.data.dt+self.main.data.t0
            tmx = nmx*self.main.data.dt+self.main.data.t0
            if x[1] > x[0]:
                if x[1] > xmx:
                    t[1] = t[0]+(t[1]-t[0])/(x[1]-x[0])*(xmx-x[0])
                    x[1] = xmx
                elif x[0] < xmn:
                    t[0] = t[0]+(t[1]-t[0])/(x[1]-x[0])*(xmn-x[0])
                    x[0] = xmn
            else:
                if x[1] < xmn:
                    t[1] = t[0]+(t[1]-t[0])/(x[1]-x[0])*(xmn-x[0])
                    x[1] = xmn
                elif x[0] > xmx:
                    t[0] = t[0]+(t[1]-t[0])/(x[1]-x[0])*(xmx-x[0])
                    x[0] = xmx
            if t[1] > t[0]:
                if t[1] > tmx:
                    x[1] = x[0]+(x[1]-x[0])/(t[1]-t[0])*(tmx-t[0])
                    t[1] = tmx
                elif t[0] < tmn:
                    x[0] = x[0]+(x[1]-x[0])/(t[1]-t[0])*(tmn-t[0])
                    t[0] = tmn
            else:
                if t[1] < tmn:
                    x[1] = x[0]+(x[1]-x[0])/(t[1]-t[0])*(tmn-t[0])
                    t[1] = tmn
                elif t[0] > tmx:
                    x[0] = x[0]+(x[1]-x[0])/(t[1]-t[0])*(tmx-t[0])
                    t[0] = tmx
            self.window.axes[self.window.fig_plotted].plot(x, t, "k")
            self.window.axes[self.window.fig_plotted].draw(renderer)
            dx = self.window.coor_x[1]-self.window.coor_x[0]
            if np.abs(dx) < np.abs(self.window.x[1]-self.window.x[0]):
                continue
            slope = (self.window.coor_y[1]-self.window.coor_y[0])/dx
            vel = np.round(1/np.abs(slope), 0)
            tint = np.round(self.window.coor_y[0]-slope*self.window.coor_x[0],
                            4)
            tt = self.window.coor_y[0]-slope*self.window.coor_x[0]
            if self.window.side > 0:
                if n_lay_r == 0:
                    tk = 0
                    depth = 0
                else:
                    depth = 0
                    if n_lay_r > 1:
                        for i in range(n_lay_r):
                            fac = self.vels_r[i]/vel
                            if fac < 1:
                                tt = tt-2*self.thick_r[i]/self.vels_r[i]\
                                     * np.sqrt(1-fac**2)
                            else:
                                tt = np.nan
                    fac = self.vels_r[n_lay_r-1]/vel
                    if fac < 1:
                        tk = np.round(tt*self.vels_r[n_lay_r-1] /
                                      (2*np.sqrt(1-fac**2)), 2)
                    else:
                        tk = np.nan
                if n_lay_r > 0:
                    depth = self.depths_r[-1]+tk
                else:
                    depth = 0
            else:
                if n_lay_l == 0:
                    tk = 0
                    depth = 0
                else:
                    depth = 0
                    if n_lay_l > 1:
                        for i in range(n_lay_l):
                            fac = self.vels_l[i]/vel
                            if fac < 1.:
                                tt = tt-2*self.thick_l[i]/self.vels_l[i]\
                                     * np.sqrt(1-fac**2)
                            else:
                                tt = np.nan
                    fac = self.vels_l[n_lay_l-1]/vel
                    if fac < 1:
                        tk = np.round(tt*self.vels_l[n_lay_l-1] /
                                      (2*np.sqrt(1-fac**2)), 2)
                    else:
                        tk = np.nan
                if n_lay_l > 0:
                    depth = self.depths_l[-1]+tk
                else:
                    depth = 0
            msgbox = QtWidgets.QMessageBox()
            msgbox.setIcon(QtWidgets.QMessageBox.Information)
            if self.window.side > 0:
                msgbox.setWindowTitle(f"Positive direction Layer {n_lay_r+1}")
            else:
                msgbox.setWindowTitle(f"Negative direction Layer {n_lay_l+1}")
            msgbox.setText("Another layer: Yes\n Finish model: No\n "
                           + "Redo this layer: Retry")
            msgbox.setInformativeText("Click 'More details' to see layer "
                                      + "parameters")
            msgbox.setDetailedText(f"Velocity: {vel:0.0f} m/s\nThickness: "
                                   + f"{tk:0.2f} m\nDepth of top: "
                                   "{depth:0.2f} m\n"
                                   f"Intercept: {tint*1000:0.1f} ms")
            msgbox.setStandardButtons(QtWidgets.QMessageBox.Yes |
                                      QtWidgets.QMessageBox.No |
                                      QtWidgets.QMessageBox.Retry)
            msgbox.setDefaultButton(QtWidgets.QMessageBox.Yes)
            msgbox.exec_()
            if msgbox.result() == QtWidgets.QMessageBox.Yes:
                if self.window.side > 0:
                    self.vels_r.append(vel)
                    self.tints_r.append(tint)
                    self.depths_r.append(depth)
                    self.thick_r.append(tk)
                    n_lay_r += 1
                else:
                    self.vels_l.append(vel)
                    self.tints_l.append(tint)
                    self.depths_l.append(depth)
                    self.thick_l.append(tk)
                    n_lay_l += 1
# If line has to be redone, remove the last two lines: the one created during
# followLine and the black line created during lin.append...
# It seems that remove cannot remove two itens at a time,
# therefore the two calls.
            elif msgbox.result() == QtWidgets.QMessageBox.Retry:
                self.window.axes[self.window.fig_plotted].lines[-1].remove()
                self.window.axes[self.window.fig_plotted].lines[-1].remove()
# Plot stored background
                figure.canvas.restore_region(background)
                figure.canvas.blit(figure.bbox)
                continue
            else:
                finish = True
                if self.window.side > 0:
                    self.vels_r.append(vel)
                    self.tints_r.append(tint)
                    self.depths_r.append(depth)
                    self.thick_r.append(tk)
                    n_lay_r += 1
                else:
                    self.vels_l.append(vel)
                    self.tints_l.append(tint)
                    self.depths_l.append(depth)
                    self.thick_l.append(tk)
                    n_lay_l += 1
                now = datetime.now()
                c_time = now.strftime("%H-%M-%S")
                today = date.today()
                d1 = today.strftime("%Y-%m-%d")
                if self.window.rg_flag:
                    trec = self.traces.receiver[self.window.actual_traces[0]]+1
                    fname = f"receiver_{trec:0>5}_{d1}_{c_time}.1D{model}"
                elif self.window.sg_flag:
                    trec = self.traces.shot[self.window.actual_traces[0]]+1
                    fname = f"shot_point_{trec:0>5}_{d1}_{c_time}.1D{model}"
                elif self.window.fg_flag:
                    trec = self.files.numbers[self.window.fig_plotted]
                    fname = f"file_{trec:0>5}_{d1}_{c_time}.1D{model}"
                with open(fname, "w") as fo:
                    if n_lay_l > 0:
                        fo.write("Negative direction:\n")
                        fo.write("vel[m/s]  depth[m]  thick[m] intercept[s]\n")
                        for i, d in enumerate(self.depths_l):
                            fo.write(f"{int(self.vels_l[i])} {d:0.2f} "
                                     + f"{self.thick_l[i]:0.2f} "
                                     + f"{self.tints_l[i]:0.4f}\n")
                        if n_lay_r > 0:
                            fo.write("\n")
                    if n_lay_r > 0:
                        fo.write("Positive direction:\n")
                        fo.write("vel[m/s]  depth[m]  thick[m] intercept[s]\n")
                        for i, d in enumerate(self.depths_r):
                            fo.write(f"{int(self.vels_r[i])} {d:0.2f} "
                                     + f"{self.thick_r[i]:0.2f} "
                                     + f"{self.tints_r[i]:0.4f}\n")
        self.window.setHelp(self.window.main_text)
        self.main.function = "main"

    def sModel(self):
        """

        Function determines line slopes and intercept times to calculate
        velocities and interface depths

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.pModel("S")

    def envel(self):
        """
        Function calculates envelope of seismic traces of actual gather and
        plots them into the actual window. The envelope is calculated using
        Scipy.signal.hilbert

        It also calculates the second derivative of the envelope and plots it
        as additional trace on top of the envelope.

        Finally the envelope and second derivative is saved into an ASCII file
           envNNNNN.asc where NNNNN is the number of the actual window

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
# Save data into temporary array
        v_data = deepcopy(self.window.v)
        v_env = np.zeros((v_data.shape[0]*2, v_data.shape[1]))
        ienv = -2
        x = np.zeros(v_data.shape[0]*2)
        traces = []
        time = self.data.t0 + self.data.dt*np.arange(self.data.nsamp)
        progressBar = QtWidgets.QProgressBar(self.window)
        self.window.mplvl.addWidget(progressBar)
        progressBar.show()
        progressBar.setValue(0)
        for i in range(v_data.shape[0]):
            ienv += 2
            v_env[ienv, :] = np.abs(scipy.signal.hilbert(self.window.v[i, :]))
# Calculate smoothened 2nd derivative of envelope:
# Calculate average slope over the l1 points before and the l2 points after
# the calculation point and make the difference (no need to divide by dt**2
# since traces are anyhow normalized before plotting)
            l1 = 30
            l2 = 10
            v_env[ienv+1, :] = self.secondDerivative(v_env[ienv, :], l1, l2)
            x[ienv] = self.window.x[i]-0.2
            x[ienv+1] = self.window.x[i]+0.2
            traces.append(ienv)
            traces.append(ienv+1)
            completed = int((i+1)/self.window.v.shape[0]*100)
            progressBar.setValue(completed)
        if self.w_env:
            self.w_env.close()
            self.w_env = None
        progressBar.setValue(0)
        self.window.mplvl.removeWidget(progressBar)
        progressBar.close()
        self.w_env = newWindow("Envelopes")
        ax = self.w_env.fig.subplots()
        if self.window.fg_flag:
            nfig = self.files.file_numbers[self.window.fig_plotted]
            text_t = f"File {nfig}: envelopes and their 2nd derivative"
        elif self.window.sg_flag:
            nfig = self.traces.shot[self.window.actual_traces[0]]+1
            text_t = f"Shot {nfig}: envelopes and their 2nd derivative"
        elif self.window.rg_flag:
            nfig = self.traces.receiver[self.window.actual_traces[0]]+1
            text_t = "Receiver " +\
                + f"{nfig}: envelopes and their 2nd derivative"
        else:
            text_t = f"distance_{self.window.actual_dist:0.0f}: "\
                + "envelopes and their 2nd derivative"
        self.window.seismogram(
            ax, time, x, v_env, fill=False, amp=self.window.amp_plt,
            traces=traces, nt_min=self.window.nt_mn, nt_max=self.window.nt_mx,
            text_x="Shot offset (m)", text_t=text_t)
        self.window.picksPlot(self.w_env.fig, ax)
        self.w_env.show()
        v_data = deepcopy(self.window.v)
        self.window.v = v_env
        file_out = f"env{self.window.fig_plotted+1:0>5d}.asc"
        self.data.saveASCII(file_out)
# Restore data array
        self.window.v = deepcopy(v_data)

    def falseColour(self, code=None):
        """
        Function plots a 2D interpolated false color plot of different
        transformed record sections. User may chose up to 3 indicators out of
        the following possibilities:

        Instant frequency (typically higher frequencies in noise than signal)
        Envelop
        2nd derivative of envelop
        2nd derivative of the data
        STA/LTA transform
        Akaike Information Criterium
        Relation of local minima and maxima (max-min_after)/(max-min_before)
        The color associated to the different indicators is defined by the
        order of clicking (first: red, second: green, third: blue)

        Parameters
        ----------
        code : int, optional
            Actually not used

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "False_Colour"
        self.window.setHelp(self.window.falseCol_text)
        self.nd_start = self.window.zooms[self.window.i_zooms][2]
        self.nd_end = self.window.zooms[self.window.i_zooms][3]
        time = self.data.t0 + self.data.dt*np.arange(self.nd_start,
                                                     self.nd_end)
        nt0 = np.where(time >= 0)[0][0]
        ndat = len(time)
        tcor_max = min(np.round(time.max()/1.5, 2), 0.3)
# Call dialog box with checkboxes for all of the mentioned indicators
        labels = ["Instant Frequency", "Envelope", "2nd Derivative Env",
                  "2nd Derivative data", "Sta/Lta", "Akaike AIC",
                  "Max-min relation", "Autocorrelation",
                  "\n2 layer model for\nminimum arrival time",
                  "Velocity layer 1 [m/s]", "Velocity layer 2 [m/s]",
                  "Interface depth [m]", "Low-cut  frequency [Hz]",
                  "High-cut frequency [Hz]", "Max lag [s]"]
        types = ["c", "c", "c", "c", "c", "c", "c", "c", "l", "e", "e", "e",
                 "e", "e", "e"]
        values = [None, None, None, None, None, None, None, None, None,
                  self.mod_v1, self.mod_v2, self.mod_h, self.fmin, self.fmax,
                  tcor_max]
        e_results, btn = self.main.dialog(labels, types, values,
                                          "Indicators to chose")
        ck_order = []
        ck_results = []
        for i, t in enumerate(types):
            if t == 'c':
                ck_order.append(int(e_results[i]))
                if ck_order[-1] > -1:
                    ck_results.append(True)
                else:
                    ck_results.append(False)
            else:
                ck_order.append(-1)
                ck_results.append(False)
        sta_lta_trigger_level = 0.2
        sta_len = 0.005
        lta_len = 0.05
        nsta = int(sta_len/self.data.dt)
        nlta = int(lta_len/self.data.dt)
        self.nd_start_slta = int(
            max(0, self.window.zooms[self.window.i_zooms][2]-nlta))
        if not btn:
            return
        self.mod_v1 = float(e_results[9])
        self.mod_v2 = float(e_results[10])
        self.mod_h = float(e_results[11])
        self.fmin = float(e_results[12])
        self.fmax = float(e_results[13])
        self.max_lag = min(
            int(float(e_results[14])/self.data.dt), int(ndat/1.5))
        if self.fmin == 0.:
            self.low_cut_flag = False
        else:
            self.low_cut_flag = True
        if np.isclose(self.fmax, 0.):
            self.high_cut_flag = False
        else:
            self.high_cut_flag = True
        if self.low_cut_flag | self.high_cut_flag:
            self.frequencyFilter(trace_filt=-1, plot_flag=False)
        v12 = self.mod_v1**2
        v22 = self.mod_v2**2
        d_pass = int(2*self.mod_h*np.sqrt((v22+v12)/(v22-v12)))

        labs = []
        nlabs = np.zeros(3, dtype=int)
        for i, r in enumerate(ck_results):
            if r:
                labs.append(labels[i])
                nlabs[ck_order[i]] = len(labs)-1
# Define data to be treated (all traces, but only data of time zoom)
        nrec = self.window.v.shape[0]
        # kurt = np.zeros_like(self.window.v)
        # skewness = np.zeros_like(self.window.v)
        c = np.zeros((ndat, nrec, 3))
        indicator = np.zeros((ndat, nrec))
        self.pk = np.zeros(nrec)
        calc = np.zeros((ndat, 3))
# Initialize progress bar
        progressBar = QtWidgets.QProgressBar(self.window)
        self.window.mplvl.addWidget(progressBar)
        progressBar.show()
        progressBar.setValue(0)
        acor_av = np.zeros(ndat)
# Start loop over traces
        for i in range(nrec):
            dat = self.window.v[i, self.nd_start:self.nd_end]
            if abs(self.window.x[i]) < d_pass:
                tmin_eval = abs(self.window.x[i])/self.mod_v1
            else:
                tmin_eval = 2*self.mod_h/self.mod_v1*np.sqrt(1-(v12/v22))\
                    + abs(self.window.x[i])/self.mod_v2
            nmin = int((tmin_eval-self.data.t0)/self.data.dt)-self.nd_start
# If trace has not data skip
            s_dat = np.std(dat)
            if np.isclose(s_dat, 0):
                calc[:, :] = 0
                continue
# If first checkbox has been checked, calculate instantaneous frequency
            if ck_results[0]:
                analyt = hilbert(dat)
                pha = np.unwrap(np.angle(analyt))
                sig = np.zeros(ndat)
                sig[1:ndat] = np.diff(pha)/(2.0*np.pi*self.data.dt)
                instFreq = medfilt(sig, 11)
                instFreq[instFreq < 0] = 1.E-3
                f_q10 = np.quantile(instFreq, 0.1)
                instFreq[instFreq < f_q10] = f_q10
                f_q90 = np.quantile(instFreq, 0.9)
                instFreq[instFreq > f_q90] = f_q90
                instFreq = instFreq.max()-instFreq
# If envelope is used, calculate it
            if ck_results[1] or ck_results[2]:
                env = obs_filter.envelope(dat)
# If first checkbox has been checked, fill the corresponding part in array c
#   with instantaneous frequency depending on the order of checing stored in
#   ch_order
            if ck_results[0]:
                instFreq[:nmin] = 0.
                c[:, i, ck_order[0]] = np.flip(instFreq)/max(instFreq)
# If second checkbox has been checked, fill the corresponding part in array c
#   with envelope depending on the order of checing stored in ch_order
            if ck_results[1]:
                env[:nmin] = 0.
                e_q90 = np.quantile(env, 0.9)
                env[env > e_q90] = e_q90
                c[:, i, ck_order[1]] = np.flip(env)/max(env)
# If third checkbox has been checked, fill the corresponding part in array c
#   with 2nd derivative of envelope depending on the order of checing stored
#   in ch_order
            if ck_results[2]:
                envDeriv = self.secondDerivative(env, 30, 8)
                envDeriv[envDeriv < 0] = 0.
                envDeriv[:nmin] = 0.
                c[:, i, ck_order[2]] = np.flip(envDeriv)/max(envDeriv)
# If fourth checkbox has been checked, fill the corresponding part in array c
#   with 2nd derivative of data depending on the order of checing stored in
#   ch_order
            if ck_results[3]:
                dataDeriv = self.secondDerivative(dat, 30, 8)
                dataDeriv[dataDeriv < 0] = 0.
                dataDeriv[:nmin] = 0.
                c[:, i, ck_order[3]] = np.flip(dataDeriv)/max(dataDeriv)
# If fifth checkbox has been checked, fill the corresponding part in array c
#   with STA-LTA transform depending on the order of checing stored in ch_order
            if ck_results[4]:
                dsta = self.window.v[i, self.nd_start_slta:self.nd_end]
                s_l = self.sta_lta(dsta, nsta, nlta)[self.nd_start
                                                     - self.nd_start_slta:]
                mn = np.nanmin(s_l)
                s_l[np.isnan(s_l)] = mn
                sl_q01 = np.quantile(s_l, 0.01)
                s_l[s_l < sl_q01] = sl_q01
                s_l[:nmin] = sl_q01
                s_l -= sl_q01
                mx = np.nanmax(s_l)
                trig = mx*sta_lta_trigger_level
                itrig = np.where(s_l >= trig)[0][0]
                s_l[itrig+1:] = trig*0.99
                c[:, i, ck_order[4]] = np.flip(s_l)/max(s_l)
# If sixth checkbox has been checked, fill the corresponding part in array c
# with negative AIC value depending on the order of checing stored in ch_order
            if ck_results[5]:
                aic = -self.akaike(dat)
                aic_q10 = np.quantile(aic, 0.1)
                aic[:nmin] = aic_q10
                aic[aic < aic_q10] = aic_q10
                aic -= aic_q10
                c[:, i, ck_order[5]] = np.flip(aic)/max(aic)
# If seventh checkbox has been checked, fill the corresponding part in array c
#   with max/min values depending on the order of checing stored in ch_order
            if ck_results[6]:
                rel = self.max_min_amplitudes(dat, 5)
                rel[:nmin] = 0.
                c[:, i, ck_order[6]] = np.flip(rel)/max(rel)
# Eighth checkbox has been checked, caculate autocorrelation
            if ck_results[7]:
                acor = sm.tsa.acf(dat, nlags=self.max_lag)
                acor_av[nt0:nt0+len(acor)] += acor
                Fourier = np.fft.fft(dat)
                if i == 0:
                    s = np.std(abs(Fourier))
                    if s > 0:
                        Fourier_av = np.copy(Fourier/s)
                    else:
                        Fourier_av = np.copy(Fourier)
                else:
                    if s > 0:
                        Fourier_av += Fourier/s
                acor = (acor-acor.min())/(acor.max()-acor.min())
                c[ndat-(nt0+len(acor)):ndat-nt0, i, ck_order[7]] = \
                    np.flip(acor)
# Advance progress bar
            if np.isclose(self.window.x[i], 0.):
                self.pk[i] = 0.
            else:
                indicator[:, i] =\
                    np.flip(c[:, i, 0]**2+c[:, i, 1]**2+c[:, i, 2]**2)
                self.pk[i] = time[np.argmax(indicator[:, i])]
            completed = int((i+1)/self.window.v.shape[0]*100)
            progressBar.setValue(completed)
# Close progress bar
        progressBar.setValue(0)
        self.window.mplvl.removeWidget(progressBar)
        progressBar.close()
# Draw false color plot
        if len(labs) == 0:
            _ = QtWidgets.QMessageBox.warning(
                None, "Warning",
                " No indicator chosen. \n Restart with SHFT-F ",
                QtWidgets.QMessageBox.Close)
            return
# Check picks for outlyers
# First in negative direction
        for _ in range(2):
            if ck_results[7]:
                continue
            n_neg = np.where(self.window.x <= 0)[0]
            if len(n_neg) > 4:
                xx = -self.window.x[n_neg]
                tt = self.pk[n_neg]
                nlin = 2
                off, slope, inter, x_rn, t_rn, t_xr, _, _ =\
                    self.window.bestLines(xx, tt, n_lines=nlin)
                print("neg off:", off)
                print("vel:", 1/slope)
                print("inter:", inter)
                d_t = t_xr-tt
                pts = np.arange(len(n_neg))
                split1 = pts[xx <= off[1]]
                split2 = pts[xx > off[1]]
                d_std = np.std(d_t[split1])
                for k in split1:
                    if abs(d_t[k]) > 1.5*d_std:
                        self.pk[n_neg[k]] = None
                d_std = np.std(d_t[split2])
                for k in split2:
                    if abs(d_t[k]) > 2*d_std:
                        self.pk[n_neg[k]] = None
# Now in positive direction
            n_pos = np.where(self.window.x >= 0)[0]
            if len(n_pos) > 4:
                xx = self.window.x[n_pos]
                tt = self.pk[n_pos]
                nlin = 2
                off, slope, inter, x_rp, t_rp, t_xr, _, _ =\
                    self.window.bestLines(xx, tt, n_lines=nlin)
                print("pos off:", off)
                print("vel:", 1/slope)
                print("inter:", inter)
                d_t = t_xr-tt
                pts = np.arange(len(n_pos))
                split1 = pts[xx <= off[1]]
                split2 = pts[xx > off[1]]
                d_std = np.std(d_t[split1])
                for k in split1:
                    if abs(d_t[k]) > 1.5*d_std:
                        self.pk[n_pos[k]] = None
                d_std = np.std(d_t[split2])
                for k in split2:
                    if abs(d_t[k]) > 2*d_std:
                        self.pk[n_pos[k]] = None
        self.w_fcol = newWindow("False colours")
        ax = self.w_fcol.fig.subplots()
        dx = (self.window.x[1]-self.window.x[0])/2.
        dy = self.data.dt/2.
        extent = [self.window.x.min()-dx, self.window.x.max()+dx,
                  time.min()-dy, time.max()+dy]
        _ = ax.imshow(c, extent=extent, aspect="auto")
        if ck_results[7]:
            fac = (self.window.x.max()-self.window.x.min())/6
            xx0 = self.window.x.mean()*2./3.
            x = xx0+fac*acor_av/len(self.window.x)
            ax.plot(x, time, "w")
            avTrace = np.real(np.fft.ifft(Fourier_av))
            xx0 *= 2.
            x = xx0+fac*avTrace/max(np.abs(avTrace))
            ax.plot(x, time, "y")
        else:
            ax.plot(self.window.x, self.pk, "m")
            try:
                ax.plot(-x_rn, t_rn, "w")
            except:
                pass
            try:
                ax.plot(x_rp, t_rp, "w")
            except:
                pass
        if len(labs) == 1:
            self.window.picksPlot(self.w_fcol.fig, ax, "y")
        elif len(labs) == 2:
            self.window.picksPlot(self.w_fcol.fig, ax, "b")
        else:
            self.window.picksPlot(self.w_fcol.fig, ax)
# Set title with indication of meaning of colors
        txt = f"red={labs[nlabs[0]]}"
        if len(labs) > 1:
            txt += f", green={labs[nlabs[1]]}"
            if len(labs) > 2:
                txt += f", blue={labs[nlabs[2]]}"
        if ck_results[7]:
            txt += ";  white: av. acor, yellow: av. trace"
        ax.set_xlabel("Distance[m]", fontsize=18)
        ax.set_ylabel("Time [s]", fontsize=18)
        if self.window.sg_flag:
            sht = self.traces.shot[self.window.actual_traces[0]]
            ax.set_title(f"shot {sht+1}: {txt}", fontsize=20)
        elif self.window.fg_flag:
            strings = self.window.plot_names[self.window.fig_plotted].split()
            nfile = int(strings[1])
            ax.set_title(f"file {nfile}: {txt}", fontsize=20)
        elif self.window.rg_flag:
            rec = self.traces.receiver[self.window.actual_traces[0]]
            ax.set_title(f"receiver {rec+1}: {txt}", fontsize=20)
# Show false colour plot and wait for key stroke
        ax.tick_params(axis='both', labelsize=18)
        self.w_fcol.show()
        self.main.function = "main"

    def secondDerivative(self, data, l_before, l_after):
        """
        Calculation of the smoothed second derivative of a data series.
        For every point, a first slope is calculated over the l_before points
        recorded before the point and a second one ove the l_after points after
        the point. The second derivative is calculated as the difference
        between the slopes, i.e. it is not normalized by the squared sampling
        rate.

        Parameters
        ----------
        data : Numpy array of one dimension; float
            Data series for which to calculate the second derivatives
        l_before : int
            Number of points over which the slope before the calculation point
            is determined
        l_after : int
            Number of points over which the slope after the calculation point
            is determined

        Returns
        -------
        deriv2 : Numpy array of the same shape and type as data
            contains the not normalized second derivatives. The first l_before
            and the last l_after points are zero.

        """
        xx1 = np.arange(l_before)
        xs1 = (l_before*(l_before-1))*0.5
        fac1 = l_before*np.dot(xx1, xx1)-xs1*xs1
        xx2 = np.arange(l_after)
        xs2 = (l_after*(l_after-1))*0.5
        fac2 = l_after*np.dot(xx2, xx2)-xs2*xs2
        ndat = len(data)
        deriv2 = np.zeros(ndat)
# Slope of l1 points before
        for k in range(l_before, ndat-l_after):
            f1 = l_before*np.dot(xx1, data[k-l_before:k])
            f2 = xs1*np.sum(data[k-l_before:k])
            slope1 = (f1-f2)/fac1
# Slope of l2 points after
            f3 = l_after*np.dot(xx2, data[k:k+l_after])
            f4 = xs2*np.sum(data[k:k+l_after])
            slope2 = (f3-f4)/fac2
            deriv2[k] = slope2-slope1
        return deriv2

    def sta_lta(self, data, nsta, nlta):
        """
        Calculate STA-LTA transform using obspy routine

        Parameters
        ----------
        data : Numpy array of one dimension; float
            Data series for which to calculate the second derivatives
        nsta : int
            Window length for STA calculation
        nlta : int
            Window length for LTA calculation

        Returns
        -------
        Array with same shape an dtype as data
            STA-LTA values

        """
        return trigger.recursive_sta_lta_py(data, nsta, nlta)

    def akaike(self, data):
        """
        Function calculates the Akaike Information Criterion for data series
        data (see Akaike, H, ( 1973), Information theory and an extension of
              the maximum likelihood principle,
              Second International Symposium on Information Theory,
              Akademiai Kiado, Budapest, 267-281.)

        Parameters
        ----------
        data : Numpy 1D array; float
            Data series to be analyzed

        Returns
        -------
        aic : Numpy array of same shape and type as data
            AIC values of data series. First 5 and las 5 samples are set to
            zero

        """
        aic = np.zeros_like(data)
        ndat = len(data)
        for i in range(5, ndat-5):
            aic[i] = i*np.log(np.var(data[:i])) +\
                (ndat-i-1)*np.log(np.var(data[i+1:]))
        return aic

    def max_min_amplitudes(self, data, half_width=10):
        """
        Function searches successive minima and maxima and calculates the
        amplitude relation beween the differences between a maximum and its
        following minimum on the one hand and between the same maximum and its
        preceeding minimum on the other hand. This calculation uses the
        observation that often this relationship is assymmetric and biased
        towards values above one for the arrival of the P-wave. The value
        of this relationship is then assigned to the samples that ly between
        the sample that has the average value of the maximum and its preceeding
        minimum and the position of the maximum. This is a bit arbitrary.
        Another possibility might be to search the highest second derivative
        between the two as starting point.

        Parameters
        ----------
        data : Numpy 1D array; float
            Data series to be analyzed
        half_width: int
            a maximum / minimum is found at position i if all values at
            positions [i-half_width,i+half_with] are larger / smaller than the
            value at i.

        Returns
        -------
        rel : Numpy array of same shape and type as data
            values of amplitude relaitonships as described above.

        """
        max_pos, max_val, min_pos, min_val = self.min_max(data, half_width)
        relation = []
        min_pos_before = []
        max_position = []
        max_value = []
        min_value = []
        i_min = 0
        n_min = len(min_pos)
        n0 = self.window.nt_0-self.nd_start
        for i, m in enumerate(max_pos):
            if min_pos[i_min] > m:
                continue
            if min_pos[-1] < m:
                break
            for j in range(i_min, n_min-1):
                if min_pos[j] < m and min_pos[j+1] > m:
                    break
            i_min = j
            if min_pos[j] < n0:
                continue
            if min_val[j] >= max_val[i] or min_val[j+1] > max_val[i]:
                continue
            relation.append((max_val[i]-min_val[j+1])/(max_val[i]-min_val[j])
                            / (m-n0))
            min_pos_before.append(min_pos[j])
            max_position.append(m)
            max_value.append(max_val[i])
            min_value.append(min_val[j])
        rel = np.zeros_like(data)
        for i, r in enumerate(relation):
            der = data[min_pos_before[i]+1:max_position[i]+1] +\
                  data[min_pos_before[i]-1:max_position[i]-1] -\
                  2*data[min_pos_before[i]:max_position[i]]
            j1 = min_pos_before[i] + np.argmax(der)
            damp = (max_value[i]-min_value[i])/2.
            for j2 in range(max_position[i], min_pos_before[i], -1):
                if data[j2] <= damp:
                    break
            k1 = min(j1, j2)
            k2 = max(j1, j2)
            kcen = int((k1+k2)/2.)
            if k2-k1 < 6:
                k1 = kcen-3
                k2 = kcen+3
            rel[k1:k2] = r
        return rel

    def spectrum(self, tr1, tr2):
        """
        Parameters
        ----------
        tr1 : int
        tr2 : int
            tr1 : Number of first trace to be filtered
            tr2: Number of last trace to be filtered
            Usually, two cases are used: tr1 = 0 and tr2 = number of traces
            to filter all traces or tr2 = tr1+1 in order to filter only trace
            number tr1

        Returns
        -------
        None.

        If a filter has already been selected, the function proposes first to
        apply the same filter again. In the dialog box, values may be changed
        manually. If the answer is "n" or the Cancel button is pressed, new
        filter coefficients are defined interactively.

        For this, function calculates the sum of spectra of traces tr1 to tr2-1
        normalized by their standard deviation and plots the absolute values of
        the corresponding average spectrum.
        It then allows the user to pick filter frequencies. For this, the
        procedure is as follows:

        The user left clicks at the low-cut frequency, keeps the mouse button
        pressed and pulls the mouse to the high-cut frquency where the mouse
        button is released.

        If the low-cut frequency has a negative value, no low-cut filter is
        applied. Also if the high-cut frequency is negative, no high-cut filter
        is applied. So, if only high-cut filter is desired, left-click at a
        negative frequency and pull the mouse until the posifive high-cut
        frequency and release it there. If only low-cut filter should be
        applied, click on a positive frequency and pull the mouse to a negative
        value and release it.

        Finally, a dialog box opens and teh frequencies may be changed manually
        and/or confirmed.

        Returns
        -------
        None.

        """
        self.main.function = "spectrum"
        if self.high_cut_flag+self.low_cut_flag:
            results, okButton = self.main.dialog(
                ["Low-cut  frequency", "High-cut frequency",
                 "Use this filter (y/n)?", "Apply to all shots"],
                ["e", "e", "e", "c"],
                [int(self.fmin), int(self.fmax), "y", None],
                "Frequency filter")

            if okButton and results[2].lower() == "y":
                if int(results[3]) < 0:
                    self.all_freq_filter = False
                else:
                    self.all_freq_filter = True
                self.fmin = float(results[0])
                if self.fmin > 0.:
                    self.low_cut_flag = True
                else:
                    self.low_cut_flag = False
                self.fmax = float(results[1])
                if self.fmax > 0.:
                    self.high_cut_flag = True
                else:
                    self.high_cut_flag = False
                print(f"Low-cut  frequency: {int(self.fmin)}\n"
                      + f"High-cut frequency: {int(self.fmax)}")
                return True
        ndat = np.size(self.window.v[0, :])
        tr_spec = self.window.actual_traces[tr1:tr2]
        F = np.zeros(ndat, dtype=np.complex128)
        for it in range(tr1, tr2):
            s = np.std(self.window.v[it, :])
            if s > 0:
                F += np.fft.fft(self.window.v[it, :])/s
        Fabs = np.abs(F[0:int(ndat/5)])
        freq = np.arange(np.size(Fabs))/(ndat*self.data.dt)
        fmax = np.max(freq)
        self.window.drawNew(False)
        self.window.axes[self.window.fig_plotted].plot(freq, Fabs)
        self.window.axes[self.window.fig_plotted].set_xlabel("Frequency[Hz]",
                                                             fontsize=18)
        self.window.axes[self.window.fig_plotted].set_ylabel(
            "Amplitude [a.u.]", fontsize=1)
        if self.window.sg_flag:
            sht = self.traces.shot[tr_spec[0]]
            self.window.axes[
                self.window.fig_plotted].set_title(f"Spectrum, shot {sht+1}",
                                                   fontsize=20)
        elif self.window.fg_flag:
            self.window.axes[self.window.fig_plotted].set_title(
                    f"Spectrum, file {self.window.actual_shot}", fontsize=20)

        elif self.window.rg_flag:
            rec = self.traces.receiver[tr_spec[0]]
            self.window.axes[self.window.fig_plotted].set_title(
                f"Spectrum, receiver {rec+1}", fontsize=20)
        else:
            dist = abs(self.traces.offset[tr_spec[0]])
            self.window.axes[self.window.fig_plotted].set_title(
                f"Spectrum, distance {dist:0.1f}", fontsize=20)
        self.window.axes[self.window.fig_plotted].tick_params(
            axis='both', labelsize=18)

        self.window.setHelp(self.window.spectrum_text)
# Pick low_cut and high-cut frequencies, first low_cut
        self.coor_x = []
        self.coor_y = []
        self.main.function = "spectrum"
        self.window.followLine(True, [])
# (xline,yline) are the frquency and amplitude coordinates of the points
# defining the line
        xline = np.array(self.window.coor_x)
        if len(xline) == 0:
            print("No filter chosen, filtering cancelled")
            return False
        print("New frequency filter defined")
        if xline[0] < 0 or xline[0] > fmax:
            self.fmin = 0
            self.low_cut_flag = False
        else:
            self.fmin = round(xline[0], 0)
            self.low_cut_flag = True
        if xline[1] < 0 or xline[1] > fmax:
            self.fmax = 0
            self.high_cut_flag = False
        else:
            self.fmax = round(xline[1], 0)
            self.high_cut_flag = True
        results, okButton = self.main.dialog(
            ["Low-cut  frequency", "High-cut frequency", "Apply to all shots"],
            ["e", "e", "c"], [int(self.fmin), int(self.fmax), None],
            "Frequency filter")

        if okButton is False:
            print("Frequency filter cancelled")
            self.high_cut_flag = False
            self.low_cut_flag = False
            return
        try:
            self.fmin = float(results[0])
        except (IndexError, ValueError):
            self.fmin = 0.
        if self.fmin > 0.:
            self.low_cut_flag = True
        else:
            self.low_cut_flag = False
        try:
            self.fmax = float(results[1])
        except (IndexError, ValueError):
            self.fmax = 0.
        if self.fmax > 0.:
            self.high_cut_flag = True
        else:
            self.high_cut_flag = False
        if int(results[2]) < 0:
            self.all_freq_filter = False
        else:
            self.all_freq_filter = True
        print(f"Low-cut  frequency: {int(self.fmin)}\n"
              + f"High-cut frequency: {int(self.fmax)}")
        return True

    def filterAll(self):
        """
        Apply frequency filter to all shots

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.frequencyFilter(-1)
        self.filtered = True
        self.filtered_all = True

    def filterTrace(self):
        """
        Chose one trace to be filtered
        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "filt_trace"
        self.finish = False
        self.picked = False
        self.window.setHelp(self.window.trace_filter_text)

        def onPress(event):
            self.window.searchTrace(event.xdata)
            self.frequencyFilter(self.window.n_tr)
            self.main.function = "main"

        self.x_coor = []
        self.y_coor = []
        self.lin, = self.window.axes[self.window.fig_plotted].\
            plot(self.x_coor, self.y_coor, animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect(
            'button_press_event', onPress)
        self.main.function = "main"

    def ffilter(self, data, fmin, fmax, dt):
        """
        Filters data of one single trace.
        Uses Obspy zero-phase filter routines with order 8.

        Needs:
            obspy.signal.filter.filter

        Parameters
        ----------
        data: Numpy 1D float array
                Contains the data to be filtered
        fmin: float
                Low-cut frequency. If 0, only lowpass filtering done
        fmax: float
                High-cut frequency. If 0, only highpass filtering done
        dt:   float
                Time step between samples [s]

        Returns
        -------
        d:    Numpy 1D float array (same shape as data)
                Array with filtered data

        """
        if fmin == 0 and fmax > 0:
            d = obs_filter.lowpass(data, fmax, 1/dt, 8, zerophase=True)
        elif fmin > 0 and fmax == 0:
            d = obs_filter.highpass(data, fmin, 1/dt, 8, zerophase=True)
        elif fmin > 0 and fmax > 0:
            d = obs_filter.bandpass(data, fmin, fmax, 1/dt, 8, zerophase=True)
        return d

    def frequencyFilter(self, trace_filt=-1, plot_flag=True):
        """
        Frequency filterig data
        Frequencies are chosen interactively by mouse click
        ObsPy filter is used with order 8

        Parameters
        ----------
        trace_filt : int, optional; default: -1
            If trace_filt is < 0, all traces on screen are filtered
            Else, only trace number trace_filt is filtered
        plot_flag : Boolean, optional; default: True
            If True, filtered data are plotted. False is usually used when a
            whole series of shot gathers is to be treated

        Returns
        -------
        None.
        """
        if trace_filt < 0 or not plot_flag:
            tr1 = 0
            tr2 = self.main.window.actual_number_traces
        else:
            tr1 = trace_filt
            tr2 = tr1+1
        tr_filt = self.main.window.actual_traces[tr1:tr2]
        flag = False
        if plot_flag:
            flag = self.spectrum(tr1, tr2)
            if flag:
                if self.all_freq_filter:
                    tr1 = 0
                    tr2 = len(self.main.traces.trace)
                    tr_filt = np.arange(tr2)
        if flag:
            if self.low_cut_flag is False and self.high_cut_flag is False:
                print("Frequency filter cancelled")
            else:
                if self.fmin == 0 and self.fmax > 0:
                    print(f"\nLow pass filter {int(self.fmax)}Hz of ",
                          f"{len(tr_filt)} traces")
                elif self.fmin > 0 and self.fmax == 0:
                    print(f"\nHigh pass filter {int(self.fmin)}Hz of ",
                          f"{len(tr_filt)} traces")
                elif self.fmin > 0 and self.fmax > 0:
                    print("\nBand pass filter "
                          + "{int(self.fmin)}-{int(self.fmax)}",
                          f"Hz of {len(tr_filt)} traces")
                else:
                    print("Error in frequencies, no filter applied")
                    return
                for i, ii in enumerate(range(tr1, tr2)):
                    if (i+1) % 500 == 0:
                        print(f"     Trace {ii+1}")
                    t = tr_filt[i]
                    fnr = self.traces.file[t]
                    tnr = self.traces.trace[t]
                    dat = self.ffilter(self.data.st[fnr][tnr].data, self.fmin,
                                       self.fmax, self.main.data.dt)
                    self.data.st[fnr][tnr].data = np.float32(dat)
                    if t in self.main.window.actual_traces:
                        iit = np.where(np.array(
                            self.main.window.actual_traces) == int(t))[0][0]
                        if self.traces.amplitudes[t] > 0:
                            self.main.window.v[iit, :] = dat
                        elif self.traces.amplitudes[t] < 0:
                            self.main.window.v[iit, :] = -dat
            self.filtered = True
            print("     All traces filtered")
        if plot_flag:
            self.window.v_set = False
            self.window.drawNew(True)
        self.window.setHelp(self.window.main_text)

    def FK_filt(self, x, dt, data, v_red, nk_el=3, direc=1, plot_flag=False):
        """
        Function does an f-k filter on data eliminating all velocities smaller
        than v_red. To do this, first, a time transform is applied using a
        reduction velocity v_red which corresponds to the velocity where 50%
        amplitude is eliminated.

        Parameters
        ----------
        x : numpy float vector
            Contains the x-coordinates of the traces to be filtered [m]
            These values are necessary for the calculation of time reduction.
            It is supposed that the traces are equidistant and dx corresponds
            to the distance between trace 1 and 2.
            All x must be >= 0 in increasing order
        dt : float
            Time step [s] between two samples
        data : numpy float 2D array [#time_samples, #traces]
            Data to be filtered.
            Before passing the data to the function, you must make sure that
            only traces having been recoreded in the same direction (positive
            or negative) are collected in "data". In addition, traces have to
            be stored in increasing x-direction. This may need some
            pre-treatment of the original data
        v_red : float
            Applied reduction velocity [m/s]. This velocity is also the one
            for which the amplitudes are reduced to 50%
            The default is 400.
        nk_el : int
            The slope to attenuate amplitudes is from -nk_el to +nk_el around
            the (vertical) frequency axis. After time reduction, the filter
            velocity is located on the vertical axis of the 2D Fourier
            transformed data. Data within the quadrants (f>0,k<0) and (f<0,k>0)
            are eliminated
            except for a stripe of +/-nk_el columns where a linear taper is
            applied.
            The default is 3.
        direc : int
            indicator whether data are in positive direction (dir=1, default)
            or negatibe one (dir = -1). Only for plotting title purposes
        plot_flag : bool, optional; default: False
            True if 2D spectrum should be plotted and velovity determined
            interactively. False mainly for filtering of many shot gathers with
            the same filter coefficients

        There is a problem which I did not understand: Positive slopes in the
        data appear on the negative side in the 2D spectrum. Therefore the
        counter-intuitive choice of the quadrants to be zeroed.

        Needs:
            colorcet as cc
            PyRefra.Seg2_Slide

        Returns
        -------
        data_filt : numpy float array like input array "data"
            Filtered data in the same trace order as input data.
            Time reduction has been undone
        v_red : Float
            Applied cut-off velocity [m/s] which may be used in further calls
        Logical value
            True if filter was applied, False if filter cancelled

        """
        from PyRefra import Seg2_Slide
        nsamp = np.size(data, 0)
        ntrace = np.size(data, 1)
        dx = x[1]-x[0]
# nf_damp is the number of coefficients to each side of a slope for which
# the velocity filter coefficients are damped
        nf_damp = 10

# The data array is expanded to twice its size (50% before and after the end)
# and filled with zeros in order to avoid edge effects and in order to allow
# for time reduction without loosing data
# The following lines prepare the expanded data array
        nsamp_exp = 2*nsamp
        ntrace_exp = 2*ntrace
        data_exp = np.copy(data)
        data_exp = np.zeros((nsamp_exp, ntrace_exp))
        nsamp1_data = int(nsamp/2)
        nsamp2_data = nsamp1_data+nsamp
        ntrace1_data = int(ntrace/2)
        ntrace2_data = ntrace1_data+ntrace

# Calculate standard deviation amplitudes of each trace for
#   trace normalization before filtering
        sigmas = np.std(data, axis=0)

# Fill array data_exp with trace-normalized data
        for i in range(ntrace):
            if sigmas[i] > 0:
                data_exp[nsamp1_data:nsamp2_data, i+ntrace1_data] =\
                    data[:, i]/sigmas[i]

# Calculate frequencies
        f0 = 1/(dt*nsamp_exp)
        freq = np.fft.fftfreq(nsamp_exp, dt)
        f_shift = np.fft.fftshift(freq)
# fmax_plot is the maximum frequency to be plotted for the f-k spectrum
        self.fmax_plot = 400
        nf_p_max = np.where(f_shift >= self.fmax_plot)[0][0]
        nf = int(nsamp_exp/2)

# Calculate wavenumbers (as inverse of wavelength, without the factor 2*pi)
        k0 = 1/(dx*ntrace_exp)
        k = np.fft.fftfreq(ntrace_exp, dx)
        self.kmax_plot = np.max(k)
# For plotting purposes, the sign of wavenumbers is changed
# No idea why, but the fft2 function returns a 2D spectrum with
# negative velocities
        k_shift = -np.fft.fftshift(k)

# Calculate and plot 2D Fourier transform
        F = np.fft.fft2(data_exp)
        Fabs = np.log10(np.abs(np.fft.fftshift(F))+1E-10)
#        Fmin = np.min(Fabs)
#        Fmax = np.max(Fabs)
        Fmin = np.quantile(Fabs, 0.05)
        Fmax = np.quantile(Fabs, 0.999)
#        cmp = plt.cm.gist_rainbow_r
        cmp = cc.cm.rainbow4
        if plot_flag:
            self.window.drawNew(False)
            c = self.window.axes[self.window.fig_plotted].\
                pcolormesh(k_shift, f_shift[nf:nf_p_max], Fabs[nf:nf_p_max, :],
                           cmap=cmp, vmin=Fmin, vmax=Fmax, shading='auto')
            self.window.axes[self.window.fig_plotted].set_xlabel(
                "Spacial frequency[1/m]", fontsize=18)
            self.window.axes[self.window.fig_plotted].set_ylabel(
                "Frequency [Hz]", fontsize=18)
            if direc > 0:
                txt = ", positive direction"
            else:
                txt = ", negative direction"
            if self.window.sg_flag:
                sht = self.traces.shot[self.window.actual_traces[0]]
                self.window.axes[self.window.fig_plotted].set_title(
                    f"Spectrum, shot {sht+1}{txt}", fontsize=20)
            elif self.window.fg_flag:
                strings =\
                    self.window.plot_names[self.window.fig_plotted].split()
                nfile = int(strings[1])
                self.window.axes[self.window.fig_plotted].set_title(
                    f"Spectrum, file {nfile}{txt}", fontsize=20)
            elif self.window.rg_flag:
                rec = self.traces.receiver[self.window.actual_traces[0]]
                self.window.axes[self.window.fig_plotted].set_title(
                    f"Spectrum, receiver {rec+1}{txt}", fontsize=20)
            _ = self.window.figs[self.window.fig_plotted].\
                colorbar(c, ax=self.window.axes[self.window.fig_plotted],
                         format='%.1f', label="log(amplitude), a.u.",
                         orientation='vertical', aspect=50)
            self.window.axes[self.window.fig_plotted].tick_params(
                axis='both', labelsize=18)
            self.window.setHelp(self.window.fk_text)
# Show f-k spectrum and wait for key stroke
            self.window.figs[self.window.fig_plotted].canvas.draw()
            self.window.figs[self.window.fig_plotted].canvas.flush_events()

            self.main.function = "vfilt"
            SL = Seg2_Slide(self, v_red)
            SL.sliderReleased = False
            while self.window.verticalSlider.isVisible():
                QtCore.QCoreApplication.processEvents()
            if np.isclose(SL.position, 0.):
                print("f-k filter cancelled")
                return data, v_red, False
            v_red = float(SL.position)
            print(f"\nCut-off velocity: {v_red:0.0f} m/s")
        elif abs(v_red) < 100:
            return data, v_red, False

# Do time reduction (t_reduced = t_original - x/v_red)
# For this, shift the data by n_red samples downward in the expanded array
        if abs(v_red) >= 100:
            for i in range(ntrace1_data, ntrace2_data):
                j = i-ntrace1_data
                t_red = np.abs(x[j])/v_red
                n_red = int(t_red/dt)
                if v_red > 0:
                    data_exp[nsamp1_data-n_red:nsamp_exp-n_red, i] =\
                             data_exp[nsamp1_data:nsamp_exp, i]
                else:
                    data_exp[nsamp1_data:nsamp_exp, i] =\
                             data_exp[nsamp1_data+n_red:nsamp_exp+n_red, i]

# Calculate and plot f-k spectrum of velocity-reduced (turned) traces
        F = np.fft.fft2(data_exp)

# Apply filter in the quadrant (f>0, k<0)
# ik designates all columns to be zeroed, ifr all lines to be zeroed
        F_filt = np.copy(F)
        ik = np.arange(len(k))[k < -nk_el*k0]
        ifr = np.arange(len(freq))[freq > nf_damp*f0]
        for i in ik:
            for j in ifr:
                F_filt[j, i] = F_filt[j, i]*0.
            if self.neg_flag:
                fr_slope = v_red*i*k0
                ifr_slope = np.arange(len(freq))[freq <= 0]
                ff = freq[freq <= 0]
                ifr_slope = ifr_slope[ff > -fr_slope+nf_damp*f0]
                ifs = np.arange(len(freq))[freq <= -fr_slope+nf_damp*f0]
                ff = freq[freq <= -fr_slope+nf_damp*f0]
                ifs = ifs[ff > -fr_slope-nf_damp*f0]
                ff = ff[ff > -fr_slope-nf_damp*f0]
                for j in ifr_slope:
                    F_filt[j, i] = F_filt[j, i]*0.
                for m, j in enumerate(ifs):
                    F_filt[j, i] = F_filt[j, i]*(max(ff)-ff[m])/(2*nf_damp)*f0

# # Now ik designates the columnes where the ramp is applied
        ikk = np.arange(len(k))
        ik0 = ikk[k >= -nk_el*k0]
        ik = ik0[k[ik0] <= nk_el*k0]
        for i in ik:
            fac = (k[i]/k0+nk_el+1)/((nk_el+1)*2)
            for j in ifr:
                F_filt[j, i] = F_filt[j, i]*fac
        if not self.neg_flag:
            iff = np.arange(len(freq))
            if0 = iff[freq >= -nf_damp*f0]
            ifr = if0[freq[if0] <= nf_damp*f0]
            ik = np.arange(len(k))[k < 0]
            for j in ifr:
                fac = 1-(freq[j]/f0+2*nf_damp+1)/(2*nf_damp+1)
                F_filt[j, ik] = F_filt[j, ik]*fac

# Apply filter in the quadrant (f<0, k>0)
# ik designates all columns to be zeroed, ifr all lines to be zeroed
        ifr = np.arange(len(freq))[freq < -5*f0]
        ik = np.arange(len(k))[k > nk_el*k0]
        for i in ik:
            for j in ifr:
                F_filt[j, i] = F_filt[j, i]*0.
            if self.neg_flag:
                fr_slope = v_red*i*k0
                ifr_slope = np.arange(len(freq))[freq >= 0]
                ff = freq[freq >= 0]
                ifr_slope = ifr_slope[ff < fr_slope-nf_damp*f0]
                ifs = np.arange(len(freq))[freq >= fr_slope-nf_damp*f0]
                ff = freq[freq >= fr_slope-nf_damp*f0]
                ifs = ifs[ff < fr_slope+nf_damp*f0]
                ff = ff[ff < fr_slope+nf_damp*f0]
                for j in ifr_slope:
                    F_filt[j, i] = F_filt[j, i]*0.
                for m, j in enumerate(ifs):
                    F_filt[j, i] = F_filt[j, i]*(ff[m]-min(ff))/(2*nf_damp)*f0
# # Now ik designates the columnes where the ramp is applied
        ik0 = ikk[k <= nk_el*k0]
        ik = ik0[k[ik0] >= -nk_el*k0]
        for i in ik:
            fac = 1-(k[i]/k0+nk_el+1)/((nk_el+1)*2)
            for j in ifr:
                F_filt[j, i] = F_filt[j, i]*fac

        if not self.neg_flag:
            if0 = iff[freq >= -nf_damp*f0]
            ifr = if0[freq[if0] <= nf_damp*f0]
            ik = np.arange(len(k))[k < 0]
            for j in ifr:
                fac = (freq[j]/f0+(2*nf_damp+1))/(2*nf_damp+1)
                F_filt[j, ik] = F_filt[j, ik]*fac

# # Plot filtered spectrum
#         if plot_flag:
#             Fabs = np.log10(np.abs(np.fft.fftshift(F_filt))+1E-10)
#             Fmin = np.min(Fabs)
#             Fmax = np.max(Fabs)
#             c = self.window.axes[self.window.fig_plotted].\
#                   pcolormesh(k_shift,f_shift[nf:nf_p_max],Fabs[nf:nf_p_max,:],\
#                   cmap=cmp,vmin=Fmin,vmax=Fmax,shading='auto')
# # Show f-k spectrum and wait for key stroke
#             self.window.figs[self.window.fig_plotted].canvas.draw()
#             self.window.figs[self.window.fig_plotted].canvas.flush_events()

# Calculate the inverse 2D Fourier transform
        d = np.fft.ifft2(F_filt)

# Undo time reduction and undo trace normalizing
        if v_red != 0:
            for i in range(ntrace1_data, ntrace2_data):
                j = i-ntrace1_data
                t_red = np.abs(x[j])/v_red
                n_red = int(t_red/dt)
                if v_red > 0:
                    d[nsamp1_data:nsamp_exp, i] =\
                        d[nsamp1_data-n_red:nsamp_exp-n_red, i]*sigmas[j]
                else:
                    d[nsamp1_data+n_red:nsamp_exp+n_red, i] =\
                        d[nsamp1_data:nsamp_exp, i]*sigmas[j]

# Final filtered data are the real part of the inverly transformed spectrum
# In addition, extract the block of real data from the expanded block.
        data_filt = np.real(d[nsamp1_data:nsamp2_data,
                              ntrace1_data:ntrace2_data])
        return data_filt, v_red, True

    def airWaveFilter(self):
        """
        Function eliminates all waves traveling with a velocity below 400 m/s
        i.e. air waves and smaller velocities.

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
# For description of v_red and nk_el, see function fk_filt
        v_red = 450.
        nk_el = 3
        self.neg_flag = False
        xx = np.array(self.main.window.x)
        vt = np.zeros_like(np.transpose(self.main.window.v))
        for i in range(np.size(self.main.window.v, 0)):
            vt[:, i] = self.main.window.v[i, :]
# Extract traces recorded in negative direction, if there are more than 5
        for idir in [-1, 1]:
            if idir < 0:
                ntrace_neg = np.size(np.where(xx < 0))
# The array is flipped in the distance-direction in order to have the trace
# nearest to the shot point on the left side
                if ntrace_neg > 5:
                    data_treat = np.flip(vt[:, xx < 0], axis=1)
# Also the x-position vectorei is flipped and the signe is inversed in order to
# simulate recording in positive direction which is needed for function fk_filt
                    x_treat = -np.flip(xx[xx < 0])
# Apply filter
                    d, _, _ = self.FK_filt(x_treat, self.data.dt, data_treat,
                                           v_red, nk_el)
# Flip filtered array back and copy filtered data back into ariginal array v
                    vt[:, xx < 0] = np.flip(d, axis=1)
# Extract traces recorded in positive direction, if there are more than 5
# Here, arrays don't have to be flipped, since they should already be in the
# correct order
            else:
                ntrace_pos = np.size(np.where(xx >= 0))
                if ntrace_pos > 5:
                    data_treat = vt[:, xx >= 0]
                    x_treat = xx[xx >= 0]
                    d, _, _ = self.FK_filt(x_treat, self.data.dt, data_treat,
                                           v_red, nk_el)
                    vt[:, xx >= 0] = d
# redraw filtered data
        self.main.window.v = np.transpose(vt)
        for it, t in enumerate(self.main.window.actual_traces):
            ifile = self.traces.shot[t]
            itrace = self.traces.trace[t]
            if self.traces.amplitudes[t] > 0:
                self.data.st[ifile][itrace].data =\
                    np.float32(self.main.window.v[it, :])
            elif self.traces.amplitudes[t] < 0:
                self.data.st[ifile][itrace].data =\
                    -np.float32(self.main.window.v[it, :])
        self.main.window.drawNew(True)
        self.main.window.setHelp(self.main.window.main_text)
        self.filtered = True

    def velocityFilter(self):
        """
        Allows definition of maximum filtered velocity

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.window.setHelp(self.window.spectrum_text)
        self.main.function = "velocity_filter"
        nk_el = 3
        xx = np.array(self.window.x)
        # vt = np.zeros_like(np.transpose(self.window.v))
        # for i in range(np.size(self.window.v,0)):
        #     vt[:,i] = self.window.v[i,:]
        vt = np.transpose(self.window.v)
# Extract traces recorded in negative direction, if there are more than 5
        ntrace_neg = np.size(np.where(xx < 0))
        ntrace_pos = np.size(np.where(xx >= 0))
        act = True
        results, okButton = self.main.dialog(
            ["Filter negative velocities ? [y/n]"], ["c"], "Velocity filter")
        if okButton is False:
            print("Velocity filter cancelled")
            self.window.drawNew(True)
            self.window.setHelp(self.window.main_text)
            self.main.function = "main"
            return
        self.neg_flag = results[0] > -1
        if ntrace_neg > 5:
            if ntrace_pos > 5:
                if ntrace_neg > ntrace_pos:
                    order = [-1, 1]
                else:
                    order = [1, -1]
            else:
                order = [-1]
        else:
            order = [1]

        for idir in order:
            if idir < 0:
                ntrace_neg = np.size(np.where(xx < 0))
# The array is flipped in the distance-direction in order to have the trace
# nearest to the shot point on the left side
                if ntrace_neg > 5:
                    data_treat = np.flip(vt[:, xx < 0], axis=1)
# Also the x-position vector is flipped and the sign is inversed in order to
# simulate recording in positive direction which is needed for function fk_filt
                    x_treat = -np.flip(xx[xx < 0])
# Apply filter
                    if idir == order[0]:
                        d, self.v_red, act = self.FK_filt(
                            x_treat, self.data.dt, data_treat, self.v_red,
                            nk_el, idir, True)
                    else:
                        if act:
                            d, _, _ = self.FK_filt(
                                x_treat, self.data.dt, data_treat, self.v_red,
                                nk_el, idir, False)
                        else:
                            d, _, _ = self.FK_filt(
                                x_treat, self.data.dt, data_treat, 0., nk_el,
                                idir, False)
# Flip filtered array back and copy filtered data back into ariginal array v
                    vt[:, xx < 0] = np.flip(d, axis=1)
# Extract traces recorded in positive direction, if there are more than 5
# Here, arrays don't have to be flipped, since they should already be in the
# correct order
            else:
                ntrace_pos = np.size(np.where(xx >= 0))
                if ntrace_pos > 5:
                    data_treat = vt[:, xx >= 0]
                    x_treat = xx[xx >= 0]
                    if idir == order[0]:
                        d, self.v_red, act = self.FK_filt(
                            x_treat, self.data.dt, data_treat, self.v_red,
                            nk_el, idir, True)
                    else:
                        if act:
                            d, _, _ = self.FK_filt(
                                x_treat, self.data.dt, data_treat, self.v_red,
                                nk_el, idir, False)
                        else:
                            d, _, _ = self.FK_filt(
                                x_treat, self.data.dt, data_treat, 0., nk_el,
                                idir, False)
                    vt[:, xx >= 0] = d

        all_flag = False
        if act:
            results, okButton = self.main.dialog(
                ["Apply to all shots ? [y/n]"], ["e"], ["n"],
                "Velocity filter")
            if okButton is False:
                print("Velocity filter cancelled")
                self.window.drawNew(True)
                self.window.setHelp(self.window.main_text)
                return
            all_flag = results[0].lower() == "y"
# copy filtered data back into stream st
        self.window.v = np.transpose(vt)
        for it, t in enumerate(self.window.actual_traces):
            ifile = self.traces.file[t]
            itrace = self.traces.trace[t]
            if self.traces.amplitudes[t] > 0:
                self.data.st[ifile][itrace].data = np.float32(
                    self.window.v[it, :])
            else:
                self.data.st[ifile][itrace].data = -np.float32(
                    self.window.v[it, :])
# If all data should be treated, do this now
        if all_flag:
            nsht = len(self.traces.sht_pt_dict)
            print(f'Number of shots to be treated: {nsht}')
            progressBar = QtWidgets.QProgressBar(self.window)
            self.window.mplvl.addWidget(progressBar)
            progressBar.show()
            progressBar.setValue(0)
# Loop over all shot points
            for jsht, isht in enumerate(self.traces.sht_pt_dict):
                ifiles = self.traces.sht_pt_dict[isht]["file"]
                itrac = self.traces.sht_pt_dict[isht]["trace"]
                irec = self.traces.sht_pt_dict[isht]["receiver"]
                itraces = []
                ntr = len(ifiles)
                vv = np.zeros((ntr, self.data.nsamp))
                for it in range(ntr):
                    itr = self.traces.sht_rec_dict[(isht, irec[it])]
                    itraces.append(itr)
                    vv[it, :] = self.data.st[ifiles[it]][itrac[it]].data *\
                        self.traces.amplitudes[itr]
                itraces = np.array(itraces, dtype=int)
                xx = np.array(self.traces.offset[itraces])
                vt = np.zeros_like(np.transpose(vv))
                for i in range(np.size(vv, 0)):
                    fac = np.std(vv[i, :])
                    if fac > 0:
                        vt[:, i] = vv[i, :]/fac
# Extract traces recorded in negative direction, if there are more than 5
                for idir in [-1, 1]:
                    if idir < 0:
                        ntrace_neg = np.size(np.where(xx < 0))
# The array is flipped in the distance-direction in order to have the trace
# nearest to the shot point on the left side
                        if ntrace_neg > 5:
                            data_treat = np.flip(vt[:, xx < 0], axis=1)
# Also the x-position vector is flipped and the sign is inversed in order to
# simulate recording in positive direction which is needed for function fk_filt
                            x_treat = -np.flip(xx[xx < 0])
# Apply filter
                            d, _, _ = self.FK_filt(
                                x_treat, self.data.dt, data_treat, self.v_red,
                                nk_el)
# Flip filtered array back and copy filtered data back into ariginal array v
                            vt[:, xx < 0] = np.flip(d, axis=1)
# Extract traces recorded in positive direction, if there are more than 5
# Here, arrays don't have to be flipped, since they should already be in the
# correct order
                    else:
                        ntrace_pos = np.size(np.where(xx >= 0))
                        if ntrace_pos > 5:
                            data_treat = vt[:, xx >= 0]
                            x_treat = xx[xx >= 0]
                            d, _, _ = self.FK_filt(
                                x_treat, self.data.dt, data_treat, self.v_red,
                                nk_el)
                            vt[:, xx >= 0] = d
# copy filtered data to stream st
                vv = np.transpose(vt)
                for it in range(ntr):
                    if self.traces.amplitudes[itraces[it]] > 0:
                        self.data.st[ifiles[it]][itrac[it]].data =\
                            np.float32(vv[it, :])
                    elif self.traces.amplitudes[itraces[it]] < 0:
                        self.data.st[ifiles[it]][itrac[it]].data =\
                            -np.float32(vv[it, :])
                completed = int((jsht+1)/nsht*100)
                print(f"Shot {jsht+1}/{nsht} filtered; completed: "
                      + "{completed}%")
                progressBar.setValue(completed)
            progressBar.setValue(0)
            self.window.mplvl.removeWidget(progressBar)
            progressBar.close()
            print("\nFiltering done\n")
# redraw filtered data
        self.window.drawNew(True)
        self.window.setHelp(self.window.main_text)
        self.data.filtered = True
        self.main.function = "main"
        return

    def v_nmo(self):
        """
        Function write file with nmo velocities for every tenths cdp in SU
        format.
        For this, it interpolates first the velocities on a grid of 0.5 m
        vertical grid-step and 10x dx_cdp. Then it calculates the two-way
        vertical travel time of a wave in each cell and finally the RMS
        velocity down to the base of every cell.
        The results are written into file v_semblance.txt, which should be
        copied to Linux for use in Seismic Unix.

        Needs:
            scipy.interpolate.LinearNDInterpolator

        Returns
        -------
        None.

        """
        coor = self.mgr.paraDomain.cellCenters().array()
        ymn = np.ceil(min(coor[:, 1]))+0.25
        ymx = float(int(max(coor[:, 1])))-0.25
        grd_x = np.unique(self.main.traces.xcdp)
        grd_y = np.arange(ymn, ymx+0.5, 0.5)
        xg, yg = np.meshgrid(grd_x, grd_y)
        interp = LinearNDInterpolator(coor[:, :2], self.mgr.model.array())
        vg = interp(xg, yg)
        x = grd_x[::10]
        y = np.flip(-grd_y)
        v = np.flip(vg[:, ::10], axis=0)
# Replace nan's in v with nearest value within column
        for i in range(len(x)):
            for j in range(len(y)):
                if not np.isnan(v[j, i]):
                    break
            if j > 0:
                v[:j, i] = v[j, i]
            if j < len(y)-1:
                for k in range(j+1, len(y)):
                    if np.isnan(v[k, i]):
                        v[k, i] = v[k-1, i]
# t contains the vertical two-way traveltime (TWTT) a wave spends in one cell
# The interpolation was done on cells with size 0.5m, therefore the TWTT is
# 2*(0.5/v) = 1/v
        t = 1./v
# t0 contains the vertical TWTT down to the base of each cell
        t0 = t*0.
        t0[0, :] = t[0, :]
        vnmo = v*0.
        vnmo[0, :] = v[0, :]**2*t[0, :]
        for i in range(1, len(y)):
            t0[i, :] = t0[i-1, :]+t[i, :]
            vnmo[i, :] = vnmo[i-1, :]+v[i, :]**2*t[i, :]
        vnmo = np.sqrt(vnmo/t0)
# Write NMO velocities into file with SU format
        with open("v_semblance.txt", "w") as fo:
            fo.write("cdp=1")
            for i in range(1, len(x)):
                fo.write(f",{int(i*10+1)}")
            fo.write(" \\")
            fo.write("\n")
            for i in range(len(x)):
                fo.write(f"tnmo={t0[0,i]:0.6f}")
                for k in range(1, len(y)):
                    fo.write(f",{t0[k,i]:0.6f}")
                fo.write(" \\")
                fo.write("\n")
                fo.write(f"vnmo={vnmo[0,i]:0.1f}")
                for k in range(1, len(y)):
                    fo.write(f",{vnmo[k,i]:0.1f}")
                fo.write(" \\")
                fo.write("\n")

    def inversion(self, code=0):
        """
        Function uses Pygimly for tomography inversion.

        For nicer plots, modify numbering format in drawDataMatrix
        The function is found in file
        Environment_path/Lib/site-packages/pygimli/viewer/mpl
        For me, Environment path is C:/Users/Hermann/anaconda3/
        There go to lines 458 and 462 starting with
        ax.set_xticklabels and
        ax.set_yticklabels
        and change the rounding to 0 ciphers instead of the default 2 ciphers

        Needs:
        pygimli.physics.TravelTimeManager
        matplotlib as mpl
        matplotlib.gridspec.GridSpec
        matplotlib.path.Path
        matplotlib.tri
        mpl_toolkits.axes_grid1.make_axes_locatable
        copy
        matplotlib.pyplot as plt
        colorcet as cc

        Returns
        -------
        None.

        """
        try:
            from pygimli.physics import TravelTimeManager
        except ModuleNotFoundError:
            _ = QtWidgets.QMessageBox.warning(
                None, "Warning", "Pygimli is not installed.\n "
                + "Tomography cannot be executed.\n",
                QtWidgets.QMessageBox.Close, QtWidgets.QMessageBox.Close)
            return

        answer = self.main.test_function()
        if not answer:
            return
        try:
            import pygimli as pg
        except ModuleNotFoundError:
            _ = QtWidgets.QMessageBox.warning(
                None, "Warning", "PyGimli is not installed\n"
                + "Tomography connot be executed\n", QtWidgets.QMessageBox.Ok,
                QtWidgets.QMessageBox.Ok)
            return

        def vel_scale(self, ncol=128, scale="specialP"):
            if "special" in scale:
                cmp = colors.LinearSegmentedColormap.from_list(
                    'velocities',
                    ['violet', 'darkgreen', 'forestgreen', 'cyan', 'blue',
                     'springgreen', 'lime', 'greenyellow', 'yellow', 'yellow',
                     'gold', 'orange', 'orangered', 'red'], N=ncol)
            elif "rainbow" in scale:
                cmp = cc.cm.rainbow4
            else:
                cmp = plt.get_cmap(scale)

            cmp.set_under('pink')
            cmp.set_over('darkred')
            cols = cmp(np.linspace(0, 1, ncol))
            return cmp, cols

# If code == 0, color scale and maximum depth for model plot are calculated
#    automatically
        self.tick_size_mod = 16
        self.tick_size_sec = 12
# Store picks in Gimli format (picks.sgt)
        if code != 67:
            self.traces.saveGimli()
# scheme contains the coordinates of shot and receiver points as  well as the
# measured travel times in seconds.
# scheme("s").array() gives the numbers of the shot points for every pick
# scheme("g").array() gives the numbers of the receiver points for every pick
# scheme("t").array() gives the measured travel times
# scheme("err").array() gives the travel time uncertainties
# scheme.sensors().array() gives a numpy array with 3 colums containing the
#       coordinates of all sensors and shots ordered by position.
            self.scheme = pg.load("picks.sgt", verbose=True)
            self.px = pg.x(self.scheme)
            self.gx = np.asarray([self.px[g] for g in self.scheme.id("g")])
            self.sx = np.asarray([self.px[s] for s in self.scheme.id("s")])
            self.sens = np.array(self.scheme.sensors().array(), dtype=float)
            self.offsets =\
                pg.physics.traveltime.shotReceiverDistances(self.scheme)
            self.zmax = np.round(np.max(self.offsets)*0.333)
            gx_min = self.gx.min()
            gx_max = self.gx.max()
            sx_min = self.sx.min()
            sx_max = self.sx.max()
            self.xax_min = min(gx_min, sx_min)
            self.xax_max = max(gx_max, sx_max)
# Call dialog window for input of a number of inversion control parameters
            results, okButton = self.main.dialog(
                    ["Maximum depth (m, positive down)",
                     "Initial smoothing parameter (<0: optimize)",
                     "Smoothing reduction per iteration",
                     "Smoothing in z direction (0..1):",
                     "Maximum iterations (0 = automatic)",
                     "Initial velocity at surface [m/s]",
                     "Initial velocity at bottom [m/s]",
                     "Minimum allowed velocity [m/s]",
                     "Maximum allowed velocity [m/s]",
                     "\nIf min or max is 0, the corresponding limit\n"
                     + " of the color scale is calculated automatically\n",
                     "Velocity color scale min [m/s]",
                     "Velocity color scale max [m/s]",
                     "Type of color scale:",
                     self.model_colors,
                     "Plot rays on final model"],
                    ["e", "e", "e", "e", "e", "e", "e", "e", "e", "l", "e",
                     "e", "l", "b", "c"],
                    [self.zmax, self.smooth, self.s_fact, self.zSmooth,
                     self.maxiter, self.vmin, self.vmax, self.vmin_limit,
                     self.vmax_limit, None, self.v_scale_min, self.v_scale_max,
                     None, None, -1], "Inversion parameters")

            if not okButton:
                print("\n Inversion cancelled")
                self.main.function = "main"
                return

            self.zmax = float(results[0])
            self.smooth = float(results[1])
            if self.smooth < 0:
                opt_flag = True
                self.smooth *= -1.
            else:
                opt_flag = False
            self.s_fact = float(results[2])
            self.zSmooth = float(results[3])
            self.maxiter = int(results[4])
# Avoid maximum number of iteration being 1, since this would not produce a
#   final model (self.endModel is empty). It seems that pygimli counts the
#   forward calculation of the starting model as first iteration.
            if self.maxiter > 0:
                self.maxiter = max(self.maxiter, 2)
            else:
                self.maxiter = 0
            self.vmin = float(results[5])
            self.vmax = float(results[6])
            self.vmin_limit = float(results[7])
            self.vmax_limit = float(results[8])
            try:
                self.v_scale_min = float(results[10])
                self.v_scale_max = float(results[11])
            except IndexError:
                self.v_scale_min = 150
                self.v_scale_max = 6000
            color_scale = self.model_colors[int(results[13])]
            if "Special" in color_scale:
                if "1500" in color_scale:
                    color_scale = "specialP"
                else:
                    color_scale = "specialS"
            if int(results[14]) > -1:
                self.rays_flag = True
            else:
                self.rays_flag = False

# Initialize PyGimli TravelTime Manager and set control parameters
            self.mgr = TravelTimeManager()
            if opt_flag:
                self.mgr.inv.inv.setOptimizeLambda(True)
# Do tomography inversion. If maxiter == 0, stop iterationautomatically if
# chi2<1 of if error does not decrease by more than 1% (dPhi=0.01) if
# maxiter>0 stop iterations latest after maxiter iterations (but earlier if
# one of the other conditions is fulfilled)
            if self.maxiter < 1:
                if self.vmin_limit == 0 and self.vmax_limit == 0:
                    self.mgr.invert(
                        self.scheme, secNodes=3, paraMaxCellSize=5.0,
                        zWeight=self.zSmooth, vTop=self.vmin,
                        vBottom=self.vmax, verbose=1, paraDepth=self.zmax,
                        dPhi=0.01, lam=self.smooth, lambdaFactor=self.s_fact,
                        maxIter=1000)
                else:
                    self.mgr.invert(
                        self.scheme, secNodes=3, paraMaxCellSize=5.0,
                        zWeight=self.zSmooth, vTop=self.vmin,
                        vBottom=self.vmax, verbose=1, paraDepth=self.zmax,
                        dPhi=0.01, lam=self.smooth,
                        limits=[self.vmin_limit, self.vmax_limit],
                        lambdaFactor=self.s_fact, maxIter=1000)
            else:
                self.mgr.invert(
                    self.scheme, secNodes=3, paraMaxCellSize=5.0,
                    zWeight=self.zSmooth, vTop=self.vmin, vBottom=self.vmax,
                    maxIter=self.maxiter, verbose=1, paraDepth=self.zmax,
                    dPhi=0.01, lam=self.smooth,
                    limits=[self.vmin_limit, self.vmax_limit],
                    lambdaFactor=self.s_fact)
# pass starting model from slowness to velocity
            self.startModel = 1./self.mgr.fop.startModel()
# get final model
            self.endModel_all = self.mgr.model.array()
# Get coverage
            self.mesh_coor_all = self.mgr.paraDomain.cellCenters().array()
# Get path where PyGimly stores its results: date-time/TravelTime/Manager
            self.path = self.mgr.saveResult()
            p = os.path.split(self.path)[:-1]
            self.p_aim = ""
            for pp in p:
                self.p_aim = os.path.join(self.p_aim, pp)
            try:
                self.cover = self.mgr.inv.cov_sum
                mesh_x = self.mesh_coor_all[:, 0]
                mesh_y = self.mesh_coor_all[:, 1]
                mesh_v = self.mesh_coor_all[:, 2]
                mesh_x = mesh_x[self.cover > 0.]
                mesh_y = mesh_y[self.cover > 0.]
                mesh_v = mesh_v[self.cover > 0.]
                self.mesh_coor = np.zeros((len(mesh_x), 3))
                self.mesh_coor[:, 0] = mesh_x
                self.mesh_coor[:, 1] = mesh_y
                self.mesh_coor[:, 2] = mesh_v
                self.cell_rays = np.array(self.mgr.inv.cell_rays)
                self.ncover = np.sum(self.cell_rays, axis=0)
                nchi = self.cell_rays.shape[0]
                with open(os.path.join(self.p_aim, "vel&cover.txt"), "w")\
                        as fo:
                    fo.write("   X     Z     V     L_cum     nRay_cum   "
                             + f"{nchi} x nRay\n")
                    for i, c in enumerate(self.cover):
                        fo.write(f"{self.mesh_coor_all[i,0]:0.3f} "
                                 + f"{self.mesh_coor_all[i,1]:0.3f} "
                                 + f"{self.endModel_all[i]:0.0f}    "
                                 + f"{c:0.2f}       {self.ncover[i]}    "
                                 + f"{' '.join(map(str,self.cell_rays[:,i]))}"
                                 + "\n")
                self.cov_txt = "Coverage (cumulated ray lengths/cell_size) "\
                    + "and rays"
                del mesh_x, mesh_y, mesh_v
                self.endModel = self.endModel_all[self.cover > 0.]
            except:
                self.cover =\
                    self.mgr.coverage()/self.mgr.mesh.cellSizes().array()
                cov_true = self.cover[self.cover > 0.]
                mesh_x = self.mesh_coor_all[:, 0]
                mesh_y = self.mesh_coor_all[:, 1]
                mesh_v = self.mesh_coor_all[:, 2]
                mesh_x = mesh_x[self.cover > 0.]
                mesh_y = mesh_y[self.cover > 0.]
                mesh_v = mesh_v[self.cover > 0.]
                self.mesh_coor = np.zeros((len(mesh_x), 3))
                self.mesh_coor[:, 0] = mesh_x
                self.mesh_coor[:, 1] = mesh_y
                self.mesh_coor[:, 2] = mesh_v
                self.endModel = self.endModel_all[self.cover > 0.]
                with open(os.path.join(self.p_aim, "vel&cover.txt"), "w")\
                        as fo:
                    fo.write("   X     Z     V     cover\n")
                    for i, c in enumerate(cov_true):
                        fo.write(f"{self.mesh_coor[i,0]:0.3f} "
                                 + f"{self.mesh_coor[i,1]:0.3f} "
                                 + f"{self.endModel[i]:0.0f}    {c:0.2f}\n")
                self.cov_txt = "log10(Coverage) and rays"
# pass pick times and calculated  from seconds to miliseconds
            self.v_nmo()
            self.dat = self.mgr.fop.data("t")*1000
            self.calc = self.mgr.inv.response.array()*1000
# Calculate minimum and maximum velocities for automatic color scale
#   First find minimum and maximum velocity of both starting and final model
#         (min_vel, max_vel)
            self.min_end = np.min(self.endModel)
            self.max_end = np.max(self.endModel)
# Now find the values of the 1% and 99% quantiles of both models and calculate
# another option for minimum and maximum of velocity color scale (min_v, max_v)
            self.q1_start = np.quantile(self.startModel, 0.01)
            self.q2_start = np.quantile(self.startModel, 0.99)
            self.q1_end = np.quantile(self.endModel, 0.01)
            self.q2_end = np.quantile(self.endModel, 0.99)
# Set maximum depth for model plotting to zmax (from dialog box)
            self.zmax_plt = self.zmax
# Store Control parameters from dialog box into file "inversion_parameters.txt"
            with open(os.path.join(self.p_aim, "inversion_parameters.txt"),
                      "w") as fo:
                fo.write(f"maximum_depth: {self.zmax:0.1f}\n")
                fo.write(f"initial_smoothing: {self.smooth:0.1f}\n")
                fo.write(f"smoothing_reduction_factor: {self.s_fact:0.3f}\n")
                fo.write(f"Z_smoothing: {self.zSmooth:0.1f}\n")
                if self.maxiter > 0:
                    fo.write(f"maximum_iterations: {self.maxiter}\n")
                fo.write(f"minimum_initial_velocity: {self.vmin:0.1f}\n")
                fo.write(f"maximum_initial_velocity: {self.vmax:0.1f}\n")
                fo.write(f"minimum_allowed_velocity: {self.vmin_limit:0.1f}\n")
                fo.write(f"maximum_allowed_velocity: {self.vmax_limit:0.1f}\n")

            print("\nInversion finished")
# If "C" was pressed, call another dialog box to change color scale
#    and/or maximum plotted depth
        else:
            res, okBut = self.main.dialog(
                ["If min or max is 0, the corresponding limit\n"
                 + " of the color scale is calculated automatically\n",
                 "Velocity color scale min [m/s]",
                 "Velocity color scale max [m/s]",
                 "Type of color scale:", self.model_colors,
                 "Maximum depth [m]", "Plot rays on final model"],
                ["l", "e", "e", "l", "b", "e", "c"],
                [None, self.v_scale_min, self.v_scale_max, None, 0,
                 self.zmax_plt, None], "Change color scale")
            if okBut is False:
                self.main.function = "main"
                return
            self.zmax_plt = float(res[5])
            try:
                self.v_scale_min = float(res[1])
                self.v_scale_max = float(res[2])
            except (IndexError, ValueError):
                self.v_scale_min = 150
                self.v_scale_max = 6000
            color_scale = self.model_colors[int(res[4])]
            if "Special" in color_scale:
                if "1500" in color_scale:
                    color_scale = "specialP"
                else:
                    color_scale = "specialS"
            if int(res[6]) > -1:
                self.rays_flag = True
            else:
                self.rays_flag = False
# If color scales are negative (default in first dialog box), they are set
#    to the quantiles rounded to the next 100 m/s
        self.main.function = "inver"
        ncol = 128
        if "special" in color_scale:
            lin_scale = False
        else:
            lin_scale = True
        if self.v_scale_min <= 0:
            self.v_scale_min = max(self.q1_start, self.q1_end)
            self.v_scale_min = np.round(self.v_scale_min/100, 0)*100
        if self.v_scale_max <= 0:
            self.v_scale_max = max(self.q2_start, self.q2_end)
            self.v_scale_max = np.round(self.v_scale_max/100, 0)*100
        if self.v_scale_max < 1600 and color_scale == "specialP":
            _ = QtWidgets.QMessageBox.warning(
                None, "Warning", "Maximum velocity < 1600\nColor scale "
                + "changed to rainbow\nOnce the plot is finished, you may "
                + "press 'C' to change color scale",
                QtWidgets.QMessageBox.Close)
            lin_scale = True
            color_scale = "rainbow"
# Define color scale and colors for values above and below extreme scale values
        cmp, cols = vel_scale(self, ncol=ncol, scale=color_scale)
        if lin_scale:
            self.levels = np.linspace(self.v_scale_min, self.v_scale_max, 128)
        else:
            if color_scale == "specialP":
                ncyan = int(ncol/16*4)
                self.levels = list(np.linspace(self.v_scale_min, 1500, ncyan))
                self.levels += list(np.linspace(1501, self.v_scale_max,
                                                ncol-ncyan))
            elif color_scale == "specialS":
                ncyan = int(ncol/16*4)
                self.levels = list(np.linspace(self.v_scale_min, 500, ncyan))
                self.levels += list(np.linspace(501, self.v_scale_max,
                                                ncol-ncyan))
            else:
                self.levels = list(np.linspace(self.v_scale_min,
                                               self.v_scale_max, 128))
            self.levels = np.array(self.levels)

# Define grid for different partial figures
        if code == 0:
            self.w_tomo = newWindow("Tomography results")
        else:
            self.figinv.clf()
        self.figinv = self.w_tomo.fig
        plt.tight_layout()
        self.gs = GridSpec(15, 13, figure=self.figinv)
# Axis for final model
        self.ax_mod = self.figinv.add_subplot(self.gs[:6, :])
# Acis for initial model
        self.ax_start = self.figinv.add_subplot(self.gs[8:11, 0:4])
# Axis for ray and coverage plot
        self.ax_rays = self.figinv.add_subplot(self.gs[12:15, 0:4])
# Axis for measured travel time plot
        self.ax_tt = self.figinv.add_subplot(self.gs[8:11, 5:8])
# Axis for difference between measured and calculated travel times
        self.ax_diff = self.figinv.add_subplot(self.gs[12:15, 5:8])
# Axis for average differences of shot gathers and receiver gathers
        self.ax_av_diff = self.figinv.add_subplot(self.gs[12:15, 9:12])
# Axis for chi2-evolution
        self.ax_chi = self.figinv.add_subplot(self.gs[8:11, 9:12])
# Define ticks for horizontal and vertical axes of model plots
# For horizontal axis, use sensor positions from pick file (picks.sgt")
# For vertical axis suppose that the topmost values is zero, which implies that
# the shot and receiver coordinates in files shots.geo and receivers.geo
# have been shifted vertically such that the topmost of both is at 0m.
# The reason for this is that PyGimly does not seem to invert for velocities
# in grid cells at negative Z coordinates, although it seems to create those
# Function Save_Gimli called at the beginning of this function takes care of
# this.
        self.ticks_x = self.window.set_ticks(np.min(self.sens),
                                             np.max(self.sens[:, 0]))
        self.ticks_y = self.window.set_ticks(-self.zmax, 0.)
# Define annotated levels of velocity color scales
        self.levs = np.linspace(self.v_scale_min, self.v_scale_max, 20,
                                endpoint=True)

# Plot starting model
        pg.viewer.showMesh(pg.Mesh(
            self.mgr.paraDomain), data=self.startModel, ax=self.ax_start,
            cMap=cmp, cMin=self.v_scale_min, cMax=self.v_scale_max,
            logScale=False, orientation="vertical", label="Velocity [m/s]",
            fitView=False)
        self.ax_start.set_xticks(self.ticks_x)
        self.ax_start.set_yticks(self.ticks_y)
        self.ax_start.set_xlim(left=self.xax_min, right=self.xax_max)
        self.ax_start.set_xlabel("Distance [m]", fontsize=self.tick_size_sec)
        self.ax_start.set_ylabel("Depth [m]", fontsize=self.tick_size_sec)
        self.ax_start.tick_params(axis='both', labelsize=self.tick_size_sec)
        self.ax_start.set_title("Starting model",
                                fontsize=self.tick_size_sec+2)
        ax_xmin, ax_xmax = self.ax_start.get_xlim()
        ax_ymin, ax_ymax = self.ax_start.get_ylim()
        xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
        ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
        txt = self.ax_start.text(xtxt, ytxt, "B", horizontalalignment="left",
                                 verticalalignment="bottom", fontsize=18)
        txt.set_bbox(dict(facecolor="white"))
# Share x and y axis parameters with ray plot and plot of final model
        self.ax_rays.sharex(self.ax_start)
        self.ax_mod.sharex(self.ax_start)
        self.ax_rays.sharey(self.ax_start)
        self.ax_mod.sharey(self.ax_start)
        print("Starting model plotted")

# Plot coverage and rays of final model
        cov_min = np.min(self.cover[self.cover > -np.inf])
        cov_max = np.max(self.cover[self.cover < np.inf])
        cmp = cc.cm.fire_r
        data = deepcopy(self.cover)
        data[np.isclose(data, 0.)] = np.nan
        pg.viewer.showMesh(pg.Mesh(
            self.mgr.paraDomain), data=data, ax=self.ax_rays, cMap=cmp,
            cMin=cov_min, cMax=cov_max, orientation="vertical",
            label="coverage [m]", fitView=False)
        rays = self.mgr.drawRayPaths(ax=self.ax_rays, color="black", lw=0.3,
                                     alpha=0.5)
        self.ax_rays.set_xticks(self.ticks_x)
        self.ax_rays.set_yticks(self.ticks_y)
        self.ax_rays.set_xlim(left=self.xax_min, right=self.xax_max)
        self.ax_rays.set_xlabel("Distance [m]", fontsize=self.tick_size_sec)
        self.ax_rays.set_ylabel("Depth [m]", fontsize=self.tick_size_sec)
        self.ax_rays.tick_params(axis='both', labelsize=self.tick_size_sec)
        self.ax_rays.set_title(self.cov_txt, fontsize=self.tick_size_sec+2)
        ax_xmin, ax_xmax = self.ax_rays.get_xlim()
        ax_ymin, ax_ymax = self.ax_rays.get_ylim()
        xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
        ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
        txt = self.ax_rays.text(xtxt, ytxt, "C", horizontalalignment="left",
                                verticalalignment="bottom", fontsize=18)
        txt.set_bbox(dict(facecolor="white"))
        print("Rays plotted")
# Interpolate coverage like final model
        self.triang = tri.Triangulation(
            self.mesh_coor_all[:, 0], self.mesh_coor_all[:, 1])
        isbad = np.isclose(self.cover, 0.)
        self.mask = np.all(np.where(isbad[self.triang.triangles], True, False),
                           axis=1)
        self.triang.set_mask(self.mask)
# Plot final model
        gci0 = self.ax_mod.tricontourf(
            self.triang, self.endModel_all, extend='both', levels=self.levels,
            colors=cols)
# self.scheme contains all shot and receiver coordinates as well as measured
# travel times, obtained at the beginning of the function from file picks.sgt
        y = pg.y(self.scheme)
# If there is topography present (not all receivers/shots at z=0), calculate
# the outline of the model (z positions of shots and receivers and the two
# corners of the base of the model) and define this outline as clipping path
        if min(y) < 0:
            x = pg.x(self.scheme)
            x = np.array(list(x)+[max(x), min(x), min(x)])
            y = np.array(list(y)+[-self.zmax_plt, -self.zmax_plt, y[0]])
            clip = np.zeros((len(x), 2))
            clip[:, 0] = x
            clip[:, 1] = y
# clip contains the coordinates of the clipping path. codes will contain the
# way to connect clipping points (move to the first point of the path, draw
# lines to all other points and finally close the path)
            codes = []
            Pat = Path
            codes += [Pat.MOVETO]
            codes += [Pat.LINETO]*(len(x)-2)
            codes += [Pat.CLOSEPOLY]
# set clipping
            clip_path = Path(clip, codes)
            del x, clip
            self.ax_mod.set_clip_on(True)
            for collection in gci0.collections:
                collection.set_clip_path(clip_path,
                                         transform=self.ax_mod.transData)
        del y
# Define positions of ticks along vertical axis such that 4 to 5 numbers are
# plotted
        dtk = round(self.zmax_plt/6., 0)
        ticks_y_mod = self.window.set_ticks(-self.zmax_plt, 0., dtick=dtk)
# Plot color bar
        if self.v_scale_max > 2000:
            ticks_vel = np.array([200, 500, 1000, 1500, 2000, 2500, 3000, 3500,
                                  4000, 4500, 5000, 5500, 6000])
        else:
            ticks_vel = np.array([200, 500, 750, 1000, 1250, 1500, 1750, 2000])
        ticks_vel = ticks_vel[ticks_vel <= self.v_scale_max]
        if ticks_vel[-1] < self.v_scale_max-100:
            ticks_vel = np.array(list(ticks_vel) + [self.v_scale_max])
        self.ax_mod.set_aspect('equal', adjustable='box', anchor='W')
        divider = make_axes_locatable(self.ax_mod)
        cax = divider.append_axes("right", size="2%", pad=0.2)
        cax2 = divider.append_axes("top", size="2%", pad="10%")
        cb = plt.colorbar(
            gci0, cmap=cmp, cax=cax, format='%.0f', label="Velocity [m/s]",
            ticks=ticks_vel, orientation='vertical', aspect=25, shrink=0.9,
            extend='both')

        cb.ax.tick_params(labelsize=14)
        if self.rays_flag:
            _ = self.mgr.drawRayPaths(ax=self.ax_mod, color="black", lw=0.3,
                                      alpha=0.5)
        ticks_x_mod = self.window.set_ticks(self.xax_min, self.xax_max,
                                            ntick=10)
        self.ax_mod.set_xticks(ticks_x_mod)
        self.ax_mod.set_yticks(ticks_y_mod)
        self.ax_mod.set_xlim(left=self.xax_min, right=self.xax_max)
        self.ax_mod.grid(which='minor', axis='both', color="gray")
        self.ax_mod.grid(which='major', axis='both', color="k")
        self.ax_mod.set_xlabel("Distance [m]", fontsize=self.tick_size_mod)
        self.ax_mod.set_ylabel("Depth [m]", fontsize=self.tick_size_mod)
        self.ax_mod.tick_params(axis='both', labelsize=self.tick_size_mod)
        self.ax_mod.set_ylim(-self.zmax_plt, 0)
        self.ax_mod.set_title(
              f"Model velocities (min:{self.min_end:0.0f}, "
              + f"max:{self.max_end:0.0f})", fontsize=self.tick_size_mod+2)
        ax_xmin, ax_xmax = self.ax_mod.get_xlim()
        ax_ymin, ax_ymax = self.ax_mod.get_ylim()
        self.ax_mod.text(
            ax_xmin, ax_ymax+(ax_ymax-ax_ymin)*0.01, self.main.dir_start,
            horizontalalignment="left", verticalalignment="bottom",
            fontsize=18)
        self.ax_mod.text(
            ax_xmax, ax_ymax+(ax_ymax-ax_ymin)*0.01, self.main.dir_end,
            horizontalalignment="right", verticalalignment="bottom",
            fontsize=18)
        xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
        ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
        txt = self.ax_mod.text(xtxt, ytxt, "A", horizontalalignment="left",
                               verticalalignment="bottom", fontsize=18)
        txt.set_bbox(dict(facecolor="white"))
        cax2.text(0.5, 0.5, self.main.title, fontsize=24, fontweight="heavy",
                  ha="center", va="bottom")
        cax2.axis('off')

        print("Final model plotted")

# Store calculated picks into file calc_picks.dat
        with open("calc_picks.dat", "w") as fo:
            spu = np.round(np.array([self.geom.sht_dict[d]["x"]
                                     for d in self.geom.sht_dict]), 2)
            rpu = np.round(np.array([self.geom.rec_dict[d]["x"]
                                     for d in self.geom.rec_dict]), 2)
            shts = list(self.geom.sht_dict.keys())
            recs = list(self.geom.rec_dict.keys())
            for i, p in enumerate(rays.get_paths()):
                ispt = self.mgr.fop.data.id("s")[i]
                irpt = self.mgr.fop.data.id("g")[i]
                cs = np.round(self.scheme.sensors().array()[ispt], 2)[0]
                cr = np.round(self.scheme.sensors().array()[irpt], 2)[0]
                try:
                    ispt = shts[np.where(spu == cs)[0][0]]
                except IndexError:
                    continue
                try:
                    irpt = recs[np.where(rpu == cr)[0][0]]
                except IndexError:
                    continue
                fo.write(f"{ispt+1} {irpt+1} {self.calc[i]:0.4f}\n")
            self.window.PlotCalculatedTimes.setEnabled(True)

# Store average differences between calculated and measured travel times for
# every shot point and every receiver point into file Differences.dat
# In the first column, the number of the shot or receiver is given, not a
#   coordinate.
            diff_abs = self.calc-self.dat
            gxu = np.unique(self.gx)
            sxu = np.unique(self.sx)
            diff_mean_shots = np.zeros(len(sxu))
            diff_mean_recs = np.zeros(len(gxu))
            for i, s in enumerate(sxu):
                diff_mean_shots[i] = np.mean(
                    diff_abs[np.where(self.sx == s)[0]])
            for i, g in enumerate(gxu):
                diff_mean_recs[i] = np.mean(
                    diff_abs[np.where(self.gx == g)[0]])
            with open(os.path.join(self.p_aim, "Differences.dat"), "w") as fo:
                fo.write("Point  shot   receiver\n")
                lsxu = len(sxu)
                lgxu = len(gxu)
                for i in range(max(lsxu, lgxu)):
                    if i < lsxu and i < lgxu:
                        fo.write(f"{i} {diff_mean_shots[i]:0.4f} "
                                 + f"{diff_mean_recs[i]:0.4f}\n")
                    elif i < lsxu:
                        fo.write(f"{i} {diff_mean_shots[i]:0.4f} 0.000\n")
                    else:
                        fo.write(f"{i} 0.0000 {diff_mean_recs[i]:0.4f}\n")

# Store ray paths for final model into file rays.dat
            with open(os.path.join(self.p_aim, "rays.dat"), "w") as fo:
                for i, p in enumerate(rays.get_paths()):
                    ispt = self.mgr.fop.data.id("s")[i]
                    irpt = self.mgr.fop.data.id("g")[i]
                    cs = self.scheme.sensors().array()[ispt]
                    cr = self.scheme.sensors().array()[irpt]
                    fo.write(f"{len(p.vertices)} : ray {i}, shot {ispt} {cs}, "
                             + f"rec {irpt} {cr}     "
                             + f"calc: {self.calc[i]:0.4f} ms, calc-meas: "
                             + f"{diff_abs[i]:0.4f} ms\n")
                    for vert in p.vertices:
                        fo.write(f"{vert[0]:0.3f} {vert[1]:0.3f}\n")
# Plot evolution of chi2
            chi_ev = np.log10(np.array(self.mgr.inv.chi2History))
            nit = len(chi_ev)
            self.ax_chi.plot(np.arange(nit)+1, chi_ev)
            self.ax_chi.set_ylabel("log10(chi2)", fontsize=self.tick_size_sec)
            self.ax_chi.set_xlabel("iteration #", fontsize=self.tick_size_sec)
            self.ax_chi.set_title(f"smoothing: ini: {int(self.smooth)}, "
                                  + f"fac: {self.s_fact:0.2f}, "
                                  + f"z: {self.zSmooth:0.2f}",
                                  fontsize=self.tick_size_sec+2)
            self.ax_chi.tick_params(axis='both', labelsize=self.tick_size_sec)
            ticks_x_chi = self.window.set_ticks(1., nit, ntick=10, dtick=1)
            self.ax_chi.set_xticks(ticks_x_chi)
            ax_xmin, ax_xmax = self.ax_chi.get_xlim()
            ax_ymin, ax_ymax = self.ax_chi.get_ylim()
            self.ax_chi.set_ylim(0., ax_ymax)
            xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
            ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
            txt = self.ax_chi.text(xtxt, ytxt, "F", horizontalalignment="left",
                                   verticalalignment="bottom", fontsize=18)
            txt.set_bbox(dict(facecolor="white"))
#            self.ax_chi.tick_params(axis='both', labelsize=18)

# Plot average traveltime differences for shots and receivers
            self.ax_av_diff.plot(sxu, diff_mean_shots, label="shots")
            self.ax_av_diff.plot(gxu, diff_mean_recs, label="receivers")
            self.ax_av_diff.set_ylabel("Time misfit [ms]",
                                       fontsize=self.tick_size_sec)
            self.ax_av_diff.set_xlabel("Position [m]",
                                       fontsize=self.tick_size_sec)
            self.ax_av_diff.set_title("Average differences calc.-meas. "
                                      + "arrival times",
                                      fontsize=self.tick_size_sec+2)
            self.ax_chi.tick_params(axis='both', labelsize=self.tick_size_sec)
            self.ax_av_diff.grid(which='major', axis='x', color='k')
            self.ax_av_diff.grid(which='major', axis='y', color='k')
            self.ax_av_diff.grid(which='minor', axis='y', color='gray')
            self.ax_av_diff.legend(bbox_to_anchor=(1, 1), loc="upper right")
            ax_xmin, ax_xmax = self.ax_av_diff.get_xlim()
            ax_ymin, ax_ymax = self.ax_av_diff.get_ylim()
            xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
            ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
            txt = self.ax_av_diff.text(xtxt, ytxt, "G",
                                       horizontalalignment="left",
                                       verticalalignment="bottom", fontsize=18)
            txt.set_bbox(dict(facecolor="white"))
#            self.ax_av_diff.tick_params(axis='both', labelsize=18)
            print("Average differences plotted")

# Plot measured travel times using the PyGimli utility
            maxval1 = np.max(self.dat)
            maxval2 = np.max(self.calc)
            maxval = max(maxval1, maxval2)
            minval = 0
            ticks = self.window.set_ticks(minval, maxval, 8, 5.)
            dmin = np.quantile(diff_abs, 0.01)
            dmax = np.quantile(diff_abs, 0.99)
            dext = max(dmax, -dmin)
            dmax = dext
            dmin = -dext
            ticks_d = self.window.set_ticks(dmin, dmax, 8, 0.5)

            gci1 = pg.viewer.mpl.dataview.drawVecMatrix(
                self.ax_tt, self.gx,  self.sx, self.dat, squeeze=True,
                logScale=False, cMin=0, cMax=maxval, aspect='equal',
                fitView=True)
            self.ax_tt.set_aspect('equal', adjustable='box', anchor='C')
            self.ax_tt.set_ylabel("Source positions [m]",
                                  fontsize=self.tick_size_sec)
            self.ax_tt.set_title("Measured travel times",
                                 fontsize=self.tick_size_sec+2)
            self.ax_tt.tick_params(axis='both', labelsize=self.tick_size_sec)
            ax_xmin, ax_xmax = self.ax_tt.get_xlim()
            ax_ymin, ax_ymax = self.ax_tt.get_ylim()
            xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
            ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
            txt = self.ax_tt.text(xtxt, ytxt, "D", horizontalalignment="left",
                                  verticalalignment="bottom", fontsize=18)
            txt.set_bbox(dict(facecolor="white"))
            _ = plt.colorbar(gci1, ax=self.ax_tt, format='%.0f',
                             label="Times (ms)", ticks=ticks,
                             orientation='vertical', aspect=20)
            print("Travel times plotted")

# Plot differences between calculated and measured travel times using the
#   Pygimli utility
            gci3 = pg.viewer.mpl.dataview.drawVecMatrix(
                self.ax_diff, self.gx, self.sx, diff_abs, squeeze=True,
                cMap="seismic", logScale=False, cMin=dmin,
                cMax=dmax, aspect='equal', fitView=True)
            self.ax_diff.set_aspect('equal', adjustable='box', anchor='C')
            self.ax_diff.set_xlabel("Receiver positions [m]",
                                    fontsize=self.tick_size_sec)
            self.ax_diff.set_ylabel("Source positions [m]",
                                    fontsize=self.tick_size_sec)
            self.ax_diff.set_title(
                f"Misfit after {self.mgr.inv.inv.iter()+1} "
                + f"iterations:\nchi2={self.mgr.inv.chi2():0.2f}; "
                + f"abs_rms={self.mgr.inv.inv.absrms()*1000:0.1f}ms",
                fontsize=self.tick_size_sec+2)
            self.ax_diff.tick_params(axis='both', labelsize=self.tick_size_sec)
            ax_xmin, ax_xmax = self.ax_diff.get_xlim()
            ax_ymin, ax_ymax = self.ax_diff.get_ylim()
            xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
            ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
            txt = self.ax_diff.text(xtxt, ytxt, "E",
                                    horizontalalignment="left",
                                    verticalalignment="bottom", fontsize=18)
            txt.set_bbox(dict(facecolor="white"))
            _ = plt.colorbar(gci3, ax=self.ax_diff, format='%.1f',
                             label="Calc - meas (ms)", ticks=ticks_d,
                             orientation='vertical', aspect=20)
            print("Misfits plotted")

# Plot title above the final model and store plot
# If plot is done with automatic scaling, the name is
#    inversion_results-auto.png
#    if not, it is inversion_results.png. So, the automatic scaling is always
#    stored, if scales are changed, only the last version is stored.
#           self.window.figs[ip].suptitle(self.main.title, fontsize="xx-large",
#                                         fontweight="heavy")
            # self.figinv.suptitle(self.main.title, fontsize="xx-large",\
            #                               fontweight="heavy")
#        self.w_tomo.showMaximized()
        self.w_tomo.show()
        if code != 67:
            self.figinv.savefig(os.path.join(self.p_aim,
                                             "inversion_results_auto.png"))
            self.figinv.savefig(os.path.join(self.p_aim,
                                             "inversion_results_auto.png"))
            with open(os.path.join(self.p_aim, "velocities.dat"), "w") as fo:
                for i, e in enumerate(self.endModel):
                    fo.write(f"{self.mesh_coor[i,0]:0.3f} "
                             + f"{self.mesh_coor[i,1]:0.3f} "
                             + f"{e:0.0f}\n")
            self.window.Change_colors.setEnabled(True)
# Interpolate model on regular quadratic grid for use with Sofi2D

            self.prepare_FWI()
# Move all files from folder TravelTimeManager to its base folder, the name of
#      which is date-hour in the format YYYYMMDD-hh.mm
# Then delete folder TravelTimeManager
# The reason is just practical, not to have to search the results deeper
#     in the path tree than necessary
            try:
                get_files = os.listdir(self.path)
                for g in get_files:
                    os.replace(os.path.join(self.path, g),
                               os.path.join(self.p_aim, g))
                os.rmdir(self.path)
            except Exception as error:
                print(f'Error: {error}\nPath {self.path} not found.\n',
                      'Tomography results stay in original folder')

        else:
            self.figinv.savefig(os.path.join(self.p_aim,
                                "inversion_results.png"))
        code = 0
        self.traces.calc_picks = False
        self.traces.readCalcPicks()

    def invCol(self):
        """

        Call function Inversion with code 67 (decimal ASCII code for letter
        "C"). This gives the possibility to change color scale and/or depth
        scale for plotting of inversion results

        This function is called when "Utilities-> Change colors tomo" or "C"
        is pressed.

        Returns
        -------
        None.

        """
        self.main.function = "main"
        self.inversion(code=67)

    def prepare_FWI(self):
        """
        Choice of output format for different Full-waveform inversions (FWI)

        Returns
        -------
        None.

        """
        res, okBut = self.main.dialog(
            ["Save model for FWI inversion with:",
             ["SOFI2D", "ShaVi", "None"]],
            ["l", "r"], [None, 2], "Choose FWI format")
        if not okBut or int(res[1]) == 2:
            print("\nNo FWI output written\n")
            return
        elif int(res[1]) == 0:
            self.prepareSofi2D()
        else:
            self.prepareShaVi()

    def prepareSOFI2D(self):
        """
        Write final model with estimated v_s and densities to binary files for
        use with SOFI2D and prepare a sample json file for SOFI2D

        Returns
        -------
        None.

        """
        x = self.mgr.paraDomain.cellCenters().array()[:, 0]
        z = self.mgr.paraDomain.cellCenters().array()[:, 1]
        v = self.mgr.model.array()
        vmin = v.min()/2.
        vmax = v.max()
# Get some control parameters
        res, okBut = self.main.dialog(
            ["Central frequency [Hz]",
             "Source signal length [number periods]",
             "Calculation time [s]",
             "Absorbing boundary [number of cells]",
             "Absorbing boundary at surface (0: free surface)",
             "Nr. of processors in X direction",
             "Nr. of processors in Z direction",
             "Cell step X for snapshots",
             "Cell step Z for snapshots"],
            ["e", "e", "e", "e", "e", "e", "e", "e", "e"],
            [50, 2, 0.1, 20, 0, 2, 2, 4, 2], "Settings for Sofi2D")
        if not okBut:
            print("\nSOFI2D output not written\n")
            return
# fc is the central frequency of the source signal
# dx is the cell size, calculated following the stability criteria given in
#    SOFI2D manual and rounded to the next lower 2.5 cm.
        fc = float(res[0])
        dx = vmin/(16*fc)
        dx = int(dx/0.025)*0.025
# ts is length of simulated source signal
        ts = float(res[1])/fc
# tmax is length of traces to be calculated
        tmax = float(res[2])
# nbound is the number of cells to be added at the left, right and lower
# edges for wave attenuation and avoiding reflections from those boundaries
# X coordinates of the shot and receiver points are increased by bound
        nbound = int(res[3])
        bound = nbound*dx
# n_surface_bound is the number of cells to be added to the model for wave
# attenuation at the upper surface. If it is 0, free surface is assumed. If it
# is >0, the shot and receiver points are placed at depth bound_s+0.05
        n_surface_bound = int(res[4])
        bound_s = n_surface_bound*dx
# proc_x*proc_z is the number of CPU processors that will be used. The whole
# domain will be split into proc_x blocks in X-direction and proc_z blocks in
# Z-direction
        proc_x = int(res[5])
        proc_z = int(res[6])
# SOFI2D writes snapshot files. In these files, every idx-th point is writte in
# X direction and every idz-th point in Z direction
        idx = int(res[7])
        idz = int(res[8])
        dt = 0.55*dx/vmax
        dt = int(dt*1E6)/1E6
        ndt = max(int(0.001/dt), 1)
        print(f"tmax: {tmax}")
        xmin = np.floor(x.min())-bound
        xmax = np.ceil(x.max())+bound
        zmin = np.floor(z.min())-bound
        zmax = np.ceil(z.max())+bound_s
        print(f"zmin: {zmin:0.3f}, zmax: {zmax:0.3f}")
        nx = int((xmax-xmin)/dx+1)
        nz = int((zmax-zmin)/dx+1)
        ix = proc_x*idx
        iz = proc_z*idz
        if nx % ix != 0:
            nx += ix - nx % ix
        if nz % iz != 0:
            nz += iz - nz % iz

# Define X and Z coordinates of interpolated mesh
        xi = np.linspace(xmin, xmax, nx)
        zi = np.linspace(zmin, zmax, nz)

# Perform linear interpolation of the data given on original positions (x,z)
# on a grid defined by (xi,zi)
        Xi, Zi = np.meshgrid(xi, zi)
        data_p = griddata(self.mgr.paraDomain.cellCenters().array()[:, :2],
                          v, (Xi, Zi), method='linear')
        dpf = data_p.flatten()
        Xif = Xi.flatten()
        Zif = Zi.flatten()
        Xi_nan = Xif[~np.isnan(dpf)]
        Zi_nan = Zif[~np.isnan(dpf)]
        centers = np.zeros((len(Xi_nan), 2))
        centers[:, 0] = Xi_nan
        centers[:, 1] = Zi_nan
        data_p_nan = dpf[~np.isnan(dpf)]
        data_p = griddata(centers, data_p_nan, (Xi, Zi), method='nearest')
        data_p = np.flip(np.transpose(data_p))
# For S-wave velocity suppose that for very small P-velocities, vs = vp/2 and
# for 6 km/s, vs = vp/sqrt(3), in between: linear function
        data_s = data_p/((np.sqrt(3)-2)/6000.*data_p+2.)
# For densities suppose that for very small P-velocities, rho = 2000 and
# for 6 km/s, rho = 2700, in between: linear function
        data_r = data_p*(1000./6000.)+1700.
        print(f"SOFI2D model: nx = {nx}, nz = {nz}, bytes = {nx*nz*4}")
        print(f"              dx: {dx:0.3f}, dt: {dt:0.6f}\n")
# Write inary files for P and S velocities and for densities
        np.array(data_p.astype('float32')).\
            tofile(os.path.join(self.p_aim, "model.vp"))
        np.array(data_s.astype('float32')).\
            tofile(os.path.join(self.p_aim, "model.vs"))
        np.array(data_r.astype('float32')).\
            tofile(os.path.join(self.p_aim, "model.rho"))
        with open(os.path.join(self.p_aim, "sofi2D.json"), "w") as fo:
            fo.write("#---------------------------------------------\n")
            fo.write("#      JSON PARAMETER FILE FOR SOFI2D\n")
            fo.write("#---------------------------------------------\n")
            fo.write("# description:\n")
            fo.write("# description/name of the model:/n#\n{\n")
            fo.write('"Domain Decomposition" : "comment",\n')
            fo.write(f'			"NPROCX" : "{proc_x}",\n')
            fo.write(f'			"NPROCY" : "{proc_z}",\n\n')
            fo.write('"FD stencil" : "comment",\n')
            fo.write('			"RSG" : "0",\n\n')
            fo.write('"FD order" : "comment",\n')
            fo.write('			"FDORDER" : "8",\n')
            fo.write('			"MAXRELERROR" : "1",\n\n')
            fo.write('"2-D Grid" : "comment",\n')
            fo.write(f'			"NX" : "{nx}",\n')
            fo.write(f'			"NY" : "{nz}",\n')
            fo.write(f'			"DH" : "{dx:0.3f}",\n\n')
            fo.write('"Time Stepping" : "comment",\n')
            fo.write(f'			"TIME" : "{tmax:0.3f}",\n')
            fo.write(f'			"DT" : "{dt:0.6f}",\n\n')
            fo.write('"Source" : "comment",\n')
            fo.write('			"SOURCE_SHAPE" : "2",\n')
            fo.write('			"SOURCE_SHAPE values: ricker=1;fumue=2;'
                     + 'from_SIGNAL_FILE=3;SIN**3=4" : "comment",\n')
            fo.write('			"SIGNAL_FILE" : "signal_mseis.tz",\n')
            fo.write('			"SOURCE_TYPE" : "3",\n')
            fo.write('			"SOURCE_TYPE values (point_source): explosive'
                     + '=1; force_in_x=2; force_in_y=3; custom_force=4" : '
                     + '"comment",\n')
            fo.write('			"SRCREC" : "1",\n')
            fo.write('			"SRCREC values :  read from SOURCE_FILE=1,'
                     + ' PLANE_WAVE=2 (internal)" : "comment"\n')
            fo.write('			"SOURCE_FILE" : "./source_pts.dat",\n')
            fo.write('			"RUN_MULTIPLE_SHOTS" : "1",\n')
            fo.write('			"PLANE_WAVE_DEPTH" : "0.0",\n')
            fo.write('			"PLANE_WAVE_ANGLE" : "0.0",\n')
            fo.write(f'			"TS" : "{ts:0.6f}",\n\n')
            fo.write('"Model" : "comment",\n')
            fo.write('			"READMOD" : "1",\n')
            fo.write('			"MFILE" : "model/model",\n')
            fo.write('			"WRITE_MODELFILES" : "0",\n\n')
            fo.write('"Q-approximation" : "comment",\n')
            fo.write('			"L" : "0",\n')
            fo.write('			"FL1" : "5.0",\n')
            fo.write('			"TAU" : "0.00001",\n\n')
            if n_surface_bound > 0:
                free = 0
            else:
                free = 1
            fo.write('"Boundary Conditions" : "comment",\n')
            fo.write(f'			"FREE_SURF" : "{free}",\n')
            fo.write('			"ABS_TYPE" : "1",\n')
            fo.write('			"ABS_TYPE values : CPML-Boundary=1; '
                     + 'Damping-Boundary=2" : "comment",\n')
            fo.write(f'			"FW" : "{nbound}",\n')
            fo.write('			"DAMPING" : "8.0",\n')
            fo.write('			"BOUNDARY" : "0",\n\n')
            fo.write('"Snapshots" : "comment",\n')
            fo.write('			"SNAP" : "1",\n')
            fo.write('			"TSNAP1" : "2e-3",\n')
            fo.write(f'			"TSNAP2" : "{np.round(tmax,3):0.3f}",\n')
            fo.write('			"TSNAPINC" : "0.002",\n')
            fo.write(f'			"IDX" : "{idx}",\n')
            fo.write(f'			"IDY" : "{idz}",\n')
            fo.write('			"SNAP_FORMAT" : "3",\n')
            fo.write('			"SNAP_FILE" : "snap/snapshots",\n\n')
            fo.write('"Receiver" : "comment",\n')
            fo.write('			"SEISMO" : "1",\n')
            fo.write('			"READREC" : "1",\n')
            fo.write('			"REC_FILE" : "receiver_pts.dat",\n')
            fo.write('			"REFRECX, REFRECY" : "0.0 , 0.0",\n')
            fo.write('			"XREC1,YREC1" : "54.0 , 2106.0",\n')
            fo.write('			"XREC2,YREC2" : "5400.0 , 2106.0",\n')
            fo.write('			"NGEOPH" : "1",\n\n')
            fo.write('"Receiver array" : "comment",\n')
            fo.write('			"REC_ARRAY" : "0",\n')
            fo.write('			"REC_ARRAY_DEPTH" : "70.0",\n')
            fo.write('			"REC_ARRAY_DIST" : "40.0",\n')
            fo.write('			"DRX" : "4",\n\n')
            fo.write('"Seismograms" : "comment",\n')
            fo.write(f'			"NDT" : "{ndt}",\n')
            fo.write('			"SEIS_FORMAT" : "1",\n')
            fo.write('			"SEIS_FILE" : "su/record",\n\n')
            fo.write('"Monitoring the simulation" : "comment",\n')
            fo.write('			"LOG_FILE" : "log/test.log",\n')
            fo.write('			"LOG" : "1",\n')
            fo.write('			"OUT_TIMESTEP_INFO" : "100",\n\n')
            fo.write('"Checkpoints" : "comment",\n')
            fo.write('			"CHECKPTREAD" : "0",\n')
            fo.write('			"CHECKPTWRITE" : "0",\n')
            fo.write('			"CHECKPT_FILE" : "tmp/checkpoint_sofi2D",\n}\n')
        nrec = np.unique(self.traces.receiver)
        xrec = np.zeros(len(nrec))
        zrec = np.zeros(len(nrec))
        for i, n in enumerate(nrec):
            if self.geom.x_dir:
                xrec[i] = self.geom.rec_dict[n]["x"]
            else:
                xrec[i] = self.geom.rec_dict[n]["y"]
            zrec[i] = self.geom.rec_dict[n]["z"]
        nshot = np.unique(self.traces.shot)
        xshot = np.zeros(len(nshot))
        zshot = np.zeros(len(nshot))
        for i, n in enumerate(nshot):
            if self.geom.x_dir:
                xshot[i] = self.geom.sht_dict[nrec[i]]["x"]
            else:
                xshot[i] = self.geom.sht_dict[nrec[i]]["y"]
            zshot[i] = self.geom.sht_dict[nrec[i]]["z"]
        xmnr = xrec.min()
        xmns = xshot.min()
        dx_stn = dx*nbound - min(xmnr, xmns)
        xrec += dx_stn
        xshot += dx_stn
        zrec += bound_s+0.05
        zshot += bound_s+0.05
        with open(os.path.join(self.p_aim, "receiver_pts.dat"), "w") as fo:
            for i, xs in enumerate(xrec):
                fo.write(f"{xs:0.3f} {zrec[i]:0.3f}\n")
        with open(os.path.join(self.p_aim, "source_pts.dat"), "w") as fo:
            for i, xs in enumerate(xshot):
                fo.write(f"{xs:0.3f} {zshot[i]:0.3f} "
                         + f"0.0 {fc:0.1f} 1.0\n")

    def prepareShaVi(self):
        """
        Write final model and data to text files for
        use with ShaVi and prepare a sample parameter file

        Returns
        -------
        None.

        """
        def ricker_wavelet(f, size, dt=1):
            """
            Create a Ricker signal starting at -1/f with a total length of size
            samples. Calculation starts at t0 = -1./f, returned as time = 0.

            Parameters
            ----------
            f : float
                Central frequency [Hz].
            size : int
                Total number of samples.
            dt : float, optional
                sampling step [s]. The default is 1.

            Returns
            -------
            t : numpy float array with length size
                time vector [s]
            y : numpy float array with length size
                ricker wavelet values

            """
            t0 = 1./f
            t = np.arange(size)*dt-t0
            y = (1.0 - 2.0*(np.pi**2)*(f**2)*(t**2)) *\
                np.exp(-(np.pi**2)*(f**2)*(t**2))
            return t+t0, y

        folder = "ShaVi"
        if not os.path.isdir(folder):
            os.makedirs(folder)
        os.chdir(folder)
        x = self.mgr.paraDomain.cellCenters().array()[:, 0]
        v = self.mgr.model.array()
# Get some control parameters
        res, okBut = self.main.dialog(
            ["Central frequency [Hz]",
             "Use every n'th shot point",
             "Reduce sampling rate by factor",
             "Maximum time to be stored [s]",
             "Grid spacing [m]",
             "Absorbing boundary [number of cells]",
             "Absorbing boundary at surface (0: free surface)",
             "Nr. of processors",
             "Nr. of iterations",
             "Objective function:",
             ["X-correlation", "Norm"]],
            ["e", "e", "e", "e", "e", "e", "e", "e", "e", "l", "r"],
            [50, 1, 4, 0.2, 1, 8, 3, 4, 50, None, 2], "Settings for ShaVi")
        if not okBut:
            print("\nShaVi output not written\n")
            return
# fc is the central frequency of the source signal
# dx is the cell size, calculated following the stability criteria given in
#    SOFI2D manual and rounded to the next lower 2.5 cm.
        fc = float(res[0])
        nd_shot = int(res[1])
        nd_time = int(res[2])
        s_inter = self.data.dt*nd_time
# tmax is length of traces to be calculated
        tmax = float(res[3])
        n_time_min = abs(int(self.data.t0/self.data.dt))
        n_time = int(tmax/self.data.dt)
        n_time_max = n_time_min+n_time
        time_samples = list(range(n_time_min, n_time_max+1, nd_time))
        nt_store = len(time_samples)
        dx = float(res[4])
# nbound is the number of cells to be added at the left, right and lower
# edges for wave attenuation and avoiding reflections from those boundaries
# X coordinates of the shot and receiver points are increased by bound
        nbound = int(res[5])
# n_surface_bound is the number of cells to be added to the model for wave
# attenuation at the upper surface. If it is 0, free surface is assumed. If it
# is >0, the shot and receiver points are placed at depth n_surface_bound
        n_surface_bound = int(res[6])
# proc is the number of CPU processors that will be used.
        proc = int(res[7])
        max_iter = int(res[8])
        object_type = int(res[10])
        x0m = dx*(nbound+5)
        depth_model = self.zmax_plt+x0m
        nz_size = int(depth_model/dx)
        x_min = x.min() - dx*(nbound+5)
        x_max = x.max() + dx*(nbound+5)
        nx_size = int((x_max-x_min)/dx)
# Write data file and store shot positions into list sh_pos and receiver
# positions into dictionary rec_pos_dir with key = conscurive number of
# stored shot and as value a list of receiver positions for each shot.
        sht_pts = list(self.traces.sht_pt_dict.keys())
        rec_pos_dir = {}
        sh_pos = []
        data = np.zeros((len(time_samples),
                         len(self.traces.sht_pt_dict[sht_pts[0]]["file"])))
        with open("data.txt", "w") as fo:
            for i, key in enumerate(sht_pts):
                if key % nd_shot != 0:
                    continue
                r = self.traces.sht_pt_dict[key]["receiver"][0]
                tr = self.traces.sht_rec_dict[(key, r)]
                sh_pos.append(self.traces.shot_pos[tr])
                rec_pos_dir[i] = []
                f_nr = []
                tf_nr = []
                r_nr = []
                for k, f in enumerate(self.traces.sht_pt_dict[key]["file"]):
                    r = self.traces.sht_pt_dict[key]["receiver"][k]
                    tr = self.traces.sht_rec_dict[(key, r)]
                    if np.isclose(self.traces.offset[tr], 0.):
                        continue
                    f = self.traces.sht_pt_dict[key]["file"][k]
                    t = self.traces.sht_pt_dict[key]["trace"][k]
                    f_nr.append(f)
                    tf_nr.append(t)
                    r_nr.append(r)
                    rec_pos_dir[i].append(self.traces.receiver_pos[tr])
                    data[:, k] = self.data.st[f][t].data[time_samples]
                index = np.argsort(np.array(r_nr))
                for j in range(data.shape[0]):
                    for k in index:
                        fo.write(f"{data[j,k]:0.4e}\n")

# Write sources and receivers per source into their files in terms of
# grid point numbers
        with open("sources.txt", "w") as fo:
            for s in sh_pos:
                fo.write(f"{int((s-x_min)/dx)}\n")
        with open("receivers.txt", "w") as fo:
            for i in range(len(rec_pos_dir[0])):
                for key in rec_pos_dir.keys():
                    ri = int((rec_pos_dir[key][i]-x_min)/dx)
                    fo.write(f"{ri} ")
                fo.write("\n")

# Write parameter file
        with open("parameter.txt", "w") as fo:
            fo.write(f"number_sample  {nt_store}\n")
            fo.write(f"sampling_intrvl  {s_inter:0.6f}\n")
            fo.write(f"depth_model  {nz_size}\n")
            fo.write(f"lateral_model  {nx_size}\n")
            fo.write(f"grid_space  {dx:0.2f}\n")
            fo.write(f"depth_source  {n_surface_bound}\n")
            fo.write(f"depth_receiver  {n_surface_bound}\n")
            fo.write(f"number_src  {len(sh_pos)}\n")
            fo.write(f"number_rcr  {len(rec_pos_dir[0])}\n")
            fo.write(f"processors  {proc}\n")
            fo.write(f"max_iteration  {max_iter}\n")
            fo.write("boundary_key  1\n")
            fo.write(f"absorb_lay  {nbound}\n")
            fo.write(f"obj_fun  {object_type}\n")
            fo.write(f"upper_avoid  {n_surface_bound}\n")
            fo.write("data_order  1\n")

# Export tomography model as starting model
        xmin = x_min
        xmax = x_max
        zmin = -3.*dx
        zmax = zmin + nz_size*dx
        print(f"zmin: {zmin:0.3f}, zmax: {zmax:0.3f}")
        nx = nx_size
        nz = nz_size

# Define X and Z coordinates of interpolated mesh
        xi = np.linspace(xmin, xmax, nx)
        zi = np.linspace(zmin, zmax, nz)

# Perform linear interpolation of the data given on original positions (x,z)
# on a grid defined by (xi,zi)
        Xi, Zi = np.meshgrid(xi, zi)
        points = self.mgr.paraDomain.cellCenters().array()[:, :2]
        points[:, 1] *= -1.
        index = np.isfinite(v)
        points = points[index, :]
        v = v[index]
        data_p = griddata(points, v, (Xi, Zi), method='linear')
# Fill nan's with neighbouring values
# First fill rows
        for i in range(data_p.shape[0]):
            ok = np.where(np.isfinite(data_p[i, :]))[0]
# If no valid data exist in the row, skip row
            if len(ok) == 0:
                continue
# If nans exist at the beginning of the row, fill them with the first valid
#    value encountered in the row
            if ok[0] > 0:
                data_p[i, :ok[0]] = data_p[i, ok[0]]
# If nans exist at the end of the row, fill them with the last valid
#    value encountered in the row
            if ok[-1] < nx_size-1:
                data_p[i, ok[-1]+1:] = data_p[i, ok[-1]]
# Fill intermediate nans with value of the left neighbour
            na = np.where(np.isnan(data_p[i, :]))[0]
            if len(na) == 0:
                continue
            for j in na:
                data_p[i, j] = data_p[i, j-1]
# Now do the same procedure for columns
        for i in range(data_p.shape[1]):
            ok = np.where(np.isfinite(data_p[:, i]))[0]
            if len(ok) == 0:
                continue
            if ok[0] > 0:
                data_p[:ok[0], i] = data_p[ok[0], i]
            if ok[-1] < nz_size-1:
                data_p[ok[-1]+1:, i] = data_p[ok[-1], i]
            na = np.where(np.isnan(data_p[:, i]))[0]
            if len(na) == 0:
                continue
            for j in na:
                data_p[j, i] = data_p[j-1, i]
        with open("initial_model.txt", "w") as fo:
            for i in range(data_p.shape[0]):
                for j in range(data_p.shape[1]):
                    fo.write(f"{data_p[i,j]:0.3f} ")
                fo.write("\n")
        w_shavi = newWindow("Shavi interpolated")
        figsha = w_shavi.fig
        _ = figsha.add_subplot(111)
        plt.imshow(data_p, extent=(Xi.min()-dx/2, Xi.max()+dx/2,
                                   -Zi.max()-dx/2, -Zi.min()+dx/2))
        plt.plot(points[:, 0], -points[:, 1], "k.", ms=1)
        plt.show()

        print(f"ShaVi model: nx = {nx}, nz = {nz}, bytes = {nx*nz*4}")
        print(f"             dx: {dx:0.3f}, dt: {s_inter:0.6f}\n")
# Create source signal
        times, data = ricker_wavelet(fc, len(time_samples), dt=s_inter)
        with open("source_signal.txt", "w") as fo:
            for d in data:
                fo.write(f"{d:0.6f}\n")
        os.chdir("..")

    def atten_amp(self):
        """
        Calculates spatial attenuation of waves using amplitude evolution

        Function searches for each trace the maximum of the envelopes of data
        plotted on the screen (you may use muting functions to focus surgically
        on certain phases). Amplitudes are multiplied by the absolute offset to
        counteract geometric spreading. Then, for each side of a shot point,
        an exponential function is fitted to the amplitude evolution, if at
        least 4 traces are available. If more than 6 traces exist on the
        corresponding side, two independent lines are fitted whose results may
        be interpreted  as attenuation near the surface and deeper down.
        A plot is presented with the fitted logarithm of the amplitudes and the
        amplitude fit itself. Instead of slope, a Q value is indicated
        (-1/slope) as well as a r2 value for the ensemble of the two lines.

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "attenuation"
        v_all = self.window.v.copy()
        x_all = self.window.x.copy()
        t_all = self.window.time.copy()
# Extract only data after trigger and exclude data for geophone at short
# distance since usually, signal is saturated
        data = v_all[:, t_all >= 0]
        t = t_all[t_all >= 0]
        data = data[np.abs(x_all) > 1.0001, :]
        x = x_all[np.abs(x_all) > 1.0001]
# Elimiate traces without data
        sdev = np.std(data, axis=1)
        data = data[sdev > 0., :]
        x = x[sdev > 0.]
# recover amplitude reduction due to geometric spreading
        data = data*t
# nx_pos and nx_neg are the indices of traces in positive and negative
#    direction respectively
        nx_pos = np.where(x > 0)[0]
        nx_neg = np.where(x < 0)[0]
# Prepare arrays for calculated amplitudes and logarithms
        amp_calc = np.zeros(len(x))
        amp_calc[:] = np.nan
        lamp_calc = np.zeros(len(x))
        lamp_calc[:] = np.nan
# Calculate envelopes of each trace and find their maxima
        for i in range(data.shape[0]):
            data[i, :] = obs_filter.envelope(data[i, :])
        amp_max = np.max(abs(data), axis=1).squeeze()
        lamp = np.log(amp_max)
# plt_txt will be used for plot title
        plt_txt = ""
# Prepare arrays for parameters of two lines on each side
        intercept_neg = np.zeros(2)
        intercept_pos = np.zeros(2)
        slope_neg = np.zeros(2)
        slope_pos = np.zeros(2)
        q_factor_neg = np.zeros(2)
        q_factor_pos = np.zeros(2)
# If there are more than 3 traces on the negative side calculate attenuation
# in negative direction
# If less than 7 traces exist, fit one single line
# Use scikit for linear fitting
        if len(nx_neg) > 3:
            if len(nx_neg) < 7:
                model = LinearRegression().fit(-x[nx_neg].reshape(-1, 1),
                                               lamp[nx_neg].reshape(-1, 1))
                intercept_neg[:] = model.intercept_
                slope_neg[:] = model.coef_.squeeze()
                q_factor_neg = -1./slope_neg
                r2_neg = model.score(-x[nx_neg].reshape(-1, 1),
                                     lamp[nx_neg].reshape(-1, 1))
                amp_calc[nx_neg] = np.exp(-x[nx_neg]*slope_neg[0]
                                          + intercept_neg[0])
                lamp_calc[nx_neg] = intercept_neg[0]-x[nx_neg]*slope_neg[0]
                plt_txt = f" neg: Q={q_factor_neg[0]:0.2f}, R2 = {r2_neg:0.3f}"
# If more than 6 traces exist, fit two line
# Use function bestLines from refraPlot.py
            else:
                _, slopes, inters, _, _, y_data, r2_neg, _ =\
                     self.window.bestLines(-x[nx_neg], lamp[nx_neg],
                                           origin=False, refra=False)
                intercept_neg = inters.copy()
                slope_neg = slopes.copy()
                q_factor_neg = -1./slope_neg
                print("negative:", slope_neg[0], intercept_neg[0], r2_neg)
                amp_calc[nx_neg] = np.exp(y_data)
                lamp_calc[nx_neg] = y_data.copy()
                plt_txt = f" neg: Q=[{q_factor_neg[0]:0.2f}, "\
                    + "{q_factor_neg[1]:0.2f}], R2 = {r2_neg:0.3f}"
# If there are more than 3 traces on the positive side calculate attenuation
# in positive direction
# If less than 7 traces exist, fit one single line
# Use scikit for linear fitting
        if len(nx_pos) > 3:
            if len(nx_pos) < 7:
                model = LinearRegression().fit(x[nx_pos].reshape(-1, 1),
                                               lamp[nx_pos].reshape(-1, 1))
                intercept_pos[:] = model.intercept_
                slope_pos[:] = model.coef_.squeeze()
                q_factor_pos = -1./slope_pos
                r2_pos = model.score(x[nx_pos].reshape(-1, 1),
                                     lamp[nx_pos].reshape(-1, 1))
                print("positive:", slope_pos[0], intercept_pos[0], r2_pos)
                amp_calc[nx_pos] = np.exp(x[nx_pos] * slope_pos[0]
                                          + intercept_pos[0])
                lamp_calc[nx_pos] = intercept_pos[0]+x[nx_pos]*slope_pos[0]
                plt_txt += f" pos: Q={q_factor_pos[0]:0.2f}, "\
                    + "R2 = {r2_pos:0.3f}"
# If more than 6 traces exist, fit two line
# Use function bestLines from refraPlot.py
            else:
                _, slopes, inters, _, _, y_data, r2_pos, _ =\
                     self.window.bestLines(x[nx_pos], lamp[nx_pos],
                                           origin=False, refra=False)
                intercept_pos = inters.copy()
                slope_pos = slopes.copy()
                q_factor_pos = -1./slope_pos
                print("positive:", slope_pos, intercept_pos, r2_pos)
                amp_calc[nx_pos] = np.exp(y_data)
                lamp_calc[nx_pos] = y_data.copy()
                plt_txt += f" pos: Q=[{q_factor_pos[0]:0.2f}, "\
                    + "{q_factor_pos[1]:0.2f}], R2 = {r2_pos:0.3f}"
# Plot results to screen
#        self.window.drawNew(False)
#        self.figatt =  self.window.figs[self.window.fig_plotted]
#        self.figatt.clear()
#        self.ax_att, self.ax_amp = self.figatt.subplots(2,1)
        self.w_amp = newWindow("Amplitudes and attenuation")
        self.ax_att, self.ax_amp = self.w_amp.fig.subplots(2, 1)
        self.ax_att.plot(x, lamp)
        self.ax_att.plot(x, lamp_calc, "r--")
        self.ax_att.set_xlabel("Offset [m]", fontsize=18)
        self.ax_att.set_ylabel("log(Max amplitude) [n.u.]", fontsize=18)
        self.ax_amp.plot(x, amp_max)
        self.ax_amp.plot(x, amp_calc, "r--")
        self.ax_amp.set_xlabel("Offset [m]", fontsize=18)
        self.ax_amp.set_ylabel("Max amplitude [n.u.]", fontsize=18)
        if self.window.fg_flag:
            text = f"file {self.files.file_numbers[self.window.fig_plotted]}:"\
                + " {plt_txt}"
            self.ax_att.set_title(f"Attenuation, {text}", fontsize=20)
        elif self.window.sg_flag:
            text = f"shot {self.traces.shot[self.window.actual_traces[0]]+1}:"\
                + " {plt_txt}"
            self.ax_att.set_title(f"Attenuation, {text}", fontsize=20)
        elif self.window.rg_flag:
            text = "receiver "\
                + f"{self.traces.receiver[self.window.actual_traces[0]]+1}: "\
                + f"{plt_txt}"
            self.ax_att.set_title(f"Attenuation, {text}", fontsize=20)
        self.ax_att.tick_params(axis='both', labelsize=18)
        self.ax_amp.tick_params(axis='both', labelsize=18)
        self.w_amp.show()
# Write results to file Q.dat. If this file exists, append a line
        if exists("Q.dat"):
            mode = "a"
        else:
            mode = "w"
        with open("Q.dat", mode) as fo:
            fo.write(f"{text}\n")
        self.main.function = "main"

    def pseudo(self):
        """
        Created on Wed Jun 28 08:23:09 2023

        @author: Hermann Zeyen

        Read pick file and plot average velocities for each pick and local
        slownesses between picks
        Needs:
        matplotlib.patches.Rectangle
        matplotlib.colors
        matplotlib.gridspec.GridSpec

        """
        answer = self.main.test_function()
        if not answer:
            return

        # Define input parameters
        # Folder where to find picks and geometry files
        # If vsp=True, z-coordinates are used for midpoint calculation, else
        #              between x or y the one having the largest extent
        res, okBut = self.main.dialog(
            ["Check if VSP configuration:\n  Z is main coordinate"],
            ["c"], [None], "Check if VSP configuration")
        if not okBut:
            print("\nPeudo velocity plot abandoned\n")
            return
        vsp = False
        if res[0] > -1:
            vsp = True
        title = self.main.title
        # Distance between cdp points
        dmid = abs(self.traces.xcdp[1]-self.traces.xcdp[0])
        # Distance between geophones or shot points (the smallest of both)
        if vsp:
            doff = self.geom.dz_geo
        else:
            doff = self.geom.dx_geo

        rec_dict = self.geom.rec_dict
        sht_dict = self.geom.sht_dict
        if vsp:
            direction = "z"
        else:
            direction = "x"

        n_sht = len(sht_dict)
        n_rec = len(rec_dict)
        off_theo = np.zeros((n_sht, n_rec))
        mid_theo = np.zeros((n_sht, n_rec))
        for i, r in enumerate(rec_dict):
            for j, s in enumerate(sht_dict):
                dx = rec_dict[r]["x"] - sht_dict[s]["x"]
                dy = rec_dict[r]["y"] - sht_dict[s]["y"]
                dz = rec_dict[r]["z"] - sht_dict[s]["z"]
                off_theo[j, i] = np.sqrt(dx*dx+dy*dy+dz*dz)
                off_theo[j, i] = np.round(off_theo[j, i]/doff, 0)*doff
                mid_theo[j, i] = (rec_dict[r][direction] +
                                  sht_dict[s][direction])*0.5
                mid_theo[j, i] = np.round(mid_theo[j, i]/dmid, 0)*dmid
        off_flat = off_theo.flatten()
        mid_flat = mid_theo.flatten()
        offsets = np.unique(off_flat)

# For each offset calculate the distance between midpoints. E.g., if shotpoints
#     are located beside every second receiver, the paired shot points have
#     distance between midpoints = 2*receiver_distance, whereas for odd shot
#     points, the distance is equal to the one of the receivers
        dx_mid = np.zeros(len(offsets))
        for i, o in enumerate(offsets):
            mid = np.unique(mid_flat[off_flat == 0.])
            if len(mid > 1):
                dx_mid[i] = abs(mid[1]-mid[0])
            else:
                dx_mid[i] = doff

        # Read pick file. pick_min_time and pick_max_time are not used
        with open("picks.dat", "r") as fh:
            lines = fh.readlines()
        pick_sht = []
        pick_rec = []
        pick_time = []
        for lin in lines:
            line = lin.split("\t")
            if len(line) < 2:
                line = lin.split(" ")
            pick_sht.append(int(line[0])-1)
            pick_rec.append(int(line[1])-1)
            pick_time.append(float(line[2]))
        pick_sht = np.array(pick_sht, dtype=int)
        pick_rec = np.array(pick_rec, dtype=int)
        pick_time = np.array(pick_time, dtype=float)

# For each pick, calculate offset between shot and receiver and the
# midpoint position
        pick_offset = np.zeros_like(pick_time)
        pick_off_rd = np.zeros_like(pick_time)
        pick_midpoint = np.zeros_like(pick_time)
        pick_vel = np.zeros_like(pick_time)
        for i, p in enumerate(pick_sht):
            dx = rec_dict[pick_rec[i]]["x"] - sht_dict[p]["x"]
            dy = rec_dict[pick_rec[i]]["y"] - sht_dict[p]["y"]
            dz = rec_dict[pick_rec[i]]["z"] - sht_dict[p]["z"]
            pick_offset[i] = np.sqrt(dx*dx+dy*dy+dz*dz)
            pick_midpoint[i] = (rec_dict[pick_rec[i]][direction]
                                + sht_dict[p][direction])/2
# Round offsets and midpoints to the values given at the beginning for doff
# and dmid
            pick_off_rd[i] = round(pick_offset[i]/doff, 0)*doff
            pick_midpoint[i] = round(pick_midpoint[i]/dmid, 0)*dmid
# Calculate average velocities for each pick. if time <= 0, exclude from
# calculation
            if pick_time[i] <= 0:
                pick_vel[i] = np.nan
            else:
                pick_vel[i] = pick_offset[i]/pick_time[i]

# Set color scale
        vmin = np.nanquantile(pick_vel, 0.01)
        vmax = np.nanquantile(pick_vel, 0.99)
        cmap = plt.get_cmap("rainbow")
        norm = colors.Normalize(vmin, vmax)
        smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

# Start plot
        self.w_pseudo = newWindow("Pseudo-velocities")
        fig_ps = self.w_pseudo.fig
        fig_ps.set_figwidth(15)
        fig_ps.set_figheight(13)
        gs = GridSpec(15, 13, figure=fig_ps)
        plt.tight_layout()
# Axis for pseudo velocities
        ax_pv = fig_ps.add_subplot(gs[:7, :])
# Axis for local slownesses
        ax_ls = fig_ps.add_subplot(gs[9:, :])
# For every offset and every pick define a rectangle with width dx_min and
#     height doff which is plotted with a color corresponding to the velocity
# For the plot, the offsets are divided by 3 so that the Y axis gives an
#     approximate idea of the depth
        ax = ax_pv
        for i, o in enumerate(offsets):
            ipk = np.where(np.isclose(pick_off_rd, o))[0]
            for j in ipk:
                ax.add_patch(Rectangle(
                    (pick_midpoint[j]-dx_mid[i]*0.5, (o-doff*0.5)/3.),
                    dx_mid[i], doff, color=cmap(norm(pick_vel[j]))))
        ax.set_xlim([pick_midpoint.min(), pick_midpoint.max()])
        ax.set_ylim([pick_offset.max()/3., pick_offset.min()/3.])
        ax.set_title(title+" average velocities", fontsize=18)
        ax.set_xlabel("midpoint position [m]", fontsize=14)
        ax.set_ylabel("offset/3 [m]", fontsize=14)
        ax.tick_params(labelsize=14)
        cb = fig_ps.colorbar(smap, ax=ax)
        cb.set_label(label="average velocity [m/s]", size=14)
        cb.ax.tick_params(labelsize=14)
        ax_xmin, ax_xmax = ax.get_xlim()
        ax_ymin, ax_ymax = ax.get_ylim()
        ax.text(ax_xmin, ax_ymax+(ax_ymax-ax_ymin)*0.01, self.main.dir_start,
                horizontalalignment="left",
                verticalalignment="bottom", fontsize=18)
        ax.text(ax_xmax, ax_ymax+(ax_ymax-ax_ymin)*0.01, self.main.dir_end,
                horizontalalignment="right",
                verticalalignment="bottom", fontsize=18)
        xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
        ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
        _ = ax.text(xtxt, ytxt, "A", horizontalalignment="left",
                    verticalalignment="bottom", fontsize=18)

# Plot local slowness
        ax2 = ax_ls
        shots = np.unique(pick_sht)
        slow = []
# Do slowness calculation like centered finite differences time and offsets
#    between points i+1 and i-1
        for i, s in enumerate(shots):
            off = pick_offset[pick_sht == s]
            mid = pick_midpoint[pick_sht == s]
            t = pick_time[pick_sht == s]
            for j in range(1, len(off)-1):
                ddx = off[j+1]-off[j-1]
                if np.isclose(ddx, 0.):
                    slow.append(np.nan)
                else:
                    slow.append((t[j+1]-t[j-1])/ddx*1000.)
                    if slow[-1] <= 0:
                        slow[-1] = np.nan
        slow = np.array(slow, dtype=float)
        smin = max(np.nanquantile(slow, 0.01), 0.)
        smax = np.nanquantile(slow, 0.99)
        cmap2 = plt.get_cmap("rainbow_r")
        norm2 = colors.LogNorm(smin, smax)
        # norm2 = colors.Normalize(smin,smax)
        sm2 = plt.cm.ScalarMappable(cmap=cmap2, norm=norm2)
        n = -1
# Loop over shot points. For every shotpoint search picks
        for i, s in enumerate(shots):
            off = pick_offset[pick_sht == s]
            mid = pick_midpoint[pick_sht == s]
            t = pick_time[pick_sht == s]
# For every offset and every pick define a rectangle with width dx_min and
#     height doff which is plotted with a color corresponding to the velocity
# For the plot, the offsets are divided by 3 so that the Y axis gives an
#     approximate idea of the depth
            for j in range(1, len(off)-1):
                n += 1
                s = slow[n]
                if np.isnan(s):
                    continue
                xx = mid[j]-doff*0.5
                yy = ((off[j+1]+off[j-1]-doff)*0.5)/3.
                ax2.add_patch(Rectangle((xx, yy), doff, doff,
                                        color=cmap2(norm2(s))))
        ax2.set_xlim([pick_midpoint.min(), pick_midpoint.max()])
        ax2.set_ylim([pick_offset.max()/3., pick_offset.min()/3.])
        ax2.set_title(title+" local slownesses", fontsize=18)
        ax2.set_xlabel("midpoint position [m]", fontsize=14)
        ax2.set_ylabel("offset/3 [m]", fontsize=14)
        ax2.tick_params(labelsize=14)
        cb2 = fig_ps.colorbar(sm2, ax=ax2, ticks=[0.2, 0.5, 1., 2., 5.],
                              format="%.1f")
        cb2.set_label(label="log10(local slowness [ms/m])", size=14)
        cb2.ax.tick_params(labelsize=14)
        ax_xmin, ax_xmax = ax2.get_xlim()
        ax_ymin, ax_ymax = ax2.get_ylim()
        ax2.text(ax_xmin, ax_ymax+(ax_ymax-ax_ymin)*0.01, self.main.dir_start,
                 horizontalalignment="left",
                 verticalalignment="bottom", fontsize=18)
        ax2.text(ax_xmax, ax_ymax+(ax_ymax-ax_ymin)*0.01, self.main.dir_end,
                 horizontalalignment="right",
                 verticalalignment="bottom", fontsize=18)
        xtxt = ax_xmin+(ax_xmax-ax_xmin)*0.02
        ytxt = ax_ymin+(ax_ymax-ax_ymin)*0.05
        _ = ax2.text(xtxt, ytxt, "B", horizontalalignment="left",
                     verticalalignment="bottom", fontsize=18)

        self.w_pseudo.show()
        # Store figure into png file
        fig_ps.savefig("pseudo_section_slowness.png")

#    atten_FFT_backup(self):

#         dt = t[1]-t[0]
#         f = np.fft.fftfreq(ndat, d=dt)
# #        omega = f*2.*np.pi
#         F = np.log(np.abs(np.fft.fft(data,axis=0))+np.nextafter(0, 1))
#         fmax = 200
#         nfmax = np.where(f>=fmax)[0][0]
#         intercept = np.zeros(nfmax)
#         slope = np.zeros(nfmax)
#         r2 = np.zeros(nfmax)
#         q_factor = np.zeros(nfmax)
#         for i in range(1,nfmax):
#             model = LinearRegression().fit(abs(x).reshape(-1, 1),
#                                            F[:,i].reshape(-1, 1))
#             intercept[i] = model.intercept_
#             slope[i] = model.coef_
#             r2[i] = model.score(abs(x).reshape(-1, 1), F[:,i].reshape(-1, 1))
#         q_factor[1:] = -1./(slope[1:]+np.nextafter(0, 1))
#         self.window.drawNew(False)
#         self.figatt =  self.window.figs[self.window.fig_plotted]
#         self.figatt.clear()
#         # self.window.drawNew(False)
#         # fig = self.window.figs[self.window.fig_plotted]
#         self.ax_att, self.ax_r2 = self.figatt.subplots(2,1)
#         self.ax_att.plot(f[1:nfmax],q_factor[1:nfmax])
#         self.ax_att.set_xlabel("Frequency [Hz]")
#         self.ax_att.set_ylabel("Quality factor Q")
#         if self.window.fg_flag:
#             self.ax_att.set_title(
#                 "Attenuation, "
#                 + f"file {self.file.file_numbers[self.fig_plotted]}")
#         elif self.window.sg_flag:
#             self.ax_att.set_title(
#                 "Attenuation, shot "
#                 + f"{self.traces.shot[self.window.actual_traces[0]]+1}")
#         elif self.window.rg_flag:
#             self.ax_att.set_title(
#                 "Attenuation, receiver "
#                 + f"{self.traces.receiver[self.window.actual_traces[0]]+1}")
#         self.ax_r2.plot(f[1:nfmax],r2[1:nfmax])
#         self.ax_r2.set_xlabel("Frequency [Hz]")
#         self.ax_r2.set_ylabel("R2 coefficient")
#         self.window.setHelp(self.window.attenuation_text)

    def checkerboard(self):
        """
        Checkerboard test
        Function uses Pygimly for forward calculation.

        User creates first a synthetic checkerboard model with or without a
        background vertical velocity gradient. Checkerboard velocity variations
        may be given in percent of the background velocity or as absolute
        difference with respect to background velocity

        Returns
        -------
        None.

        """
        try:
            from pygimli.physics import TravelTimeManager
        except ModuleNotFoundError:
            _ = QtWidgets.QMessageBox.warning(
                None, "Warning", "Pygimli is not installed.\n "
                + "Tomography cannot be executed.\n",
                QtWidgets.QMessageBox.Close, QtWidgets.QMessageBox.Close)
            return

        answer = self.main.test_function()
        if not answer:
            return
        try:
            import pygimli as pg
            import pygimli.meshtools as mt
            import pygimli.physics.traveltime as tt
        except ModuleNotFoundError:
            _ = QtWidgets.QMessageBox.warning(
                None, "Warning", "PyGimli is not installed\n"
                + "Tomography connot be executed\n",
                QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.Ok)
            return
# check whether file "picks.sgt" exists
        try:
            scheme = pg.load("picks.sgt", verbose=True)
            positions = np.array(scheme.sensorPositions())
        except Exception as error:
            _ = QtWidgets.QMessageBox.warning(
                None, 'Warning', 'File "picks.sgt" does not exist\n'
                + 'Execute "Picking -> save Gimli format" and come back\n',
                QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.Ok)
            return
# Call dialog window for input of a number of inversion control parameters
        self.zmax_check = np.round((positions[:, 0].max() -
                                    positions[:, 0].min())*0.333)
        self.xmin_check = positions[:, 0].min()
        self.xmax_check = positions[:, 0].max()
        results, okButton = self.main.dialog(
            ["Maximum depth (m, positive down)",
             "Horizontal size of blocks [m]",
             "Vertical size of blocks [m]",
             "Starting position X [% block size]",
             "Starting position Z [% block size]",
             "Initial velocity at surface [m/s]",
             "Initial velocity at bottom [m/s]",
             "Velocity difference",
             "Noise level [s]"],
            ["e", "e", "e", "e", "e", "e", "e", "e", "e"],
            [self.zmax_check, (self.xmax_check-self.xmin_check)/5.,
             self.zmax_check/5., 0., 0., self.vmin, self.vmax, 0.1, 0.001],
            "Checkerboard parameters")

        if not okButton:
            print("\n Checkerboard test cancelled")
            self.main.function = "main"
            return
        self.zmax_check = float(results[0])
        self.hsize_check = float(results[1])
        self.vsize_check = float(results[2])
        self.startx_check = float(results[3])/100.*self.hsize_check
        self.startz_check = float(results[4])/100.*self.vsize_check
        self.v_surf_check = float(results[5])
        self.v_bott_check = float(results[6])
        self.v_grad = (self.v_bott_check-self.v_surf_check)/self.zmax_check
        self.diff = float(results[7])
        self.noise = float(results[8])

# Initialize PyGimli TravelTime Manager
        self.mgr_check = TravelTimeManager()
# Calculate edge positions of checkerboard fields
        if self.startx_check > 0:
            xck = [self.xmin]
        else:
            xck = []
        xck += list(np.arange(self.xmin_check+self.startx_check,
                              self.xmax_check, self.hsize_check))
        if xck[-1] < self.xmax_check:
            xck += [self.xmax_check]
        xck = np.array(xck)
        if self.startz_check > 0:
            zck = [0.]
        else:
            zck = []
        zck += list(np.arange(self.startz_check, self.zmax_check,
                              self.vsize_check))
        if zck[-1] < self.zmax_check:
            zck += [self.zmax_check]
        zck = -np.array(zck)
        layers = []
        for iz, z1 in enumerate(zck[:-1]):
            z2 = zck[iz+1]
            for ix, x1 in enumerate(xck[:-1]):
                x2 = xck[ix+1]
                m = (iz+ix) % 2 + 1
                layers.append(mt.createPolygon(
                    [[x1, z1], [x2, z1], [x2, z2], [x1, z2]],
                    isClosed=True, marker=m, area=m))
        geom = layers[0]
        for i in layers[1:]:
            geom += i
        self.mesh_check = mt.createMesh(geom, Quality=34.3, area=2,
                                        smooth=[1, 10])
        vel = []
        vp = np.array(self.mesh_check.cellMarkers())
        for node in self.mesh_check.nodes():
            vel.append(self.v_surf_check-node.y()*self.v_grad)
        self.v_check = pg.meshtools.nodeDataToCellData(self.mesh_check,
                                                       np.array(vel))
        if self.diff < 1:
            for i, vv in enumerate(vp):
                if vv == 1:
                    self.v_check[i] *= 1.-self.diff
                else:
                    self.v_check[i] *= 1.+self.diff
        else:
            for i, vv in enumerate(vp):
                if vv == 1:
                    self.v_check[i] -= self.diff
                else:
                    self.v_check[i] += self.diff
        data = tt.simulate(slowness=1.0/self.v_check, scheme=scheme,
                           mesh=self.mesh_check, verbose=True,
                           noiseAbs=self.noise, seed=1137, secNodes=5)
# Create window for plotting of checkerboard results
        self.w_check = newWindow("Checkerboard results")
        self.figck = self.w_check.fig
        plt.tight_layout()
        self.gs = GridSpec(15, 15, figure=self.figck)
# Axis for forward model
        self.ax_forward = self.figck.add_subplot(self.gs[:6, :6])
# Axis for inverse model
        self.ax_inverse = self.figck.add_subplot(self.gs[8:, :6])
# Axis for travel times
        self.ax_times = self.figck.add_subplot(self.gs[:, 8:])
        pg.viewer.show(self.mesh_check, self.v_check, colorBar=True,
                       logScale=True, label='Synthetic model [m/s]',
                       ax=self.ax_forward, xlabel="Distance [m]",
                       ylabel="depth [m]", showBoundary=True)
        # pg.viewer.mpl.drawSensors(self.ax_times, scheme.sensors(), diam=1.0,
        #                           facecolor='white', edgecolor='black')

        results, okButton = self.main.dialog(
            ["Initial smoothing parameter (<0: optimize)",
             "Smoothing reduction per iteration",
             "Smoothing in z direction (0..1):",
             "Maximum iterations (0 = automatic)",
             "Initial velocity at surface [m/s]",
             "Initial velocity at bottom [m/s]",
             "Minimum allowed velocity [m/s]",
             "Maximum allowed velocity [m/s]",],
            ["e", "e", "e", "e", "e", "e", "e", "e"],
            [self.smooth, self.s_fact, self.zSmooth, self.maxiter,
             self.vmin, self.vmax, self.vmin_limit, self.vmax_limit],
            "Inversion parameters for checkerboard")

        smooth = float(results[0])
        if self.smooth < 0:
            opt_flag = True
            smooth *= -1.
        else:
            opt_flag = False
        s_fact = float(results[1])
        zSmooth = float(results[2])
        maxiter = int(results[3])
        vini_surf = float(results[4])
        vini_bott = float(results[5])
        min_vel = float(results[6])
        max_vel = float(results[7])
        mgr = tt.TravelTimeManager(data)
        mgr.data["t"] = abs(mgr.data["t"])
        if self.noise == 0.:
            mgr.data["err"] = 0.001
        if opt_flag:
            self.mgr.inv.inv.setOptimizeLambda(True)
# Do tomography inversion. If maxiter == 0, stop iterationautomatically if
#    chi2<1 of if error does not decrease by more than 1% (dPhi=0.01)
#    if maxiter>0 stop iterations latest after maxiter iterations (but earlier
#    if one of the other conditions is fulfilled)
        if maxiter < 1:
            if min_vel == 0 and max_vel == 0:
                vest = mgr.invert(data=data, secNodes=3, paraMaxCellSize=5.0,
                                  zWeight=zSmooth, vTop=vini_surf,
                                  vBottom=vini_bott, verbose=1,
                                  paraDepth=self.zmax_check, dPhi=0.01,
                                  lam=smooth, lambdaFactor=s_fact,
                                  maxIter=1000)
            else:
                vest = mgr.invert(data=data, secNodes=3, paraMaxCellSize=5.0,
                                  zWeight=zSmooth, vTop=vini_surf,
                                  vBottom=vini_bott, verbose=1,
                                  paraDepth=self.zmax_check, dPhi=0.01,
                                  lam=smooth, limits=[min_vel, max_vel],
                                  lambdaFactor=s_fact, maxIter=1000,
                                  robustData=True, blockyModel=True)
        else:
            vest = mgr.invert(data=data, secNodes=3, paraMaxCellSize=15.0,
                              zWeight=zSmooth, vTop=vini_surf,
                              vBottom=vini_bott, maxIter=self.maxiter,
                              verbose=1, paraDepth=self.zmax_check,
                              dPhi=0.01, lam=smooth, limits=[min_vel, max_vel],
                              lambdaFactor=s_fact)

        # vest = mgr.invert(data=data, secNodes=3, paraMaxCellSize=5.
        #                   maxIter=10,
        #                   verbose=True, paraDepth=self.zmax_check,zWeight=0.)
        # vest = mgr.invert(data=data, secNodes=3, paraMaxCellSize=5.,
        #                   maxIter=10,
        #                   verbose=True, paraDepth=self.zmax_check,zWeight=0.)
        # np.testing.assert_array_less(mgr.inv.inv.chi2(), 0.8)
        mgr.showResult(cMin=min(vest), cMax=max(vest), logScale=True,
                       label='Inverse model [m/s]', ax=self.ax_inverse,
                       xlabel="Distance [m]", ylabel="depth [m]",
                       showBoundary=True)
        for iz, z1 in enumerate(zck[:-1]):
            z2 = zck[iz+1]
            for ix, x1 in enumerate(xck[:-1]):
                x2 = xck[ix+1]
                m = (iz+ix) % 2 + 1
                self.ax_inverse.plot([x1, x2, x2, x1, x1],
                                     [z1, z1, z2, z2, z1], "k")
        # ax2, _ = pg.show(geom, ax=self.ax_inverse, fillRegion=False,
        #                  regionMarker=False)
        # rays = mgr.drawRayPaths(ax=self.ax_inverse, color="k", lw=0.3,
        #                         alpha=0.5)
        mgr.showFit(firstPicks=True, ax=self.ax_times)
        self.w_check.show()
