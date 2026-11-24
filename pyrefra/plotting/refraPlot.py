# -*- coding: utf-8 -*-
"""
@author: Hermann Zeyen, University Paris-Saclay, France

Contains the following Class:
    Window
        Data window of PyRefra. Window itself only sets up widgets and the
        data-dependent state (__init__). The remaining state is defined in
        view_state.ViewState and is reachable as window.<name>; the methods are
        grouped in the following mixin classes:

    canvas_manager.CanvasMixin: addMPL, rmMPL, addFig, drawNew, setHelp,
        savePlot; also class newWindow (re-exported here)
    zoom.ZoomMixin: addZoom, changeZoom, zooming, zoomOut, zoomIn, zoomIni,
        nextBlock
    render.RenderMixin: seismogram, agcCalc, set_ticks
    gathers.GatherMixin: changeFig, toggleNumber, chooseComponent, plotRG,
        plotReceiver, plotDG, plotDistance, plotSG, plotShot, plotFG, plotFile
    interaction.InteractionMixin: followRect, followLine, follow_air
    picking_display.PickDisplayMixin: picksPlot, pickPlot, plotPickSection,
        plotCalcPicks, searchTrace, searchPick, searchNearPicks
    picking_manual.ManualPickMixin: pickMove, movePick, uncertainty,
        changeUnc, shiftPickTrace, erase_Picks, pickManual
    picking_auto.AutoPickMixin: findNearest2ndDerivative, corrPick,
        search_corr_pick, Sta_Lta, ampPick, ampPicks, bestLines, L1_regres
    trace_processing.ProcessingMixin: original, originalScreen, changeAmp,
        tNorm, tGain, dGain, AGC, phasePlot, phaseCalc, traceMute, muteAir,
        muteBefore, muteAfter, muteTime, traceSign, changeSign, animateLine
"""

import os
import numpy as np
from PyQt5 import QtWidgets, QtGui
from PyQt5.uic import loadUiType
from matplotlib.colors import TABLEAU_COLORS
from matplotlib.figure import Figure
from .canvas_manager import CanvasMixin, newWindow  # noqa: F401
from .gathers import GatherMixin
from .interaction import InteractionMixin
from .picking_auto import AutoPickMixin
from .picking_display import PickDisplayMixin
from .picking_manual import ManualPickMixin
from .render import RenderMixin
from .trace_processing import ProcessingMixin
from .view_state import ViewState, install_state_access
from .zoom import ZoomMixin

# Ui_MainWindow, QMainWindow = loadUiType("refraWindow.ui")
Ui_MainWindow, QMainWindow = loadUiType(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "refraWindow.ui"))


class Window(CanvasMixin, ZoomMixin, RenderMixin, GatherMixin,
             InteractionMixin, PickDisplayMixin, ManualPickMixin,
             AutoPickMixin, ProcessingMixin, QMainWindow, Ui_MainWindow):
    """
    Data window of PyRefra: plotting, zooming and picking.
    The methods are implemented in the mixin classes listed above.
    """

    def __init__(self, main):
        super(Window, self).__init__()
        self.view_state = ViewState()
# Display of trace numbers at the cursor (toggled with key X)
        self.hover_on = False
        self.hover_ann = None
        self.hover_bg = None
# Set up main window based on file window.ui created with QT Designer
        self.setupUi(self)
# main points to class Main (file PyRefra.py)
        self.main = main
# files points to class Files (file refraData.py)
        self.files = main.files
# data points to class Data (file refraData.py)
        self.data = main.data
# traces points to class Traces (file refraData.py)
        self.traces = main.traces
# geom points to class Geometry (file refraData.py)
        self.geom = main.geo
# Title of main window
        self.setWindowTitle("Data window")
# With several screens, place the window on the physically biggest one (size
# in mm, not pixels); the later showMaximized() then maximizes it there
        screens = QtWidgets.QApplication.screens()
        if len(screens) > 1:
            biggest = max(screens, key=lambda s: s.physicalSize().width()
                          * s.physicalSize().height())
            self.setGeometry(biggest.availableGeometry())
# Get program icon
        # self.setWindowIcon(QtGui.QIcon(main.top+"/PyRefra_Logo.png"))
        self.setWindowIcon(QtGui.QIcon(main.top+"/PyRefra_Logo.ico"))
# create a first figure in central widget
        self.fig = Figure()
        self.addMPL(self.fig)
# Initialize vertical slider as invisible (used in refraData.Utilities.FK_filt)
        self.verticalSlider.setVisible(False)
# Figs contains a list of all plots (all gathers)
        self.figs.append(self.fig)
# placeholder fo progressbar used for possibly slow processes
        self.progressBar = None
        self.line = None
        self.ax = None
        self.axl = None
        if len(self.geom.types) > 0:
            self.phaseAngles.setEnabled(True)
        self.general_sign = self.data.general_sign
# unc_base is the number of samples assigned by default (automatic picking or
# double click during manual picking) to picking uncertainties of 0.5 ms
        self.unc_base = int(0.0005/self.data.dt)
        self.lin = None
        self.tab_colors = list(TABLEAU_COLORS.keys())
# Next line not yet used!
# All potential triggers for Sta-Lta later then t = offset/v_min_trigger are
#     ignored
        self.v_min_trigger = 100
# all potential triggers for Sta-Lta earlier then t = offset/v_max_trigger are
#     ignored
        self.v_max_trigger = 2000
        self.w_anim = None
        self.w_picks = None
# Set initial zoom limits such that all data may be plotted on the screen
        self.nt_mx = self.data.nsamp
        self.time_plt_min = self.data.t0
        self.time_plt_max = self.data.t0+self.data.dt*(self.nt_mx-1)
        self.nt_0 = int(-self.data.t0/self.data.dt)
        self.dxmax = self.traces.off_max*1
        self.dxmin = self.traces.off_min*1
        self.rxmax = self.traces.off_max*1
        self.rxmin = self.traces.off_min*1
        self.sxmax = self.traces.off_max*1
        self.sxmin = self.traces.off_min*1
        self.time = self.data.t0 + self.data.dt*np.arange(self.data.nsamp)
        self.addZoom(self.sxmin, self.sxmax, 0, self.nt_mx-1)
# Set texts printed at bottom of sceen as help for different actions
        self.main_text = "Main module: on right window choose data set to be "\
            + "plotted, Numpad+/Numpad- to change amplitudes; "\
            + "right&left arrow: shift window right&left; "\
            + "X: show/hide trace numbers at cursor"
        self.pick_amp_text = "Automatic picking by amplitude analysis"
        self.pick_manual_text = "Manual picking module: left click for pick, "\
            + "second left click for uncertainty (double click: "\
            + "2 samples uncertainty), central click to erase pick, "\
            + "right click to finish; "\
            + "SHFT+left click: accept picks along green line; "\
            + "blue dashed line: air-wave arrivals; "\
            + "green dashed line: nearby earlier picks; "\
            + "if green line exists: SHFT+left click: accept picks along the line"
        self.pick_move_text = "Move picks module: left click to choose pick, "\
            + "up/down or Numpad8/Numpad2 key to move pick, with "\
            + "SHFT 10 samples, with CTRL 1 ms; Right/Left key to "\
            + "chose next trace to the right/left; right click to "\
            + "finish;    dashed line corresponds to air-wave arrivals"
        self.uncertainty_text = "Change uncertainty module: left click to "\
            + "choose pick, up/down or Numpad8/Numpad2 key to change "\
            + "uncertainty, with SHFT 10 samples, with CTRL  ms; Right/Left "\
            + "key to chose next trace to the right/left; right click to "\
            + "finish"
        self.zoom_text = "Zooming module: left click, keep pressed. On "\
            + "release zoomed window is shown"
        self.p_model_text = "P-wave model module: left click and keep pressed"\
            + " for each phase. Release to accept. Right click to finish"
        self.tau_text = "Tau-Pi module: to return to seismograms click on "\
            + "file in right window"
        self.attenuation_text = "Attenuation module: to return to seismograms"\
            + "click on file in right window"
        self.trace_mute_text = "Click mouse on all traces which should be "\
            + "muted or recovered, finish with right click.  Traces marked "\
            + "at bottom with asterisc: red=eliminate, green=recover"
        self.trace_sign_text = "Click mouse on all traces which should "\
            + "inverse polarity, finish with right click. Traces marked at "\
            + "bottom with asterisc"
        self.mute_before_text = "Mute before module: left click and keep "\
            + "pressed for each partial line. Release to accept. Right click "\
            + "to finish (point is used). Data are muted before the line drawn"
        self.spectrum_text = "Mark low-cut with left mouse, keep pressed and "\
            + "release at high-cut. Point on negative horizontal axis means "\
            + "that corresponding filter is not applied"
        self.fk_text = "Pull slider left of screen to desired velocity and "\
            + "release mouse to accept velocity, or simply click on slider "\
            + "if default velocity is ok."
        self.trace_filter_text = "Click on trace to be filtered"
        self.inver_text = "Inversion results: to return to seismograms click "\
            + "on file in right window; Press C to change color scale and "\
            + "maximum plotted depth"
        self.cpick_text = "Correlation picking: Do one first manual pick on "\
            + "good trace; Right click = cancel"
        self.falseCol_text = "Chose up to 3 indicators for false colour plot"
        self.animate_text = "Animation of spatial wave image. To return to "\
            + "seismograms click on file in right window"

        self.help = QtWidgets.QLabel(self)
        self.setHelp(self.main_text)


install_state_access(Window)

