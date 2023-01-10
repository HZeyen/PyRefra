# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 18:30:59 2019
last modified on Tue Dec 06, 2022
@author: Hermann

Contains the following Class:
    Window

Contains the following functions:
    addMPL
    rmMPL
    addFig
    changeFig
    SavePlot
    original
    originalScreen
    drawNew
    setHelp
    changeAmp
    tNorm
    tGain
    dGain
    AGC
    traceMute
        on Press
    muteAir
    muteBefore
    muteAfter
    muteTime
    traceSign
        onPress
    changeSign
    zooming
    zoomOut
    zoomIn
    zoomIni
    agcCalc
    seismogram
    plotRG
    plotReceiver
    plotDG
    plotDistance
    plotSG
    plotShot
    plotFG
    plotFile
    set_ticks
    nextBlock
    picksPlot
    pickPlot
    plotCalcPicks
    searchTrace
    searchPick
    pickMove
        onPress
    movePick
    uncertainty
        onPress
    changeUnc
    shiftPickTrace
    searchNearPicks
    erase_Picks
    pickManual
        onPress
    corrPick
        onPress
    Sta_Lta
    ampPick
    ampPicks
    followRect
        onPress
        onRelease
        onMove
    followLine
        onPress
        onRelease
        onMove
    bestLines
    L1_regres
        cost_function
    animateLine
"""

import numpy as np
from matplotlib.figure import Figure
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.uic import loadUiType
import copy
from matplotlib.backends.backend_qt5agg import(
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)
Ui_MainWindow, QMainWindow = loadUiType("refraWindow.ui")

class Window(QMainWindow, Ui_MainWindow):
    def __init__(self, main, files, data, traces, geom):
        super(Window, self).__init__()
        self.setupUi(self) #Set up main window based on file window.ui created with QT Designer
        self.setWindowTitle("Data window")
        self.setWindowIcon(QtGui.QIcon("C:\Documents\Orsay\logoIDES.jpg"))
        self.fig = Figure() #create a first figure in central widget
        self.addMPL(self.fig)
        self.verticalSlider.setVisible(False)
        self.figs=[]
        self.figs.append(self.fig)
        self.fig_plotted = 0
        self.amp_plt = 1
        self.gain = "tnorm"
        self.general_sign = 1
        self.n_zooms = 0
        self.i_zooms = -1
        self.zooms = []
        self.v_set = False
        self.v_min_trigger = 100 # not yet used! all potential triggers for Sta-Lta later then t = offset/v_min_trigger are ignored
        self.v_max_trigger = 2000 # all potential triggers for Sta-Lta earlier then t = offset/v_max_trigger are ignored
        self.main = main
        self.files = files
        self.data = data
        self.traces = traces
        self.geom = geom
# Set initial zoom limits such that all data may be plotted on the screen
        self.nt_mn = 0
        self.nt_mx = self.data.nsamp
        self.time_plt_min = self.data.t0
        self.time_plt_max = self.data.t0+self.data.dt*(self.nt_mx-1)
        self.nt_0 = np.int(-self.data.t0/self.data.dt)
        self.dxmax = self.traces.off_max*1
        self.dxmin = self.traces.off_min*1
        self.rxmax = self.traces.off_max*1
        self.rxmin = self.traces.off_min*1
        self.sxmax = self.traces.off_max*1
        self.sxmin = self.traces.off_min*1
        self.time = self.data.t0 + self.data.dt*np.arange(self.data.nsamp)
        self.addZoom(self.sxmin, self.sxmax, 0, self.nt_mx-1)
# Set texts printed at bottom of sceen as help for different actions
        self.main_text = "Main module: on right window choose data set to be "+\
                         "plotted, Numpad+/Numpad- to change amplitudes; "+\
                         "right&left arrow: shift window right&left"
        self.pick_amp_text = "Automatic picking by amplitude analysis"
        self.pick_manual_text = "Manual picking module: left click for pick, "+\
                         "second left click for uncertainty (double click: "+\
                         "2 samples uncertainty), central click to erase pick, "+\
                         "right click to finish;    dashed line corresponds "+\
                         "to air-wave arrivals"
        self.pick_move_text = "Move picks module: left click to choose pick, "+\
                         "up/down or Numpad8/Numpad2 key to move pick, with "+\
                         "SHFT 10 samples, with CTRL 1 ms; Right/Left key to "+\
                         "chose next trace to the right/left; right click to "+\
                         "finish;    dashed line corresponds to air-wave "+\
                         "arrivals"
        self.uncertainty_text = "Change uncertainty module: left click to "+\
                         "choose pick, up/down or Numpad8/Numpad2 key to "+\
                         "change uncertainty, with SHFT 10 samples, with CTRL "+\
                         "1 ms; Right/Left key to chose next trace to the "+\
                         "right/left; right click to finish"
        self.zoom_text = "Zooming module: left click, keep pressed. On release "+\
                         "zoomed window is shown"
        self.p_model_text = "P-wave model module: left click and keep pressed "+\
                         "for each phase. Release to accept. Right click to "+\
                         "finish"
        self.tau_text = "Tau-Pi module: to return to seismograms click on file "+\
                        "in right window"
        self.attenuation_text = "Attenuation module: to return to seismograms "+\
                        "click on file in right window"
        self.trace_mute_text = "Click mouse on all traces which should be muted"+\
                        " or recovered, finish with right click.  Traces "+\
                        "marked at bottom with asterisc: red=eliminate, "+\
                        "green=recover"
        self.trace_sign_text = "Click mouse on all traces which should inverse "+\
                        "polarity, finish with right click. Traces marked at "+\
                        "bottom with asterisc"
        self.mute_before_text = "Mute before module: left click and keep "+\
                        "pressed for each partial line. Release to accept. "+\
                        "Right click to finish (point is used). Data are "+\
                        "muted before the line drawn"
        self.spectrum_text = "Mark low-cut with left mouse, keep pressed and "+\
                        "release at high-cut. Point on negative horizontal "+\
                        "axis means that corresponding filter is not applied"
        self.fk_text = "Pull slider left of screen to desired velocity and "+\
                       "release mouse to accept velocity, or simply press "+\
                       "ENTER if default velocity is ok."
        self.trace_filter_text = "Click on trace to be filtered"
        self.inver_text = "Inversion results: to return to seismograms click "+\
                       "on file in right window; Press C to change color "+\
                       "scale and maximum plotted depth"
        self.cpick_text = "Correlation picking: Do one first manual pick on "+\
                       "good trace"
        self.falseCol_text = "Chose up to 3 indicators for false colour plot"
        self.animate_text = "Animation of spatial wave image. To return to "+\
                       "seismograms click on file in right window"

        self.help = QtWidgets.QLabel(self)
        self.setHelp(self.main_text)


    def addMPL(self, fig):
        """
        Add widget to actual Figure
        Input: fig is the actually used Figure (matplotlib.figure.Figure object)
        """
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,self)
        self.mplvl.addWidget(self.toolbar)

    def rmMPL(self, ):
        """
        Remove widget from actual Figure
        """
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()

    def addFig(self, name, fig):
        """
        add data set to list of plottable data sets (right side of screen)
        input:
            name: name of the data set to be added (text)
            fig: Actual figure (matplotlib.figure.Figure object)
        """
        self.fig_dict[name] = fig
        self.mplfigs.addItem(name)

    def changeFig(self, item):
        """
        After having clicked on a plottable data set on the right side of the screen
            this function is called to change the plot
        input:
            item: index of the data set to be added,
                    automatically assigned in line "self.mplfigs.itemClicked.connect"
        """
        self.v_set = False
        text = item.text()
        self.fig_plotted = self.plot_names.index(text) #keep name of plotted data set
        self.drawNew(False) # Draw only a frame
# Plot record section, depending on which type of plot is wanted
#    (shot section, file section, record section, distance section)
        if self.sg_flag:
# Keep zoom in time, but plot all traces
            self.zooms[self.i_zooms][0] = self.traces.off_min
            self.zooms[self.i_zooms][1] = self.traces.off_max
            self.fig_plotted_S = self.fig_plotted
            self.plotShot(self.axes[self.fig_plotted],self.fig_plotted)
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
        elif self.fg_flag:
# Keep zoom in time, but plot all traces
            self.zooms[self.i_zooms][0] = self.traces.off_min
            self.zooms[self.i_zooms][1] = self.traces.off_max
            self.fig_plotted_F = self.fig_plotted
            self.plotFile(self.axes[self.fig_plotted],self.fig_plotted)
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
        elif self.rg_flag:
# Keep zoom in time, but plot all traces. Offset of shots with respect to
#      receivers has reverse sign than receiver vs. offset, therefore the
#      negative sign and inversion of min and max
            self.zooms[self.i_zooms][0] = -self.traces.off_max
            self.zooms[self.i_zooms][1] = -self.traces.off_min
            self.fig_plotted_R = self.fig_plotted
            self.plotReceiver(self.axes[self.fig_plotted],self.receivers[self.fig_plotted])
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
        else:
# Keep zoom in time, but plot all traces. Distances are plotted with respect to
#      common depth point coordinate
            self.fig_plotted_D = self.fig_plotted
            self.zooms[self.i_zooms][0] = self.traces.xcdp.min()
            self.zooms[self.i_zooms][1] = self.traces.xcdp.max()
            self.plotDistance(self.axes[self.fig_plotted],self.distances[self.fig_plotted])
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
        if self.i_zooms == 0:
            self.i_zooms = 1
            try:
                self.zooms[1] =  self.zooms[0]
                self.zooms[self.i_zooms][0] = self.traces.off_min
                self.zooms[self.i_zooms][1] = self.traces.off_max
            except:
                self.zooms.append(self.zooms[0])
                self.zooms[1][0] = self.traces.off_min
                self.zooms[1][1] = self.traces.off_max
        else:
            self.zooms[self.i_zooms][0] = self.traces.off_min
            self.zooms[self.i_zooms][1] = self.traces.off_max
        self.picksPlot() # Plot picks if available
        self.setHelp(self.main_text) # plot help line for self module
        self.v_set = True
        self.main.function = "main"

    def savePlot(self):
        """
        Function saves plot inside actual window into a png file
        The file name is prefix_date_time.png
        Prefix depends on the actual image type (inversion result, tau-p,
            data gather type).

        Returns
        -------
        None.

        """
        from datetime import datetime,date
        now = datetime.now()
        c_time = now.strftime("%H-%M-%S")
        today = date.today()
        d1 = today.strftime("%Y-%m-%d")
        if self.main.function == "inver":
            fname = f"Inversion_{d1}_{c_time}.png"
            self.figs[-1].savefig(fname)
        elif self.main.function == "tau_p":
            tit = self.axes[self.fig_plotted].get_title()
            tits = tit.split()
            fname = f"Tau-P_{tits[1]}-{tits[2]}_{d1}_{c_time}.png"
            self.figs[self.fig_plotted].savefig(fname)
        else:
            if self.rg_flag:
                fname = f"receiver_{self.actual_receivers[0]+1:0>5}_{d1}_{c_time}.png"
            if self.fg_flag:
                fil_num = int(self.plot_names[self.fig_plotted_F].split()[1])
                fname = f"file_{fil_num:0>5}_{d1}_{c_time}.png"
            elif self.sg_flag:
                fname = f"shot_point_{self.actual_shot+1:0>5}_{d1}_{c_time}.png"
            elif self.dg_flag:
                fname = f"distance_{int(round(self.actual_dist,0)):0>5}_{d1}_{c_time}.png"
            else:
                fname = f"file_{self.actual_shot:0>5}_{d1}_{c_time}.png"
            self.figs[self.fig_plotted].savefig(fname)


    def original(self):
        self.main.function = "plot_original"
        self.data.st = copy.deepcopy(self.data.st_ori)
        self.v_set = False
        self.drawNew(True)
        self.v_set = True
        self.setHelp(self.main_text) # Change help text to self module
        self.filtered = False

    def originalScreen(self):
        for t in self.actual_traces:
            nf = self.traces.file[t]
            nt = self.traces.trace[t]
            self.data.st[nf][nt].data = self.data.st_ori[nf][nt].data
        self.v_set = False
        self.drawNew(True)
        self.v_set = True
        self.setHelp(self.main_text) # Change help text to self module

    def drawNew(self, plot_flag):
        """
        Draw a new data set
        input:
            plot_flag (Boolean): If True plot seismograms
                                 if False only plot frame
        """
        self.rmMPL() # remove actual window
        self.figs[self.fig_plotted].clear() # remove actual Figure
        self.figs[self.fig_plotted] = Figure() # create new Figure
        self.axes[self.fig_plotted] = self.figs[self.fig_plotted].add_subplot() #add axes to Figure
        self.addMPL(self.figs[self.fig_plotted]) #Add Figure to central widget
        if plot_flag:
# Plot record section, depending on which type of plot is wanted (shot section, record section, distance section)
            if self.dg_flag:
                self.plotDistance(self.axes[self.fig_plotted],\
                                  self.distances[self.fig_plotted])
                if self.PlotCalculatedTimes.isChecked():
                    self.plotCalcPicks()
            elif self.rg_flag:
                self.plotReceiver(self.axes[self.fig_plotted],\
                                  self.receivers[self.fig_plotted])
                if self.PlotCalculatedTimes.isChecked():
                    self.plotCalcPicks()
            elif self.fg_flag:
                self.plotFile(self.axes[self.fig_plotted],self.fig_plotted)
                if self.PlotCalculatedTimes.isChecked():
                    self.plotCalcPicks()
            else:
                self.plotShot(self.axes[self.fig_plotted],self.fig_plotted)
                if self.PlotCalculatedTimes.isChecked():
                    self.plotCalcPicks()

            self.picksPlot() # Plot picks if available
        self.main.function = "main"

    def setHelp(self,text):
        """
        Set help text at bottom of screen
        Input: text is the text to be printed (defined in __init__)
        Text is written in QLabel widget. In order to keep the widget at the bottom of the screen,
            the existing one is first erased and then reopened
        This call has to be done after any modification of the graphics window
        """
        self.help.close()
        self.help = QtWidgets.QLabel(self)
        self.help.setMaximumHeight(15)
        self.mplvl.addWidget(self.help)
        self.help.show()
        self.help.setText(text)
        self.figs[self.fig_plotted].canvas.flush_events()

    def changeAmp(self,sens):
        """
        Change plotting amplitude of data set
        The value is selftained for the whole run until called
        again, also when changing data set.
        input:
           sens (int of float): If > 0 increase amplitude by a factor 2
                                if < 0 decrease amplitude by a factor 1/sqrt(2)
        """
        self.main.function = "change_amp"
        if sens>0:
            self.amp_plt *= 2
        else:
            self.amp_plt /= np.sqrt(2)
        print("Amp:",self.amp_plt)
        self.drawNew(True)
        self.setHelp(self.main_text) # Change help text to self module

    def tNorm(self):
        """
        Function set gain as trace normalized
        """
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
        self.setHelp(self.main_text) # Change help text to self module

    def tGain(self):
        """
        Function set gain as time**time_gain
        """
        results, okButton = self.main.dialog(\
                                    ["Exponent for time gain"],["e"],["2"],\
                                    "Time gain")
        if okButton == False:
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
            self.setHelp(self.main_text) # Change help text to self module

    def dGain(self):
        """
        Function set gain as distance**dist_gain
        """
        results, okButton = self.main.dialog(\
                                            ["Exponent for distance gain"],\
                                             ["e"],["2"],"Distance gain")
        if okButton == False:
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
            self.setHelp(self.main_text) # Change help text to self module

    def AGC(self):
        """
        Function set gain as AGC
        """
        results, okButton = self.main.dialog(\
                                            ["AGC window length [ms]"],\
                                            ["e"],["100"],"AGC gain")
        if okButton == False:
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
            self.setHelp(self.main_text) # Change help text to self module

    def agcCalc(self,data,n_agc_win):
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
            i1 = max(j-n_agc_win,0)
            i2 = min(j+n_agc_win,nsamp)
            m = np.max(x[i1:i2])
            if m > 0:
                data[j] = data[j]/m
        return data

    def traceMute(self,QselfWindow):
        """
        Eliminate one trace which is chosen by a mouse click, i.e. set
        amplitudes to zero. The choice is selftained for the whole run
        until called again, also when changing data set.
        """
        self.main.function = "trace_null"
        self.finish = False
        self.picked = False
        self.setHelp(self.trace_mute_text)
        def onPress(event):
            if event.button == 1:
                self.searchTrace(event.xdata)
                trace = self.itrace
                nf = self.traces.file[trace]
                nt = self.traces.trace[trace]
                if np.isclose(self.traces.amplitudes[trace], 0.):
                    self.traces.amplitudes[trace] = self.general_sign
                    try:
                        fac = self.data.receiver_corr_dict[self.p_r]["amp"]
                        if fac != 0:
                            self.traces.amplitudes[trace] *= fac
                    except:
                        pass
                    self.v[self.n_tr,:] = self.data.st[nf][nt].data*self.general_sign
                    self.axes[self.fig_plotted].plot(self.x[self.n_tr],\
                                                self.time[self.nt_mn+5],"*g")
                    print("trace ", self.n_tr, " recovered")
                else:
                    self.traces.amplitudes[trace] = 0
                    self.v[self.n_tr,:] *= 0
                    self.axes[self.fig_plotted].plot(self.x[self.n_tr],\
                                                self.time[self.nt_mn+5],"*r")
                    print("trace ", self.n_tr, " muted")
                self.canvas.draw()
            elif event.button == 3:
                print("finish trace_mute")
                self.lin.set_animated(False)
                self.finish = True
                self.lin.figure.canvas.mpl_disconnect(self.cidpress)
                self.drawNew(True)
                self.setHelp(self.main_text) # Change help text to self module
        self.x_coor=[]
        self.y_coor=[]
        self.lin, = self.axes[self.fig_plotted].plot(self.x_coor, self.y_coor,\
                                                     animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect('button_press_event',\
                                                           onPress)

    def muteAir(self,QselfWindow):
        """
        Routine zeroes data around the air wave arrival for all traces

        Parameters
        ----------
        QselfWindow

        Returns
        -------
        modified traces

        """
        self.main.function = "Mute_air"
        results, okButton = self.main.dialog(\
                                    ["Air velocity [m/s]",\
                                     "Mute-time [ms]",\
                                     "Width sin**2 taper [ms]"],\
                                    ["e","e","e"],[340,30,5],"Mute air wave")
        if okButton == False:
            print("Air mute cancelled")
            return
        air_speed = np.float(results[0])
        mute = np.float(results[1])
        sin_width = np.float(results[2])
        nmute = np.int(mute/self.data.dt/1000)
        nsin = np.int(sin_width/self.data.dt/1000)
        nfac = np.int(nmute+2*nsin)
        factor = np.zeros(nfac)
        factor[0:nsin] = (np.cos(0.5*np.pi/nsin*np.arange(nsin)))**2
        factor[nfac-1:nfac-nsin-1:-1] = factor[0:nsin]
        for i in range(self.traces.number_of_traces):
            nf = self.traces.file[i]
            nt = self.traces.trace[i]
            tini = np.abs(self.traces.offset[i])/air_speed
            nini = np.int((tini-self.data.t0)/self.data.dt-nsin)
            nend = nini+nfac
            nini_fac = 0
            nend_fac = nfac
            if nini<0:
                nini_fac = -nini
                nini = 0
            elif nend > self.data.nsamp:
                nend_fac = nfac-(nend-self.data.nsamp)
                nend = self.data.nsamp
            self.data.st[nf][nt].data[nini:nend]*=factor[nini_fac:nend_fac]
        self.v_set=False
        self.drawNew(True)
        self.setHelp(self.main_text)

    def muteBefore(self):
        self.muteTime(-1)

    def muteAfter(self):
        self.muteTime(1)

    def muteTime(self,sign):
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
        """
        self.main.function = "Mute_before"
        self.setHelp(self.mute_before_text)
# Define the line before which the samples are zeroed
        self.followLine(False,1,1)
# (xline,yline) are the distance and time coordinates of the points defining the line
        xline = np.array(self.coor_x)
        yline = np.array(self.coor_y)
# Sort the points in increasing X order
        line_index = np.argsort(xline)
        ntraces = np.arange(len(self.x))
        for i in range(len(line_index)-1):
# Treat one segment after the other
            i0 = line_index[i]
            i1 = line_index[i+1]
# If too a point is double clicked, skip the second occurrence
            if i0>0:
                if xline[i0]==xline[line_index[i-1]]:
                    continue
# Search traces lying inside a line segment
            npos = ntraces[self.x >= xline[i0]]
# Not sure why the folloing line does not work, replaced by the following for loop
#                ntreat = npos[self.x[npos] <= xline[i1]]
            nt = []
            for j in npos:
                if self.x[j]>xline[i1]:
                    break
                else:
                    nt.append(j)
# ntreat contains the numbers of all traces found within the segment
            ntreat = np.array(nt,dtype='int')
            for j in ntreat:
# Treat the corresponding traces
# First find time where the line crosses the trace and the sample number
                nt_cut = yline[i0]+(self.x[j]-xline[i0])/(xline[i1]-xline[i0])*\
                    (yline[i1]-yline[i0])
                ncut = np.int((nt_cut-self.data.t0)/self.data.dt)
# Find the sample after and including the line cut position which has a (relative)
# minimum amplitude
                test = False
                k = ncut
                while test==False:
                    if np.abs(self.v[j,k])<=np.abs(self.v[j,k-1]):
                        if np.abs(self.v[j,k])<=np.abs(self.v[j,k+1]):
                            test = True
                            break
                    k += 1
                kpos = k-1
# Find the sample before the line cut position which has a (relative)
# minimum amplitude
                test = False
                k = ncut-1
                while test==False:
                    if np.abs(self.v[j,k])<=np.abs(self.v[j,k-1]):
                        if np.abs(self.v[j,k])<=np.abs(self.v[j,k+1]):
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
                if sign<0:
                    self.v[j,:ncut] = 0
                else:
                    self.v[j,ncut:] = 0
                trace = self.actual_traces[j]
                nf = self.traces.file[trace]
                nt = self.traces.trace[trace]
                self.data.st[nf][nt].data = self.v[j,:]
# Plot the treated record section
        self.drawNew(True)
        self.setHelp(self.main_text)


    def traceSign(self,QselfWindow):
        """
        Change polarity of one trace which is chosen by a mouse click
        (the value is selftained for the whole run until called again,
          also when changing data set)
        """
        self.main.function = "trace_sign"
        self.finish = False
        self.picked = False
        self.setHelp(self.trace_sign_text)
        def onPress(event):
            if event.button == 1:
                self.searchTrace(event.xdata)
                self.traces.amplitudes[self.n_tr] *= -1
                self.v[self.n_tr,:] *= -1
                print("trace ", self.n_tr, " changed sign")
                self.setHelp(self.trace_sign_text) # Rewrite help text
                self.picked = True
                self.axes[self.fig_plotted].plot(self.x[self.n_tr],self.data.time[self.nt_mn+5],"*r")
                self.canvas.draw()
            elif event.button == 3:
                print("finish trace_sign")
                self.lin.set_animated(False)
                self.finish = True
                self.lin.figure.canvas.mpl_disconnect(self.cidpress)
                self.drawNew(True)
                self.setHelp(self.main_text) # Change help text to self module
        self.x_coor=[]
        self.y_coor=[]
        self.lin, = self.axes[self.fig_plotted].plot(self.x_coor, self.y_coor, animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect('button_press_event', onPress)

    def changeSign(self):
        """
        Change polarity of all traces
        (the sign is selftained for the whole run until called again,
          also when changing data set)
        """
        self.main.function = "change_sign"
        self.traces.amplitudes *= -1
        self.data.general_sign *= -1
        self.v *= -1
        self.drawNew(True)
        self.setHelp(self.main_text) # Change help text to self module

    def addZoom(self, xmin, xmax, nt_min, nt_max):
        self.n_zooms += 1
        self.i_zooms += 1
        self.zooms.append([xmin, xmax, nt_min, nt_max])

    def changeZoom(self, i, xmin, xmax, nt_min, nt_max):
        self.i_zooms = i
        self.zooms[i][:] =[xmin, xmax, nt_min, nt_max]

    def zooming(self):
        """
        Zoom by drawing a rectangle
        Input: non
        Output: Bool
                False if rectangle too small (usually zero size by inadvertancy)
        """
        self.main.function = "zoom"
        zm = False
        self.setHelp(self.zoom_text) # Change help text
        self.followRect() # Pull a rectangle over the window to define area to be plotted
        if self.coor_x[0] is None or self.coor_x[1] is None or\
            self.coor_y[0] is None or self.coor_y[1] is None:
            print("\nZoom coordinates not correctly detected. Zoom not executed")
            self.setHelp(self.main_text) # Change help text to self module
            self.drawNew(True)
            self.setHelp(self.main_text) # Change help text to self module
            return zm
        xmi = min(self.coor_x)
        xma = max(self.coor_x)
        tmi = min(self.coor_y)
        tma = max(self.coor_y)
        if tmi>=self.time_plt_max or tma<=self.time_plt_min:
            print("\nTime window outside actual zoom. Zoom not executed")
            self.setHelp(self.main_text) # Change help text to self module
            self.drawNew(True)
            self.setHelp(self.main_text) # Change help text to self module
            return zm
        if xmi>=self.x.max() or xma<=self.x.min():
            print("\nX window outside actual zoom. Zoom not executed")
            self.setHelp(self.main_text) # Change help text to self module
            self.drawNew(True)
            self.setHelp(self.main_text) # Change help text to self module
            return zm
        if tma-tmi < 10*self.data.dt:
            print("\nTime axis too much zoomed. Zoom not executed")
            self.setHelp(self.main_text) # Change help text to self module
            self.drawNew(True)
            self.setHelp(self.main_text) # Change help text to self module
            return zm
        if xma-xmi < 2*self.geom.dx_geo:
            print("\nX axis too much zoomed. Zoom not executed")
            self.setHelp(self.main_text) # Change help text to self module
            self.drawNew(True)
            self.setHelp(self.main_text) # Change help text to self module
            return zm
        self.time_plt_min = tmi # Define start time of zoomed plot
        self.time_plt_max = tma # Define end time of zoomed plot
        self.nt_mn = max(int((self.time_plt_min-self.data.t0)/self.data.dt),0) # Define start sample of zoomed plot
        self.nt_mx = min(int((self.time_plt_max-self.data.t0)/self.data.dt),\
                         max(self.traces.nsample_trace)) # Define end sample of zoomed plot
        self.rmMPL() # remove actual window
        self.figs[self.fig_plotted].clear() # Remove actual Figure
        self.figs[self.fig_plotted] = Figure() # create new Figure
        self.axes[self.fig_plotted] = self.figs[self.fig_plotted].add_subplot(111) #add axes to Figure
        self.addMPL(self.figs[self.fig_plotted]) #Add Figure to central widget
# Add zoom limits to list of existing zooms to be used by "zoom_out"
#     actual values depend on type of sections to be plotted
# Then plot section with actual zoom parameters
        if self.dg_flag:
            self.dxmin = xmi
            self.dxmax = xma
            self.i_zooms += 1
            if self.i_zooms >= self.n_zooms:
                self.zooms.append([self.dxmin, self.dxmax, self.nt_mn, self.nt_mx])
                self.n_zooms += 1
            else:
                self.zooms[self.i_zooms][:] = [self.dxmin, self.dxmax, self.nt_mn, self.nt_mx]
            self.plotDistance(self.axes[self.fig_plotted],self.distances[self.fig_plotted])
        elif self.rg_flag:
            self.rxmin = xmi
            self.rxmax = xma
            self.i_zooms += 1
            if self.i_zooms >= self.n_zooms:
                self.zooms.append([self.rxmin, self.rxmax, self.nt_mn, self.nt_mx])
                self.n_zooms += 1
            else:
                self.zooms[self.i_zooms][:] = [self.rxmin, self.rxmax, self.nt_mn, self.nt_mx]
            self.plotReceiver(self.axes[self.fig_plotted],self.receivers[self.fig_plotted])
        elif self.fg_flag:
            self.sxmin = xmi
            self.sxmax = xma
            self.i_zooms += 1
            if self.i_zooms >= self.n_zooms:
                self.zooms.append([self.sxmin, self.sxmax, self.nt_mn, self.nt_mx])
                self.n_zooms += 1
            else:
                self.zooms[self.i_zooms][:] = [self.sxmin, self.sxmax, self.nt_mn, self.nt_mx]
            self.plotFile(self.axes[self.fig_plotted],self.fig_plotted)
        else:
            self.sxmin = xmi
            self.sxmax = xma
            self.i_zooms += 1
            if self.i_zooms >= self.n_zooms:
                self.zooms.append([self.sxmin, self.sxmax, self.nt_mn, self.nt_mx])
                self.n_zooms += 1
            else:
                self.zooms[self.i_zooms][:] = [self.sxmin, self.sxmax, self.nt_mn, self.nt_mx]
            self.plotShot(self.axes[self.fig_plotted],self.fig_plotted)
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
        self.picksPlot() # Plot picks if available
        self.setHelp(self.main_text) # Change help text to self module
        self.main.function = "main"
        zm = True
        return zm

    def zoomOut(self):
        """
        If several zooms have been used, zoom out one step
        """
        if self.i_zooms == 0:
            print("No zoom out possible, already at initial zoom")
            return False
        self.i_zooms = max(self.i_zooms-1,0)
        self.nt_mn = self.zooms[self.i_zooms][2] # Set actual zoom parameters to earlier window
        self.nt_mx = self.zooms[self.i_zooms][3]
        self.time_plt_min = self.data.t0+self.nt_mn*self.data.dt
        self.time_plt_max = self.data.t0+self.nt_mx*self.data.dt
        if self.i_zooms == 0:
            self.zooms[0][0] = min(self.x)-1
            self.zooms[0][1] = max(self.x)+1
        self.drawNew(True) # draw data section with actual zoom paramters
        self.setHelp(self.main_text) # Change help text to self module

    def zoomIn(self):
        """
        If several zooms have been used, zoom in one step
        """
        if self.i_zooms == self.n_zooms:
            print("No tighter zoom available, use Z to define new zoom")
            return False
        self.i_zooms = min(self.i_zooms+1,self.n_zooms-1)
        self.nt_mn = self.zooms[self.i_zooms][2] # Set actual zoom parameters to earlier window
        self.nt_mx = self.zooms[self.i_zooms][3]
        self.time_plt_min = self.data.t0+self.nt_mn*self.data.dt
        self.time_plt_max = self.data.t0+self.nt_mx*self.data.dt
        self.drawNew(True) # draw data section with actual zoom paramters
        self.setHelp(self.main_text) # Change help text to self module

    def zoomIni(self):
        """
        Go back to initial zoom (all data are plotted)
        """
        self.i_zooms = 0
        self.nt_mn = self.zooms[0][2]
        self.nt_mx = self.zooms[0][3]
        self.time_plt_min = self.data.t0
        self.time_plt_max = self.data.t0+self.nt_mx*self.data.dt
        self.zooms[0][0] = min(self.x)-1
        self.zooms[0][1] = max(self.x)+1
        self.drawNew(True) # draw full data section
        self.setHelp(self.main_text) # Change help text to self module

    def seismogram(self,ax,time,x_pos,data,nt_min=0,nt_max=0,traces=None,text_x="Distance [m]",
             text_y="Time [s]",text_t="",amp=0.667,fill=False,trace_norm=True):
        """
        plots a seismogram section

        Obligatory input:
        time (array[nsamples,]) time of each sample
        x_pos (array[ntraces,]) position of geophones
        data (array[nsamples,ntraces]) measured amplitudes

        Optional input:
        nt_min (=0): number of first sample to be plotted (counting starts at 0!)
        nt_max (=0): number of last sample to be plotted. If 0: plot all samples
        traces (array[ntraces] = None): list of traces to be plotted. If None: plot all traces
        text_x (string = "Distance [m]"): Text to be written on x-axis
        text_y (string = "Time [s]"): Text to be written on y-axis
        amp (=0.667): Amplitude of maximum value of a trace in units of x-axis
        fill (bool = False): Fill positive part of traces
        trace_norm (bool = True): Normalize each trace. If False, global normalization
        """
        try:
            self.actual_axis = ax
            t = time[nt_min:nt_max]
            ntrace = np.alen(x_pos)
            if ntrace == 0:
                _ = QtWidgets.QMessageBox.warning(None,"Warning",\
                    f"No data in {text_t}\nChoose other gather",
                    QtWidgets.QMessageBox.Close)
                return
            nsamp = np.alen(t)
            if nt_max == 0:
                nt_max = nsamp
            if traces == None or len(traces) == 0:
                traces = np.arange(ntrace,dtype='int')
            else:
                traces = np.array(traces,dtype='int')
            xmin = x_pos[traces].min()
            xmax = x_pos[traces].max()
            self.x_zoom_min = xmin
            self.x_zoom_max = xmax
            self.t_zoom_min = nt_min*self.data.dt+self.data.t0
            self.t_zoom_max = nt_max*self.data.dt+self.data.t0
            xmin -= 1
            xmax += 1
            d = np.copy(data[:,nt_min:nt_max])
            dmax = np.zeros(ntrace)
            if np.std(d) > 0.:
                if self.gain == "tnorm" or (self.gain=="dist" and self.dg_flag):
                    dmax = np.abs(d).max(axis=1)
                elif self.gain == "time":
                    fac = np.zeros(len(t))
                    fac[t>0] = t[t>0]**self.time_gain
                    for i in traces:
                        d[i,:] = d[i,:]*fac
                    dmax = np.abs(d).max(axis=0)
                elif self.gain == "dist":
                    fac = np.abs(x_pos)**self.dist_gain
                    for i in traces:
                        d[i,:] = d[i,:]*fac[i]
                    dmax = np.ones(ntrace)*np.max(np.abs(d))
                elif self.gain == "agc":
                    for i in traces:
                        d[i,:] = self.agcCalc(d[i,:],self.n_agc_win)
                    dmax = np.abs(d).max(axis=0)
            if np.isclose(dmax.max(), 0.):
                _ = QtWidgets.QMessageBox.warning(None,"Warning",\
                    "All traces are empty\n"+\
                    "Chose other gather",\
                    QtWidgets.QMessageBox.Close)
                return
            for i in traces:
                if np.isclose(dmax[i], 0.):
                    continue
                x = np.array(x_pos[i]+d[i,:]/dmax[i]*amp)
                x = np.clip(x,xmin,xmax)
                ax.plot(x,time[nt_min:nt_max])
                if fill:
                    nulls = np.ones_like(t)*x_pos[i]
                    ax.fill_betweenx(t,nulls,x,where=nulls<=x)
            ax.grid(b="on",which="major",axis="y")
            ax.set_xlabel(text_x)
            ax.set_ylabel(text_y)
            ax.set_title(text_t)
            return None
        except:
            print("Error in seismogram")
            return None

    def plotRG(self):
        """
        Prepare all plots for the different existing receiver gathers

        Returns
        -------
        None.

        """
# Set flags correctly
        self.main.function = "plot_RG"
        self.mplfigs.clear()
        self.v_set = False
        self.rg_flag = True
        self.sg_flag = False
        self.dg_flag = False
        self.fg_flag = False
        self.figs = []
        self.axes = []
        self.plot_names = []
        self.keys_held = set()
# Find all unique receiver numbers. For every single one, a plot gather is prepared
        self.receivers = np.sort(np.array(list(self.traces.rec_pt_dict.keys()),\
                         dtype = int))
# Loop over all found receivers
        for i in range(len(self.receivers)):
            self.figs.append(Figure())
            self.axes.append(self.figs[-1].subplots())
            if i==0:
# plot data of the first receiver to the screen
                self.zooms[self.i_zooms][0] = self.traces.off_min
                self.zooms[self.i_zooms][1] = self.traces.off_max
                self.plotReceiver(self.axes[-1],self.receivers[0])
# Add the figure to the list of figures in the right-side tool bar
            self.addFig(f"Receiver {self.receivers[i]+1}", self.figs[-1])
            self.plot_names.append(f"Receiver {self.receivers[i]+1}")
# Erase actual plot in central widget and plot the new one
        self.rmMPL()
        self.addMPL(self.figs[0])
        self.fig_plotted = 0
        self.fig_plotted_R = 0
# Plot picks
        self.picksPlot()
        if self.PlotCalculatedTimes.isChecked():
            self.plotCalcPicks()
        self.main.function = "main"

    def plotReceiver(self,ax,irec):
        self.time = self.data.t0+np.arange(int(self.data.nsamp))*self.data.dt
        self.x = []
        self.tr = []
        self.stdev = []
# Search traces having been recorded by receiver irec
        ntraces = len(np.where(self.traces.receiver == irec)[0])
        self.traces.plotted[:] = False
        if self.v_set != True:
            self.v = np.zeros((ntraces,self.data.nsamp))
        self.v_norm = np.zeros((ntraces,self.data.nsamp))
        self.actual_traces = []
        self.actual_number_traces = -1
# Loop over all traces of the receiver gather
        for i,isht in enumerate(self.traces.rec_pt_dict[irec]["shot"]):
            ifile = self.traces.rec_pt_dict[irec]["file"][i]
            itrace = self.traces.rec_pt_dict[irec]["trace"][i]
            ntr = self.traces.sht_rec_dict[(isht,irec)]
            self.actual_traces.append(ntr)
            self.traces.plotted[ntr] = True
            xx = -self.traces.offset[ntr]
            self.x.append(xx)
# If a trace exists already at the actual position, shift the first one by -0.2m,
#    the second one by +0.2m so they may be distinguished on the screen
            if i > 0:
                if self.x[-1] == self.x[-2]:
                    self.x[-1] += 0.2
                    self.x[-2] -= 0.2
            self.actual_number_traces += 1
# Copy data into array self.v
# v_set = True means that this work has already been done earlier
            if self.v_set != True:
                self.v[self.actual_number_traces,:] = \
                    self.data.st[ifile][itrace].data*\
                    self.traces.amplitudes[ntr]
                if self.nt_0>0:
                    self.v[self.actual_number_traces,:] -= \
                        np.mean(self.v[self.actual_number_traces,:self.nt_0])
# Normalize data by the standard deviation of each trace
            s = np.std(self.v[self.actual_number_traces,:])
            self.stdev.append(s)
            if s > 0:
                self.v_norm[self.actual_number_traces,:] = self.v[self.actual_number_traces,:]/s
            else:
                self.v_norm[self.actual_number_traces,:] = 0.
# If trace position is within the zoom limits, include it in list of traces to
#    be plotted (self.tr)
            if (self.x[-1] >= self.zooms[self.i_zooms][0]-0.2 and\
                self.x[-1] <= self.zooms[self.i_zooms][1]+0.2)\
                or self.i_zooms == 0:
                self.tr.append(i)
        self.x = np.array(self.x)
        self.actual_number_traces += 1
        self.indices = np.argsort(self.x)
        self.stdev = np.array(self.stdev)
# Seismogram does the plotting
        self.seismogram(ax,self.time,self.x,self.v,fill=True,\
                        amp=self.amp_plt,traces=self.tr,\
                        nt_min=self.nt_mn,nt_max=self.nt_mx,\
                        text_x="Shot offset (m)",
                        text_t=f"Receiver {irec+1}")
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
        self.main.function = "plot_DG"
        self.v_set = False
        self.mplfigs.clear()
        self.dg_flag = True
        self.sg_flag = False
        self.rg_flag = False
        self.fg_flag = False
        self.figs = []
        self.axes = []
        self.plot_names = []
        self.keys_held = set()

# Find all unique distance. For every single one, a plot gather is prepared
        self.distances = np.sort(np.unique(np.abs(self.traces.off_round)))
# Since typîcally, the preparation of distance gathers is rather time consuming,
#   activate a progress bar
        self.progressBar = QtWidgets.QProgressBar(self)
        self.mplvl.addWidget(self.progressBar)
        self.progressBar.show()
        self.progressBar.setValue(0)
# Loop over all found offsets
        for i in range(len(self.distances)):
# prepare a figure for each offset
            self.figs.append(Figure())
            self.axes.append(self.figs[-1].subplots())
            if i==0:
# plot the smallest found offset to the screen
                self.zooms[self.i_zooms][0] = np.min(self.traces.xcdp)
                self.zooms[self.i_zooms][1] = np.max(self.traces.xcdp)
                self.plotDistance(self.axes[-1],self.distances[i])
# Add the figure to the list of figures in the right-side tool bar
            self.addFig(f"Distance {self.distances[i]}", self.figs[-1])
            self.plot_names.append(f"Distance {self.distances[i]}")
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
        self.picksPlot()
        if self.PlotCalculatedTimes.isChecked():
            self.Plot_calc_picks()
        self.main.function = "main"

    def plotDistance(self,ax,dist):
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
        self.tr = []
        self.stdev=[]
# Loop over all traces of the distance gather
        for i in range(nrecplt):
            self.x.append(xx[i])
            self.actual_traces.append(traces[i])
            self.traces.plotted[traces[i]] = True
            self.actual_number_traces += 1
# If a trace exists already at the actual position (reciprocal positions of
#    shot and receiver), shift the first one by -0.2m, the second one by +0.2m
#    so they may be distinguished on the screen
            if i > 0:
                if self.x[-1] == self.x[-2]:
                    self.x[-1] += 0.2
                    self.x[-2] -= 0.2
# If trace position is within the zoom limits, include it in list of traces to
#    be plotted (self.tr)
            if (self.x[-1] >= self.zooms[self.i_zooms][0]-0.2 and\
                self.x[-1] <= self.zooms[self.i_zooms][1]+0.2)\
                or self.i_zooms == 0:
                    self.tr.append(i)
        self.x = np.array(self.x)
        ns = self.data.nsamp
# Copy data into array self.v
# v_set = True means that this work has already been done earlier
        if self.v_set != True:
            self.v = np.zeros((nrecplt,ns))
        self.v_norm = np.zeros((nrecplt,ns))
        for i in range(nrecplt):
            ifile = self.traces.file[traces[i]]
            itrace = self.traces.trace[traces[i]]
            if self.v_set != True:
                self.v[i,:] = self.data.st[ifile][itrace].data*\
                              self.traces.amplitudes[traces[i]]
                if self.nt_0>0:
                    self.v[i,:] -= np.mean(self.v[i,:self.nt_0])
# Normalize data by the standard deviation of each trace
            s = np.std(self.v[i,:])
            self.stdev.append(s)
            if s > 0:
                self.v_norm[i,:] = self.v[i,:]/s
            else:
                self.v_norm[i,:] = 0.
        self.indices = np.argsort(self.x)
        self.stdev = np.array(self.stdev)
# Seismogram does the plotting
        self.seismogram(ax,self.time,self.x,self.v,fill=True,\
                        amp=self.amp_plt,traces=self.tr,\
                        nt_min=self.nt_mn,nt_max=self.nt_mx,\
                        text_x="Midpoint position (m)",\
                        text_t=f"Distance {dist} m")
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
# Set flags correctly
        self.v_set = False
        self.fg_flag = False
        self.sg_flag = True
        self.rg_flag = False
        self.dg_flag = False
        self.mplfigs.clear()
        self.figs = []
        self.axes = []
        self.plot_names = []
        self.sh_plot = []
        self.stdev=[]
        self.keys_held = set()
# Loop over all shot points
        for i,sh in enumerate(self.traces.sht_pt_dict):
            self.sh_plot.append(sh)
            # print ("shot:", i, sh)
            self.figs.append(Figure())
            self.axes.append(self.figs[-1].subplots())
            self.traces.sht_pt_dict[sh]["axes"].append(i)
            if i == 0:
# Set zoom limits in x such that all traces are plotted to the screen
                self.zooms[self.i_zooms][0] = self.traces.off_min
                self.zooms[self.i_zooms][1] = self.traces.off_max
# Plot data of first shot point to the screen
                self.plotShot(self.axes[-1],i)
# Add the figure to the list of figures in the right-side tool bar
            self.addFig(f"Shot {sh+1}", self.figs[-1])
            self.plot_names.append(f"Shot {sh+1}")
# Erase actual plot in central widget and plot the new one
        self.rmMPL()
        self.addMPL(self.figs[isht])
        self.fig_plotted = 0
        self.fig_plotted_S = 0
# Plot picks
        self.picksPlot()
        if self.PlotCalculatedTimes.isChecked():
            self.plotCalcPicks()
        self.main.function = "main"

    def plotShot(self,ax,sh0):
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
        ntraces = len(np.where(self.traces.shot==sh)[0])
        if ntraces == 0:
            return
        self.x = []
        self.tr = []
        self.actual_traces = []
        self.traces.plotted[:] = False
        self.actual_number_traces = 0
        self.stdev = []
        if self.v_set != True:
            self.v = np.zeros((ntraces,ns))
        self.v_norm = np.zeros((ntraces,ns))
        self.actual_shot = sh
# Loop over all traces of the shot point gather
        for i,nt in enumerate(self.traces.sht_pt_dict[sh]["trace"]):
            ifile = self.traces.sht_pt_dict[sh]["file"][i]
            irec = self.traces.sht_pt_dict[sh]["receiver"][i]
            ntr = self.traces.sht_rec_dict[(sh,irec)]
            self.actual_traces.append(ntr)
            self.traces.plotted[ntr] = True
            self.actual_number_traces += 1
            xx = self.traces.offset[ntr]
            self.x.append(xx)
# If a trace exists already at the actual position, shift the first one by -0.2m,
#    the second one by +0.2m so they may be distinguished on the screen
            if i > 0:
                if self.x[-1] == self.x[-2]:
                    self.x[-1] += 0.2
                    self.x[-2] -= 0.2
# Copy data into array self.v
# v_set = True means that this work has already been done earlier
            if self.v_set != True:
                self.v[i,:] = self.data.st[ifile][nt].data*\
                              self.traces.amplitudes[ntr]
                if self.nt_0>0:
                    self.v[i,:] -= np.mean(self.v[i,:self.nt_0])
# Normalize data by the standard deviation of each trace
            s = np.std(self.v[i,:])
            self.stdev.append(s)
            if s > 0:
                self.v_norm[i,:] = self.v[i,:]/s
            else:
                self.v_norm[i,:] = 0.
# If trace position is within the zoom limits, include it in list of traces to
#    be plotted (self.tr)
            if (self.x[-1] >= self.zooms[self.i_zooms][0]-0.2 and\
                self.x[-1] <= self.zooms[self.i_zooms][1]+0.2)\
                or self.i_zooms == 0:
                self.tr.append(i)
        self.x = np.array(self.x)
        self.indices = np.argsort(self.x)
        self.stdev = np.array(self.stdev)
# Seismogram does the plotting
        self.seismogram(ax, self.data.time, self.x, self.v, fill=True,\
                        amp=self.amp_plt, traces=self.tr,\
                        nt_min=self.nt_mn, nt_max=self.nt_mx,\
                        text_x="Offset (m)", text_t=f"Shot point {sh+1}")
# The following lines are only for testing purpose, to show potential picks from
# false colour plots
        # try:
        #     self.axes[self.fig_plotted].plot(self.x,self.main.utilities.pk,"m")
        #     self.axes[self.fig_plotted].set_ylim([self.time_plt_min, self.time_plt_max])
        #     renderer = self.figs[self.fig_plotted].canvas.renderer
        #     self.axes[self.fig_plotted].draw(renderer)
        #     self.canvas.blit(self.axes[self.fig_plotted].bbox)
        # except:
        #     pass
# Change_colors is only meant for plot of tomography result, so deactivate it
        self.Change_colors.setEnabled(False)
        self.main.function = "main"

    def plotFG(self,isht=0):
        """
        Prepare all plots for the different files

        Returns
        -------
        None.

        """
# Set flags correctly
        self.v_set = False
        self.sg_flag = False
        self.rg_flag = False
        self.dg_flag = False
        self.fg_flag = True
        self.mplfigs.clear()
        self.figs = []
        self.axes = []
        self.plot_names = []
        self.stdev=[]
        self.keys_held = set()
# Loop over all recorded files
        for i in range(self.files.file_count):
            self.figs.append(Figure())
            self.axes.append(self.figs[-1].subplots())
            if i==0:
# Set zoom limits in x such that all traces are plotted to the screen
                self.zooms[self.i_zooms][0] = self.zooms[0][0]
                self.zooms[self.i_zooms][1] = self.zooms[0][1]
# Plot data of first file to the screen
                self.plotFile(self.axes[-1],i)
# Add the figure to the list of figures in the right-side tool bar
            self.addFig(f"File {self.files.numbers[i]}", self.figs[-1])
            self.plot_names.append(f"File {self.files.numbers[i]}")
# Erase actual plot in central widget and plot the new one
        self.rmMPL()
        self.addMPL(self.figs[isht])
        self.fig_plotted = 0
        self.fig_plotted_F = 0
# Plot picks
        self.picksPlot()
        if self.PlotCalculatedTimes.isChecked():
            self.plotCalcPicks()
        self.main.function = "main"

    def plotFile(self,ax,isht):
        self.x = []
        self.tr = []
        self.actual_traces = []
        self.traces.plotted[:] = False
        self.actual_number_traces = 0
        self.stdev = []
# Set a few plotting parameters
        ntr = len(self.data.st[isht])
        nsamp = len(self.data.st[isht][0].data)
        sht = int(self.data.st[isht][0].stats.seg2['SOURCE_STATION_NUMBER'])

        if self.v_set != True:
            self.v = np.zeros((ntr,nsamp))
        self.v_norm = np.zeros((ntr,nsamp))
        self.actual_shot = self.files.numbers[isht]
# Loop over all traces of the file
        for i,t in enumerate(self.files.file_dict[isht]["traces"]):
            self.actual_traces.append(t)
            self.traces.plotted[t] = True
            self.actual_number_traces += 1
            xx = self.traces.offset[t]
            self.x.append(xx)
# Copy data into array self.v
# v_set = True means that this work has already been done earlier
            if self.v_set != True:
                self.v[i,:] = self.data.st[isht][i].data*\
                              self.traces.amplitudes[t]
                if self.nt_0>0:
                    self.v[i,:] -= np.mean(self.v[i,:self.nt_0])
# Normalize data by the standard deviation of each trace
            s = np.std(self.v[i,:])
            self.stdev.append(s)
            if s > 0:
                self.v_norm[i,:] = self.v[i,:]/s
            else:
                self.v_norm[i,:] = 0.
# If trace position is within the zoom limits, include it in list of traces to
#    be plotted (self.tr)
            if (xx>=self.zooms[self.i_zooms][0] and xx<=self.zooms[self.i_zooms][1])\
                or self.i_zooms == 0:
                self.tr.append(i)
        self.x = np.array(self.x)
        self.indices = np.argsort(self.x)
        self.stdev = np.array(self.stdev)
# Seismogram does the plotting
        self.seismogram(ax, self.data.time, self.x, self.v, fill=True,\
                        amp=self.amp_plt, traces=self.tr,\
                        nt_min=self.nt_mn, nt_max=self.nt_mx,\
                        text_x="Offset (m)",\
                        text_t=f"File {self.files.numbers[isht]}, "+\
                               f"shot point {sht}")
# Change_colors is only meant for plot of tomography result, so deactivate it
        self.Change_colors.setEnabled(False)
        self.main.function = "main"

    def set_ticks(self,minval,maxval,ntick=5,dtick=5.):
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
        d_ticks = (maxval-minval)/ntick
        mult = round(d_ticks/dtick)
        mult = max(mult, 1)
        d_ticks = dtick*mult
        min_tick = np.ceil(minval/d_ticks)*d_ticks
        max_tick = int(maxval/d_ticks)*d_ticks
        ticks = np.arange(min_tick,max_tick+d_ticks/2,d_ticks)
        return ticks

    def nextBlock(self,direction):
        if self.i_zooms == 0:
            return
        d = int(direction)
        dwin = self.zooms[self.i_zooms][1]-self.zooms[self.i_zooms][0]
        if d > 0:
            xma = max(self.x)
            if self.zooms[self.i_zooms][1] >= xma:
                return
            zm1 = min(self.zooms[self.i_zooms][1]+dwin-1,xma+0.5)
            zm0 = zm1 - dwin
        else:
            xmi = min(self.x)
            if self.zooms[self.i_zooms][0] <= xmi:
                return
            zm0 = max(self.zooms[self.i_zooms][0]-dwin+1,xmi-0.5)
            zm1 = zm0 + dwin
        self.zooms[self.i_zooms][0] = zm0
        self.zooms[self.i_zooms][1] = zm1
        self.drawNew(True) # draw data section with actual zoom paramters
        self.setHelp(self.main_text) # Change help text to self module
        return True

    def picksPlot(self, col="r"):
        """
        Plot all measured picks into the actual gather
        Uses function pickPlot to plot each pick found

        Parameters:
        ----------
        col: str
            color to be used for pcik plot. Default : red ("r")

        Returns
        -------
        None.

        """
        if np.sum(np.array(self.traces.npick, dtype=int)[self.actual_traces]) > 0:
            for i in range(self.actual_number_traces):
                ntr = self.actual_traces[i]
                if self.traces.npick[ntr]>0:
                    for j in range(self.traces.npick[ntr]):
                        self.pickPlot(i,ntr,j,col)
            renderer = self.figs[self.fig_plotted].canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.PlotPicks.setEnabled(True)
            self.MovePicks.setEnabled(True)
            self.Tomography.setEnabled(True)
        else:
            pass

    def pickPlot(self,i,ntr,j,col):
        """
        Plot one single pick

        Parameters
        ----------
        i : TYPE int
            number of trace in gather
        ntr: TYPE int
            number of trace in list of all traces
        j : TYPE int
            number of pick for that trace
        col : TYPE string
            Color to be used for plot

        Returns
        -------
        None.

        """
        t = self.traces.pick_times[ntr][j]
        if  t >= self.time_plt_min and t <= self.time_plt_max and\
            self.x[i]>=self.zooms[self.i_zooms][0] and\
            self.x[i]<=self.zooms[self.i_zooms][1]:
            x_c = []
            x_c.append(self.x[i]-0.3)
            x_c.append(self.x[i]+0.3)
            y_c = [t,t]
            self.line, = self.axes[self.fig_plotted].plot(x_c, y_c,col)
# If plot pick is called before the uncertainty of a travel time has been
#    picked, i.e. before the picking has been finished, add standard uncertainty
            tmn = self.traces.pick_times_min[ntr][j]
            try:
                tmx = self.traces.pick_times_max[ntr][j]
            except:
                self.traces.pick_times_max[ntr][j] = t-2*self.data.dt
                tmn = self.traces.pick_times_max[ntr][j]
            x_c = [self.x[i],self.x[i]]
            y_c = [max(tmn,self.time_plt_min),min(tmx,self.time_plt_max)]
            self.line, = self.axes[self.fig_plotted].plot(x_c, y_c,col)

    def plotCalcPicks(self):
        """
        Plot calculated travel times, see function Traces.readCalcPicks for
        further explanations

        Returns
        -------
        None.

        """
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
# Next lines are important for file gather to avoid picks plotted outside of the
# plotted traces
        xc = np.array(xc)
        tc = np.array(tc)
        tc = tc[(xc<=max(self.x)) & (xc>=min(self.x))]
        xc = xc[(xc<=max(self.x)) & (xc>=min(self.x))]
        nsort = np.argsort(xc)
        xx = xc[nsort]
        yy = tc[nsort]
# The following try-except construction avoids problems in case no calculated
# travel times are found for a gather. This is mainly the case for distance
# gathers at zero distance, since there, no times are calculated.
        l = len(xx)
        try:
            if l < 1:
                print("\nNo calculated travel times available")
            else:
                xmn = self.zooms[self.i_zooms][0]
                xmx = self.zooms[self.i_zooms][1]
                xxx = xx[xx>=xmn]
                yyy = yy[xx>=xmn]
                xx = xxx[xxx<=xmx]
                yy = yyy[xxx<=xmx]
                self.line, = self.axes[self.fig_plotted].plot(xx,yy,"g")
                self.axes[self.fig_plotted].set_ylim([self.time_plt_min, self.time_plt_max])
                renderer = self.figs[self.fig_plotted].canvas.renderer
                self.axes[self.fig_plotted].draw(renderer)
                self.canvas.blit(self.axes[self.fig_plotted].bbox)
        except:
            if l < 1:
                print("\nNo calculated travel times available")
            else:
                print(f"\nProblems plotting {l} calculated travel times")
            pass

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
        dmin = 1e10
        self.n_tr = 0
        for i in range(self.actual_number_traces):
            dist = np.abs(self.x[i]-x)
            if dist<dmin:
                dmin = dist
                self.n_tr = i
        self.itrace = self.actual_traces[self.n_tr]
        self.p_s = self.traces.shot[self.itrace]
        self.p_r = self.traces.receiver[self.itrace]

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
        except:
            itr = -1
            for i,it in enumerate(self.actual_traces):
                if self.traces.npick[it]>0:
                    itr = it
                    break
            if itr < 0:
                _ = QtWidgets.QMessageBox.warning(None, "Warning",
                        "No pick exists for this gather.\n",\
                        QtWidgets.QMessageBox.Close)
                return False
        self.searchTrace(x)
# If there are several picks for the trace, search pick nearest to time of mouse position
        if self.traces.npick[self.itrace]>0:
            dmin = np.abs(self.traces.pick_times[self.itrace][0]-t)
            self.ipk = 0
            if self.traces.npick[self.itrace]>1:
                for i in range(1,self.traces.npick[self.itrace]):
                    dist = np.abs(self.traces.pick_times[self.itrace][i]-t)
                    if dist < dmin:
                        dmin = dist
                        self.ipk = i
        else:
# If clicked trace has no pick, stay on last trace
            self.itrace = itr
            self.ipk = 0
        return True

    def pickMove(self):
        """
        Routine waits for left mouse click, searches pick nearest to
        mouse position and marks this pick in green color. Then use
        keyboard buttons up and down, with or without SHFT to move
        selected pick (see routine movePick).
        Right mouse click finishes routine

        Returns
        -------
        None.

        """
        global figure,r_flag
        self.setHelp(self.pick_move_text)
        self.main.function = "pick_move"
        def onPress(event):
# If left öouse button has been clicked, search nearest pick
            if event.button == 1:
                self.keys_held = set()
                result = self.searchPick(event.xdata,event.ydata)
                if not result:
                    self.keys_held = set()
                    self.drawNew(True)
                    self.picksPlot()
                    self.end = True
                    return
                self.update()
# Save position of chosen pick in memory
                time_back = self.traces.pick_times[self.itrace][self.ipk]
                # print("Move pick nr %d at shot %d, receiver %d" %\
                #       (self.ipk, sht+1, stn+1))
# Re-plot pick in white color (i.e. background color)
                self.pickPlot(self.n_tr, self.itrace, self.ipk, "w")
# Erase pick temporally to redraw canvas without this pick
                self.traces.pick_times[self.itrace][self.ipk] = -1
                self.picksPlot()
                self.back = figure.canvas.copy_from_bbox(figure.bbox)
                self.axes[self.fig_plotted].figure.canvas.draw()
                self.traces.pick_times[self.itrace][self.ipk] = time_back
                self.x_c = []
                self.x_c.append(self.x[self.n_tr]-0.3)
                self.x_c.append(self.x[self.n_tr]+0.3)
                self.y_c = [time_back,time_back]
                self.lin, = self.axes[self.fig_plotted].plot(self.x_c, self.y_c,"g")
                self.axes[self.fig_plotted].draw_artist(self.lin)
                self.canvas.blit(self.axes[self.fig_plotted].bbox)
            elif event.button == 3:
# If right mouse has been clicked, plot all picks, including the moved
# one in red and leave routine
                self.keys_held = set()
                self.drawNew(True)
                self.picksPlot()
                self.end = True
        self.end = False
        self.coor_x = []
        self.coor_y = []
        figure = self.figs[self.fig_plotted]

# Plot lines in positive and negative direction corresponding to air wave arrivals
        xx = np.arange(self.x_zoom_min,self.x_zoom_max)
        tt = np.abs(xx)/345.
        xxx = xx[tt>=self.t_zoom_min]
        ttt = tt[tt>=self.t_zoom_min]
        plt_air = True
        if self.dg_flag:
            dist = self.distances[self.fig_plotted_D]
            x_air = xx
            t_air = np.ones_like(x_air)*dist/345.
            if t_air[0] > self.t_zoom_max or t_air[0]<self.t_zoom_min:
                plt_air = False
        else:
            x_air = xxx[ttt<=self.t_zoom_max]
            t_air = ttt[ttt<=self.t_zoom_max]
            if len(t_air) < 2:
                plt_air = False
        if plt_air:
            self.axes[self.fig_plotted].plot(x_air,t_air,ls="--")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)

        self.lin, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y, animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect('button_press_event',\
                            onPress)
        while (self.end != True):
            QtCore.QCoreApplication.processEvents()
        self.setHelp(self.main_text)

    def movePick(self, sign):
        """
        Routine moves selected pick (from routine Pick_move) by one or
        ten samples.

        Parameters
        ----------
        sign : float
            Direction of movement and step length. Comes from
            Interactive_routines.py -> on_key_down
            Value is +1  if upward key has been pressed
                         (pick time increses by one sample)
            Value is -1  if downward key has been pressed
                         (pick time decreses by one sample)
            Value is +10 if SHFT+upward key has been pressed
                         (pick time increses by ten samples)
            Value is +10 if SHFT+downward key has been pressed
                         (pick time decreses by ten samples)
            Value is +/- samples in 1 ms if CTRL has been pressed

        Returns
        -------
        None.

        """
        figure.canvas.restore_region(self.back)
        ds = sign*self.data.dt
        self.y_c[0] += ds
        self.y_c[1] += ds
        self.traces.pick_times_min[self.itrace][self.ipk] += ds
        self.traces.pick_times_max[self.itrace][self.ipk] += ds
        self.traces.pick_times[self.itrace][self.ipk] += ds
        self.lin.set_data(self.x_c, self.y_c)
        self.axes[self.fig_plotted].draw_artist(self.lin)
        self.canvas.blit(self.axes[self.fig_plotted].bbox)

    def uncertainty(self):
        """
        Routine waits for left mouse click, searches pick nearest to
        mouse position and marks this pick in green color. Then use
        keyboard buttons up and down, with or without SHFT to change
        uncerteinty of selected pick (see routine Unc_change).
        Right mouse click finishes routine

        Returns
        -------
        None.

        """
        global figure,r_flag
        self.setHelp(self.uncertainty_text)
        self.main.function = "uncertainty"
        def onPress(event):
# If left öouse button has been clicked, search nearest pick
            if event.button == 1:
                self.keys_held = set()
                self.searchPick(event.xdata,event.ydata)
                self.update()
# Save position of chosen pick in memory
                trace = self.itrace
                time_back = self.traces.pick_times[self.itrace][self.ipk]
# Erase pick temporally to redraw canvas without this pick
                self.traces.pick_times[self.itrace][self.ipk] = -1
                self.picksPlot()
                self.back = figure.canvas.copy_from_bbox(figure.bbox)
                self.axes[self.fig_plotted].figure.canvas.draw()
                self.traces.pick_times[self.itrace][self.ipk] = time_back
                self.x_c = [self.x[self.n_tr],self.x[self.n_tr]]
                self.y_c = [self.traces.pick_times_max[trace][self.ipk],\
                            self.traces.pick_times_min[trace][self.ipk]]
                self.lin, = self.axes[self.fig_plotted].plot(self.x_c, self.y_c,"b")
                self.axes[self.fig_plotted].draw_artist(self.lin)
                self.canvas.blit(self.axes[self.fig_plotted].bbox)
            elif event.button == 3:
# If right mouse has been clicked, plot all picks, including the moved
# one in red and leave routine
                self.keys_held = set()
                self.drawNew(True)
                self.picksPlot()
                self.end = True
        self.end = False
        self.coor_x = []
        self.coor_y = []
        figure = self.figs[self.fig_plotted]
        self.lin, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y,\
                                                     animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect('button_press_event',\
                                                           onPress)
        while (self.end != True):
            QtCore.QCoreApplication.processEvents()
        self.setHelp(self.main_text)

    def changeUnc(self,sign):
        """
        Routine increases or decreases uncertainty of pick by one or
        ten samples or one millisecond.

        Parameters
        ----------
        sign : float
            Direction of movement and step length. Comes from
            Interactive_routines.py -> on_key_down
            Value is +1  if upward key has been pressed
                         (uncertainty increses by one sample)
            Value is -1  if downward key has been pressed
                         (uncertainty decreses by one sample)
            Value is +10 if SHFT+upward key has been pressed
                         (uncertainty increses by ten samples)
            Value is +10 if SHFT+downward key has been pressed
                         (uncertainty decreses by ten samples)
            Value is +/- samples in 1 ms if CTRL has been pressed

        Returns
        -------
        None.

        """
        figure.canvas.restore_region(self.back)
        ds = sign*self.data.dt
        self.y_c[0] += ds
        self.y_c[1] -= ds
        self.traces.pick_times_min[self.itrace][self.ipk] -= ds
        self.traces.pick_times_max[self.itrace][self.ipk] += ds
        self.lin.set_data(self.x_c, self.y_c)
        self.axes[self.fig_plotted].draw_artist(self.lin)
        self.canvas.blit(self.axes[self.fig_plotted].bbox)

    def shiftPickTrace(self,direction):
        """
        Go to next trace to the left or to the right, depending on "direction".
        Function is called from self.Interactive_routines.On_key_down.

        Parameters
        ----------
        direction : int
            If positive, go one trace to the right, if negative go one trace to
            the left.

        Returns
        -------
        None.

        """
        global figure, r_flag
        self.d_time = 0
        self.keys_held = set()
        self.pickPlot(self.n_tr, self.itrace, self.ipk, "r")
        nt = self.n_tr
        while True:
            self.n_tr += int(direction)
# Avoid going left of first trace
            if self.n_tr < 0:
                self.n_tr = 0
# Avoid going right of last trace
            if self.n_tr == len(self.stdev):
                self.n_tr = len(self.stdev)-1
# If standard deviation of samples is zero (no data) skip trace
            if self.stdev[self.n_tr] > 0:
                break
        while True:
            self.itrace = self.actual_traces[self.n_tr]
            if self.traces.npick[self.itrace] > 0:
                break
            self.n_tr += int(direction)
            if ((self.n_tr==len(self.actual_traces) and int(direction)>0) or\
                (self.n_tr<0 and int(direction)<0)):
                self.n_tr = nt
                self.itrace = self.actual_traces[self.n_tr]
                break
        self.ipk = 0
        self.update()
# Save position of chosen pick in memory
        trace = self.itrace
        time_back = self.traces.pick_times[trace][self.ipk]
# Re-plot pick in white color (i.e. background color)
        self.pickPlot(self.n_tr, self.itrace, self.ipk, "w")
# Erase pick temporally to redraw canvas without this pick
        self.traces.pick_times[trace][self.ipk] = -1
        self.picksPlot()
        self.back = figure.canvas.copy_from_bbox(figure.bbox)
        self.axes[self.fig_plotted].figure.canvas.draw()
        self.traces.pick_times[trace][self.ipk] = time_back
        if self.main.function == "pick_move":
# If move picks is active, plot only horizontal line of pick position
            self.x_c = []
            self.x_c.append(self.x[self.indices[self.n_tr]]-0.3)
            self.x_c.append(self.x[self.indices[self.n_tr]]+0.3)
            self.y_c = [self.traces.pick_times[trace][self.ipk],\
                        self.traces.pick_times[trace][self.ipk]]
        else:
# If change uncertainty is active, plot vertical line of pick pick, i.e. uncertainty
            self.x_c = [self.x[self.indices[self.n_tr]],\
                        self.x[self.indices[self.n_tr]]]
            self.y_c = [self.traces.pick_times_max[trace][self.ipk],\
                        self.traces.pick_times_min[trace][self.ipk]]
        self.lin, = self.axes[self.fig_plotted].plot(self.x_c, self.y_c,"b")
        self.axes[self.fig_plotted].draw_artist(self.lin)
        self.canvas.blit(self.axes[self.fig_plotted].bbox)

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
            the picked traveltime measured for the offset abs(x_pk_near) in nearby
            gather
        i_tr_near : numpy int
            1D array with maximum length = number of traces on screen. Contains
            the picked traveltime measured for the offset abs(x_pk_near) in nearby
            gather

        If in other gathers no pick has been measured for the offset
        corresponding to a trace on the screen, the arrays do not contain any
        value for this trace.
        """
        import numpy as np
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
            for s in similar[0]:
                if self.traces.npick[s] < 1:
                    continue
                if s == trace:
                    continue
                sim_t.append(s)
                sim_s.append(self.traces.shot[s])
                sim_r.append(self.traces.receiver[s])
                sim_cdp.append(self.traces.xcdp[s])
            if len(sim_t) < 1:
                continue
            sim_t = np.array(sim_t, dtype = int)
            sim_s = np.array(sim_s, dtype = int)
            sim_r = np.array(sim_r, dtype = int)
            sim_cdp = np.array(sim_cdp, dtype = float)
            if self.sg_flag or self.fg_flag:
#                itr = np.argsort(abs(sim_x-self.geom.sht_dict[p_s]["x"]))
                itr = np.argsort(abs(sim_cdp-self.traces.xcdp[trace]))
                t_pk = None
                for i in range(min(len(itr),6)):
                    try:
                        t_pk = self.traces.pick_times[sim_t[itr[i]]][0]
                        d = abs(self.traces.xcdp[sim_t[itr[i]]]-self.traces.xcdp[trace])
                        break
                    except:
                        continue
                if not t_pk or d>10:
                    continue
                t_pk_near.append(t_pk)
                x_pk_near.append(self.x[n_tr])
                i_tr_near.append(n_tr)
            elif self.rg_flag:
                itr = np.argsort(np.abs(sim_r-p_r))
                t_pk =  None
                for i in range(min(len(itr),6)):
                    try:
                        t_pk = self.traces.pick_times[sim_t[i]][0]
                        break
                    except:
                        continue
                if not t_pk:
                    continue
                t_pk_near.append(t_pk)
                x_pk_near.append(self.x[n_tr])
                i_tr_near.append(n_tr)
        return np.array(x_pk_near),np.array(t_pk_near),np.array(i_tr_near,dtype=int)

    def erase_Picks(self):
        """
        Erases all picks of the actually shown gather

        Returns
        -------
        None.

        """
        for i in range(self.actual_number_traces):
            ipt = self.actual_traces[i]
            self.traces.npick[ipt] = 0
        self.drawNew(True)
        self.setHelp(self.main_text)

    def pickManual(self):
        """
        Manual picking

        Press left mouse button to define the position of a pick. Then click again left mouse button
            to define the uncertainty range around the pick.
        Define as many picks a necessary (maximum 5 picks on every trace)

        Click central mouse button (wheel) to erase a pick

        Click right mouse button to finish manual picking
        """
        global figure
        self.setHelp(self.pick_manual_text)
        self.main.function = "mpick"
        self.press = 0
        self.finish = False
        if self.data.filtered:
            umin = 4.*self.data.dt
        else:
            umin = 2.*self.data.dt
        def onPress(event):
            global figure, tpk, npi
            if event.button == 1:
# Left mouse button pressed, define a new pick
                if self.press == 0:
# First click at the position of the pick
                    self.press = 1
                    tpk = event.ydata
                    self.searchTrace(event.xdata)
                    trace = self.itrace
                    sht = self.p_s
                    stn = self.p_r
                    npi = self.traces.npick[trace]
                    self.traces.npick[trace] += 1
                    print(f"Pick nr {self.traces.npick[trace]} at shot "+\
                          f"{sht+1}, receiver {stn+1}")
                    try:
                        self.traces.pick_times[trace][npi] = tpk
                    except:
                        self.traces.pick_times[trace].append(tpk)
                    self.canvas = self.lin.figure.canvas
                    x_coor = [self.x[self.n_tr]-0.3, self.x[self.n_tr]+0.3]
                    y_coor = [tpk, tpk]
                    self.lin, = self.axes[self.fig_plotted].plot(x_coor,y_coor,"r")
                    self.axes[self.fig_plotted].draw_artist(self.lin)
                    self.canvas.blit(self.axes[self.fig_plotted].bbox)
                else:
                    self.press = 0
# Second, press nearby to define the uncertainty range
                    trace = self.itrace
                    sht = self.p_s
                    stn = self.p_r
                    npi = self.traces.npick[trace]-1
                    tpk = self.traces.pick_times[trace][npi]
                    unc = max(np.abs(tpk - event.ydata),umin)
                    try:
                        self.traces.pick_times_min[trace][npi] = tpk-unc
                        self.traces.pick_times_max[trace][npi] = tpk+unc
                    except:
                        self.traces.pick_times_min[trace].append(tpk-unc)
                        self.traces.pick_times_max[trace].append(tpk+unc)
                    x_coor = [self.x[self.n_tr], self.x[self.n_tr]]
                    y_coor = [self.traces.pick_times_min[trace][npi],\
                              self.traces.pick_times_max[trace][npi]]
                    self.lin, = self.axes[self.fig_plotted].plot(x_coor,y_coor,"r")
                    self.axes[self.fig_plotted].draw_artist(self.lin)
                    self.canvas.blit(self.axes[self.fig_plotted].bbox)
            elif event.button == 2:
                self.searchPick(event.xdata, event.ydata)
                trace = self.itrace
                sht = self.p_s
                stn = self.p_r
                npi = self.traces.npick[trace]
                print(f"Erase pick nr {self.ipk+1} at shot {sht+1},"+\
                      f"receiver {stn*1}")
                xx = self.x[self.n_tr]
                x_c = [xx-0.3,xx+0.3,xx,xx,xx]
                yy = self.traces.pick_times[trace][self.ipk]
                yy1 = self.traces.pick_times_min[trace][self.ipk]
                yy2 = self.traces.pick_times_max[trace][self.ipk]
                y_c = [yy,yy,yy,yy1,yy2]
# TO DO: I did not yet find a way to simply erase the drawing of a pick (except
#        redrawing). Therefore, for the moment, the pick is only hidden, plotting
#        it in white
                self.lin, = self.axes[self.fig_plotted].plot(x_c,y_c,"w")
                self.axes[self.fig_plotted].draw_artist(self.lin)
                self.canvas.blit(self.axes[self.fig_plotted].bbox)
                if npi>1:
                    if self.ipk<npi-1:
                        for i in range(self.ipk+1,npi):
                            self.traces.pick_times[trace][i-1]=\
                                self.traces.pick_times[trace][i]
                            self.traces.pick_times_min[trace][i-1] =\
                                self.traces.pick_times_min[trace][i]
                            self.traces.pick_times_max[trace][i-1] =\
                                self.traces.pick_times_max[trace][i]
                    self.traces.npick[trace] -= 1
                elif npi == 0:
                    return
                else:
                    self.traces.npick[trace] = 0
            else:
                if self.press == 1:
                    npi = self.traces.npick[self.itrace]-1
                    tpk = self.traces.pick_times[self.itrace][npi]
                    try:
                        self.traces.pick_times_min[self.itrace][npi] = tpk-umin
                        self.traces.pick_times_max[self.itrace][npi] = tpk+umin
                    except:
                        self.traces.pick_times_min[self.itrace].append(tpk-umin)
                        self.traces.pick_times_max[self.itrace].append(tpk+umin)
                self.lin.set_animated(False)
                self.finish = True
                self.lin.figure.canvas.mpl_disconnect(self.cidpress)
                self.drawNew(True)
                self.setHelp(self.main_text) # Change help text to main module
                self.traces.storePicks()
                self.end = True

        self.end = False
        figure = self.figs[self.fig_plotted]
        xpk=[]
        tpk=[]
# Plot lines in positive and negative direction corresponding to air wave arrivals
        xx = np.arange(self.x_zoom_min,self.x_zoom_max)
        tt = np.abs(xx)/345.
        xxx = xx[tt>=self.t_zoom_min]
        ttt = tt[tt>=self.t_zoom_min]
        plt_air = True
        if self.dg_flag:
            dist = self.distances[self.fig_plotted_D]
            x_air = xx
            t_air = np.ones_like(x_air)*dist/345.
            if t_air[0] > self.t_zoom_max or t_air[0]<self.t_zoom_min:
                plt_air = False
        else:
            x_air = xxx[ttt<=self.t_zoom_max]
            t_air = ttt[ttt<=self.t_zoom_max]
            if len(t_air) < 2:
                plt_air = False
        if plt_air:
            self.axes[self.fig_plotted].plot(x_air, t_air, c="b", ls="--")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)

# Search picks done earlier in nearby gathers an plot a line connecting them
        x_pk_help, t_pk_help, i_pk_help = self.searchNearPicks()
        x1 = x_pk_help[(x_pk_help >= self.x_zoom_min) & \
                       (x_pk_help <= self.x_zoom_max)]
        t1 = t_pk_help[(x_pk_help >= self.x_zoom_min) &\
                       (x_pk_help <= self.x_zoom_max)]
        x_pk = x1[x1 <= self.x_zoom_max]
        t_pk = t1[x1 <= self.x_zoom_max]
        if len(x_pk):
            self.axes[self.fig_plotted].plot(x_pk, t_pk, c="g", ls="--")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)

# Prepare self.lin to plot picks into.
        self.lin, = self.axes[self.fig_plotted].plot(xpk, tpk, animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect('button_press_event',\
                                                onPress)
        while (self.end != True):
            QtCore.QCoreApplication.processEvents()
            self.Tomography.setEnabled(True)
        self.PlotPicks.setEnabled(True)
        self.MovePicks.setEnabled(True)
        self.setHelp(self.main_text)

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
        the reference trace are not the full trace but only a length of +/-20 ms
        (as for the reference trace) around a guessed time. There are two
        possible ways of guessing this position:
            If there are already picks from other shots nearby, search for the
            nearest pick having the same offset (supposing that the arrival times
            change only slowly from one shot point to the next)
            If there are no picks available, take the same time as found for the
            last pick (initially, the pick of the reference trace, then the one
            from the last trace treated). For small offsets (actually up to <5m,
            very low apparent velocities are normal, which is why there, the
            times are change by +/-5ms depending on whether the offset increases
            or decreases. If this should be changed, search
            the following lines:
                off_lim = 5
                dt_pick = 0.005
            and replace off_lim by the distance limit (in meters) up to which
            very low velocities are supposed and dt_pick by the time shift per
            trace (in seconds) for thes short offsets.
        The picks are only searched for the traces visible in the actual zoom.
        This may be convenient if data are noisy such that the user may first
        pick near-offset traces, then filter the data and pick far-offsets.
        Uncertainty limits are fixed to 2 samples for raw data and to 4 samples
        for filtered data. If you want to change this do it in the if structure:
            if self.filtered:
                umin = 4*self.dt
            else:
                umin = 2*self.dt

        Returns
        -------
        None.

        """
        self.main.function = "cpick"
        global figure,r_flag
        self.setHelp(self.cpick_text)

        def onPress(event):
            """
            Come here if a mouse button has been pressed
            """
            global figure, r_flag
            self.d_time = 0
# Set some default values explained further up in the explanation of the function
            off_lim = 5
            dt_pick = 0.005
            t_half = 0.01
            if self.data.filtered:
                umin = 4*self.data.dt
            else:
                umin = 2*self.data.dt
# If the left mouse button was pressed (manual pick set), come here
            if event.button == 1:
# Search the trace where the pick has been defined
                self.searchTrace(event.xdata)
                trace = self.itrace
                sht = self.traces.shot[trace]
                stn = self.traces.receiver[trace]
                print(f"shot {sht}, receiver {stn}")
# Set the reference time the the pick position
                time_ref = event.ydata
                self.traces.npick[trace] = 1
                try:
                    self.traces.pick_times[trace][0] = time_ref
                    self.traces.pick_times_min[trace][0] = \
                        self.traces.pick_times[trace][0]-umin
                    self.traces.pick_times_max[trace][0] = \
                        self.traces.pick_times[trace][0]+umin
                except:
                    self.traces.pick_times[trace].append(time_ref)
                    self.traces.pick_times_min[trace].append(\
                        self.traces.pick_times[trace][0]-umin)
                    self.traces.pick_times_max[trace].append(\
                        self.traces.pick_times[trace][0]+umin)
                n_pick_ref = int(np.round((time_ref-self.data.t0)/self.data.dt,0))
                trace_ref = self.n_tr
# Calculate the width of the zone to find maxima in the correlation function
# For this, search the first maximum behind the manual pick. The distance
# between the pick and this first maximum is taken as halfwidth of the
# maximum to be searched (the full width correspond thus approximately to
# half a wavelenght of the signal)
                mapo,mava,mip,miva = \
                    self.main.utilities.min_max(self.v[trace_ref,n_pick_ref:],\
                                                half_width=10)
                half_width = max(mapo[0],10)
                print(f"Half-width: {half_width}")
                n_ref1 = 0
                nd = int(np.round(t_half/self.data.dt,0))
                n_ref1 = max(0,n_pick_ref-nd)
                n_max = self.v.shape[1]
                n_ref2 = min(n_pick_ref+nd,n_max)
# tr_ref will contain the "source signal"
                tr_ref = self.v[trace_ref,n_ref1:n_ref2]
                tr_ref = (tr_ref-np.mean(tr_ref))/np.std(tr_ref)
                n_ref = len(tr_ref)
                off0 = self.traces.offset[trace]
                off0_a = abs(off0)
                if trace_ref > 0:
                    n_act = n_pick_ref
# Loop over the traces left of the reference trace
                    for i in range(trace_ref-1,-1,-1):
                        if i<np.array(self.tr,dtype=int).min():
                            break
                        trace = self.actual_traces[i]
                        off = self.traces.offset[trace]
                        off_a = abs(off)
                        if off_a < 0.01:
                            self.traces.npick[trace] = 1
                            try:
                                self.traces.pick_times[trace][0] = 0.
                                self.traces.pick_times_min[trace][0] = -2*self.data.dt
                                self.traces.pick_times_max[trace][0] = 2*self.data.dt
                            except:
                                self.traces.pick_times[trace].append(0.)
                                self.traces.pick_times_min[trace].append(\
                                                            -2*self.data.dt)
                                self.traces.pick_times_max[trace].append(\
                                                            2*self.data.dt)
                            continue
                        nt_help = np.where(self.i_tr_help == i)[0]
                        help_flag = False
# Set signal tr2 to be correlated with the reference signal
                        if nt_help > 0:
                            n_act = int((self.t_pk_help[nt_help[0]]-self.data.t0)\
                                        /self.data.dt)
                            help_flag = True
                        elif off_a > off0_a and off_a<off_lim:
                            n_act += int(dt_pick/self.data.dt)
                        elif off_a < off0_a and off_a<off_lim:
                            n_act -= int(dt_pick/self.data.dt)
                        off0 = off
                        n_sam1 = max(0,n_act-nd)
                        n_sam2 = min(n_act+nd,n_max)
                        tr2 = self.v[i,n_sam1:n_sam2]
                        s = np.std(tr2)
                        if np.isclose(s, 0.):
                            continue
                        tr2 = (tr2-np.mean(tr2))/s
# Do correlation
                        cor = np.correlate(tr2,tr_ref,'full')
# Find maximum of correlation function
                        n_disp = int(len(cor)/2)
                        if help_flag:
# If nearby picks exist, search maximum nearest to that time
                            max_pos,max_val,min_pos,min_val = \
                                self.main.utilities.min_max(cor,\
                                          half_width=half_width)
                            if len(max_pos) > 0:
                                i_max_pos = max_pos[np.argmin(np.abs(max_pos-n_disp))]
                                dm = i_max_pos-n_disp+n_act-n_pick_ref
                            else:
                                i_max_pos = np.argmax(cor)
                                dm = i_max_pos-n_ref+1+n_sam1-n_ref1
                        else:
# If not, use absolute maximum
                            i_max_pos = np.argmax(cor)
                            dm = i_max_pos-n_ref+1+n_sam1-n_ref1
                        n_act = n_pick_ref+dm
                        self.traces.pick_times[trace].\
                            append(time_ref+dm*self.data.dt)
                        self.traces.npick[trace] += 1
                        self.traces.pick_times_min[trace].\
                            append(self.traces.pick_times[trace][-1]-umin)
                        self.traces.pick_times_max[trace].\
                            append(self.traces.pick_times[trace][-1]+umin)
# Now treat traces on the right side of the reference trace
                if trace_ref < self.actual_number_traces:
                    n_act = n_pick_ref
                    for i in range(trace_ref+1,self.actual_number_traces):
                        if i>np.array(self.tr,dtype=int).max():
                            break
                        trace = self.actual_traces[i]
                        off = self.traces.offset[trace]
                        off_a = abs(off)
                        if off_a < 0.01:
                            self.traces.npick[trace] = 1
                            try:
                                self.traces.pick_times[trace][0] = 0.
                                self.traces.pick_times_min[trace][0] = -2*self.data.dt
                                self.traces.pick_times_max[trace][0] = 2*self.data.dt
                            except:
                                self.traces.pick_times[trace].append(0.)
                                self.traces.pick_times_min[trace].append(\
                                                            -2*self.data.dt)
                                self.traces.pick_times_max[trace].append(\
                                                            2*self.data.dt)
                            continue
# Set signal tr2 to be correlated with the reference signal
                        nt_help = np.where(self.i_tr_help == i)[0]
                        help_flag = False
                        if nt_help > 0:
                            n_act = int((self.t_pk_help[nt_help[0]]-self.data.t0)\
                                        /self.data.dt)
                            help_flag = True
                        elif off_a > off0_a and off_a<off_lim:
                            n_act += int(dt_pick/self.dt)
                        elif off_a < off0_a and off_a<off_lim:
                            n_act -= int(dt_pick/self.data.dt)
                        off0 = off
                        n_sam1 = max(0,n_act-nd)
                        n_sam2 = min(n_act+nd,n_max)
                        tr2 = self.v[i,n_sam1:n_sam2]
                        s = np.std(tr2)
                        if np.isclose(s, 0.):
                            continue
                        tr2 = (tr2-np.mean(tr2))/s

# Do correlation
                        cor = np.correlate(tr2,tr_ref,'full')
# Find maximum of correlation function
                        n_disp = int(len(cor)/2)
                        if help_flag:
# If nearby picks exist, search maximum nearest to that time
                            max_pos,max_val,min_pos,min_val = \
                                self.main.utilities.min_max(cor,\
                                          half_width=half_width)
                            if len(max_pos > 0):
                                i_max_pos = max_pos[np.argmin(np.abs(max_pos-n_disp))]
                            else:
                                i_max_pos = np.argmax(cor)
                            dm = i_max_pos-n_disp+n_act-n_pick_ref
                        else:
# If not, use absolute maximum
                            i_max_pos = np.argmax(cor)
                            dm = i_max_pos-n_ref+1+n_sam1-n_ref1
                        n_act = n_pick_ref+dm
                        self.traces.pick_times[trace].\
                            append(time_ref+dm*self.data.dt)
                        self.traces.npick[trace] += 1
                        self.traces.pick_times_min[trace].\
                            append(self.traces.pick_times[trace][-1]-umin)
                        self.traces.pick_times_max[trace].\
                            append(self.traces.pick_times[trace][-1]+umin)
# When all traces inside the zoom have been treated, leave the function
                self.end = True

# Start function

# Search picks done earlier in nearby gathers an plot a line connecting them
        self.x_pk_help,self.t_pk_help,self.i_tr_help = self.searchNearPicks()
# Plot found picks as dotted green line
        figure = self.figs[self.fig_plotted]
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
        self.line, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y,\
                                                      animated=True)
        self.cidpress = self.line.figure.canvas.mpl_connect('button_press_event',\
                                                            onPress)

# Wait for mouse event
        while (self.end != True):
            QtCore.QCoreApplication.processEvents()

# Redraw record section with new picks
        self.PlotPicks.setEnabled(True)
        self.Tomography.setEnabled(True)
        self.MovePicks.setEnabled(True)
        self.drawNew(True)
        self.picksPlot()
        if len(x_pk):
            self.axes[self.fig_plotted].plot(x_pk, t_pk, c="g", ls="--")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)
        self.setHelp(self.main_text)

    def Sta_Lta(self):
        import obspy.signal.trigger as trigger

        self.main.function = "slpick"
        visible_flag = True
        figure = self.figs[self.fig_plotted] # set figure (simply shorter variable name)
        results, okButton = self.main.dialog(\
                                             ["Sta length [ms]",\
                                              "Lta length [ms]",\
                                              "Only visible traces (y/n)"],\
                                             ["e","e","e"],[4,50,"y"],\
                                              "Sta-Lta parameter input")
        trig = []

        if okButton == False:
            print("Picking cancelled")
        else:
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
                n2 = np.size(self.v,0)
            for i in range(n1,n2):
                trace = self.actual_traces[i]
                npi = self.traces.npick[trace]
# If a pick exists for this trace, skip search (the idea is that for some
# traces, picks may have been placed before manually or calling Sta_Lta with other
# settings
                if npi > 0:
                    continue
# Check whether the trace has data
                if np.isclose(np.max(np.abs(self.v[i,:])), 0.):
                    continue
                t_trig_min = abs(self.traces.offset[trace])/self.v_max_trigger
                n_trig_min = int((t_trig_min - self.data.t0)/self.data.dt)
                trg = trigger.recursive_sta_lta(self.v[i,:], nsta, nlta)
                trg[0:n_trig_min] = 0
                indices = np.arange(len(trg))
                j50 = indices[trg>max(trg)*0.5][0]
                j25 = indices[trg>max(trg)*0.45][0]
                j75 = indices[trg>max(trg)*0.55][0]
                dminus = j25-j50
                dplus = j75-j50
# Search next maximum and go back to near zero crossing
                alim = self.v[i,j50]*0.1
                for k in range(j50,0,-1):
                    if self.v[i,k] < alim:
                        j50 = k
                        j25 = k+dminus
                        j75 = k+dplus
                        break
                t50 = self.data.time[j50]
                t25 = min(self.data.time[j25],t50-2*self.data.dt)
                t75 = max(self.data.time[j75],t50+2*self.data.dt)
                self.traces.npick[trace] += 1
                self.traces.pick_times[trace].append(t50)
                self.traces.pick_times_max[trace].append(t75)
                self.traces.pick_times_min[trace].append(t25)
                trig.append(t50)
            self.background = figure.canvas.copy_from_bbox(figure.bbox)
            self.picksPlot()
            self.Tomography.setEnabled(True)
            self.PlotPicks.setEnabled(True)
            self.MovePicks.setEnabled(True)

    def ampPick(self,data,istart=0,iend=0,half_width=7):
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
        The pick is located in a first try near the highest value of dmax, located
        later than istart, precisely at the position of this maximum minus
        half the distance between the corrresponding maximum position and the following
        minimum position. If the maximum value of dmin is earlier than that of dmax
        and its value is larger than the one of dmax, then the pick is switched to
        the position of maximum dmin, shifted again by half the distance between this
        minimum and the following maximum.
        Then as a last test, the function searches for a (near) zero-crossing between
        the chosen pick position and the corresponding maximum (minimum) of dmax (dmin)
        If there is one, it will be the final pick, if not the chosen pick is
        maintained.

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
            The default is 7. This values has been found to give the best results.

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
        if istart>0:
            zero_pos = istart
    # Search relative maxima and minima from beginning of data set to iend
        max_pos,max_val,min_pos,min_val =\
            self.main.utilities.min_max(data[:ie],half_width)
    # Prepare vectors for amplitude parameters
    #   sums_pos will contain the sums for relative maxima. i.e.
    #       maximum - (minimum earler+minimum later)
    #   sums_neg will contain the sums for relative minima. i.e.
    #       (maximum earler+maximum later)-minimum
        lx = len(max_pos)
        ln = len(min_pos)
        sums_pos = np.zeros(lx)
        sums_neg = np.zeros(ln)
    # Calculate the sums
        if min_pos[0]<max_pos[0]:
    # First relative minimum comes earlier than first relative maximum
    #   In this case, calculate sums_pos from the first maximum on. For later
    #   search of the pick for minima, the following maximum will have
    #   an index higher by 1 than the value of the positionof the mminimum
    #   (iadd_pos)
            ipos = 0
            iadd_pos = 1
            while ipos<lx-1 and ipos<ln-1:
                sums_pos[ipos] = np.abs(max_val[ipos]-(min_val[ipos]+\
                                        min_val[ipos+1]))
                if sums_pos[ipos]==0:
                    sums_pos[ipos] = np.nan
                ipos += 1
            ipos = 1
            iadd_neg = 0
            while ipos<lx-1 and ipos<ln-1:
                sums_neg[ipos] = np.abs(-min_val[ipos]+(max_val[ipos-1]+\
                                        max_val[ipos]))
                if sums_neg[ipos]==0:
                    sums_neg[ipos] = np.nan
                ipos += 1
        else:
    # First relative minimum comes later than first relative maximum
    #   In this case, calculate sums_pos from the second maximum on. For later
    #   search of the pick for maxima, the following minimum will have
    #   an index higher by 1 than the value of the position of the maximum
    #   (iadd_neg)
            ipos = 1
            iadd_pos = 0
            while ipos<lx and ipos<ln:
                sums_pos[ipos] = np.abs(max_val[ipos]-(min_val[ipos-1]+\
                                        min_val[ipos]))
                if sums_pos[ipos]==0:
                    sums_pos[ipos] = np.nan
                ipos += 1
            ipos = 0
            iadd_neg = 1
            while ipos<lx-1 and ipos<ln-1:
                sums_neg[ipos] = np.abs(-min_val[ipos]+(max_val[ipos]+\
                                        max_val[ipos+1]))
                if sums_neg[ipos]==0:
                    sums_neg[ipos] = np.nan
                ipos += 1
    # Calculate amplitude relations for maxima (dsums_pos) and minima (dsums_neg)
        dsums_pos = np.zeros_like(sums_pos)
        dsums_neg = np.zeros_like(sums_neg)
    # For possible plotting purposes, the amplitude coefficients are
    #   normalized to the data maximum. This may be eliminated without
    #   changing the result
        max_where = np.where(max_pos>zero_pos)
        min_where = np.where(min_pos>zero_pos)
        if len(sums_pos)<4 or len(sums_neg)<4 or np.size(max_where)==0 or\
            np.size(min_where)==0:
            return np.nan,np.nan
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
        if iadd_pos>0:
            nmx_pos = np.argmax(dsums_pos[ini_pos:])+ini_pos
        else:
            nmx_pos = np.argmax(dsums_pos[ini_pos:])+ini_pos
    # Search the pick position corresponding to the maximum value of the
    #   amplitude relations for the minima.
        nmx_neg = ini_neg
        if iadd_neg>0:
            nmx_neg = np.argmax(dsums_neg[ini_neg:])+ini_neg
        else:
            nmx_neg = np.argmax(dsums_neg[ini_neg:])+ini_neg
    # The position of the relations for maxima is shifted by half
    #   the distance between this maximum relation and the next relative minimum
    #   (i.e., in principle, a quarter of a wavelength)
        pos = np.int(max_pos[nmx_pos]-(min_pos[nmx_pos+iadd_pos]-max_pos[nmx_pos])/2)
    # The position of the relations for minima is shifted by half the
    #   distance between this maximum relation and the next relative maximum
    #   (i.e., in principle, a quarter of a wavelength)
        neg = np.int(min_pos[nmx_neg]-(max_pos[nmx_neg+iadd_neg]-min_pos[nmx_neg])/2)
        tpos = pos*self.data.dt+self.data.t0
        tneg = neg*self.data.dt+self.data.t0
    # Now check which of the two pick positions to choose:
    #   normally, the one of the maxima is used, supposing that most traces
    #   have a positive first arrival (it will halp the program if the input
    #   data are signed in this sense).
    #   If, however, the pick calculated from the minima is earlier and has
    #   in addition a higher relationship indicator take that position.
        if tneg<tpos and dsums_neg[nmx_neg]>dsums_pos[nmx_pos]:
            tpk = tneg
            npk1 = neg
            npk2 = np.int(min_pos[nmx_neg])
    # Check wether there is a zero crossing or near zero crossing between
    #   the shifted position of the minima pick and the corresponding minimum.
    #   If one is found, this will be the final pick position.
            a_comp = min_val[nmx_neg]/100
            for i in range(npk2,npk1-1,-1):
                if np.abs(data[i])<np.abs(a_comp) or a_comp*data[i]<=0:
                    break
        else:
    # Similar procededure as above if the pick is obtained from relations of
    #   maxima
            tpk = tpos
            npk1 = pos
            npk2 = max_pos[nmx_pos]
            a_comp = max_val[nmx_pos]/100
            for i in range(npk2,npk1-1,-1):
                if np.abs(data[i])<np.abs(a_comp) or a_comp*data[i]<=0:
                    break
    # Calculate and return final pick time
        tpk = (i-1)*self.data.dt+self.data.t0
        amp = data[i]
        return tpk,amp

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
            The default is 5. This values has been found to give the best results.
        nlines : Int, optional
            Number of lines to be fitted to first-guess picks. The first line
            is a direct wave, i.e., it passes through the origin. the other
            lines are refracted waves. The slope has to increase for every line.
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
        self.main.function = "apick"
        figure = self.figs[self.fig_plotted] # set figure (simply shorter variable name)
        self.keys_held = set()
        results, okButton = self.main.dialog(\
                                        ["Half_width for max",\
                                         "Number of lines fitted",\
                                         "Initial distance [m]",\
                                         "signal length [ms]",\
                                         "Maximum velocity [m/s]"],\
                                         ["e","e","e","e","e"],\
                                         [5,3,15,50,self.v_max_trigger],\
                                         "Amp_Pick parameter input")

        if okButton == False:
            pass
        else:
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
            n_end = np.zeros(ntrac,dtype=np.int)
            n_len =  np.int(sig_len/self.data.dt)
            nt0 = np.int(-self.data.t0/self.data.dt)
            if visible_flag:
                n1 = self.tr[0]
                n2 = self.tr[-1]+1
            else:
                n1 = 0
                n2 = ntrac
            for j in range(n1,n2):
                i = j
                itry = 0
                data = self.v[i,:]
                if np.std(self.v[i,:])>0:
                    if self.v_max_trigger>0:
                        n_start = np.int((np.abs(self.x[self.indices[i]])/\
                                self.v_max_trigger-self.data.t0)/self.data.dt)
                    else:
                        n_start = nt0
                    n_end[j] = n_start+n_len
                    ne = n_end[j]
                    dne = n_end[j]-n_start
                    while np.isnan(t_pick[j]):
                        itry += 1
                        t_pick[j],a_pick[j] = self.ampPick(data,n_start,ne,\
                            half_width=half_width)
                        ne = min(ne+dne,len(data))
            for iter in range(2):
                x = np.array(self.x)[self.indices]
                best_off,best_slope,best_intercept,x_regres,y_regres,\
                    y_points_regres,r2,medians_pos =\
                    self.bestLines(x,t_pick,nlines,dist1,1)
                dist = 1.5*np.nanstd(y_points_regres-t_pick)
                for i in range(ntrac):
                    if i<n1 or i > n2:
                        continue
                    if np.abs(t_pick[i]-y_points_regres[i])>dist:
                        nt_reg = int((y_points_regres[i]-self.data.t0)/self.data.dt)
                        max_pos,max_val,min_pos,min_val = \
                        self.main.utilities.min_max(data[:n_end[i]],5)
                        nmx_pos = len(max_pos)
                        nmn_pos = len(min_pos)
                        dpos = np.zeros(nmx_pos)
                        if max_pos[0]>min_pos[0]:
                            for k in range(nmx_pos):
                                if k>nmn_pos-2:
                                    break
                                dpos = np.int((min_pos[k+1]-max_pos[k])/2)
                        else:
                            for k in range(1,nmx_pos):
                                if k>nmn_pos-1:
                                    break
                                dpos = np.int((min_pos[k]-max_pos[k])/2)
                        mxp = max_pos-dpos
                        k = np.argmin(np.abs(mxp-nt_reg))
                        nt_pik = mxp[k]
                        for kk in range(max_pos[k],mxp[k],-1):
                            if data[kk]<max_val[k]/100 or data[kk]*max_val[k]<0:
                                nt_pik = kk
                                break
                        t_pick[i] = nt_pik*self.data.dt+self.data.t0
                        a_pick[i] = data[nt_pik]
            if self.data.filtered:
                unc = 4*self.data.dt
            else:
                unc = 2*self.data.dt
            for j in range(ntrac):
                if j<n1 or j > n2:
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
            self.picksPlot()
            self.Tomography.setEnabled(True)
            self.PlotPicks.setEnabled(True)
            self.MovePicks.setEnabled(True)



    def followRect(self):
        """
        Pull rectangle across plot
        """
        global figure
        def onPress(event):
            global figure
            if event.button == 1: # when left mouse botton is pressed initialize rectangle
                if len(self.coor_x) == 0: # in the beginning store clicked point as first corner of rectangle
                    self.coor_x.append(event.xdata)
                    self.coor_y.append(event.ydata)
                    self.background = figure.canvas.copy_from_bbox(figure.bbox) # store background figure
                    self.coor_x.append(event.xdata)
                    self.coor_y.append(event.ydata)
                self.canvas = self.line.figure.canvas #define canvas and new axes
                self.axl = self.line.axes
                self.line.set_data(self.coor_x,self.coor_y)
                self.axl.draw_artist(self.line) # Draw point
                self.cidmotion = self.line.figure.canvas.mpl_connect\
                    ('motion_notify_event', onMotion) # set action on mouse motion
                self.cidrelease = self.line.figure.canvas.mpl_connect\
                    ('button_release_event', onRelease) # set action on mouse release

        def onRelease(event):
            global figure
            self.line.figure.canvas.mpl_disconnect(self.cidpress) # free actions
            self.line.figure.canvas.mpl_disconnect(self.cidmotion)
            self.line.figure.canvas.mpl_disconnect(self.cidrelease)
            self.line.set_animated(False)
            if len(self.coor_x) > 0:
                figure.canvas.restore_region(self.background) # plot background image without rectangle
            self.background = None # remove stored background
            self.released = True # ste flag to finish module

        def onMotion(event):
            global figure
            self.coor_x[-1] = event.xdata # replace last coordinate point by actual mouse position
            self.coor_y[-1] = event.ydata
            line_coor_x = [self.coor_x[0], self.coor_x[1], self.coor_x[1],
                           self.coor_x[0], self.coor_x[0]] # set x and y coordiantes of rectangle plot
            line_coor_y = [self.coor_y[0], self.coor_y[0], self.coor_y[1], self.coor_y[1], self.coor_y[0]]
            self.line.set_data(line_coor_x, line_coor_y)
            self.canvas = self.line.figure.canvas
            self.axl = self.line.axes # define new axes values
            figure.canvas.restore_region(self.background) # plot background
            self.axl.draw_artist(self.line) # plot new rectangle
            self.canvas.blit(self.axl.bbox)

        self.released = False # initialize flags and coordinate vector
        self.coor_x = []
        self.coor_y = []
        figure = self.figs[self.fig_plotted] # set actual figure (simply shorter variable name)
        self.line, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y,\
                               animated=True) # define line value
        self.cidpress = self.line.figure.canvas.mpl_connect\
            ('button_press_event', onPress) # set action on mouse press
        while (self.released != True):
            QtCore.QCoreApplication.processEvents() # cycle event detection as long as released flag is not set

    def followLine(self, release_flag = False, nleft=1, nright=1):
        """
        Pull line across plot
        Input:
                release_flag (boolean): If True, end of line when left button released
                                        if False, end of line triggered by pressing right button
                nleft (int): if 0 start line for negative direction at origin
                             if not, start line at the position of first click
                nright (int): if 0 start line for positive direction at origin
                             if not, start line at the position of first click
        """
        global figure,r_flag
        def onPress(event):
            global figure, r_flag
            self.line_click = False
            if event.button == 1: # left mouse button is pressed
                if event.xdata==None or event.ydata==None:
                    return False
                if len(self.coor_x) == 0:
                    if (event.xdata<0 and nleft==0) or\
                       (event.xdata>=0 and nright==0):
                        self.start = [0, 0]
                        self.coor_x.append(0)
                        self.coor_y.append(0)
                    else:
                        self.start = [event.xdata, event.ydata]
                        self.coor_x.append(event.xdata)
                        self.coor_y.append(event.ydata)
                    if event.xdata < 0:
                        self.side = -1
                    else:
                        self.side = 1
                    self.background = figure.canvas.copy_from_bbox(figure.bbox)
                self.coor_x.append(event.xdata) # set starting point initially also as end point
                self.coor_y.append(event.ydata)
                self.canvas = self.line.figure.canvas
                self.axl = self.line.axes
                self.line.set_data(self.coor_x,self.coor_y)
                self.axl.draw_artist(self.line)
                self.cidmotion = self.line.figure.canvas.mpl_connect\
                    ('motion_notify_event', onMotion) # set action on mouse motion
                if release_flag:
                    self.cidrelease = self.line.figure.canvas.mpl_connect\
                        ('button_release_event', onRelease) # set action on mouse release
            else: # if right button is pressed, finish
                self.line.figure.canvas.mpl_disconnect(self.cidpress)
                try:
                    self.line.figure.canvas.mpl_disconnect(self.cidmotion)
                    self.line_click = True
                except:
                    self.line_click = False
                self.line.set_animated(False)
                if self.line_click:
                    figure.canvas.restore_region(self.background)
                self.background = None
                self.released = True
                return True

        def onRelease(event):
# If line finishes when button is released do this here
            global figure
            self.line.figure.canvas.mpl_disconnect(self.cidpress)
            self.line.figure.canvas.mpl_disconnect(self.cidmotion)
            self.line.figure.canvas.mpl_disconnect(self.cidrelease)
            self.line.set_animated(False)
            if len(self.coor_x) > 0:
                figure.canvas.restore_region(self.background)
            self.background = None
            self.released = True
            return False

        def onMotion(event):
            global figure
            if event.xdata==None or event.ydata==None:
                return False
            if self.main.function == "P_model":
                if (self.side<0 and event.xdata>0) or\
                   (self.side>0 and event.xdata<0):
                    event.xdata = 0.
            self.coor_x[-1] = event.xdata # set second point of line as actual mouse position
            self.coor_y[-1] = event.ydata
            if len(self.start)>0 and len(self.coor_x)==2 and \
                self.main.function != "spectrum": # Check if starting point is outside window
                zxmin = self.zooms[self.i_zooms][0]
                zxmax = self.zooms[self.i_zooms][1]
                if self.start[0]<zxmin and zxmin>=0:
                    self.coor_x[0] = zxmin
                    self.coor_y[0] = self.start[1]+(self.coor_y[-1]-self.start[1])\
                        *(zxmin-self.start[0])/(self.coor_x[-1]-self.start[0])
                    if self.coor_y[0]<self.time_plt_min:
                        self.coor_y[0]=self.time_plt_min
                        self.coor_x[0] = self.start[0]+(self.coor_x[-1]-\
                            self.start[0])*(self.time_plt_min-self.start[1])/\
                            (self.coor_y[-1]-self.start[1])
                elif self.start[0]>zxmax and zxmax<=0:
                    self.coor_x[0] = zxmax
                    self.coor_y[0] = self.start[1]+(self.coor_y[-1]-self.start[1])\
                        *(zxmax-self.start[0])/(self.coor_x[-1]-self.start[0])
                    if self.coor_y[0]<self.time_plt_min:
                        self.coor_y[0]=self.time_plt_min
                        self.coor_x[0] = self.start[0]+\
                            (self.coor_x[-1]-self.start[0])*\
                            (self.time_plt_min-self.start[1])/\
                            (self.coor_y[-1]-self.start[1])
            self.line.set_data(self.coor_x, self.coor_y) # Draw new line
            self.canvas = self.line.figure.canvas
            self.axl = self.line.axes
            figure.canvas.restore_region(self.background)
            self.axl.draw_artist(self.line)
            self.canvas.blit(self.axl.bbox)

        r_flag = release_flag # set flags and initialize coordinates
        self.released = False
        self.start = []
        self.coor_x = []
        self.coor_y = []
        figure = self.figs[self.fig_plotted] # set figure (simply shorter variable name)
        self.line, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y,\
                                                      animated=True)
        self.cidpress = self.line.figure.canvas.mpl_connect\
                                ('button_press_event', onPress)
        while (self.released != True): # as long as release flag is not set listen to events
            QtCore.QCoreApplication.processEvents()


    def bestLines(self,x_points,y_points,n_lines=2,xoff1=0,Lnorm=2,origin=True,
                  refra=True):
        """
        Function to calculate best fitting series of n_lines lines through points
        (x_points, y_points).
        First line is supposed to pass through the origin.
        Slopes have to decrease monotonously
        The function calculates first the unique absolute positions of x_points, which
        means that the lines are supposed to be symmetric around the origin.
        The function then searches the optimum points where to break slopes
        If more than 2 lines have to be fitted, it is possible to give a maximum
        distance for the points to be fitted for the two first lines.

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
            offsets of line breaks.
            first line goes from best_off[0] until best_off[1] inclusive
            second line goes from best_off[1] exclusive to best_off[2] inclusive etc.
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

        import scipy.stats
# check best regression lines crossing the obtained picks
    #   The direct wave and i_refra_max refractions are fitted
        if Lnorm<1 or Lnorm>2:
            print("Lnorm is not correctly defined, must be 1 or 2")
            print("Default value of 2 is used")
            Lnorm = 2
        n_points = len(x_points)
        best_slope = np.zeros(n_lines,dtype='float')
        best_intercept = np.zeros(n_lines,dtype='float')
        best_i = np.zeros(n_lines+1,dtype='int')
        best_off = np.zeros(n_lines+1,dtype='float')
    # x_off contains the ordered offsets, i_off the corresponding trace number
    # ind_c and ind_l contain the positions of the traces in the data stream st
        x_abs = np.abs(x_points)
        i_off = np.argsort(x_abs)
        x_off = x_abs[i_off]
        x_unique = np.unique(x_off)
        n_unique = len(x_unique)
    # t_pick_test contains the picks ordered by offset
        y_test = y_points[i_off]
        y_regres = np.zeros(n_unique)
        x_regres = np.zeros(n_unique)
        medians = np.zeros(n_unique)

    # First check for the direct wave, for this take the 15 traces with smalles offset
    #       and split picks into two groups (direct wave and first refracted wave)
        if n_lines > 2 and xoff1>0:
            max_off = xoff1
        else:
            max_off = x_unique[n_unique-1]
    # xpk will contain the offsets of existing pickks, ypk the corresponding times
    #   (t_env_der_pick may contain NaNs). i_tr is the trace number of each pick
        xpk = []
        ypk = []
        i_tr = []
    # npts contains the number of picks until (and including) each position
        npts = np.zeros(n_unique+2,dtype='int')
        n_off = 0
        i=0
        npts[0] = -1
        for n_off in range(n_unique):
            if x_unique[n_off] <= max_off:
                i_off_test = n_off
            while x_off[i]==x_unique[n_off]:
                if np.isnan(y_test[i]) == False:
                    xpk.append(x_off[i])
                    ypk.append(y_test[i])
                    i_tr.append(i)
                    npts[n_off] += 1
                i += 1
                if i>=n_points:
                    break
            npts[n_off+1] = npts[n_off]
    # create numpy arrays from lists
        xpk = np.array(xpk,dtype='float')
        ypk = np.array(ypk,dtype='float')
        yy = np.zeros_like(xpk,dtype='float')
        sig_best = 1E20
        i_best = -1
        i_best_sig = 2

    # If n_lines==1 calculate the best fitting line through all points
        if n_lines == 1 or n_unique<4:
            if origin:
                slope1 = np.dot(xpk,ypk)/np.dot(xpk,xpk)
                intercept1 = 0.
                if Lnorm==1:
                    slope1 = self.L1_regres(xpk,ypk,[slope1])
            else:
                slope1,intercept1,r1,p1,stderror1 = scipy.stats.linregress(xpk,ypk)
                if Lnorm==1:
                    par = self.L1_regres(xpk,ypk,[slope1,intercept1])
                    slope1 = par[0]
                    intercept1 = par[1]
            y_regres = x_off*slope1+intercept1
            x_regres = x_off
            y_points_regres = np.abs(x_points)*slope1+intercept1
            best_off[1] = np.max(x_off)
        else:
    # Start loop over offsets to see for which point a line split gives the best fit
            for i in range(2,max(i_off_test-4,3)):
                ne = npts[i_off_test]
    #            print ("direct:", i, npts[i]+1,ne)
    # slope1 corresponds to the slope of the direct wave for the first i picks
                if origin:
                    slope1 = np.dot(xpk[0:npts[i]+1],ypk[0:npts[i]+1])/\
                             np.dot(xpk[0:npts[i]+1],xpk[0:npts[i]+1])
                    if Lnorm==1:
                        slope1 = self.L1_regres(xpk[0:npts[i]+1],\
                                                    ypk[0:npts[i]+1],[slope1])
                    intercept1 =0.
                else:
                    slope1,intercept1,r1,p1,stderror1 = \
                        scipy.stats.linregress(xpk[0:npts[i]+1],ypk[0:npts[i]+1])
                    if Lnorm==1:
                        par = self.L1_regres(xpk[npts[i+1]:ne],ypk[npts[i+1]:ne],\
                                             [slope1,intercept1])
                        slope1 = par[0]
                        intercept1 = par[1]
    # slope2 corresponds to the slope of the refacted wave for picks i+1 until 15
                slope2,intercept2,r2,p2,stderror2 = \
                    scipy.stats.linregress(xpk[npts[i+1]:ne],ypk[npts[i+1]:ne])
                if Lnorm==1:
                    par = self.L1_regres(xpk[npts[i+1]:ne],ypk[npts[i+1]:ne],\
                                         [slope2,intercept2])
                    slope2 = par[0]
                    intercept2 = par[1]
    # yy contians the theoretical arrival times for the pick positions
                yy[0:npts[i]+1] = slope1*xpk[0:npts[i]+1]+intercept1
                yy[npts[i+1]:ne] = slope2*xpk[npts[i+1]:ne]+intercept2
    # Calculate misfit
                diff = np.abs(ypk[:ne]-yy[:ne])
                sig = np.sqrt(np.dot(diff,diff)/len(diff))
    # If misfit is smallest, save parameters
                if sig<sig_best and ((slope1>0 and slope2>0 and slope2<slope1\
                    and intercept2>0) or not refra):
                    i_best = i
                    s1_best = slope1
                    s2_best = slope2
                    int1_best = intercept1
                    int2_best = intercept2
                    sig_best = sig
    # Store best fitting parameters and synthetic times (in regres_pick)
                if i_best < 0:
                    i_best = i
                    s1_best = slope1
                    s2_best = slope2
                    int1_best = intercept1
                    int2_best = intercept2
                    sig_best = sig
            best_i[1] = i_best
            best_slope[0] = s1_best
            best_intercept[0] = int1_best
            best_off[1] = x_unique[i_best]
            for i in range(best_i[1]+1):
                y_regres[i] = x_unique[i]*s1_best+int1_best
                x_regres[i] = x_unique[i]
            if n_lines>2:
                i_best = -1
    # Now calculate best fitting refractions
                for i_slope in range(1,n_lines-1):
                    sig_best = 1E20
                    n1 = npts[best_i[i_slope]]+1
                    for n_off in range(best_i[i_slope]+3,n_unique-3):
                        n2 = npts[n_off]+1
                        slope1,intercept1,r1,p1,stderror1 = scipy.stats.linregress(xpk[n1:n2],ypk[n1:n2])
                        if Lnorm==1:
                            par = self.L1_regres(xpk[n1:n2],ypk[n1:n2],[slope1,intercept1])
                            slope1 = par[0]
                            intercept1 = par[1]
                        slope2,intercept2,r2,p2,stderror2 = scipy.stats.linregress(xpk[n2:],ypk[n2:])
                        if Lnorm==1:
                            par = self.L1_regres(xpk[n2:],ypk[n2:],[slope1,intercept1])
                            slope2 = par[0]
                            intercept2 = par[1]
                        yy[n1:n2] = slope1*xpk[n1:n2]+intercept1
                        yy[n2:] = slope2*xpk[n2:]+intercept2
                        diff = np.abs(ypk[n1:]-yy[n1:])
                        sig = np.sqrt(np.dot(diff,diff)/len(diff))
                        if sig < sig_best and ((slope1>0 and slope2>0 and slope2<slope1\
                            and intercept1>0 and intercept2>0) or not refra):
                            i_best = n_off
                            s1_best = slope1
                            s2_best = slope2
                            int1_best = intercept1
                            int2_best = intercept2
                            sig_best = sig
                        if i_best < 0:
                            i_best = i
                            s1_best = slope1
                            s2_best = slope2
                            int1_best = intercept1
                            int2_best = intercept2
                            sig_best = sig
                    best_i[i_slope+1] = i_best
                    best_off[i_slope+1] = x_unique[i_best]
                    best_slope[i_slope] = s1_best
                    best_slope[i_slope+1] = s2_best
                    best_intercept[i_slope] = int1_best
                    best_intercept[i_slope+1] = int2_best
                    for i in range(best_i[i_slope],best_i[i_slope+1]+1):
                        y_regres[i] = x_unique[i]*s1_best+int1_best
                best_i[n_lines] = n_unique
                best_slope[n_lines-1] = x_unique[n_unique-1]
            elif n_lines == 2:
                best_slope[1] = s2_best
                best_intercept[1] = int2_best
            for i in range(best_i[n_lines-1],n_unique):
                y_regres[i] = x_unique[i]*s2_best+int2_best
                x_regres[i] = x_unique[i]
            y_points_regres = np.empty(n_points)
            y_points_regres.fill(np.nan)
            for i in range(n_points):
                no = np.where(x_unique==x_abs[i])
                y_points_regres[i] = y_regres[no]
            medians[0] = np.median(ypk[0:npts[1]+1])
            for i in range(1,n_unique):
                try:
                    if npts[i+1] > npts[i]:
                        medians[i] = np.median(ypk[npts[i]+1:npts[i+1]+1])
                except:
                    pass
        del x_abs
        r = np.dot((x_regres-x_regres.mean()),(y_regres-y_regres.mean()))/\
            (len(x_regres)*np.std(x_regres)*np.std(y_regres))
        r2 = r*r
        return best_off,best_slope,best_intercept,x_regres,y_regres,\
               y_points_regres,r2,medians

    def L1_regres(self,x,y,par_ini=[0,0]):
        """
        Function calculates linear regression line over points[x,y] using an L1 norm

        Parameters
        ----------
        x : numpy.float 1D array
            Contains X coordinates of points to be fitted
        y : numpy.float 1D array
            Contains X coordinates of points to be fitted
        par_ini : float 1D list with one or two values optional, default: [0,0]
            Initial guess (usually obtained from a linear regression fit)
            if len(par_ini)==1, it is supposed that the given value is an initial guess
            of the slope and that the regression line passes through the origin
            if len(par_ini)==2, first parameter is slope, second intercept and
            both parameters are optimized
        Returns
        -------
        param : numpy.float 1D array , same size as par_ini
            optimized parameters
            if len(par_ini)==1, param contains only one value, the slope
            if len(par_ini)==2, param contains [slope,intercept]

        """
        import numpy as np
        from scipy.optimize import minimize

        def cost_function(params,X,y):
            return np.sum(np.abs(y -X.dot(params)))

        if len(par_ini)==1:
            X = np.ones((len(x),1))
        else:
            X = np.ones((len(x),2))
            X[:,0] = x
        param = minimize(cost_function, par_ini, args=(X,y))
        return param.x

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
        import matplotlib.pyplot as plt
        self.main.function = "Animated"
# chose data to be plotted and set clipping to 99% of extreme amplitudes
        nd_start = self.zooms[self.i_zooms][2]
        nd_end = self.zooms[self.i_zooms][3]
        clip_level = 0.99
        vt = self.v[:,nd_start:nd_end].transpose()
        for i in range(vt.shape[0]):
            vt[i,1:-1] = (vt[i,:-2]+vt[i,1:-1]+vt[i,2:])/3
        t = self.data.t0+self.data.dt*np.arange(nd_start,nd_end)
        vmax = np.quantile(vt,clip_level)
        vt[np.abs(vt)>vmax] = vmax*np.sign(vt[np.abs(vt)>vmax])
# Define new figure
        self.drawNew(False)
        fig = self.figs[self.fig_plotted]
        ax = self.axes[self.fig_plotted]
        ax.set_xlim(self.x.min(),self.x.max())
        ax.set_ylim(vt.min(),vt.max())
        self.setHelp(self.animate_text)
# Set matplotlib animated mode
        plt.ion()
        line = ax.plot(self.x,vt[0,:])[0]
        fig.canvas.draw()
        fig.canvas.flush_events()
        for i in range(len(t)):
            line.set_ydata(vt[i,:])
            ax.set_title(f"Time {t[i]*1000:0.2f}ms")
            fig.canvas.draw()
            fig.canvas.flush_events()
            # plt.pause(0.002)
        plt.ioff()



