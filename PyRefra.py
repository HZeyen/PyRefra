# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 15:33:26 2019
last modified on Tue Dec 05, 2023

@author: Hermann Zeyen, University Paris-Saclay, France

Before running the program the first time, change the folders around line 95:
    sys_path is the path to the python scripts
    dir_0 is the path to the data files

Needed Python packages:
    obspy
    pygimli
        In obspy and pygimli, some changes should be done.
            See installation manual for details.
    PyQt5
    matplotlib
    numpy
    os
    sys
    copy
    struct
    signal
    scikit-learn
    datetime
    statsmodels

Needed private packages:

    refraData.py
    refraPlot.py
    refraWindow.ui


Contains the following classes:
    Main
        With the following functions:
            __init__
            eventFilter
            fileOpen
            Handler
            dialog
            closeApp

    Dialog
        With the following functions:
            __init__
            checked
            on_YesButton_clicked
            on_CancelButton_clicked

    Seg2_Slide
        With the following functions:
            __init__
            updateLabel

"""

import sys
import os

#The next two lines may have to be modified:
#  sys_path is the folder where all python program files are located
#  dir0 is the folder where the data are located

# Paths on HZ computer
sys_path = r"E:/Sources_2010/Python_programs"
dir0 = r"E:/Seg2Dat/Fontaines-Salees/2021/2021-10-17_Profil5"
#dir0 = r"E:/Seg2Dat/Fontaines-Salees/2023/Line2"
#dir0 = r"E:/Seg2Dat/Siscarb/Falaise"
#dir0 = r"E:/Seg2Dat/Barbeau"
#dir0 = r"E:/Seg2Dat/Erreurs"
#dir0 = r"E:/Seg2Dat/Beaufremont"

# Example of paths for programs and data on department desktop
# sys_path = r"C:/Users/Utilisateur/Desktop/Geophysique/Sismique"
# dir0 = r"C:/Users/Utilisateur/Desktop/Geophysique/Sismique/data"

# Example of paths for LIONEL
#sys_path = r"C:/Users/marsj/Desktop/Leo/Sismique"
#dir0 = r"C:/Users/marsj/Desktop/Leo/Line2"

#Example of paths for Linux
# sys_path = r"/home/zeyen/src/Python"
# dir0 = r"/home/zeyen/Seismics/Fontaines_salees/2021-10-17_Profil5"

if sys_path not in sys.path:
    sys.path.append(sys_path)

import numpy as np
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QSlider, QLabel, QRadioButton, QButtonGroup)
from signal import signal, SIGINT

import refraData as rD
import refraPlot as rP

class Main(QtWidgets.QWidget):
    def __init__(self, dir0 = None):
        super(Main, self).__init__()

        self.function = "main"
        try:
            os.chdir(dir0)
            self.dir0 = dir0
            self.dir_flag = True
            self.sys_path = sys_path
        except:
            self.dir0 = None
            self.sys_path = sys_path

# Input data
        self.fileOpen()
# Initialize main frame widget
        self.window = rP.Window(self, self.files, self.data, self.traces, self.geo)
        if len(self.geo.types) > 1:
            self.window.component.setEnabled(True)
# Check whether measured picks exist and if so, activate tomography button
        sp = sum(self.traces.npick)
        if sp > 0:
            self.window.Tomography.setEnabled(True)
            print(f"\n{sp} measured picks read")
# Initialize utility routines
        self.utilities = rD.Utilities(self, self.files, self.data,\
                                      self.traces, self.geo, self.window)
        QtWidgets.qApp.installEventFilter(self)

# Define actions for Menu buttons
# Actions for Menu File
        self.window.openFile.triggered.connect(self.fileOpen) # In PyRefra.py
        self.window.Save_SEGY.triggered.connect(self.data.saveSEGY) # In refraData.py
        self.window.Save_SU.triggered.connect(self.data.saveSU) # In refraData.py
        self.window.Save_SEG2.triggered.connect(self.data.saveSEG2) # In refraData.py
        self.window.Save_Bin.triggered.connect(self.data.saveBinary) # In refraData.py
        self.window.Save_ASCII.triggered.connect(self.data.saveASCII) # In refraData.py
        self.window.Save_headers.triggered.connect(self.data.saveHeader) # In refraData.py
        self.window.Save_plot.triggered.connect(self.window.savePlot) # In refraPlot.py
        self.window.quitAction.triggered.connect(self.closeApp) # In PyRefra.py
# Actions for menu Display
        self.window.originalDataScreen.triggered.connect(self.window.originalScreen) # In refraPlot.py
        self.window.originalDataAll.triggered.connect(self.window.original) # In refraPlot.py
        self.window.shotGather.triggered.connect(self.window.plotSG) # In refraPlot.py
        self.window.fileGather.triggered.connect(self.window.plotFG) # refraPlot.py
        self.window.receiverGather.triggered.connect(self.window.plotRG) # refraPlot.py
        self.window.distanceGather.triggered.connect(self.window.plotDG) # refraPlot.py
        self.window.component.triggered.connect(self.window.chooseComponent) # refraPlot.py
        self.window.phaseAngles.triggered.connect(self.window.phasePlot) # In refraPlot.py
        self.window.zoom.triggered.connect(self.window.zooming) # In refraPlot.py
        self.window.zoom_Out.triggered.connect(self.window.zoomOut) # In refraPlot.py
        self.window.zoom_In.triggered.connect(self.window.zoomIn) # In refraPlot.py
        self.window.zoom_Initial.triggered.connect(self.window.zoomIni) # In refraPlot.py
        self.window.t_Norm.triggered.connect(self.window.tNorm) # In refraPlot.py
        self.window.t_Gain.triggered.connect(self.window.tGain) # In refraPlot.py
        self.window.d_Gain.triggered.connect(self.window.dGain) # In refraPlot.py
        self.window.Agc.triggered.connect(self.window.AGC) # In refraPlot.py
# Actions for menu Utilities
        self.window.PModel.triggered.connect(self.utilities.pModel) # In refraData.py
        self.window.SModel.triggered.connect(self.utilities.sModel) # In refraData.py
        self.window.Tomography.triggered.connect(self.utilities.inversion) # In refraData.py
        self.window.Envelopes.triggered.connect(self.utilities.envel) # In refraData.py
        self.window.TauP.triggered.connect(self.utilities.tauP) # In refraData.py
        self.window.falseCol.triggered.connect(self.utilities.falseColour) # In refraData.py
        self.window.TraceSign.triggered.connect(self.window.traceSign) # In refraPlot.py
        self.window.ChangeSign.triggered.connect(self.window.changeSign) # In refraPlot.py
        self.window.Change_colors.triggered.connect(self.utilities.invCol) # In refraData.py
        self.window.Animation.triggered.connect(self.window.animateLine) # In PyRefra.py
        self.window.Attenuation.triggered.connect(self.utilities.atten_amp) # In refraData.py
        self.window.Pseudo_velocity.triggered.connect(self.utilities.pseudo) # In refraData.py
# Actions for menu Picking
        self.window.ManualPicks.triggered.connect(self.window.pickManual) # In refraPlot.py
        self.window.AmpPicks.triggered.connect(self.window.ampPicks) # In refraPlot.py
        self.window.StaLta.triggered.connect(self.window.Sta_Lta) # In refraPlot.py
        self.window.CorrelationPicks.triggered.connect(self.window.corrPick) # In refraPlot.py
        self.window.MovePicks.triggered.connect(self.window.pickMove) # In refraPlot.py
        self.window.Uncertainty_change.triggered.connect(self.window.uncertainty) # In refraPlot.py
        self.window.erasePicks.triggered.connect(self.window.erase_Picks) # In refraPlot.py
        self.window.plotAllPicks.triggered.connect(self.window.plotPickSection) # In refraPlot.py
        self.window.PlotCalculatedTimes.triggered.connect(self.window.plotCalcPicks) # In refraPlot.py
        if self.traces.calc_picks:
            self.window.PlotCalculatedTimes.setEnabled(True)
        self.window.StorePicks.triggered.connect(self.traces.storePicks) # In refraData.py
        self.window.StoreGimli.triggered.connect(self.traces.saveGimli) # In refraData.py
# Actions for menu Filter
        self.window.FrequencyFilter.triggered.connect(self.utilities.filterAll) # In refraData.py
        self.window.FiltTrace.triggered.connect(self.utilities.filterTrace) # In refraData.py
        self.window.Air_wave_fk.triggered.connect(self.utilities.airWaveFilter) # In refraData.py
        self.window.Vel_filter.triggered.connect(self.utilities.velocityFilter) # In refraData.py
# Actions for menu Mute
        self.window.MuteTrace.triggered.connect(self.window.traceMute) # In refraPlot.py
        self.window.MuteAir.triggered.connect(self.window.muteAir) # In refraPlot.py
        self.window.MuteBefore.triggered.connect(self.window.muteBefore) # In refraPlot.py
        self.window.MuteAfter.triggered.connect(self.window.muteAfter) # In refraPlot.py

        self.window.fig_dict = {} #Prepare dictionary for data sets to be plotted

        self.window.mplfigs.itemClicked.connect(self.window.changeFig) # Allow for changing of data set to be plotted
        # self.window.Envelopes.setEnabled(True)

# Intercept CTL-C to exit in a controlled way
        signal(SIGINT, self.Handler)
# Plot data of first shot point
        self.window.plotSG()

    def eventFilter(self, obj, event):
        """
        Function interprets keyboard events. It reacts only on key release !

        Parameters
        ----------
        obj : TYPE Object
            The object that contains the window where the cursor should be so
            the keystroke is captured, i.e. Main.window

        event : Qt event
            Any QT event, although only key-release events are treated.

        Returns
        -------
        None.

        """
        from PyQt5.QtGui import QWindow
        if (event.type() == QtCore.QEvent.KeyRelease and obj is self.window):
# In main function, only "+" and "-" keys have a meaning, changing amplitude
            if self.function == "main":
                if event.key() == 43:
                    self.window.changeAmp(1)
                elif event.key() == 45:
                    self.window.changeAmp(-1)
                elif event.key() == 16777234:
#     If Left arrow is pressed, choose next block to be plotted to the left
                    self.window.nextBlock(-1)
                elif event.key() == 16777236:
#     If Right arrow is pressed, choose next block to be plotted to the right
                    self.window.nextBlock(1)
# in pickMove function, the keyboard arrows are checkes with possible SHFT and CTRL
            elif self.function == "pick_move":
                if event.key() == 16777234:
#     If Left arrow is pressed, choose next trace to the left
                    self.window.shiftPickTrace(-1)
                elif event.key() == 16777236:
#     If Right arrow is pressed, choose next trace to the right
                    self.window.shiftPickTrace(1)
#     If Up key is pressed, move pick upwards
                elif event.key() == 16777235:
                    if event.modifiers() & Qt.ShiftModifier:
#       With shift by 10 samples
                        self.window.movePick(10)
                    elif event.modifiers() & Qt.ControlModifier:
#       With CTRL by 1 milisecond
                        nshift = 0.001/self.data.dt
                        self.window.movePick(nshift)
                    else:
#       Without modifier by 1 sample
                        self.window.movePick(1)
#     If Down key is pressed, move pick upwards
                elif event.key() == 16777237:
                    if event.modifiers() & Qt.ShiftModifier:
#       With shift by 10 samples
                        self.window.movePick(-10)
                    elif event.modifiers() & Qt.ControlModifier:
#       With CTRL by 1 milisecond
                        nshift = 0.001/self.data.dt
                        self.window.movePick(-nshift)
                    else:
#       Without modifier by 1 sample
                        self.window.movePick(-1)
# in pickMove function, the keyboard arrows are checkes with possible SHFT and CTRL
            elif self.function == "uncertainty":
                if event.key() == 16777234:
#     If Left arrow is pressed, choose next trace to the left
                    self.window.shiftPickTrace(-1)
                elif event.key() == 16777236:
#     If Right arrow is pressed, choose next trace to the right
                    self.window.shiftPickTrace(1)
#     If Up key is pressed, move pick upwards
                elif event.key() == 16777235:
                    if event.modifiers() & Qt.ShiftModifier:
#       With shift by 10 samples
                        self.window.changeUnc(10)
                    elif event.modifiers() & Qt.ControlModifier:
#       With CTRL by 1 milisecond
                        nshift = 0.001/self.data.dt
                        self.window.changeUnc(nshift)
                    else:
#       Without modifier by 1 sample
                        self.window.changeUnc(1)
#     If Down key is pressed, move pick upwards
                elif event.key() == 16777237:
                    if event.modifiers() & Qt.ShiftModifier:
#       With shift by 10 samples
                        self.window.changeUnc(-10)
                    elif event.modifiers() & Qt.ControlModifier:
#       With CTRL by 1 milisecond
                        nshift = 0.001/self.data.dt
                        self.window.changeUnc(-nshift)
                    else:
#       Without modifier by 1 sample
                        self.window.changeUnc(-1)
            elif self.function == "vfilt":
                if self.window.verticalSlider.isVisible() and event.key() == 16777220:
                    self.window.verticalSlider.setVisible(False)
# Check whether the letter C has been types in tomography result window
#       No idea why, but after typing C the first time, the event gets associated to
#       2 types of objects, rP.newWindow and PyQt5.QtGui.QWindow and the dialogue
#       window appears again after plotting the tomography results with the new
#       settings. Only when the object type is QWindow, the event should be
#       accepted.
# This works as long as only one key is pressed (so, do not press SHFT+C), just "c"
        elif event.type() == QtCore.QEvent.KeyRelease and self.function == "inver":
            if type(obj) is QWindow:
                if event.key() == 67 or event.key() == 16777248:
                    self.utilities.invCol()
        return QtWidgets.QWidget.eventFilter(self, obj, event)

    def fileOpen(self):
        """
        Gets the name of all data files to be treated, opens them dans stores
        the data in list st (one Stream per data file in class Data),
        distributes trace pointers into class Traces and checks whether measured
        and calculated picks are available in file picks.dat and picks_calc.dat
        respectively. If so, picks are distributed into the container for
        traces. In addition, the function reads geometry information from files
        shots.geo and receivers.geo and stores them in class Geometry.

        Returns
        -------
        None.

        """
        self.direction_start = np.array(["S","SSW","SW","WSW","W","WNW","NW","NNW",\
                          "N","NNE","NE","ENE","E","ESE","SE","SSE"], dtype=str)
        self.direction_end = np.array(["N","NNE","NE","ENE","E","ESE","SE","SSE",\
                        "S","SSW","SW","WSW","W","WNW","NW","NNW"], dtype=str)
        if os.path.exists("PyRefra.config"):
            with open("PyRefra.config","r") as fo:
                self.dir_flag = True
                self.title = fo.readline().split("\n")[0]
                print(f"Title: {self.title}")
                self.dir_start = fo.readline().split("\n")[0]
                self.config = True
            if len(np.where(self.direction_start==self.dir_start)[0]) > 0:
                self.dir_end = self.direction_end[\
                        np.where(self.direction_start==self.dir_start)[0][0]]
                print(f"Profile direction: {self.dir_start}-{self.dir_end}")
            else:
                self.dir_start = ""
                self.dir_end = ""
                print(f"Profile direction in file PyRefra incorrect: *{self.dir_start}*")
        else:
            print("File PyRefra.config not found")
            results, okButton = self.dialog(\
                    ["General title (name of profile...)",\
                     "Geographic direction of profile start",\
                     self.direction_start],
                    ["e","l","b"],\
                    [self.dir0,None,None],"PyRefra configuration")
            if okButton:
                self.title = results[0]
                self.dir_start = self.direction_start[int(results[2])]
                self.dir_end = self.direction_end[int(results[2])]
                self.config = True
                with open("PyRefra.config","w") as fo:
                    fo.write(f"{self.title}\n")
                    fo.write(self.dir_start)
                print(f"{self.title}\n"+\
                      "Profile direction: {self.dir_start}-{self.dir_end}")
            else:
                self.title = self.dir0
                self.dir_start = ""
                self.dir_end = ""
                self.config = False
        self.files = rD.Files(self.dir0)
        self.function = "file_open"
        self.files.get_files()
        self.geo = rD.Geometry()
        self.geo.readGeom()
        self.data = rD.Data(self)
        self.data.readData(self.files)
        self.traces = rD.Traces(self, self.data, self.geo)
        print("")

    def Handler(self,signal_received,frame):
        """
        Handles CTRL-C key stroke

        Parameters
        ----------
        signal_received : CTRL-C keystroke
        frame : No idea

        Returns
        -------
        None.

        """
        self.close_app()

    def dialog(self, labels, types, values, title="Title"):
        """
        Wrapper for class Dialog. Two buttons are shown for finishing: Ok and
        Cancel

        Parameters
        ----------
        labels : list of strings
            Explanatory text put beside the fields foe data entrance.
            If values==None, label is interpreted as information/Warning text
            For a series of radiobuttons or Comboboxes, all corresponding labels
            are given as a list within the list.
        types : list of strings (length as labels).
            Possible values:
                "c": to define a check box
                "e": to define an editable box (lineEdit)
                "l": To define a simple label, comment text
                'r': to define a series of radiobuttons (the corresponding list of
                     labels is considered as one single label for the numbering of
                     labels, types and values)
                "b": to define a combobox (dropdown menu)
        values : list of texts, numbers or Nones (length as labels)
            Initial values given to editable fields. Should be None for check boxes
                should be the number (natural numbering, starting at 1, not at 0)
                of the radiobutton to be activated by default.
                For labels, it may be "b" (bold text), "i" (italic) or anything
                else, including None for standard text. Not used for combo box
        title : str, default: "Title"
            Title of the dialog box.

        Returns
        -------
        results : list of str
            Response of each data entrance field. Should be transformed to the
            needed data format (int, float...) in the calling function
    		For radiobuttons, the returned value indicates the number of the
    		active button (counting starts at 0). For checkboxes, the returned
    		value is -1 if the box was not checked and the position at which
    		the box was checked (starting at 0) if it has been checked
        Dbutton: bool
            If True, "Apply" button has been pressed to finish dialog, if False
            "Cancel" button has been pressed.

        """
        D = Dialog(self, labels, types, values, title)
        D.Dfinish = False
        while (D.Dfinish != True):
            QtCore.QCoreApplication.processEvents()

        results = []
        iline = 0
        l = len(values)
        if l > 0:
            for it,t in enumerate(types):
                if t.lower() == "e":
                    results.append(D.dlines[iline].text())
                    iline += 1
                elif t.lower() == "r":
                    results.append(None)
                    for i in range(len(labels[it])):
                        iline += 1
                        if D.rbtn[it][i].isChecked():
                            results[-1] = i
#                            break
                elif t.lower() == 'c':
                    results.append(D.ck_order[it]-1)
                    iline += 1
                elif t.lower() == "b":
                    results.append(D.combo[it].currentIndex())
                    iline += 1
                else:
                    results.append(None)
                    iline += 1
        return results, D.Dbutton

    def closeApp(self):
        """
        Finishes application:
            Stores picks into file pick.dat
            Stores information concerning possibly modified traces (muted or
                sign-inversed)
            Deletes unneeded folder if tomography was calculated
            Closes window
        """
        import os
        choice = QtWidgets.QMessageBox.question(None, "Confirm", "Are you sure?",
                           QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            if max(self.traces.amplitudes) != self.window.general_sign:
                if os.path.exists("receivers_modify.dat"):
                    os.remove("receivers_modify.dat")
                with open("receivers_modify.dat","w") as fh:
                    for i in range(self.traces.number_of_traces):
                        if self.traces.amplitudes[i] != self.window.general_sign:
                            a = self.traces.amplitudes[i]*\
                                self.window.general_sign
                            rec = self.traces.receiver[i]
                            fh.write(f"{rec} {a:0.0f}\n")
            print("application finished")
            if self.utilities.w_tau:
                self.utilities.w_tau.close()
            if self.utilities.w_env:
                self.utilities.w_env.close()
            if self.utilities.w_tomo:
                self.utilities.w_tomo.close()
            if self.utilities.w_fcol:
                self.utilities.w_fcol.close()
            if self.utilities.w_amp:
                self.utilities.w_amp.close()
            if self.window.w_anim:
                self.window.w_anim.close()
            if self.window.w_picks:
                self.window.w_picks.close()
            if np.amax(self.traces.npick) > 0:
                self.traces.storePicks()
            try:
                os.rmdir(self.utilities.path)
            except:
                pass
            self.window.close()
            QtWidgets.QApplication.quit()
            return True
        else:
            pass


class Dialog(QtWidgets.QWidget):
    """
    Created on Tue Sep 13 12:50:26 2022
    @author: Hermann
    General purpose dialog box

    Input:
        parent (Class object): Certain values are needed from the parent class
                               usually, the call will thus be:
                               Dialog(self,labels,types,values,title)
        labels (list, string): labels explaining the different fields
                               For a series of radiobuttons, labels[i] is itself
                               a list
        types (list,string): In the order of "labels", indicates the type
                            of field the label belongs to. If labels[i] is itself
                            a list, this is considered as one single label.
                            May be
                            'l' for label (no dialog entry, only information)
                            'c' for checkbox
                            'e' for LineEdit (to enter values)
                            'r' for radiobutton
                            'b' for combobox (dropdown menu)
                            may be capital or not
        values (list, string, float or int): initial values for LineEdit fields,
                            number of radiobutton activated by default.
                            for labels, it may be "b" (bold), "i" (italic), or
                            anything else (usually None) for standard text
                            Ignored for checkbox.
                            Optional, default value: None
        title (string): Title of the dialog window
                            Optional, default value: "Title"
    """
    def __init__(self, parent, labels, types, values = None, title="Title"):
        super(Dialog, self).__init__()
        self.parent = parent
        nlab = len(labels)
        self.labels = labels
        self.ck_order = np.zeros(nlab, dtype = int)
        self.n_checked = 0
        self.Dfinish = False
        self.Dbutton = False
        self.dlabels = []
        self.dlines = []
        self.ckState = []
        self.ckb = []
        self.rbtn = []
        self.btngroup = []
        self.combo = []
        if parent.function == "False_Colour":
            self.label = QtWidgets.QLabel("Check up to 3 items; chose red item")
        for i in range(nlab):
            self.ckState.append(False)
        self.YesBtn = QtWidgets.QPushButton('Ok',self)
        self.YesBtn.move(10,20*(nlab+3))
        self.CancelBtn = QtWidgets.QPushButton('Cancel',self)
        self.CancelBtn.move(150,20*(nlab+3))
        self.mainLayout = QtWidgets.QGridLayout()
        self.setLayout(self.mainLayout)
        il_add = 0
        if parent.function == "False_Colour":
            self.mainLayout.addWidget(self.label, 0, 0, 1, 2)
            il_add = 1
        ilin = 0
        for i in range(len(labels)):
            il = ilin + il_add
            if types[i].lower() == 'l':
                if values[i]:
                    if values[i].lower() == 'b':
                        self.dlabels.append(QtWidgets.QLabel("<b>"+labels[i]+"</b>"))
                    elif values[i].lower() == 'i':
                        self.dlabels.append(QtWidgets.QLabel("<i>"+labels[i]+"</i>"))
                    else:
                        self.dlabels.append(QtWidgets.QLabel(labels[i]))
                else:
                    self.dlabels.append(QtWidgets.QLabel(labels[i]))
                self.dlines.append(None)
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.mainLayout.addWidget(self.dlabels[ilin], il, 0, 1, 2)
                ilin += 1
            elif types[i].lower() == 'e':
                self.dlabels.append(QtWidgets.QLabel(labels[i]))
                self.dlines.append(QtWidgets.QLineEdit())
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.mainLayout.addWidget(self.dlabels[ilin], il, 0, 1, 1)
                self.mainLayout.addWidget(self.dlines[ilin], il, 1, 1, 1)
                try:
                    s = str(values[i])
                    self.dlines[-1].setText(s)
                except:
                    pass
                ilin += 1
            elif types[i].lower() == 'r':
                self.ckb.append(None)
                self.combo.append(None)
                self.rbtn.append([])
                self.btngroup.append(QButtonGroup())
                rck = int(values[i])-1
                if rck<0 or rck>=len(labels[i]):
                    rck = 0
                for ir,l in enumerate(labels[i]):
                    self.dlabels.append(None)
                    self.dlines.append(None)
                    self.rbtn[i].append(QRadioButton(l))
                    self.btngroup[-1].addButton(self.rbtn[i][-1])
                    self.mainLayout.addWidget(self.rbtn[i][-1], il, 0, 1, 2)
                    if ir == rck:
                        self.rbtn[i][-1].setChecked(True)
                    else:
                        self.rbtn[i][-1].setChecked(False)
                    il += 1
                    ilin += 1
            elif types[i].lower() == 'c':
                self.dlabels.append(None)
                self.dlines.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.ckb.append(QtWidgets.QCheckBox(self))
                self.ckb[i].setText(self.labels[i])
                self.mainLayout.addWidget(self.ckb[i], il, 0, 1, 2)
                self.ckb[i].stateChanged.connect(self.checked)
                ilin += 1
            elif types[i].lower() == 'b':
                self.dlabels.append(None)
                self.dlines.append(None)
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(QtWidgets.QComboBox())
                for il,l in enumerate(labels[i]):
                    self.combo[i].addItem(l)
                ilin += 1
                self.mainLayout.addWidget(self.combo[i], ilin, 0, 1, 1)
                il_add += 1
        ilin += 2
        il = ilin + il_add
        self.mainLayout.addWidget(self.YesBtn, il, 0)
        self.mainLayout.addWidget(self.CancelBtn, il, 1)
        self.YesBtn.setDefault(True)
        self.YesBtn.clicked.connect(self.on_YesButton_clicked)
        self.CancelBtn.clicked.connect(self.on_CancelButton_clicked)

        self.setWindowTitle(title)
        self.show()

    def checked(self, checked):
        """
        Actions executed if a check box has changed its state.
        if a box has been checked, the function searches the one which was
        checked using self.ckState as indicator (this variable contains the
        state of all check boxes before the click) Its click-order is stored
        and self.ckState is changed. In addition, if the calling function is
        "False_Color", the explanation text is changed and the color to be
        used for the following clicked item is indicated. If in this case
        the third box is clicked, a message appears that no other box can be
        checked unless one is unchecked. If nevertheless a fourth box is
        checked, the corresponding box is automatically unchecked.
        If a box is unchecked, this is stored in self.ckState and the colors
        indicated are changed if necessary and if function is "False_Color".

        Parameters
        ----------
        checked : QtWidget.QCheckBox.checkState
            state of a checkbox after clicking into it

        Returns
        -------
        None.

        """
        cols = [" (red)", " (green)", " (blue)"]
# If a nex box is checked, search the one which has been checked
        if checked == Qt.Checked:
# If alread 3 boxes were checked, undo the checking
            if self.parent.function == "False_Colour":
                if self.n_checked >= 3:
                    for j in range(len(self.ck_order)):
                        if self.ck_order[j] == 0:
                            self.ckb[j].setChecked(False)
                            self.ckState[j] = False
                    return
# Else do necessary changes
            for i,ck in enumerate(self.ckb):
# If self.ckb has None value, the corresponding entry is not a checkbox
                if not ck:
                    continue
# If self.ckb.checkState is checked after click, set ckState to True and do
#    changes
                if ck.checkState()==Qt.Checked:
                    self.ckState[i] = True
# If checkbox nr i was not checked, increase the number of checked boxes
#    (n_checked) and store the order of checkin in self.ck_order
                    if self.ck_order[i] == 0:
                        self.n_checked += 1
                        self.ck_order[i] = self.n_checked
# if checkboxes are called from function "False_Color", change the label,
#    indicating the clor which will be used for the corresponding indicator.
                        if self.parent.function == "False_Colour":
                            self.ckb[i].setText(self.labels[i]+\
                                        cols[self.ck_order[i]-1])
# If the third box has been checked, give a warning, if not write the color
#    which will be used for the following check
                            if self.n_checked == 3:
                                self.label.setText(\
                                    "Before checking another item unckeck one item")
                            else:
                                self.label.setText("Check up to 3 items;"+\
                                    f"Chose {cols[self.n_checked]} item")
                        else:
                            self.ckb[i].setText(f"{self.labels[i]} ("+\
                                                f"{self.ck_order[i]})")
                        break
# If self.ckb.checkState is still unchecked, set ckState to Falsee
                else:
                    self.ckState[i] = False
# If click has unchecked a checkbox, do necessary changes
        else:
            for i,ck in enumerate(self.ckb):
# If self.ckb has None value, the corresponding entry is not a checkbox
                if not ck:
                    continue
# If self.ckb.checkState is still checked, set ckState to True
                if ck.checkState()==Qt.Checked:
                    self.ckState[i] = True
# If checkbox is no longer checked but it was (self?ckState), the unchecked box
#    is found
                else:
                    if self.ckState[i] == True:
                        self.ckState[i] = False
                        n = self.ck_order[i]
# reset ck_order to 0 (indicating also unchecked box)
                        self.ck_order[i] = 0
# Reset label to initial value (changes only for function "False_Color")
                        self.ckb[i].setText(self.labels[i])
# For all boxes that were checked later than the unchecked one, reduce their
#    checking order by 1 and, if function is False_Color, indicate the new
#    color uses for plotting
                        for j in range(len(self.ck_order)):
                            if self.ck_order[j] > n:
                                self.ck_order[j] -= 1
                                if self.parent.function == "False_Colour":
                                    self.ckb[j].setText(self.labels[j]+\
                                                   cols[self.ck_order[j]-1])
                                else:
                                    self.ckb[j].setText(f"{self.labels[j]} ("+\
                                                       f"{self.ck_order[j]})")
                        self.n_checked -= 1
                        if self.parent.function == "False_Colour":
                            self.label.setText("Check up to 3 items\n"+\
                                f"Chose {cols[self.n_checked]} item")
                        break
        self.show()

    def on_YesButton_clicked(self):
        n_checked = 0
        for ck in self.ckb:
            if not ck:
                continue
            if ck.checkState() == Qt.Checked:
                n_checked += 1
        self.Dfinish = True
        self.Dbutton = True
        self.close()

    def on_CancelButton_clicked(self):
        self.Dfinish = True
        self.Dbutton = False
        self.close()



class Seg2_Slide(QtWidgets.QSlider):
    """
    Created on Mon Dec 28 12:50:26 2020
    @author: Hermann
    Funtion calls a slider to define cut-off velocity in an f-k spectrum
    A line indicating the actaul velocity with the value is plotted on the spectrum,
    following the slider.

    Input:
        parent class, should be self in the call
        Initial cut-off velocity [m/s]
    """
    def __init__(self, parent,pos_ini):
        super(Seg2_Slide, self).__init__()
#        vmin = -2000
        vmin = 0
        vmax = 2000
# In this application, the sliding bar is only called from the FK-filter module
#    which is located in the Utility class. Therefore, "parent" is Main.utility.
        self.main = parent
        self.fig = self.main.window.fig_plotted
        self.figure = self.main.window.figs[self.fig]
        self.background = self.figure.canvas.copy_from_bbox(self.figure.bbox) # copy background picture
        self.ax = self.main.window.axes[self.fig]
        self.setMouseTracking(True)
        self.position = pos_ini
# fmax_plot and kmax_plot are maximum frequency and spatial frequency plotted
#    in the window. xmax and ymax, as well as x and y are used for plotting the
#    f-k line of the actual cutting velocity
        self.xmax = self.main.kmax_plot
        self.ymax = self.main.fmax_plot
        y1 = self.xmax*abs(self.position)
        if y1 <= self.ymax:
            if self.position > 0:
                x = [0,self.xmax]
            else:
                x = [0,-self.xmax]
            y = [0,y1]
        else:
            x = [0,self.ymax/self.position]
            y = [0,self.ymax]
# xt and yt are the positions where to plot the actual velocity value beside the
#    line [x,y]
        if self.position > 0:
            xt = x[1]-self.xmax*0.05
        else:
            xt = x[1]+self.xmax*0.25
        yt = y[1]-self.ymax*0.01
# Plot initial velocity and activate animated line in f-k plot
        self.line, = self.ax.plot(x,y)
        self.text = self.ax.text(xt,yt,f"{self.position:0.0f} m/s",\
                    ha='right',va='top',fontsize='xx-large')
        self.canvas = self.line.figure.canvas
        self.axl = self.line.axes
        self.axl.draw_artist(self.line)
        self.ax.draw_artist(self.text)
# Activate slider
        self.sld = QSlider(Qt.Horizontal, self)
        self.sld = self.main.window.verticalSlider
        self.sld.setVisible(True)
# Set the maximum available cutting velocity to arbitrary 2000 m/s
        self.sld.setRange(vmin, vmax)
        self.sld.setValue(int(pos_ini))
        self.sld.setTickPosition(QSlider.TicksBothSides)
        self.sld.setTickInterval(100)
        self.sld.setSingleStep(10)
        self.sld.setPageStep(100)
        self.label = QLabel(f"{int(vmin)}", self)
        self.label.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
        self.sld.valueChanged.connect(self.updateLabel)
        self.sld.sliderReleased.connect(self.released)

    def updateLabel(self):
        self.label.setText(str(self.sld.value()))
        self.position = self.sld.value()
        y1 = self.xmax*abs(self.position)
        if y1 <= self.ymax:
            if self.position > 0:
                x = [0,self.xmax]
            else:
                x = [0,-self.xmax]
            y = [0,y1]
        else:
            x = [0,self.ymax/self.position]
            y = [0,self.ymax]
        if self.position > 0:
            xt = x[1]-self.xmax*0.05
        else:
            xt = x[1]+self.xmax*0.25
        yt = y[1]-self.ymax*0.01
        self.text.set_position([xt,yt])
        self.text.set_text(f"{self.position:0.0f} m/s")
        self.axl = self.line.axes
        self.line.set_data(x,y)
        self.figure.canvas.restore_region(self.background)
        self.axl.draw_artist(self.line)
        self.ax.draw_artist(self.text)
        self.canvas.blit(self.axl.bbox)

    def released(self):
        self.sld.setVisible(False)


def my_exception_hook(exctype, value, tracebk):
    print(exctype, value, tracebk)
    sys._excepthook(exctype, value, tracebk)
    sys.exit(1)

if __name__ == "__main__":

    try:
        app = QtWidgets.QApplication(sys.argv)
        sys._excepthook = sys.excepthook
        sys.excepthook = my_exception_hook
        main = Main(dir0)
        main.window.showMaximized()
        sys.exit(app.exec_())
    except:
#        sys.exit()
        pass

