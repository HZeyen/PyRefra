# -*- coding: utf-8 -*-
"""
last modified on Nov 25, 2024

@author: Hermann Zeyen, University Paris-Saclay, France
         hermann.zeyen@universite-paris-saclay.fr

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
    colorcet

Needed private packages:

    refraData.py
    refraPlot.py
    refraWindow.ui


Contains the following classes:
    Main
        With the following functions:
            __init__
            test_function
            event_filter
            file_open
            handler
            dialog
            close_app

    Dialog
        With the following functions:
            __init__
            checked
            on_yes_button_clicked
            on_cancel_button_clicked

    Seg2_Slide
        With the following functions:
            __init__
            updateLabel

"""

import os
from signal import signal, SIGINT
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QSlider, QLabel, QRadioButton, QButtonGroup)
from PyQt5.QtGui import QWindow
import numpy as np
from .data.files import Files
from .data.data import Data
from .data.geometry import Geometry
from .data.traces import Traces
from .utilities.utilities import Utilities
from .plotting.refraPlot import Window

# # The next two lines may have to be modified:
# #  SYS_PATH is the folder where all python program files are located
# #  DIR0 is the folder where the data are located

# # Example of paths for programs and data on HZ desktop
# SYS_PATH = r"E:/Sources_2010/Python_programs"
# # DIR0 = r"E:/Seg2Dat/Fontaines-Salees/2021/2021-10-17_Profil5"
# # DIR0 = r"E:/Seg2Dat/Brigaud/Beaufremont"
# # DIR0 = r"E:/Seg2Dat/Feroes/Eidi_21_07_23"
# # DIR0 = r"E:/Seg2Dat/Erreurs/Chidozie"
# # DIR0 = r"E:/Seg2Dat/Campus_haut"
# # DIR0 = r"E:/Seg2Dat/Rouvres2"
# DIR0 = r"E:\Seg2Dat\Alina\L3_2024"

# # Example of paths for Linux
# # sys_path = r"/home/zeyen/src/Python"
# # DIR0 = r"/home/zeyen/Seismics/Fontaines_salees/2021-10-17_Profil5"

# if SYS_PATH not in sys.path:
#     sys.path.append(SYS_PATH)


class Main(QtWidgets.QWidget):
    """
    Main class for PyRefra

    """

    def __init__(self, top, dir0=None):
        """
        Initialization of pyrefra Main class

        Parameters
        ----------
        top : str
            Base folder for program pyrefra
        dir0 : str, optional; default: None
            Data folder if given

        Returns
        -------
        None.

        """
        super(Main, self).__init__()

        self.function = "main"
        self.top = top
        try:
            os.chdir(dir0)
            self.dir0 = dir0
            self.dir_flag = True
            # self.sys_path = SYS_PATH
        except NameError:
            self.dir0 = None
            # self.sys_path = SYS_PATH
            pass

# Input data
        self.file_open()
# Initialize main frame widget
        self.window = Window(self, self.files, self.data, self.traces,
                             self.geo)
        if len(self.geo.types) > 1:
            self.window.component.setEnabled(True)
# Check whether measured picks exist and if so, activate tomography button
        sp = sum(self.traces.npick)
        if sp > 0:
            self.window.Tomography.setEnabled(True)
            print(f"\n{sp} measured picks read")
# Initialize utility routines
        self.utilities = Utilities(self, self.files, self.data, self.traces,
                                   self.geo, self.window)
        QtWidgets.qApp.installEventFilter(self)

# Define actions for Menu buttons
# Actions for Menu File
        self.window.openFile.triggered.connect(self.file_open)
        self.window.Save_SEGY.triggered.connect(self.data.saveSEGY)
        self.window.Save_SU.triggered.connect(self.data.saveSU)
        self.window.Save_SEG2.triggered.connect(self.data.saveSEG2)
        self.window.Save_Bin.triggered.connect(self.data.saveBinary)
        self.window.Save_ASCII.triggered.connect(self.data.saveASCII)
        self.window.Save_headers.triggered.connect(self.data.saveHeader)
        self.window.Save_plot.triggered.connect(self.window.savePlot)
        self.window.quitAction.triggered.connect(self.close_app)
# Actions for menu Display
        self.window.originalDataScreen.triggered.connect(
            self.window.originalScreen)
        self.window.originalDataAll.triggered.connect(self.window.original)
        self.window.shotGather.triggered.connect(self.window.plotSG)
        self.window.fileGather.triggered.connect(self.window.plotFG)
        self.window.receiverGather.triggered.connect(self.window.plotRG)
        self.window.distanceGather.triggered.connect(self.window.plotDG)
        self.window.component.triggered.connect(self.window.chooseComponent)
        self.window.phaseAngles.triggered.connect(self.window.phasePlot)
        self.window.zoom.triggered.connect(self.window.zooming)
        self.window.zoom_Out.triggered.connect(self.window.zoomOut)
        self.window.zoom_In.triggered.connect(self.window.zoomIn)
        self.window.zoom_Initial.triggered.connect(self.window.zoomIni)
        self.window.t_Norm.triggered.connect(self.window.tNorm)
        self.window.t_Gain.triggered.connect(self.window.tGain)
        self.window.d_Gain.triggered.connect(self.window.dGain)
        self.window.Agc.triggered.connect(self.window.AGC)
# Actions for menu Utilities
        self.window.PModel.triggered.connect(self.utilities.pModel)
        self.window.SModel.triggered.connect(self.utilities.sModel)
        self.window.Tomography.triggered.connect(self.utilities.inversion)
        self.window.Checker.triggered.connect(self.utilities.checkerboard)
        self.window.Envelopes.triggered.connect(self.utilities.envel)
        self.window.TauP.triggered.connect(self.utilities.tauP)
        self.window.falseCol.triggered.connect(self.utilities.falseColour)
        self.window.TraceSign.triggered.connect(self.window.traceSign)
        self.window.ChangeSign.triggered.connect(self.window.changeSign)
        self.window.Change_colors.triggered.connect(self.utilities.invCol)
        self.window.Animation.triggered.connect(self.window.animateLine)
        self.window.Attenuation.triggered.connect(self.utilities.atten_amp)
        self.window.Pseudo_velocity.triggered.connect(self.utilities.pseudo)
# Actions for menu Picking
        self.window.ManualPicks.triggered.connect(self.window.pickManual)
        self.window.AmpPicks.triggered.connect(self.window.ampPicks)
        self.window.StaLta.triggered.connect(self.window.Sta_Lta)
        self.window.CorrelationPicks.triggered.connect(self.window.corrPick)
        self.window.MovePicks.triggered.connect(self.window.pickMove)
        self.window.Uncertainty_change.triggered.connect(
            self.window.uncertainty)
        self.window.erasePicks.triggered.connect(self.window.erase_Picks)
        self.window.plotAllPicks.triggered.connect(self.window.plotPickSection)
        self.window.PlotCalculatedTimes.triggered.connect(
            self.window.plotCalcPicks)
        if self.traces.calc_picks:
            self.window.PlotCalculatedTimes.setEnabled(True)
        self.window.StorePicks.triggered.connect(self.traces.storePicks)
        self.window.StoreGimli.triggered.connect(self.traces.saveGimli)
# Actions for menu Filter
        self.window.FrequencyFilter.triggered.connect(self.utilities.filterAll)
        self.window.FiltTrace.triggered.connect(self.utilities.filterTrace)
        self.window.Air_wave_fk.triggered.connect(self.utilities.airWaveFilter)
        self.window.Vel_filter.triggered.connect(self.utilities.velocityFilter)
# Actions for menu Mute
        self.window.MuteTrace.triggered.connect(self.window.traceMute)
        self.window.MuteAir.triggered.connect(self.window.muteAir)
        self.window.MuteBefore.triggered.connect(self.window.muteBefore)
        self.window.MuteAfter.triggered.connect(self.window.muteAfter)
# Prepare dictionary for data sets to be plotted
        self.window.fig_dict = {}

# Allow for changing of data set to be plotted
        self.window.mplfigs.itemClicked.connect(self.window.changeFig)
        # self.window.Envelopes.setEnabled(True)

# Intercept CTL-C to exit in a controlled way
        signal(SIGINT, self.handler)
# Plot data of first shot point
        self.window.plotSG()

    def test_function(self):
        """
        Test whether a call is made from Main or from within another function.
        The test is done on the value of self.function, which should be "main".
        If the call comes from a function that has not been correctly finished,
        the program may return to tha function to finish it correctly or reset
        the value of self.function to "main" and leave the calling function in
        an uncontrolled way.

        If calling function is a picking function, it is not allowed to leave
        the function in an uncontrolled way, since this will dammage pick
        arrays.

        Returns
        -------
        bool
            True if call is done from Main or if self.function is reset to
            "main", False else.

        """
        if self.function == "inver":
            self.function = "main"
            return True
        if not self.function == "main":
            if "pick" in self.function:
                answer = QtWidgets.QMessageBox.warning(
                    None, "Warning",
                    f"You did not finish the action {self.function}\n"
                    + "Close and finish function correctly\n",
                    QtWidgets.QMessageBox.Close, QtWidgets.QMessageBox.Close)
            else:
                answer = QtWidgets.QMessageBox.warning(
                    None, "Warning",
                    f"You did not finish the action {self.function}\n"
                    + "Close and finish function correctly or\n"
                    + "Ignore and leave function uncontrolled\n",
                    QtWidgets.QMessageBox.Ignore | QtWidgets.QMessageBox.Close,
                    QtWidgets.QMessageBox.Close)
            if answer == QtWidgets.QMessageBox.Close:
                return False
            else:
                self.function = "main"
        return True

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
# In main function, only "+" and "-" keys have a meaning, changing amplitude
        if (event.type() == QtCore.QEvent.KeyRelease and obj is self.window):
            if self.function == "main":
                if event.key() == 43:
                    self.window.changeAmp(1)
                elif event.key() == 45:
                    self.window.changeAmp(-1)
# If Left arrow is pressed, choose next block to be plotted to the left
                elif event.key() == 16777234:
                    self.window.nextBlock(-1)
# If Right arrow is pressed, choose next block to be plotted to the right
                elif event.key() == 16777236:
                    self.window.nextBlock(1)
# in pickMove function, the keyboard arrows are checkes with possible SHFT
# and CTRL
# If Left arrow is pressed, choose next trace to the left
            elif self.function == "pick_move":
                if event.key() == 16777234:
                    self.window.shiftPickTrace(-1)
# If Right arrow is pressed, choose next trace to the right
# If Up key is pressed, move pick upwards
                elif event.key() == 16777236:
                    self.window.shiftPickTrace(1)
# With shift by 10 samples
                elif event.key() == 16777235:
                    if event.modifiers() & Qt.ShiftModifier:
                        self.window.movePick(10)
# With CTRL by 1 milisecond
                    elif event.modifiers() & Qt.ControlModifier:
                        nshift = 0.001/self.data.dt
                        self.window.movePick(nshift)
#  Without modifier by 1 sample
                    else:
                        self.window.movePick(1)
# If Down key is pressed, move pick upwards
# With shift by 10 samples
                elif event.key() == 16777237:
                    if event.modifiers() & Qt.ShiftModifier:
                        self.window.movePick(-10)
# With CTRL by 1 milisecond
                    elif event.modifiers() & Qt.ControlModifier:
                        nshift = 0.001/self.data.dt
                        self.window.movePick(-nshift)
# Without modifier by 1 sample
                    else:
                        self.window.movePick(-1)
# in pickMove function, the keyboard arrows are checkes with possible SHFT
# and CTRL
# If Left arrow is pressed, choose next trace to the left
            elif self.function == "change_pick_uncertainty":
                if event.key() == 16777234:
                    self.window.shiftPickTrace(-1)
# If Right arrow is pressed, choose next trace to the right
                elif event.key() == 16777236:
                    self.window.shiftPickTrace(1)
# If Up key is pressed, move pick upwards
# With shift by 10 samples
                elif event.key() == 16777235:
                    if event.modifiers() & Qt.ShiftModifier:
                        self.window.changeUnc(10)
# With CTRL by 1 milisecond
                    elif event.modifiers() & Qt.ControlModifier:
                        nshift = 0.001/self.data.dt
                        self.window.changeUnc(nshift)
# Without modifier by 1 sample
                    else:
                        self.window.changeUnc(1)
# If Down key is pressed, move pick upwards
# With shift by 10 samples
                elif event.key() == 16777237:
                    if event.modifiers() & Qt.ShiftModifier:
                        self.window.changeUnc(-10)
# With CTRL by 1 milisecond
                    elif event.modifiers() & Qt.ControlModifier:
                        nshift = 0.001/self.data.dt
                        self.window.changeUnc(-nshift)
# Without modifier by 1 sample
                    else:
                        self.window.changeUnc(-1)
            elif self.function == "vfilt":
                if self.window.verticalSlider.isVisible() and \
                        event.key() == 16777220:
                    self.window.verticalSlider.setVisible(False)
# Check whether the letter C has been types in tomography result window
#       No idea why, but after typing C the first time, the event gets
#       associated to 2 types of objects, rP.newWindow and PyQt5.QtGui.QWindow
#       and the dialogue window appears again after plotting the tomography
#       results with the new settings. Only when the object type is QWindow,
#       the event should be accepted.
# This works as long as only one key is pressed (so, do not press SHFT+C),
# just "c"
# if event.key() == 67 or event.key() == 16777248:
        elif event.type() == QtCore.QEvent.KeyRelease and\
                self.function == "inver":
            if isinstance(obj, QWindow):
                if event.key() == 67:
                    self.utilities.invCol()
        return QtWidgets.QWidget.eventFilter(self, obj, event)

    def file_open(self):
        """
        Gets the name of all data files to be treated, opens them dans stores
        the data in list st (one Stream per data file in class Data),
        distributes trace pointers into class Traces and checks whether
        measured and calculated picks are available in file picks.dat and
        picks_calc.dat respectively. If so, picks are distributed into the
        container for traces. In addition, the function reads geometry
        information from files shots.geo and receivers.geo and stores them in
        class Geometry.

        Returns
        -------
        None.

        """
        self.direction_start = np.array(["S", "SSW", "SW", "WSW", "W", "WNW",
                                         "NW", "NNW", "N", "NNE", "NE", "ENE",
                                         "E", "ESE", "SE", "SSE"], dtype=str)
        self.direction_end = np.array(["N", "NNE", "NE", "ENE", "E", "ESE",
                                       "SE", "SSE", "S", "SSW", "SW", "WSW",
                                       "W", "WNW", "NW", "NNW"], dtype=str)
        if os.path.exists("PyRefra.config"):
            with open("PyRefra.config", "r") as fo:
                self.dir_flag = True
                self.title = fo.readline().split("\n")[0]
                print(f"Title: {self.title}")
                self.dir_start = fo.readline().split("\n")[0]
                self.config = True
            if len(np.where(self.direction_start == self.dir_start)[0]) > 0:
                self.dir_end = self.direction_end[
                        np.where(self.direction_start == self.dir_start)[0][0]]
                print(f"Profile direction: {self.dir_start}-{self.dir_end}")
            else:
                self.dir_start = ""
                self.dir_end = ""
                print("Profile direction in file PyRefra incorrect: "
                      + f"*{self.dir_start}*")
        else:
            print("File PyRefra.config not found")
            results, ok_button = self.dialog(
                    ["General title (name of profile...)",
                     "Geographic direction of profile start",
                     self.direction_start], ["e", "l", "b"],
                    [self.dir0, None, None], "PyRefra configuration")
            if ok_button:
                self.title = results[0]
                self.dir_start = self.direction_start[int(results[2])]
                self.dir_end = self.direction_end[int(results[2])]
                self.config = True
                with open("PyRefra.config", "w") as fo:
                    fo.write(f"{self.title}\n")
                    fo.write(self.dir_start)
                print(f"{self.title}\n"
                      + f"Profile direction: {self.dir_start}-{self.dir_end}")
            else:
                self.title = self.dir0
                self.dir_start = ""
                self.dir_end = ""
                self.config = False
        self.files = Files(self.dir0)
        self.function = "file_open"
        self.files.get_files()
        self.geo = Geometry()
        self.geo.readGeom()
        self.data = Data(self)
        self.data.readData(self.files)
        self.traces = Traces(self, self.data, self.geo)
        print("")
        self.function = "main"

    def handler(self, signal_received, frame):
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
            For a series of radiobuttons or Comboboxes, all corresponding
            labels are given as a list within the list.
        types : list of strings (length as labels).
            Possible values:
                "c": to define a check box
                "e": to define an editable box (lineEdit)
                "l": To define a simple label, comment text
                'r': to define a series of radiobuttons (the corresponding list
                     of labels is considered as one single label for the
                     numbering of labels, types and values)
                "b": to define a combobox (dropdown menu)
        values : list of texts, numbers or Nones (length as labels)
            Initial values given to editable fields. Should be None for check
                boxes should be the number (natural numbering, starting at 1,
                not at 0) of the radiobutton to be activated by default.
                For labels, it may be "b" (bold text), "i" (italic) or anything
                else, including None for standard text. Not used for combo box
        title : str, default: "Title"
            Title of the dialog box.

        Returns
        -------
        results : list of str
            Response of each data entrance field. Should be transformed to the
            needed data format (int, float...) in the calling function.
            For radiobuttons, the returned value indicates the number of the
            active button (counting starts at 0). For checkboxes, the returned
            value is -1 if the box was not checked and the position at which
            the box was checked (starting at 0) if it has been checked
        dbutton: bool
            If True, "Apply" button has been pressed to finish dialog, if False
            "Cancel" button has been pressed.

        """
        d = Dialog(self, labels, types, values, title)
        d.dfinish = False
        while not d.dfinish:
            QtCore.QCoreApplication.processEvents()

        results = []
        iline = 0
        n = len(values)
        if n > 0:
            for it, t in enumerate(types):
                if t.lower() == "e":
                    results.append(d.dlines[iline].text())
                    iline += 1
                elif t.lower() == "r":
                    results.append(None)
                    for i in range(len(labels[it])):
                        iline += 1
                        if d.rbtn[it][i].isChecked():
                            results[-1] = i
#                            break
                elif t.lower() == 'c':
                    results.append(d.ck_order[it]-1)
                    iline += 1
                elif t.lower() == "b":
                    results.append(d.combo[it].currentIndex())
                    iline += 1
                else:
                    results.append(None)
                    iline += 1
        return results, d.dbutton

    def close_app(self):
        """
        Finishes application:
            Stores picks into file pick.dat
            Stores information concerning possibly modified traces (muted or
                sign-inversed)
            Deletes unneeded folder if tomography was calculated
            Closes window
        """
        answer = self.test_function()
        if not answer:
            return
        choice = QtWidgets.QMessageBox.question(
            None, "Confirm", "Are you sure?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            if max(self.traces.amplitudes) != self.window.general_sign:
                if os.path.exists("receivers_modify.dat"):
                    os.remove("receivers_modify.dat")
                with open("receivers_modify.dat", "w") as fh:
                    for i in range(self.traces.number_of_traces):
                        if self.traces.amplitudes[i] !=\
                                self.window.general_sign:
                            a = self.traces.amplitudes[i]\
                                * self.window.general_sign
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
        return False


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
                               For a series of radiobuttons, labels[i] is
                               itself a list
        types (list,string): In the order of "labels", indicates the type
                            of field the label belongs to. If labels[i] is
                            itself a list, this is considered as one single
                            label.
                            May be
                            'l' for label (no dialog entry, only information)
                            'c' for checkbox
                            'e' for LineEdit (to enter values)
                            'r' for radiobutton
                            'b' for combobox (dropdown menu)
                            may be capital or not
        values (list, string, float or int): initial values for LineEdit
                           fields, number of radiobutton activated by default.
                            for labels, it may be "b" (bold), "i" (italic), or
                            anything else (usually None) for standard text
                            Ignored for checkbox.
                            Optional, default value: None
        title (string): Title of the dialog window
                            Optional, default value: "Title"
    """

    def __init__(self, parent, labels, types, values=None, title="Title"):
        super(Dialog, self).__init__()
        self.parent = parent
        nlab = len(labels)
        self.labels = labels
        self.ck_order = np.zeros(nlab, dtype=int)
        self.n_checked = 0
        self.dfinish = False
        self.dbutton = False
        self.dlabels = []
        self.dlines = []
        self.ck_state = []
        self.ckb = []
        self.rbtn = []
        self.btngroup = []
        self.combo = []
        if parent.function == "False_Colour":
            self.label = QtWidgets.QLabel(
                "Check up to 3 items; chose red item")
        for i in range(nlab):
            self.ck_state.append(False)
        self.yes_btn = QtWidgets.QPushButton('Ok', self)
        self.yes_btn.move(10, 20*(nlab+3))
        self.cancel_btn = QtWidgets.QPushButton('Cancel', self)
        self.cancel_btn.move(150, 20*(nlab+3))
        self.main_layout = QtWidgets.QGridLayout()
        self.setLayout(self.main_layout)
        il_add = 0
        if parent.function == "False_Colour":
            self.main_layout.addWidget(self.label, 0, 0, 1, 2)
            il_add = 1
        ilin = 0
        for i, lab in enumerate(labels):
            il = ilin + il_add
            if types[i].lower() == 'l':
                if values[i]:
                    if values[i].lower() == 'b':
                        self.dlabels.append(QtWidgets.QLabel("<b>"+lab+"</b>"))
                    elif values[i].lower() == 'i':
                        self.dlabels.append(QtWidgets.QLabel("<i>"+lab+"</i>"))
                    else:
                        self.dlabels.append(QtWidgets.QLabel(lab))
                else:
                    self.dlabels.append(QtWidgets.QLabel(lab))
                self.dlines.append(None)
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.main_layout.addWidget(self.dlabels[ilin], il, 0, 1, 2)
                ilin += 1
            elif types[i].lower() == 'e':
                self.dlabels.append(QtWidgets.QLabel(lab))
                self.dlines.append(QtWidgets.QLineEdit())
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.main_layout.addWidget(self.dlabels[ilin], il, 0, 1, 1)
                self.main_layout.addWidget(self.dlines[ilin], il, 1, 1, 1)
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
                if rck < 0 or rck >= len(lab):
                    rck = 0
                for ir, l in enumerate(lab):
                    self.dlabels.append(None)
                    self.dlines.append(None)
                    self.rbtn[i].append(QRadioButton(l))
                    self.btngroup[-1].addButton(self.rbtn[i][-1])
                    self.main_layout.addWidget(self.rbtn[i][-1], il, 0, 1, 2)
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
                self.ckb[i].setText(lab)
                self.main_layout.addWidget(self.ckb[i], il, 0, 1, 2)
                self.ckb[i].stateChanged.connect(self.checked)
                ilin += 1
            elif types[i].lower() == 'b':
                self.dlabels.append(None)
                self.dlines.append(None)
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(QtWidgets.QComboBox())
                for il, l in enumerate(lab):
                    self.combo[i].addItem(l)
                ilin += 1
                self.main_layout.addWidget(self.combo[i], ilin, 0, 1, 1)
                il_add += 1
        ilin += 2
        il = ilin + il_add
        self.main_layout.addWidget(self.yes_btn, il, 0)
        self.main_layout.addWidget(self.cancel_btn, il, 1)
        self.yes_btn.setDefault(True)
        self.yes_btn.clicked.connect(self.on_yes_button_clicked)
        self.cancel_btn.clicked.connect(self.on_cancel_button_clicked)

        self.setWindowTitle(title)
        self.show()

    def checked(self, checked):
        """
        Actions executed if a check box has changed its state.
        if a box has been checked, the function searches the one which was
        checked using self.ck_state as indicator (this variable contains the
        state of all check boxes before the click) Its click-order is stored
        and self.ck_state is changed. In addition, if the calling function is
        "False_Color", the explanation text is changed and the color to be
        used for the following clicked item is indicated. If in this case
        the third box is clicked, a message appears that no other box can be
        checked unless one is unchecked. If nevertheless a fourth box is
        checked, the corresponding box is automatically unchecked.
        If a box is unchecked, this is stored in self.ck_state and the colors
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
# If alread 3 boxes were checked, undo the checking
        if checked == Qt.Checked:
            if self.parent.function == "False_Colour" and self.n_checked >= 3:
                for j, co in enumerate(self.ck_order):
                    if co == 0:
                        self.ckb[j].setChecked(False)
                        self.ck_state[j] = False
                return
# Else do necessary changes
# If self.ckb has None value, the corresponding entry is not a checkbox
            for i, ck in enumerate(self.ckb):
                if not ck:
                    continue
# If self.ckb.checkState is checked after click, set ck_state to True and do
# changes
                if ck.checkState() == Qt.Checked:
                    self.ck_state[i] = True
# If checkbox nr i was not checked, increase the number of checked boxes
# (n_checked) and store the order of checkin in self.ck_order
                    if self.ck_order[i] == 0:
                        self.n_checked += 1
                        self.ck_order[i] = self.n_checked
# if checkboxes are called from function "False_Color", change the label,
# indicating the clor which will be used for the corresponding indicator.
                        if self.parent.function == "False_Colour":
                            self.ckb[i].setText(self.labels[i]
                                                + cols[self.ck_order[i]-1])
# If the third box has been checked, give a warning, if not write the color
# which will be used for the following check
                            if self.n_checked == 3:
                                self.label.setText(
                                    "Before checking another item unckeck "
                                    + "one item")
                            else:
                                self.label.setText(
                                    "Check up to 3 items; Chose "
                                    + f"{cols[self.n_checked]} item")
                        else:
                            self.ckb[i].setText(f"{self.labels[i]} ("
                                                + f"{self.ck_order[i]})")
                        break
# If self.ckb.checkState is still unchecked, set ck_state to Falsee
                else:
                    self.ck_state[i] = False
# If click has unchecked a checkbox, do necessary changes
# If self.ckb has None value, the corresponding entry is not a checkbox
        else:
            for i, ck in enumerate(self.ckb):
                if not ck:
                    continue
# If self.ckb.checkState is still checked, set ck_state to True
                if ck.checkState() == Qt.Checked:
                    self.ck_state[i] = True
# If checkbox is no longer checked but it was (self.ck_state), the unchecked
# box is found
                else:
                    if self.ck_state[i]:
                        self.ck_state[i] = False
                        n = self.ck_order[i]
# reset ck_order to 0 (indicating also unchecked box)
                        self.ck_order[i] = 0
# Reset label to initial value (changes only for function "False_Color")
                        self.ckb[i].setText(self.labels[i])
# For all boxes that were checked later than the unchecked one, reduce their
# checking order by 1 and, if function is False_Color, indicate the new
# color uses for plotting
                        for j, co in enumerate(self.ck_order):
                            if co > n:
                                co -= 1
                                if self.parent.function == "False_Colour":
                                    self.ckb[j].setText(self.labels[j]
                                                        + cols[co-1])
                                else:
                                    self.ckb[j].setText(
                                        f"{self.labels[j]} ("
                                        + f"{self.ck_order[j]})")
                        self.n_checked -= 1
                        if self.parent.function == "False_Colour":
                            self.label.setText(
                                "Check up to 3 items\n"
                                + f"Chose {cols[self.n_checked]} item")
                        break
        self.show()

    def on_yes_button_clicked(self):
        """
        Finish Dialogue pressing acceptation button
        """
        n_checked = 0
        for ck in self.ckb:
            if not ck:
                continue
            if ck.checkState() == Qt.Checked:
                n_checked += 1
        self.dfinish = True
        self.dbutton = True
        self.close()

    def on_cancel_button_clicked(self):
        """
        Finish Dialogue pressing cancel button
        """
        self.dfinish = True
        self.dbutton = False
        self.close()


class Seg2_Slide(QtWidgets.QSlider):
    """
    Created on Mon Dec 28 12:50:26 2020
    @author: Hermann
    Funtion calls a slider to define cut-off velocity in an f-k spectrum
    A line indicating the actaul velocity with the value is plotted on the
    spectrum, following the slider.

    Input:
        parent class, should be self in the call
        Initial cut-off velocity [m/s]
    """

    def __init__(self, parent, pos_ini):
        super(Seg2_Slide, self).__init__()
#        vmin = -2000
        vmin = 0
        vmax = 2000
# In this application, the sliding bar is only called from the FK-filter module
# which is located in the Utility class. Therefore, "parent" is Main.utility.
        self.main = parent
        self.fig = self.main.window.fig_plotted
        self.figure = self.main.window.figs[self.fig]
# copy background picture
        self.background = self.figure.canvas.copy_from_bbox(self.figure.bbox)
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
                x = [0, self.xmax]
            else:
                x = [0, -self.xmax]
            y = [0, y1]
        else:
            x = [0, self.ymax/self.position]
            y = [0, self.ymax]
# xt and yt are the positions where to plot the actual velocity value beside
# the line [x,y]
        if self.position > 0:
            xt = x[1]-self.xmax*0.05
        else:
            xt = x[1]+self.xmax*0.25
        yt = y[1]-self.ymax*0.01
# Plot initial velocity and activate animated line in f-k plot
        self.line, = self.ax.plot(x, y)
        self.text = self.ax.text(xt, yt, f"{self.position:0.0f} m/s",
                                 ha='right', va='top', fontsize='xx-large')
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
        """
        Rewrite label written along the velocity line

        Returns
        -------
        None.

        """
        self.label.setText(str(self.sld.value()))
        self.position = self.sld.value()
        y1 = self.xmax*abs(self.position)
        if y1 <= self.ymax:
            if self.position > 0:
                x = [0, self.xmax]
            else:
                x = [0, self.xmax]
            y = [0, y1]
        else:
            x = [0, self.ymax/self.position]
            y = [0, self.ymax]
        if self.position > 0:
            xt = x[1]-self.xmax*0.05
        else:
            xt = x[1]+self.xmax*0.25
        yt = y[1]-self.ymax*0.01
        self.text.set_position([xt, yt])
        self.text.set_text(f"{self.position:0.0f} m/s")
        self.axl = self.line.axes
        self.line.set_data(x, y)
        self.figure.canvas.restore_region(self.background)
        self.axl.draw_artist(self.line)
        self.ax.draw_artist(self.text)
        self.canvas.blit(self.axl.bbox)

    def released(self):
        """
        Close slider when mouse button released
        """
        self.sld.setVisible(False)
