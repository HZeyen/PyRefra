# -*- coding: utf-8 -*-
"""
Matplotlib/Qt widgets shared by the plotting modules.

Contains the following Class:
    newWindow
"""

from PyQt5.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from datetime import datetime, date
from PyQt5 import QtWidgets


class newWindow(QWidget):
    """
    This "window" is a QWidget. If it has no parent, it
    will appear as a free-floating window as we want.
    """

    def __init__(self, title):
        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.resize(1800, 1200)
        self.setWindowTitle(title)
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.layout.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout.addWidget(self.toolbar)


class CanvasMixin:
    """Canvas, figure list, help line and plot saving of the data window."""

    def addMPL(self, fig):
        """
        Add widget to actual Figure
        Parameters
        ----------
        fig : matplotlib.figure.Figure object
            Is the actually used Figure
        """
        self.canvas = FigureCanvas(fig)
        self.hover_ann = None
        self.hover_bg = None
        self.canvas.mpl_connect("draw_event", self.hoverSaveBackground)
        self.canvas.mpl_connect("motion_notify_event", self.hoverMotion)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.mplvl.addWidget(self.toolbar)
# Canvas and toolbar are added after an existing help label: move it back to
# the bottom
        help_label = getattr(self, "help", None)
        if help_label is not None:
            self.mplvl.removeWidget(help_label)
            self.mplvl.addWidget(help_label)

    def rmMPL(self):
        """
        Remove widget from actual Figure
        """
        if self.canvas:
            self.mplvl.removeWidget(self.canvas)
            self.canvas.close()
        if self.toolbar:
            self.mplvl.removeWidget(self.toolbar)
            self.toolbar.close()

    def addFig(self, name, fig):
        """
        add data set to list of plottable data sets (right side of screen)

        Parameters
        ----------
        name : str
            Name of the data set to be added (text)
        fig : matplotlib.figure.Figure object
            Actual figure
        """
        self.fig_dict[name] = fig
        self.mplfigs.addItem(name)

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
        answer = self.main.test_function()
        if not answer:
            return
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
        if self.main.function == "attenuation":
            tit = self.axes[self.fig_plotted].get_title()
            tits = tit.split()
            fname = f"Attenuation_{tits[1]}{tits[2][:-1]}_{d1}_{c_time}.png"
        else:
            if self.rg_flag:
                fname = "receiver_"\
                    + f"{self.actual_receivers[0]+1:0>5}_{d1}_{c_time}.png"
            if self.fg_flag:
                fil_num = int(self.plot_names[self.fig_plotted_F].split()[1])
                fname = f"file_{fil_num:0>5}_{d1}_{c_time}.png"
            elif self.sg_flag:
                fname = "shot_point_"\
                    + f"{self.actual_shot+1:0>5}_{d1}_{c_time}.png"
            elif self.dg_flag:
                fname = "distance_"\
                    + f"{int(round(self.actual_dist,0)):0>5}_{d1}_{c_time}.png"
            else:
                fname = f"file_{self.actual_shot:0>5}_{d1}_{c_time}.png"
            self.figs[self.fig_plotted].savefig(fname)

    def drawNew(self, plot_flag):
        """
        Draw a new data set
        input:
            plot_flag (Boolean): If True plot seismograms
                                 if False only plot frame
        """
# remove actual window
        self.rmMPL()
# remove actual Figure
        self.figs[self.fig_plotted].clear()
# create new Figure
        self.figs[self.fig_plotted] = Figure()
# add axes to Figure
        self.axes[self.fig_plotted] = self.figs[self.fig_plotted].add_subplot()
# Add Figure to central widget
        self.addMPL(self.figs[self.fig_plotted])
# Plot record section, depending on which type of plot is wanted
# (shot section, record section, distance section)
        if plot_flag:
            if self.dg_flag:
                self.plotDistance(self.axes[self.fig_plotted],
                                  self.distances[self.fig_plotted])
                if self.PlotCalculatedTimes.isChecked():
                    self.plotCalcPicks()
            elif self.rg_flag:
                self.plotReceiver(self.axes[self.fig_plotted],
                                  self.receivers[self.fig_plotted])
                if self.PlotCalculatedTimes.isChecked():
                    self.plotCalcPicks()
            elif self.fg_flag:
                self.plotFile(self.axes[self.fig_plotted], self.fig_plotted)
                if self.PlotCalculatedTimes.isChecked():
                    self.plotCalcPicks()
            else:
                self.plotShot(self.axes[self.fig_plotted], self.fig_plotted)
                if self.PlotCalculatedTimes.isChecked():
                    self.plotCalcPicks()

# Plot picks if available
            self.picksPlot(self.figs[self.fig_plotted],
                           self.axes[self.fig_plotted])
        self.main.function = "main"

    def setHelp(self, text):
        """
        Set help text at bottom of screen

        Parameters
        ----------
        Input : str
            text is the text to be printed (defined in __init__).
            Text is written in QLabel widget. In order to keep the widget at
            the bottom of the screen, the existing one is first erased and then
            reopened

        This call has to be done after any modification of the graphics window
        """
        self.mplvl.removeWidget(self.help)
        self.help.close()
        self.help = QtWidgets.QLabel(self)
        self.help.setMaximumHeight(15)
        self.mplvl.addWidget(self.help)
        self.help.show()
        self.help.setText(text)
        self.figs[self.fig_plotted].canvas.flush_events()
