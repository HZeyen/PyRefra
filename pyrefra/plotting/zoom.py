# -*- coding: utf-8 -*-
"""
Zoom handling of the data window.

Contains the following Class:
    ZoomMixin
"""

from matplotlib.figure import Figure


class ZoomMixin:
    """Zoom handling of the data window."""

    def addZoom(self, xmin, xmax, nt_min, nt_max):
        """
        add a new zoom to the list of zooms

        Parameters
        ----------
        xmin : float
            Minimum offset to be plotted.
        xmax : float
            Maximum offset to be plotted.
        nt_min : int
            First sample to be plotted.
        nt_max : int
            Last sample to be plotted.

        Returns
        -------
        None.

        """
        self.n_zooms += 1
        self.i_zooms += 1
        self.zooms.append([xmin, xmax, nt_min, nt_max])

    def changeZoom(self, i, xmin, xmax, nt_min, nt_max):
        """
        change to another existing zoom

        Parameters
        ----------
        xmin : float
            Minimum offset to be plotted.
        xmax : float
            Maximum offset to be plotted.
        nt_min : int
            First sample to be plotted.
        nt_max : int
            Last sample to be plotted.

        Returns
        -------
        None.

        """
        self.i_zooms = i
        self.zooms[i][:] = [xmin, xmax, nt_min, nt_max]

    def zooming(self):
        """
        Zoom by drawing a rectangle
        Input: none
        Output: Bool
            False if rectangle too small (usually zero size by inadvertancy)
        """
        answer = self.main.test_function()
        if not answer:
            return None

        self.main.function = "zoom"
        zm = False
# Change help text
        self.setHelp(self.zoom_text)
# Pull a rectangle over the window to define area to be plotted
        self.followRect()
        if len(self.coor_x) == 0:
            self.coor_x.append(None)
            self.coor_y.append(None)
        if self.coor_x[0] is None or self.coor_x[1] is None or\
                self.coor_y[0] is None or self.coor_y[1] is None:
            print("\nZoom coordinates not correctly detected. "
                  + "Zoom not executed")
# Change help text to self module
            self.setHelp(self.main_text)
            self.drawNew(True)
# Change help text to self module
            self.setHelp(self.main_text)
            return zm
        xmi = min(self.coor_x)
        xma = max(self.coor_x)
        tmi = min(self.coor_y)
        tma = max(self.coor_y)
        if tmi >= self.time_plt_max or tma <= self.time_plt_min:
            print("\nTime window outside actual zoom. Zoom not executed")
# Change help text to self module
            self.setHelp(self.main_text)
            self.drawNew(True)
# Change help text to self module
            self.setHelp(self.main_text)
            return zm
        if xmi >= self.x.max() or xma <= self.x.min():
            print("\nX window outside actual zoom. Zoom not executed")
# Change help text to self module
            self.setHelp(self.main_text)
            self.drawNew(True)
# Change help text to self module
            self.setHelp(self.main_text)
            return zm
        if tma-tmi < 10*self.data.dt:
            print("\nTime axis too much zoomed. Zoom not executed")
# Change help text to self module
            self.setHelp(self.main_text)
            self.drawNew(True)
# Change help text to self module
            self.setHelp(self.main_text)
            return zm
        if xma-xmi < 2*self.geom.dx_geo:
            print("\nX axis too much zoomed. Zoom not executed")
# Change help text to self module
            self.setHelp(self.main_text)
            self.drawNew(True)
# Change help text to self module
            self.setHelp(self.main_text)
            return zm
# Define start and end time of zoomed plot
        self.time_plt_min = tmi
        self.time_plt_max = tma
# Define start sample of zoomed plot
        self.nt_mn = max(int((self.time_plt_min-self.data.t0)/self.data.dt), 0)
# Define end sample of zoomed plot
        self.nt_mx = min(int((self.time_plt_max-self.data.t0)/self.data.dt),
                         max(self.traces.nsample_trace))
# remove actual window
        self.rmMPL()
# Remove actual Figure
        self.figs[self.fig_plotted].clear()
# create new Figure
        self.figs[self.fig_plotted] = Figure()
# add axes to Figure
        self.axes[self.fig_plotted] =\
            self.figs[self.fig_plotted].add_subplot(111)
# Add Figure to central widget
        self.addMPL(self.figs[self.fig_plotted])
# Add zoom limits to list of existing zooms to be used by "zoom_out"
#     actual values depend on type of sections to be plotted
# Then plot section with actual zoom parameters
        if self.dg_flag:
            self.dxmin = xmi
            self.dxmax = xma
            self.i_zooms += 1
            if self.i_zooms >= self.n_zooms:
                self.zooms.append([self.dxmin, self.dxmax, self.nt_mn,
                                   self.nt_mx])
                self.n_zooms += 1
            else:
                self.zooms[self.i_zooms][:] = [self.dxmin, self.dxmax,
                                               self.nt_mn, self.nt_mx]
            self.plotDistance(self.axes[self.fig_plotted],
                              self.distances[self.fig_plotted])
        elif self.rg_flag:
            self.rxmin = xmi
            self.rxmax = xma
            self.i_zooms += 1
            if self.i_zooms >= self.n_zooms:
                self.zooms.append([self.rxmin, self.rxmax, self.nt_mn,
                                   self.nt_mx])
                self.n_zooms += 1
            else:
                self.zooms[self.i_zooms][:] = [self.rxmin, self.rxmax,
                                               self.nt_mn, self.nt_mx]
            self.plotReceiver(self.axes[self.fig_plotted],
                              self.receivers[self.fig_plotted])
        elif self.fg_flag:
            self.sxmin = xmi
            self.sxmax = xma
            self.i_zooms += 1
            if self.i_zooms >= self.n_zooms:
                self.zooms.append([self.sxmin, self.sxmax, self.nt_mn,
                                   self.nt_mx])
                self.n_zooms += 1
            else:
                self.zooms[self.i_zooms][:] = [self.sxmin, self.sxmax,
                                               self.nt_mn, self.nt_mx]
            self.plotFile(self.axes[self.fig_plotted], self.fig_plotted)
        else:
            if self.phase_plot:
                self.v_set = False
            self.sxmin = xmi
            self.sxmax = xma
            self.i_zooms += 1
            if self.i_zooms >= self.n_zooms:
                self.zooms.append([self.sxmin, self.sxmax, self.nt_mn,
                                   self.nt_mx])
                self.n_zooms += 1
            else:
                self.zooms[self.i_zooms][:] = [self.sxmin, self.sxmax,
                                               self.nt_mn, self.nt_mx]
            self.plotShot(self.axes[self.fig_plotted], self.fig_plotted)
            if self.PlotCalculatedTimes.isChecked():
                self.plotCalcPicks()
# Plot picks if available
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
# Change help text to self module
        self.setHelp(self.main_text)
        self.main.function = "main"
        zm = True
        return zm

    def zoomOut(self):
        """
        If several zooms have been used, zoom out one step
        """
        answer = self.main.test_function()
        if not answer:
            return False

        if self.i_zooms == 0:
            print("No zoom out possible, already at initial zoom")
            return False
        self.i_zooms = max(self.i_zooms-1, 0)
# Set actual zoom parameters to earlier window
        self.nt_mn = self.zooms[self.i_zooms][2]
        self.nt_mx = self.zooms[self.i_zooms][3]
        self.time_plt_min = self.data.t0+self.nt_mn*self.data.dt
        self.time_plt_max = self.data.t0+self.nt_mx*self.data.dt
        if self.i_zooms == 0:
            self.zooms[0][0] = min(self.x)-1
            self.zooms[0][1] = max(self.x)+1
# draw data section with actual zoom paramters
        self.drawNew(True)
# Change help text to self module
        self.setHelp(self.main_text)
        return True

    def zoomIn(self):
        """
        If several zooms have been used, zoom in one step
        """
        answer = self.main.test_function()
        if not answer:
            return False

        if self.i_zooms == self.n_zooms:
            print("No tighter zoom available, use Z to define new zoom")
            return False
        self.i_zooms = min(self.i_zooms+1, self.n_zooms-1)
# Set actual zoom parameters to earlier window
        self.nt_mn = self.zooms[self.i_zooms][2]
        self.nt_mx = self.zooms[self.i_zooms][3]
        self.time_plt_min = self.data.t0+self.nt_mn*self.data.dt
        self.time_plt_max = self.data.t0+self.nt_mx*self.data.dt
# draw data section with actual zoom paramters
        self.drawNew(True)
# Change help text to self module
        self.setHelp(self.main_text)
        return True

    def zoomIni(self):
        """
        Go back to initial zoom (all data are plotted)
        """
        answer = self.main.test_function()
        if not answer:
            return

        self.i_zooms = 0
        self.nt_mn = self.zooms[0][2]
        self.nt_mx = self.zooms[0][3]
        self.time_plt_min = self.data.t0
        self.time_plt_max = self.data.t0+self.nt_mx*self.data.dt
        self.zooms[0][0] = min(self.x)-1
        self.zooms[0][1] = max(self.x)+1
# draw full data section
        self.drawNew(True)
    # Change help text to self module
        self.setHelp(self.main_text)

    def nextBlock(self, direction):
        """
        if zoom on x-axis is active, show next traces to teh right of left

        Parameters
        ----------
        direction : int
            if > 0 show traces to the right, else to the left

        Returns
        -------
        bool
            DESCRIPTION.

        """
        if self.i_zooms == 0:
            return
        d = int(direction)
        dwin = self.zooms[self.i_zooms][1]-self.zooms[self.i_zooms][0]
        if d > 0:
            xma = max(self.x)
            if self.zooms[self.i_zooms][1] >= xma:
                return
            zm1 = min(self.zooms[self.i_zooms][1]+dwin, xma+0.2)
            zm0 = zm1 - dwin
        else:
            xmi = min(self.x)
            if self.zooms[self.i_zooms][0] <= xmi:
                return
            zm0 = max(self.zooms[self.i_zooms][0]-dwin, xmi-0.2)
            zm1 = zm0 + dwin
        self.zooms[self.i_zooms][0] = zm0
        self.zooms[self.i_zooms][1] = zm1
# draw data section with actual zoom paramters
        self.drawNew(True)
# Change help text to self module
        self.setHelp(self.main_text)
