# -*- coding: utf-8 -*-
"""
Interactive editing of picks.

Contains the following Class:
    ManualPickMixin
"""

import numpy as np
from PyQt5 import QtWidgets, QtCore


class ManualPickMixin:
    """Interactive editing of picks."""

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
        answer = self.main.test_function()
        if not answer:
            return
        self.setHelp(self.pick_move_text)
        self.main.function = "pick_move"
        self.end = False
        self.coor_x = []
        self.coor_y = []
        figure = self.figs[self.fig_plotted]

        def onPress(event):
            """
            Action if mouse button was pressed
            """
# If left mouse button has been clicked, search nearest pick
            if event.button == 1:
                self.keys_held = set()
                result = self.searchPick(event.xdata, event.ydata)
                if not result:
                    self.keys_held = set()
                    self.drawNew(True)
                    self.picksPlot(self.figs[self.fig_plotted],
                                   self.axes[self.fig_plotted])
                    self.end = True
                    return
                self.update()
# Save position of chosen pick in memory
                time_back = self.traces.pick_times[self.itrace][self.ipk]
# Re-plot pick in white color (i.e. background color)
                # plt_tr = np.where(np.array(self.tr) == self.n_tr)[0][0]
                # plt_col = self.line_colors[plt_tr]
                self.pickPlot(self.n_tr, self.itrace, self.ipk, "w", v=False)
# Erase pick temporally to redraw canvas without this pick
                self.traces.pick_times[self.itrace][self.ipk] = -1
                self.picksPlot(self.figs[self.fig_plotted],
                               self.axes[self.fig_plotted])
                self.back = figure.canvas.copy_from_bbox(figure.bbox)
                self.axes[self.fig_plotted].figure.canvas.draw()
                self.traces.pick_times[self.itrace][self.ipk] = time_back
                self.x_c = []
                self.x_c.append(self.x[self.n_tr]-0.3)
                self.x_c.append(self.x[self.n_tr]+0.3)
                self.y_c = [time_back, time_back]
                self.lin, = self.axes[self.fig_plotted].plot(self.x_c,
                                                             self.y_c, "k")
                self.axes[self.fig_plotted].draw_artist(self.lin)
                self.canvas.blit(self.axes[self.fig_plotted].bbox)
# If right mouse has been clicked, plot all picks, including the moved
# one in red and leave routine
            elif event.button == 3:
                self.keys_held = set()
                self.drawNew(True)
                self.picksPlot(self.figs[self.fig_plotted],
                               self.axes[self.fig_plotted])
                self.end = True

# Plot lines in positive and negative direction corresponding to air wave
# arrivals
#         xx = np.arange(self.x_zoom_min, self.x_zoom_max)
        xx = self.x[(self.x >= self.x_zoom_min) & (self.x <= self.x_zoom_max)]
        index = np.argsort(xx)
        xx = xx[index]
        xxx = self.x_off[(self.x >= self.x_zoom_min) & (self.x <= self.x_zoom_max)]
        xxx = xxx[index]
        tt = np.abs(xxx)/345.
        xxx = xx[tt >= self.t_zoom_min]
        ttt = tt[tt >= self.t_zoom_min]
        plt_air = True
# Define a first pick just in case user does not choose a first one by mouse
# click and tries to change trace using keyboarch arrows
        _ = self.searchPick(self.x_zoom_min, self.t_zoom_min)
        if self.dg_flag:
            dist = self.distances[self.fig_plotted_D]
            x_air = xx
            t_air = np.ones_like(x_air)*dist/345.
            if t_air[0] > self.t_zoom_max or t_air[0] < self.t_zoom_min:
                plt_air = False
        else:
            index = np.where(ttt <= self.t_zoom_max)[0]
            if len(index) > 1:
                ind = [0]
                for ii in range(1, len(index)):
                    if index[ii]-index[ii-1] > 1:
                        ind.append(ii)
                ind.append(len(index))
            else:
                plt_air = False
            # x_air = xxx[ttt <= self.t_zoom_max]
            # t_air = ttt[ttt <= self.t_zoom_max]
            # if len(t_air) < 2:
            #     plt_air = False
        if plt_air:
            for i in range(len(ind)-1):
                # self.axes[self.fig_plotted].plot(x_air, t_air, ls="--", c="b")
                i0 = index[ind[i]]
                i1 = index[ind[i+1]-1]+1
                self.axes[self.fig_plotted].plot(xxx[i0:i1], ttt[i0:i1], ls="--", c="b")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)

        self.lin, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y,
                                                     animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect(
            'button_press_event', onPress)
        while self.end is False:
            QtCore.QCoreApplication.processEvents()
        self.setHelp(self.main_text)
        self.traces.storePicks()
        self.main.function = "main"

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
        figure = self.figs[self.fig_plotted]
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
        answer = self.main.test_function()
        if not answer:
            return
        self.setHelp(self.uncertainty_text)
        self.main.function = "change_pick_uncertainty"
        self.end = False
        self.coor_x = []
        self.coor_y = []
        figure = self.figs[self.fig_plotted]

        def onPress(event):
            """
            Action when mouse button was pressed
            """
# If left öouse button has been clicked, search nearest pick
            if event.button == 1:
                self.keys_held = set()
                self.searchPick(event.xdata, event.ydata)
                self.update()
# Save position of chosen pick in memory
                trace = self.itrace
                time_back = self.traces.pick_times[self.itrace][self.ipk]
# Erase pick temporally to redraw canvas without this pick
                self.traces.pick_times[self.itrace][self.ipk] = -1
                self.picksPlot(self.figs[self.fig_plotted],
                               self.axes[self.fig_plotted])
                self.back = figure.canvas.copy_from_bbox(figure.bbox)
                self.axes[self.fig_plotted].figure.canvas.draw()
                self.traces.pick_times[self.itrace][self.ipk] = time_back
                self.x_c = [self.x[self.n_tr], self.x[self.n_tr]]
                self.y_c = [self.traces.pick_times_max[trace][self.ipk],
                            self.traces.pick_times_min[trace][self.ipk]]
                self.lin, = self.axes[self.fig_plotted].plot(self.x_c,
                                                             self.y_c, "b")
                self.axes[self.fig_plotted].draw_artist(self.lin)
                self.canvas.blit(self.axes[self.fig_plotted].bbox)
# If right mouse has been clicked, plot all picks, including the moved
# one in red and leave routine
            elif event.button == 3:
                self.keys_held = set()
                self.drawNew(True)
                self.picksPlot(self.figs[self.fig_plotted],
                               self.axes[self.fig_plotted])
                self.end = True

        self.lin, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y,
                                                     animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect(
            'button_press_event', onPress)
        while self.end is False:
            QtCore.QCoreApplication.processEvents()
        self.setHelp(self.main_text)
        self.traces.storePicks()
        self.main.function = "main"

    def changeUnc(self, sign):
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
        figure = self.figs[self.fig_plotted]
        figure.canvas.restore_region(self.back)
        ds = sign*self.data.dt
        self.y_c[0] += ds
        self.y_c[1] -= ds
        self.traces.pick_times_min[self.itrace][self.ipk] -= ds
        self.traces.pick_times_max[self.itrace][self.ipk] += ds
        self.lin.set_data(self.x_c, self.y_c)
        self.axes[self.fig_plotted].draw_artist(self.lin)
        self.canvas.blit(self.axes[self.fig_plotted].bbox)

    def shiftPickTrace(self, direction):
        """
        Go to next trace to the left or to the right, depending on "direction".
        Function is called from self.main.event_filter.

        Parameters
        ----------
        direction : int
            If positive, go one trace to the right, if negative go one trace to
            the left.

        Returns
        -------
        None.

        """
        figure = self.figs[self.fig_plotted]
        self.d_time = 0
        self.keys_held = set()
        self.pickPlot(self.n_tr, self.itrace, self.ipk, "r")
        nt = np.where(self.tr_plot == self.n_tr)[0][0]
        nt1 = nt
        while True:
            nt1 += int(direction)
            if nt1 == -1:
                break
            if nt1 == len(self.tr_plot):
                break
            self.itrace = self.actual_traces[self.tr_plot[nt1]]
# If next trace doesn't have data or no pick exists, continue search
            if self.traces.npick[self.itrace] > 0 and self.stdev[self.tr_plot[nt1]]>0.:
                break
        if ((nt1 == len(self.tr_plot) and int(direction) > 0)
                or (nt1 < 0 and int(direction) < 0)):
            self.n_tr = self.tr_plot[nt]
            self.itrace = self.actual_traces[self.tr_plot[nt]]
            nt1 = nt
        self.n_tr = self.tr_plot[nt1]
        self.ipk = 0
        self.update()
# Save position of chosen pick in memory
        trace = self.itrace
        time_back = self.traces.pick_times[trace][self.ipk]
# Re-plot pick in white color (i.e. background color)
        self.pickPlot(self.n_tr, self.itrace, self.ipk, "w", v=False)
# Erase pick temporally to redraw canvas without this pick
        self.traces.pick_times[trace][self.ipk] = -1
        self.picksPlot(self.figs[self.fig_plotted],
                       self.axes[self.fig_plotted])
        self.back = figure.canvas.copy_from_bbox(figure.bbox)
        self.axes[self.fig_plotted].figure.canvas.draw()
        self.traces.pick_times[trace][self.ipk] = time_back
# If move picks is active, plot only horizontal line of pick position
        if self.main.function == "pick_move":
            self.x_c = []
            self.x_c.append(self.x_plot[nt1]-0.3)
            self.x_c.append(self.x_plot[nt1]+0.3)
            self.y_c = [self.traces.pick_times[trace][self.ipk],
                        self.traces.pick_times[trace][self.ipk]]
# If change uncertainty is active, plot vertical line of pick,
# i.e. uncertainty
        else:
            self.pickPlot(self.n_tr, self.itrace, self.ipk, "k", v=False)
            self.x_c = [self.x_plot[nt1],
                        self.x_plot[nt1]]
            self.y_c = [self.traces.pick_times_max[trace][self.ipk],
                        self.traces.pick_times_min[trace][self.ipk]]
        self.lin, = self.axes[self.fig_plotted].plot(self.x_c, self.y_c, "k")
        self.axes[self.fig_plotted].draw_artist(self.lin)
        self.canvas.blit(self.axes[self.fig_plotted].bbox)

    def erase_Picks(self):
        """
        Erases all picks of the actually shown gather

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        for i in range(self.actual_number_traces):
            ipt = self.actual_traces[i]
            self.traces.npick[ipt] = 0
            self.traces.pick_times[ipt] = []
            self.traces.pick_times_min[ipt] = []
            self.traces.pick_times_max[ipt] = []
        self.drawNew(True)
        self.setHelp(self.main_text)
        self.traces.storePicks()

    def pickManual(self):
        """
        Manual picking

        Press left mouse button to define the position of a pick. Then click
        again left mouse button to define the uncertainty range around the
        pick. Define as many picks a necessary (maximum 5 picks on every trace)

        Click central mouse button (wheel) to erase a pick

        Click right mouse button to finish manual picking
        """
        answer = self.main.test_function()
        if not answer:
            return
        self.setHelp(self.pick_manual_text)
        self.main.function = "manual_pick"
        self.press = 0
        self.finish = False
        if self.main.utilities.filtered:
            umin = 2*self.unc_base*self.data.dt
        else:
            umin = self.unc_base*self.data.dt
        self.end = False
        figure = self.figs[self.fig_plotted]
        xpk = []
        tpk = []
        npi = 0

        def onPress(event):
            nonlocal tpk, npi
            if not event.xdata:
                return
            if event.button == 1:
                modifiers = QtWidgets.QApplication.keyboardModifiers()
                if modifiers == QtCore.Qt.ShiftModifier:
                    print("Shift modifier")
                    for i, _ in enumerate(self.i_pk):
                        _ = self.searchTrace(self.x_pk[i])
                        trace = self.itrace
                        sht = self.p_s
                        stn = self.p_r
                        npi = self.traces.npick[trace]
                        self.traces.npick[trace] += 1
                        print(f"Pick nr {self.traces.npick[trace]} at shot "
                              + f"{sht+1}, receiver {stn+1}")
                        try:
                            self.traces.pick_times[trace][npi] = self.t_pk[i]
                            self.traces.pick_times_min[trace][npi] =\
                                self.t_pk[i]-umin
                            self.traces.pick_times_max[trace][npi] =\
                                self.t_pk[i]+umin
                        except (KeyError, IndexError, NameError):
                            self.traces.pick_times[trace].append(self.t_pk[i])
                            self.traces.pick_times_min[trace].append(
                                self.t_pk[i]-umin)
                            self.traces.pick_times_max[trace].append(
                                self.t_pk[i]+umin)
                    self.drawNew(True)
                    self.traces.storePicks()
                    self.end = True
                    return
# Left mouse button pressed, define a new pick
# First click at the position of the pick
                if self.press == 0:
                    tpk = event.ydata
                    ret = self.searchTrace(event.xdata)
                    if not ret:
                        print("Mouse outside axis")
                        return
                    self.press = 1
                    trace = self.itrace
                    sht = self.p_s
                    stn = self.p_r
                    npi = self.traces.npick[trace]
                    self.traces.npick[trace] += 1
                    print(f"Pick nr {self.traces.npick[trace]} at shot "
                          + f"{sht+1}, receiver {stn+1}")
                    try:
                        self.traces.pick_times[trace][npi] = tpk
                    except (KeyError, IndexError, NameError):
                        self.traces.pick_times[trace].append(tpk)
                    self.canvas = self.lin.figure.canvas
                    x_coor = [self.x[self.n_tr]-0.3, self.x[self.n_tr]+0.3]
                    y_coor = [tpk, tpk]
                    self.lin, = self.axes[self.fig_plotted].plot(
                        x_coor, y_coor, "r")
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
                    unc = max(np.abs(tpk-event.ydata), umin)
                    try:
                        self.traces.pick_times_min[trace][npi] = tpk-unc
                        self.traces.pick_times_max[trace][npi] = tpk+unc
                    except (KeyError, IndexError, NameError):
                        self.traces.pick_times_min[trace].append(tpk-unc)
                        self.traces.pick_times_max[trace].append(tpk+unc)
                    x_coor = [self.x[self.n_tr], self.x[self.n_tr]]
                    y_coor = [self.traces.pick_times_min[trace][npi],
                              self.traces.pick_times_max[trace][npi]]
                    self.lin, = self.axes[self.fig_plotted].plot(x_coor,
                                                                 y_coor, "r")
                    self.axes[self.fig_plotted].draw_artist(self.lin)
                    self.canvas.blit(self.axes[self.fig_plotted].bbox)
            elif event.button == 2:
                self.searchPick(event.xdata, event.ydata)
                trace = self.itrace
                sht = self.p_s
                stn = self.p_r
                npi = self.traces.npick[trace]
                if npi < 1:
                    self.main.function = "main"
                    return
                print(f"Erase pick nr {self.ipk+1} at shot {sht+1},"
                      + f"receiver {stn+1}")
                xx = self.x[self.n_tr]
                x_c = [xx-0.3, xx+0.3, xx, xx, xx]
                yy = self.traces.pick_times[trace][self.ipk]
                try:
                    yy1 = self.traces.pick_times_min[trace][self.ipk]
                    yy2 = self.traces.pick_times_max[trace][self.ipk]
                    y_c = [yy, yy, yy, yy1, yy2]
                except IndexError:
                    y_c = [yy, yy, yy, yy, yy]
# TO DO: I did not yet find a way to simply erase the drawing of a pick (except
#        redrawing). Therefore, for the moment, the pick is only hidden,
#        plotting it in white
                self.lin, = self.axes[self.fig_plotted].plot(x_c, y_c, "w")
                self.axes[self.fig_plotted].draw_artist(self.lin)
                self.canvas.blit(self.axes[self.fig_plotted].bbox)
                if npi > 1:
                    if self.ipk < npi-1:
                        for i in range(self.ipk+1, npi):
                            self.traces.pick_times[trace][i-1] =\
                                self.traces.pick_times[trace][i]
                            self.traces.pick_times_min[trace][i-1] =\
                                self.traces.pick_times_min[trace][i]
                            self.traces.pick_times_max[trace][i-1] =\
                                self.traces.pick_times_max[trace][i]
                    self.traces.npick[trace] -= 1
                elif npi == 0:
                    self.main.function = "main"
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
                    except (KeyError, IndexError, NameError):
                        self.traces.pick_times_min[self.itrace].\
                            append(tpk-umin)
                        self.traces.pick_times_max[self.itrace].\
                            append(tpk+umin)
                self.lin.set_animated(False)
                self.finish = True
                self.lin.figure.canvas.mpl_disconnect(self.cidpress)
                self.drawNew(True)
# Change help text to main module
                self.setHelp(self.main_text)
                self.traces.storePicks()
                self.end = True

# Plot lines in positive and negative direction corresponding to air wave
# arrivals
        xx = self.x[(self.x >= self.x_zoom_min) & (self.x <= self.x_zoom_max)]
        index = np.argsort(xx)
        xx = xx[index]
        xxx = self.x_off[(self.x >= self.x_zoom_min) & (self.x <= self.x_zoom_max)]
        xxx = xxx[index]
        tt = np.abs(xxx)/345.
        xxx = xx[tt >= self.t_zoom_min]
        ttt = tt[tt >= self.t_zoom_min]
        plt_air = True
        if self.dg_flag:
            dist = self.distances[self.fig_plotted_D]
            x_air = xx
            t_air = np.ones_like(x_air)*dist/345.
            if t_air[0] > self.t_zoom_max or t_air[0] < self.t_zoom_min:
                plt_air = False
        else:
            index = np.where(ttt <= self.t_zoom_max)[0]
            if len(index) > 1:
                ind = [0]
                for ii in range(1, len(index)):
                    if index[ii]-index[ii-1] > 1:
                        ind.append(ii)
                ind.append(len(index))
            else:
                plt_air = False
            # x_air = xxx[ttt <= self.t_zoom_max]
            # t_air = ttt[ttt <= self.t_zoom_max]
            # if len(t_air) < 2:
            #     plt_air = False
        if plt_air:
            # self.axes[self.fig_plotted].plot(x_air, t_air, c="b", ls="--")
            for i in range(len(ind)-1):
                i0 = index[ind[i]]
                i1 = index[ind[i+1]-1]+1
                self.axes[self.fig_plotted].plot(xxx[i0:i1], ttt[i0:i1], c="b", ls="--")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)

# Search picks done earlier in nearby gathers and plot a line connecting them
        x_pk_help, t_pk_help, i_pk_help = self.searchNearPicks()
        x1 = x_pk_help[(x_pk_help >= self.x_zoom_min) &
                       (x_pk_help <= self.x_zoom_max)]
        t1 = t_pk_help[(x_pk_help >= self.x_zoom_min) &
                       (x_pk_help <= self.x_zoom_max)]
        pk1 = i_pk_help[(x_pk_help >= self.x_zoom_min) &
                        (x_pk_help <= self.x_zoom_max)]
        self.x_pk = x1[x1 <= self.x_zoom_max]
        self.t_pk = t1[x1 <= self.x_zoom_max]
        self.i_pk = pk1[x1 <= self.x_zoom_max]
        index = np.where((self.t_pk >= self.t_zoom_min) &
                         (self.t_pk<=self.t_zoom_max))[0]
        self.x_pk = self.x_pk[index]
        self.t_pk = self.t_pk[index]
        self.i_pk = self.i_pk[index]
        if len(self.x_pk):
            self.axes[self.fig_plotted].plot(self.x_pk, self.t_pk, c="g",
                                             ls="--")
            renderer = figure.canvas.renderer
            self.axes[self.fig_plotted].draw(renderer)
            self.canvas.blit(self.axes[self.fig_plotted].bbox)

# Prepare self.lin to plot picks into.
        self.lin, = self.axes[self.fig_plotted].plot(xpk, tpk, animated=True)
        self.cidpress = self.lin.figure.canvas.mpl_connect(
            'button_press_event', onPress)
        while self.end is False:
            QtCore.QCoreApplication.processEvents()
            self.Tomography.setEnabled(True)
        self.PlotPicks.setEnabled(True)
        self.MovePicks.setEnabled(True)
        self.setHelp(self.main_text)
        self.main.function = "main"
