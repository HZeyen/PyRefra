# -*- coding: utf-8 -*-
"""
Mouse-driven rectangle and line tools of the data window.

Contains the following Class:
    InteractionMixin
"""

import numpy as np
from PyQt5 import QtCore


class InteractionMixin:
    """Mouse-driven rectangle and line tools of the data window."""

    def toggleHover(self):
        """Switch the display of trace numbers at the cursor on or off."""
        self.hover_on = not self.hover_on
        if not self.hover_on:
            self.hoverHide()
        print("Trace number display "
              + ("on" if self.hover_on else "off"))

    def hoverHide(self):
        """Remove the box with the trace numbers from the canvas."""
        ann = self.hover_ann
        if ann is None or not ann.get_visible():
            return
        ann.set_visible(False)
        if self.hover_bg is not None and ann.figure is not None:
            canvas = ann.figure.canvas
            canvas.restore_region(self.hover_bg)
            canvas.blit(ann.figure.bbox)

    def hoverSaveBackground(self, event):
        """Store the clean canvas after each full redraw (draw_event)."""
        self.hover_bg = event.canvas.copy_from_bbox(event.canvas.figure.bbox)

    def hoverText(self, x):
        """Text with the shot/receiver number of the trace nearest to x."""
        n = min(len(self.x), len(self.actual_traces))
        if n == 0:
            return None
        i = int(np.argmin(np.abs(np.asarray(self.x[:n], dtype=float) - x)))
        trace = self.actual_traces[i]
        shot = int(self.traces.shot[trace])+1
        rec = int(self.traces.receiver[trace])+1
        if self.sg_flag or self.fg_flag:
            return f"Receiver {rec}"
        if self.rg_flag:
            return f"Shot {shot}"
        return f"Shot {shot}, receiver {rec}"

    def hoverMotion(self, event):
        """Show a box with trace numbers near the cursor (motion_notify)."""
        if not self.hover_on or self.main.function != "main":
            return
        try:
            ax = self.axes[self.fig_plotted]
        except (IndexError, TypeError):
            return
        if (event.inaxes is not ax or event.xdata is None
                or self.hover_bg is None):
            self.hoverHide()
            return
        text = self.hoverText(event.xdata)
        if text is None:
            self.hoverHide()
            return
        ann = self.hover_ann
        if ann is None or ann.axes is not ax:
            ann = ax.annotate(
                "", xy=(0, 0), xytext=(15, 15), textcoords="offset pixels",
                bbox=dict(boxstyle="round", fc="lightyellow", ec="k"),
                animated=True, zorder=10)
            self.hover_ann = ann
        ann.xy = (event.xdata, event.ydata)
        ann.set_text(text)
        ann.set_visible(True)
        canvas = event.canvas
        canvas.restore_region(self.hover_bg)
        ax.draw_artist(ann)
        canvas.blit(canvas.figure.bbox)

    def followRect(self):
        """
        Pull rectangle across plot
        """
# set actual figure (simply shorter variable name)
        figure = self.figs[self.fig_plotted]

        def onPress(event):
            """
            Activity when mouse button is pressed
            """
# when left mouse botton is pressed initialize rectangle
# in the beginning store clicked point as first corner of rectangle
            if event.button == 1:
                if len(self.coor_x) == 0:
                    self.coor_x.append(event.xdata)
                    self.coor_y.append(event.ydata)
# store background figure
                    self.background = figure.canvas.copy_from_bbox(figure.bbox)
                    self.coor_x.append(event.xdata)
                    self.coor_y.append(event.ydata)
# define canvas and new axes
                self.canvas = self.line.figure.canvas
                self.axl = self.line.axes
                self.line.set_data(self.coor_x, self.coor_y)
# Draw point
                self.axl.draw_artist(self.line)
# set action on mouse motion
                self.cidmotion = self.line.figure.canvas.mpl_connect(
                    'motion_notify_event', onMotion)
# set action on mouse release
                self.cidrelease = self.line.figure.canvas.mpl_connect(
                    'button_release_event', onRelease)

        def onRelease(event):
            """
            Action taken when mouse button is released
            """
# free actions
            self.line.figure.canvas.mpl_disconnect(self.cidpress)
            self.line.figure.canvas.mpl_disconnect(self.cidmotion)
            self.line.figure.canvas.mpl_disconnect(self.cidrelease)
            self.line.set_animated(False)
# plot background image without rectangle
            if len(self.coor_x) > 0:
                figure.canvas.restore_region(self.background)
# remove stored background
            self.background = None
# set flag to finish module
            self.released = True

        def onMotion(event):
            """
            Action taken when mouse is moved
            """
# replace last coordinate point by actual mouse position
            self.coor_x[-1] = event.xdata
            self.coor_y[-1] = event.ydata
# set x and y coordiantes of rectangle plot
            line_coor_x = [self.coor_x[0], self.coor_x[1], self.coor_x[1],
                           self.coor_x[0], self.coor_x[0]]
            line_coor_y = [self.coor_y[0], self.coor_y[0], self.coor_y[1],
                           self.coor_y[1], self.coor_y[0]]
            self.line.set_data(line_coor_x, line_coor_y)
            self.canvas = self.line.figure.canvas
# define new axes values
            self.axl = self.line.axes
# plot background
            figure.canvas.restore_region(self.background)
# plot new rectangle
            self.axl.draw_artist(self.line)
            self.canvas.blit(self.axl.bbox)
# initialize flags and coordinate vector
        self.released = False
        self.coor_x = []
        self.coor_y = []
# define line value
        self.line, = self.axes[self.fig_plotted].plot(self.coor_x, self.coor_y,
                                                      animated=True)
# set action on mouse press
        self.cidpress = self.line.figure.canvas.mpl_connect(
            'button_press_event', onPress)
# cycle event detection as long as released flag is not set
        while self.released is False:
            QtCore.QCoreApplication.processEvents()

    def followLine(self, release_flag=False, nleft=1, nright=1, fig=None, ax=None):
        """
        Pull line across plot

        Parameters
        ----------
        release_flag : bool, optional; default: False
            If True, end of line when left button released; if False, end of
            line triggered by pressing right button
        nleft : int, optional; default: 1
            if 0 start line for negative direction at origin. If not, start
            line at the position of first click
        nright : int, optional; default: 1
            if 0 start line for positive direction at origin. If not, start
            line at the position of first click
        fig : matplotlib Figure, optional; default None
            Figure within which the line is drawn. If None, the main window
            is used. Else, a figure must have been created and its address is
            given here.
        ax : Matplotlib Axis, optional; Default: None
            Axis where to plot the line. Similar usabe as fig.

        Returns
        -------
        x_coor : 1D numpy float array
            Contains the x X-coordinates of all points forming the line
        """
# set figure (simply shorter variable name)
        if fig:
            figure = fig
            axis = ax
        else:
            figure = self.figs[self.fig_plotted]
            axis = self.axes[self.fig_plotted]

        def onPress(event):
            """
            action taken when mouse button is pressed
            """
            self.line_click = False
# left mouse button is pressed
            if event.button == 1:
                if event.xdata is None or event.ydata is None:
                    return
                if len(self.coor_x) == 0:
                    if (event.xdata < 0 and nleft == 0) or\
                       (event.xdata >= 0 and nright == 0):
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
# set starting point initially also as end point
                self.coor_x.append(event.xdata)
                self.coor_y.append(event.ydata)
                self.canvas = self.line.figure.canvas
                self.axl = self.line.axes
                self.line.set_data(self.coor_x, self.coor_y)
                self.axl.draw_artist(self.line)
# set action on mouse motion
                self.cidmotion = self.line.figure.canvas.mpl_connect(
                    'motion_notify_event', onMotion)
# set action on mouse release
                if release_flag:
                    self.cidrelease = self.line.figure.canvas.mpl_connect(
                        'button_release_event', onRelease)
# if right button is pressed, finish
            else:
                self.line.figure.canvas.mpl_disconnect(self.cidpress)
                try:
                    self.line.figure.canvas.mpl_disconnect(self.cidmotion)
                    self.line_click = True
                except (KeyError, IndexError, NameError):
                    self.line_click = False
                self.line.set_animated(False)
                if self.line_click:
                    figure.canvas.restore_region(self.background)
                self.background = None
                self.released = True

        def onRelease(event):
            """
            Action taken when mouse button is released
            """
# If line finishes when button is released do this here
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
            """
            Action taken when mouse is moved
            """
            if event.xdata is None or event.ydata is None:
                return False
            if self.main.function == "P_model":
                if (self.side < 0 < event.xdata) or\
                        (event.xdata < 0 < self.side):
                    event.xdata = 0.
# set second point of line as actual mouse position
            self.coor_x[-1] = event.xdata
            self.coor_y[-1] = event.ydata
# Check if starting point is outside window
            if len(self.start) > 0 and len(self.coor_x) == 2 and \
                    self.main.function not in ["spectrum", "air_wave"]:
                zxmin = self.zooms[self.i_zooms][0]
                zxmax = self.zooms[self.i_zooms][1]
                if self.start[0] < zxmin and zxmin >= 0:
                    self.coor_x[0] = zxmin
                    self.coor_y[0] = self.start[1]+(
                        self.coor_y[-1]-self.start[1])\
                        * (zxmin-self.start[0])/(self.coor_x[-1]-self.start[0])
                    if self.coor_y[0] < self.time_plt_min:
                        self.coor_y[0] = self.time_plt_min
                        self.coor_x[0] = self.start[0]+(
                            self.coor_x[-1]-self.start[0]) *\
                            (self.time_plt_min-self.start[1]) /\
                            (self.coor_y[-1]-self.start[1])
                elif self.start[0] > zxmax and zxmax <= 0:
                    self.coor_x[0] = zxmax
                    self.coor_y[0] = self.start[1]+(
                        self.coor_y[-1]-self.start[1])\
                        * (zxmax-self.start[0])/(self.coor_x[-1]-self.start[0])
                    if self.coor_y[0] < self.time_plt_min:
                        self.coor_y[0] = self.time_plt_min
                        self.coor_x[0] = self.start[0] +\
                            (self.coor_x[-1]-self.start[0]) *\
                            (self.time_plt_min-self.start[1]) /\
                            (self.coor_y[-1]-self.start[1])
# Draw new line
            self.line.set_data(self.coor_x, self.coor_y)
            self.canvas = self.line.figure.canvas
            self.axl = self.line.axes
            figure.canvas.restore_region(self.background)
            self.axl.draw_artist(self.line)
            self.canvas.blit(self.axl.bbox)
            return True

#        r_flag = release_flag # set flags and initialize coordinates
        self.released = False
        self.start = []
        self.coor_x = []
        self.coor_y = []
        self.line, = axis.plot(self.coor_x, self.coor_y, animated=True)
        self.cidpress = self.line.figure.canvas.mpl_connect(
            'button_press_event', onPress)
# as long as release flag is not set listen to events
        while self.released is False:
            QtCore.QCoreApplication.processEvents()
        return np.array(self.coor_x), np.array(self.coor_y)

    def follow_air(self, fig=None, ax=None):
        """
        Plot line across plot with slope corresponding the the inverse of the
        air wave velocity. Attached to the mouse position, the actually used
        velocity and the corresponding time break is plotted ina text box.

        Parameters
        ----------
        fig : matplotlib Figure, optional; default None
            Figure within which the line is drawn. If None, the main window
            is used. Else, a figure must have been created and its address is
            given here.
        ax : Matplotlib Axis, optional; Default: None
            Axis where to plot the line. Similar usabe as fig.

        Returns
        -------
        time_break : float
            Origine time of the line, corresponding to the time break
        vair : float
            Air wave velocity used for the calculation
        save_flag : bool
            If True, function has been finished pressing the left mouse button
            If False, function has been finished pressing the right button
        """
# set figure (simply shorter variable name)
        if fig:
            figure = fig
            axis = ax
        else:
            figure = self.figs[self.fig_plotted]
            axis = self.axes[self.fig_plotted]
        self.vair = 340.
        self.clicked = False
        save_flag = True

        def onPress(event):
            """
            action taken when mouse button is pressed
            """
# left mouse button is pressed
            nonlocal save_flag
            if event.button == 1 or event.button == 3:
                if event.xdata is None or event.ydata is None:
                    return
                self.time_break = event.ydata-abs(event.xdata)/self.vair
                self.clicked = True
                self.line.figure.canvas.mpl_disconnect(self.cidpress)
                self.line.figure.canvas.mpl_disconnect(self.cidmotion)
                self.line.figure.canvas.mpl_disconnect(self.cidrelease)
                self.line.set_animated(False)
                figure.canvas.restore_region(self.background)
                self.background = None
                if event.button == 1:
                    save_flag = True
                else:
                    save_flag = False
# if right button is pressed, finish
            else:
                self.vair += 3.
                if self.vair > 350.:
                    self.vair = 331.
                print(f"Air wave velocity: {self.vair} m/s")

        def onMotion(event):
            """
            Action taken when mouse is moved
            """
            if event.xdata is None or event.ydata is None:
                return False
            xline = []
            tline = []
            nlines= 0
            vair = self.vair
            neg_flag = True
            if self.x_zoom_min < 0.:
                nlines += 1
                tmax = (abs(self.x_zoom_min)-abs(event.xdata))/vair+event.ydata
                if tmax < self.t_zoom_min:
                    neg_flag = False
                elif tmax > self.t_zoom_max:
                    x1 = -((self.t_zoom_max-event.ydata)*vair+abs(event.xdata))
                    xline.append(x1)
                    tline.append(self.t_zoom_max)
                else:
                    xline.append(self.x_zoom_min)
                    tline.append(tmax)
                if neg_flag:
                    if self.x_zoom_max <= 0.:
                        x1 = self.x_zoom_max
                    else:
                        x1 = 0.
                    tmin = (abs(x1)-abs(event.xdata))/vair+event.ydata
                    if tmin > self.t_zoom_max:
                        return False
                    if tmin < self.t_zoom_min:
                        x1 = -((self.t_zoom_min-event.ydata)*vair+abs(event.xdata))
                        xline.append(x1)
                        tline.append(self.t_zoom_min)
                        if x1 < self.x_zoom_min:
                            xline = []
                            tline = []
                    else:
                        xline.append(x1)
                        tline.append(tmin)
            if self.x_zoom_max > 0:
                nlines += 1
                # pos_flag = True
                if self.x_zoom_min > 0.:
                    x1 = self.x_zoom_min
                else:
                    x1 = 0.
                tmin = (abs(x1)-abs(event.xdata))/vair+event.ydata
                if tmin > self.t_zoom_max:
                    # pos_flag = False
                    pass
                elif tmin < self.t_zoom_min:
                    x1 = (self.t_zoom_min-event.ydata)*vair+abs(event.xdata)
                    if x1 <= self.x_zoom_max:
                        xline.append(x1)
                        tline.append(self.t_zoom_min)
                else:
                    xline.append(x1)
                    tline.append(tmin)
                if x1 <= self.x_zoom_max:
                    tmax = (self.x_zoom_max-abs(event.xdata))/vair+event.ydata
                    if tmax < self.t_zoom_min:
                        return False
                    if tmax > self.t_zoom_max:
                        x1 = (self.t_zoom_max-event.ydata)*vair+abs(event.xdata)
                        xline.append(x1)
                        tline.append(self.t_zoom_max)
                    else:
                        xline.append(self.x_zoom_max)
                        tline.append(tmax)
# Draw new line
            if len(tline) > 0:
                t0 = (tline[-2]-(tline[-1]-tline[-2])/(xline[-1]-xline[-2])*
                      xline[-2])*1000.
            else:
                t0 = None
            self.line.set_data(xline[:2], tline[:2])
            if len(xline) > 2:
                self.line2.set_data(xline[2:], tline[2:])
            self.txt.set_position((event.xdata, event.ydata))
            if t0 is None:
                self.txt.set_text("")
            else:
                self.txt.set_text(f"v={self.vair}\nt0={t0:0.2f}")
            if event.xdata < 0.:
                self.txt.set_horizontalalignment("left")
            else:
                self.txt.set_horizontalalignment("right")
            self.canvas = self.line.figure.canvas
            self.axl = self.line.axes
            figure.canvas.restore_region(self.background)
            self.axl.draw_artist(self.line)
            if len(xline) > 2:
                self.axl.draw_artist(self.line2)
            self.axl.draw_artist(self.txt)
            self.canvas.blit(self.axl.bbox)
            return True

#        r_flag = release_flag # set flags and initialize coordinates
        xline = []
        tline = []
        self.background = figure.canvas.copy_from_bbox(figure.bbox)
        self.line, = axis.plot(xline, tline, "k", animated=True)
        self.line2, = axis.plot(xline, tline, "k", animated=True)
        self.txt = axis.text(self.x_zoom_min, self.t_zoom_min, "",
                             animated=True, verticalalignment='bottom',
                             bbox=dict(boxstyle="round,pad=0.1", fc="yellow",
                                       alpha=0.6, ec="black"))
        self.cidpress = self.line.figure.canvas.mpl_connect(
            'button_press_event', onPress)
        self.cidmotion = self.line.figure.canvas.mpl_connect(
            'motion_notify_event', onMotion)
# as long as release flag is not set listen to events
        while self.clicked is False:
            QtCore.QCoreApplication.processEvents()
        return self.time_break, self.vair, save_flag
