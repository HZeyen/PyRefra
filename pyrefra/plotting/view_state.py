# -*- coding: utf-8 -*-
"""
State of the data window, grouped by topic.

Window keeps one ViewState (Window.view_state). Every field of the groups is
also reachable as an attribute of Window itself (window.v, window.zooms, ...)
through the properties installed by install_state_access(), so code outside
the plotting package does not have to know about the grouping.

The type of the gather on the screen is stored once, in
GatherState.gather_type. The legacy flags sg_flag, fg_flag, rg_flag and
dg_flag are derived from it.

Contains the following Classes:
    GatherType
    FigureState
    GatherState
    DisplayState
    ZoomState
    PickState
    MouseState
    ViewState

Contains the following function:
    install_state_access
"""

import enum
from dataclasses import dataclass, field, fields
from typing import Any, Optional


class GatherType(enum.Enum):
    """Kind of gather shown in the data window."""
    SHOT = "shot"
    FILE = "file"
    RECEIVER = "receiver"
    DISTANCE = "distance"


@dataclass
class FigureState:
    """Figures of all gathers and the one actually on the screen."""
    figs: list = field(default_factory=list)
    axes: Any = None
    actual_axis: Any = None
    plot_names: list = field(default_factory=list)
# fig_plotted is the number of the figure actually on the screen
#   (figs[fig_plotted])
    fig_plotted: int = 0
# only fig_plotted_D is used, the others are for bookkeeping of the figure
#   plotted for the different gathers
    fig_plotted_S: Optional[int] = None
    fig_plotted_F: Optional[int] = None
    fig_plotted_R: Optional[int] = None
    fig_plotted_D: Optional[int] = None


@dataclass
class GatherState:
    """Traces and sample arrays of the gather on the screen."""
    gather_type: Optional[GatherType] = GatherType.SHOT
    actual_traces: list = field(default_factory=list)
    actual_number_traces: int = 0
    actual_shot: int = 0
    actual_receivers: list = field(default_factory=list)
    actual_dist: float = 0.
    receivers: list = field(default_factory=list)
    distances: list = field(default_factory=list)
    sh_plot: list = field(default_factory=list)
# x is the plotted trace position (offset or station number, depending on
#   DisplayState.order_number), x_off the true signed offset, x_ori x without
#   the shift applied to overlapping traces
    x: Any = field(default_factory=list)
    x_off: Any = field(default_factory=list)
    x_ori: Any = field(default_factory=list)
    tr: Any = field(default_factory=list)
    x_plot: Any = field(default_factory=list)
    tr_plot: Any = field(default_factory=list)
    ind_plot: Any = field(default_factory=list)
    indices: Any = field(default_factory=list)
    stdev: Any = field(default_factory=list)
# v and v_norm contain the (normalized) data, one row per trace
    v: Any = field(default_factory=list)
    v_norm: Any = field(default_factory=list)
# v_set is True if v has already been filled for the actual gather
    v_set: bool = False


@dataclass
class DisplayState:
    """Gain, amplitude and component settings of the plot."""
    plotComponent: str = "All"
    order_number: bool = False
    x_text: str = "Offset [m]"
    amp_plt: float = 1
    gain: str = "tnorm"
    trace_norm: bool = True
    time_gain: Any = False
    dist_gain: Any = False
    phase_plot: bool = False
    agc_window: float = 0.1
    n_agc_win: int = 0
    line_colors: list = field(default_factory=list)


@dataclass
class ZoomState:
    """Zoom stack and limits of the plotted window."""
    n_zooms: int = 0
    i_zooms: int = -1
    zooms: list = field(default_factory=list)
    nt_mn: int = 0
    nt_mx: int = 0
    nt_0: int = 0
    time_plt_min: float = 0.
    time_plt_max: float = 0.
    x_zoom_min: float = 0.
    x_zoom_max: float = 0.
    t_zoom_min: float = 0.
    t_zoom_max: float = 0.
    dxmin: float = 0.
    dxmax: float = 0.
    rxmin: float = 0.
    rxmax: float = 0.
    sxmin: float = 0.
    sxmax: float = 0.
# shift is the distance by which two traces at the same position are separated
    shift: float = 0.1


@dataclass
class PickState:
    """Selected trace and pick, and helper picks of nearby gathers."""
    n_tr: int = 0
    itrace: int = 0
    p_s: int = 0
    p_r: int = 0
    ipk: int = 0
    x_pk: Any = field(default_factory=list)
    t_pk: Any = field(default_factory=list)
    i_pk: Any = field(default_factory=list)
    i_tr: Any = field(default_factory=list)
    x_pk_help: Any = field(default_factory=list)
    t_pk_help: Any = field(default_factory=list)
    i_tr_help: Any = field(default_factory=list)
    x_c: list = field(default_factory=list)
    y_c: list = field(default_factory=list)
    d_time: float = 0.
    press: int = 0
    end: bool = False
# finish indicates end of treatment, picked the existence of measured picks
    finish: bool = False
    picked: bool = False


@dataclass
class MouseState:
    """Event connections and coordinates of the mouse tools."""
    x_coor: list = field(default_factory=list)
    y_coor: list = field(default_factory=list)
    coor_x: list = field(default_factory=list)
    coor_y: list = field(default_factory=list)
    cidpress: Any = None
    cidmotion: Any = None
    cidrelease: Any = None
    back: Any = None
    background: Any = None
    line_click: Any = None
    released: bool = False
    start: list = field(default_factory=lambda: [0, 0])
    side: int = 1
    keys_held: set = field(default_factory=set)


@dataclass
class ViewState:
    figure: FigureState = field(default_factory=FigureState)
    gather: GatherState = field(default_factory=GatherState)
    display: DisplayState = field(default_factory=DisplayState)
    zoom: ZoomState = field(default_factory=ZoomState)
    pick: PickState = field(default_factory=PickState)
    mouse: MouseState = field(default_factory=MouseState)


def _forward(group, name):
    def getter(self):
        return getattr(getattr(self.view_state, group), name)

    def setter(self, value):
        setattr(getattr(self.view_state, group), name, value)

    return property(getter, setter, doc=f"ViewState.{group}.{name}")


def _gather_flag(kind):
    def getter(self):
        return self.view_state.gather.gather_type is kind

    def setter(self, value):
        gather = self.view_state.gather
        if value:
            gather.gather_type = kind
        elif gather.gather_type is kind:
            gather.gather_type = None

    return property(getter, setter, doc=f"True for {kind} gathers")


def install_state_access(window_cls):
    """
    Make all fields of ViewState available as attributes of window_cls and
    derive the flags sg_flag, fg_flag, rg_flag and dg_flag from gather_type.
    """
    owner = {}
    for group, state_cls in (
            ("figure", FigureState), ("gather", GatherState),
            ("display", DisplayState), ("zoom", ZoomState),
            ("pick", PickState), ("mouse", MouseState)):
        for fld in fields(state_cls):
            if fld.name in owner:
                raise ValueError(f"{fld.name} is in {owner[fld.name]} and "
                                 f"{group}")
            owner[fld.name] = group
            setattr(window_cls, fld.name, _forward(group, fld.name))
    for flag, kind in (("sg_flag", GatherType.SHOT),
                       ("fg_flag", GatherType.FILE),
                       ("rg_flag", GatherType.RECEIVER),
                       ("dg_flag", GatherType.DISTANCE)):
        setattr(window_cls, flag, _gather_flag(kind))
