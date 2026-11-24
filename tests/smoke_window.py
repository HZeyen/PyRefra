# -*- coding: utf-8 -*-
"""
Offscreen smoke test for pyrefra.plotting.refraPlot.Window.

Loads the data set in Example/ (copied to a temporary folder) and exercises the
non-interactive plotting, zoom, gain and mute paths. Run with

    python tests/smoke_window.py

Exit code is 0 if all checks pass.
"""

import os
import shutil
import sys
import tempfile
import traceback
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from PyQt5 import QtWidgets  # noqa: E402
from pyrefra import Pyrefra  # noqa: E402


def _forbid_modal_boxes():
    """Any unexpected message box would block the test, so make it fail."""
    def boom(*args, **kwargs):
        raise RuntimeError(f"Unexpected message box: {args[1:3]}")
    for name in ("warning", "critical", "information", "question"):
        setattr(QtWidgets.QMessageBox, name, staticmethod(boom))


def _load_main(work):
    files = sorted(str(p) for p in Path(work).glob("*.seg2"))
    QtWidgets.QFileDialog.getOpenFileNames = staticmethod(
        lambda *args, **kwargs: (files, ""))
    return Pyrefra.Main(str(ROOT), str(work))


def _visit_gathers(window):
    """Select first, middle and last entry of the right-hand gather list."""
    n = window.mplfigs.count()
    for i in sorted({0, n // 2, n - 1}):
        window.changeFig(window.mplfigs.item(i))


def _add_picks(main):
    traces = main.traces
    for t in range(min(10, traces.number_of_traces)):
        traces.npick[t] = 1
        traces.pick_times[t] = [0.02]
        traces.pick_times_min[t] = [0.019]
        traces.pick_times_max[t] = [0.021]


def _with_dialog(main, answers, call):
    """Run call() with Main.dialog replaced by a fixed answer."""
    original = main.dialog
    main.dialog = lambda *args, **kwargs: (answers, True)
    try:
        call()
    finally:
        main.dialog = original


def _check_regression(window):
    import numpy as np
    x = np.arange(-40., 41.)
    ax = np.abs(x)
    t = np.where(ax <= 20., ax/400., 0.03667+ax/1500.)
    off, slope, inter, *_ = window.bestLines(x, t, n_lines=2)
    assert np.allclose(slope, [1/400., 1/1500.], rtol=0.15), slope
    assert np.allclose(window.set_ticks(0., 20., 4, 5.), [0., 5., 10., 15., 20.])
    window.L1_regres(ax, t, [1/400.])


def _check_state(window):
    from dataclasses import fields
    from pyrefra.plotting import view_state as vs
    from pyrefra.plotting.refraPlot import Ui_MainWindow, QMainWindow

    ui = Ui_MainWindow()
    ui.setupUi(QMainWindow())
    names = {f.name for cls in (vs.FigureState, vs.GatherState,
                                vs.DisplayState, vs.ZoomState, vs.PickState,
                                vs.MouseState) for f in fields(cls)}
    assert not names & set(vars(ui)), names & set(vars(ui))

    kinds = {"sg_flag": vs.GatherType.SHOT, "fg_flag": vs.GatherType.FILE,
             "rg_flag": vs.GatherType.RECEIVER,
             "dg_flag": vs.GatherType.DISTANCE}
    for plot, kind in ((window.plotSG, vs.GatherType.SHOT),
                       (window.plotFG, vs.GatherType.FILE),
                       (window.plotRG, vs.GatherType.RECEIVER),
                       (window.plotDG, vs.GatherType.DISTANCE)):
        plot()
        assert window.gather_type is kind
        for flag, flag_kind in kinds.items():
            assert getattr(window, flag) == (flag_kind is kind), flag
    assert window.zooms is window.view_state.zoom.zooms
    window.v = window.v
    window.sg_flag = True
    assert window.gather_type is vs.GatherType.SHOT
    window.sg_flag = False
    assert window.gather_type is None
    window.plotSG()


def _check_all_shown(window):
    """Every trace and pick of a gather must stay visible when the gather or
    the x-axis ordering changes."""
    def check(tag):
        ax = window.axes[window.fig_plotted]
        assert len(window.tr) == len(window.actual_traces), \
            (tag, len(window.tr), len(window.actual_traces))
        shown = sum(1 for line in ax.lines if len(line.get_xdata()) == 2)
        total = sum(window.traces.npick[t] for t in window.actual_traces)
        assert shown == 2*total, (tag, shown, 2*total)

    for plot in (window.plotSG, window.plotFG, window.plotRG,
                 window.plotDG):
        for _ in range(2):
            plot()
            n = window.mplfigs.count()
            for i in sorted({0, n // 2, n - 1}):
                window.changeFig(window.mplfigs.item(i))
                check((plot.__name__, i, window.order_number, "changeFig"))
                window.toggleNumber()
                check((plot.__name__, i, window.order_number, "toggle"))
                window.changeFig(window.mplfigs.item(i))
                check((plot.__name__, i, window.order_number, "changeFig"))
    if window.order_number:
        window.toggleNumber()
    window.plotSG()


def run_checks(main):
    window = main.window
    checks = []

    def gather(plot, order_number=False):
        def run():
            if window.order_number != order_number:
                window.toggleNumber()
            plot()
            _visit_gathers(window)
            window.zoomIni()
        return run

    for order_number in (False, True):
        tag = "number" if order_number else "offset"
        checks += [
            (f"shot gathers ({tag})", gather(window.plotSG, order_number)),
            (f"file gathers ({tag})", gather(window.plotFG, order_number)),
            (f"receiver gathers ({tag})", gather(window.plotRG, order_number)),
        ]
    checks += [
        ("distance gathers", gather(window.plotDG, False)),
        ("picks", lambda: (_add_picks(main), window.plotSG(),
                           _visit_gathers(window))),
        ("pick section by position", lambda: _with_dialog(
            main, [None, 0, None, 0, 0], window.plotPickSection)),
        ("pick section by offset and direction", lambda: _with_dialog(
            main, [None, 1, None, 1, -1], window.plotPickSection)),
        ("gain and amplitude", lambda: (window.plotSG(), window.tNorm(),
                                        window.changeAmp(1),
                                        window.changeAmp(-1))),
        ("zoom out/in", lambda: (window.zoomOut(), window.zoomIn())),
        ("restore original data", lambda: (window.original(),
                                           window.originalScreen())),
        ("mute air wave", lambda: _with_dialog(
            main, ["340", "10", "3"], lambda: window.muteAir(None))),
        ("regression and tick helpers", lambda: _check_regression(window)),
        ("view state and gather flags", lambda: _check_state(window)),
        ("all traces and picks stay visible",
         lambda: _check_all_shown(window)),
    ]

    failures = []
    for name, check in checks:
        try:
            check()
            print(f"ok    {name}")
        except Exception:
            failures.append(name)
            print(f"FAIL  {name}")
            traceback.print_exc()
    return failures


def main():
    _forbid_modal_boxes()
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(
        sys.argv)
    cwd = os.getcwd()
    work = tempfile.mkdtemp(prefix="pyrefra_smoke_")
    try:
        shutil.copytree(ROOT / "Example", work, dirs_exist_ok=True)
        failures = run_checks(_load_main(work))
    finally:
        os.chdir(cwd)
        shutil.rmtree(work, ignore_errors=True)
    del app
    print(f"\n{len(failures)} failing check(s)")
    return 1 if failures else 0


def test_window_smoke():
    assert main() == 0


if __name__ == "__main__":
    sys.exit(main())
