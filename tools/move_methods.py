# -*- coding: utf-8 -*-
"""
One-off refactoring helper: cut named methods out of class Window in
pyrefra/plotting/refraPlot.py and paste them, verbatim, into mixin classes in
other modules of the same package.

Method text is copied unchanged, so behaviour is preserved; only the imports
needed by the moved code are generated. Run from the repository root:

    python tools/move_methods.py --dry-run
    python tools/move_methods.py

Imports of the new mixins into refraPlot.py and the base-class list of Window
are not touched and have to be edited by hand afterwards.
"""

import argparse
import ast
from pathlib import Path

PKG = Path(__file__).resolve().parents[1] / "pyrefra" / "plotting"
SOURCE = PKG / "refraPlot.py"
SOURCE_CLASS = "Window"

# target module -> (mixin class name, description, methods to move)
PLAN = {
    "canvas_manager": (
        "CanvasMixin",
        "Canvas, figure list, help line and plot saving of the data window.",
        ["addMPL", "rmMPL", "addFig", "drawNew", "setHelp", "savePlot"]),
    "zoom": (
        "ZoomMixin",
        "Zoom handling of the data window.",
        ["addZoom", "changeZoom", "zooming", "zoomOut", "zoomIn", "zoomIni",
         "nextBlock"]),
    "render": (
        "RenderMixin",
        "Seismogram drawing of the data window.",
        ["seismogram", "agcCalc", "set_ticks"]),
    "gathers": (
        "GatherMixin",
        "Construction and plotting of shot, file, receiver and distance "
        "gathers.",
        ["changeFig", "toggleNumber", "chooseComponent", "plotRG",
         "plotReceiver", "plotDG", "plotDistance", "plotSG", "plotShot",
         "plotFG", "plotFile"]),
    "interaction": (
        "InteractionMixin",
        "Mouse-driven rectangle and line tools of the data window.",
        ["followRect", "followLine", "follow_air"]),
    "picking_display": (
        "PickDisplayMixin",
        "Display and lookup of picks in the data window.",
        ["picksPlot", "pickPlot", "plotPickSection", "plotCalcPicks",
         "searchTrace", "searchPick", "searchNearPicks"]),
    "picking_manual": (
        "ManualPickMixin",
        "Interactive editing of picks.",
        ["pickMove", "movePick", "uncertainty", "changeUnc",
         "shiftPickTrace", "erase_Picks", "pickManual"]),
    "picking_auto": (
        "AutoPickMixin",
        "Semi-automatic and automatic picking and line fitting.",
        ["findNearest2ndDerivative", "corrPick", "search_corr_pick",
         "Sta_Lta", "ampPick", "ampPicks", "bestLines", "L1_regres"]),
    "trace_processing": (
        "ProcessingMixin",
        "Gain, muting, polarity and display of original data.",
        ["original", "originalScreen", "changeAmp", "tNorm", "tGain",
         "dGain", "AGC", "phasePlot", "phaseCalc", "traceMute", "muteAir",
         "muteBefore", "muteAfter", "muteTime", "traceSign", "changeSign",
         "animateLine"]),
}


def collect_imports(tree):
    """Return list of (bound_name, module or None, name, asname)."""
    items = []
    for node in tree.body:
        if isinstance(node, ast.Import):
            for a in node.names:
                bound = a.asname or a.name.split(".")[0]
                items.append((bound, None, a.name, a.asname))
        elif isinstance(node, ast.ImportFrom):
            mod = "." * node.level + (node.module or "")
            for a in node.names:
                items.append((a.asname or a.name, mod, a.name, a.asname))
    return items


def render_imports(items):
    lines = []
    from_groups = {}
    for _, mod, name, asname in items:
        alias = f"{name} as {asname}" if asname else name
        if mod is None:
            lines.append(f"import {alias}")
        else:
            from_groups.setdefault(mod, []).append(alias)
    for mod, aliases in from_groups.items():
        lines.append(f"from {mod} import {', '.join(aliases)}")
    return lines


def used_names(nodes):
    names = set()
    for node in nodes:
        for sub in ast.walk(node):
            if isinstance(sub, ast.Name):
                names.add(sub.id)
    return names


def method_span(lines, node):
    start = min([node.lineno] + [d.lineno for d in node.decorator_list])
    while start > 1 and lines[start-2].strip().startswith("#"):
        start -= 1
    return start, node.end_lineno


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    text = SOURCE.read_text(encoding="utf-8")
    lines = text.split("\n")
    tree = ast.parse(text)
    cls = next(n for n in tree.body
               if isinstance(n, ast.ClassDef) and n.name == SOURCE_CLASS)
    methods = {n.name: n for n in cls.body if isinstance(n, ast.FunctionDef)}
    src_imports = collect_imports(tree)

    wanted = [m for _, _, names in PLAN.values() for m in names]
    missing = [m for m in wanted if m not in methods]
    dupes = {m for m in wanted if wanted.count(m) > 1}
    if missing or dupes:
        raise SystemExit(f"missing: {missing}, duplicated: {sorted(dupes)}")

    remove = []
    outputs = {}
    for module, (cls_name, doc, names) in PLAN.items():
        nodes = sorted((methods[n] for n in names), key=lambda n: n.lineno)
        chunks = []
        for node in nodes:
            start, end = method_span(lines, node)
            chunks.append("\n".join(lines[start-1:end]))
            remove.append((start, end))
        target = PKG / f"{module}.py"
        existing = target.read_text(encoding="utf-8") if target.exists() \
            else ""
        have = {i[0] for i in collect_imports(ast.parse(existing))} \
            if existing else set()
        used = used_names(nodes)
        need = [i for i in src_imports
                if i[0] in used and i[0] not in have
                and i[1] != f".{module}"]
        body = (f"class {cls_name}:\n    \"\"\"{doc}\"\"\"\n\n"
                + "\n\n".join(chunks) + "\n")
        outputs[module] = (existing, need, body)
        print(f"{module}: {cls_name} <- {len(names)} methods, "
              f"{sum(len(c.split(chr(10))) for c in chunks)} lines, "
              f"imports: {[i[0] for i in need]}")

    if args.dry_run:
        return

    for module, (existing, need, body) in outputs.items():
        target = PKG / f"{module}.py"
        imp = render_imports(need)
        if existing:
            etree = ast.parse(existing)
            elines = existing.split("\n")
            last = max((n.end_lineno for n in etree.body
                        if isinstance(n, (ast.Import, ast.ImportFrom))),
                       default=None)
            if last is None:
                raise SystemExit(f"{module}: no import block to extend")
            elines[last:last] = imp
            new_text = "\n".join(elines).rstrip("\n") + "\n\n\n" + body
        else:
            cls_name, doc, names = PLAN[module]
            head = (f"# -*- coding: utf-8 -*-\n\"\"\"\n{doc}\n\nContains the "
                    f"following Class:\n    {cls_name}\n\"\"\"\n\n")
            new_text = head + "\n".join(imp) + "\n\n\n" + body
        target.write_text(new_text, encoding="utf-8")

    drop = set()
    for start, end in remove:
        drop.update(range(start, end + 1))
        if end < len(lines) and not lines[end].strip():
            drop.add(end + 1)
    kept = [ln for i, ln in enumerate(lines, 1) if i not in drop]
    SOURCE.write_text("\n".join(kept), encoding="utf-8")
    ast.parse(SOURCE.read_text(encoding="utf-8"))
    print("done")


if __name__ == "__main__":
    main()
