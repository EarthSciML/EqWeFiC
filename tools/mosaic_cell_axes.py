"""Put single-bin MOSAIC test data on the cell axes [lev, bin].

Since 2026-11-12 the MOSAIC thermodynamics components carry the grid-cell axes
(lev, bin) ahead of their own axes (NLEV = NBIN = 1 by default), and the fill
scripts still build each test for ONE bin of one cell.  `cell_axes` reshapes such a
test in place against the shapes the document declares: an array override short
of its variable's rank by the leading cell axes is nested once per missing axis
(a scalar override already broadcasts), and an assertion on a variable with cell
axes gets coordinate 1 on each of them unless it names one.  It is idempotent, and
only numbers move: nothing about a component's equations is generated here.

    from mosaic_cell_axes import cell_axes
    cell_axes(doc, model_name, tests, doc_path)
"""
import json
import os

CELL = ("lev", "bin")


def _shapes_of(model):
    return {n: list(v.get("shape") or []) for n, v in model.get("variables", {}).items()}


def doc_shapes(doc, model_name, doc_path):
    """{name as a test of `model_name` spells it: full shape}."""
    out = {}
    models = doc.get("models", {})
    for mn, m in models.items():
        if "ref" in m:
            sub = json.load(open(os.path.join(os.path.dirname(os.path.abspath(doc_path)), m["ref"])))
            sm = next(iter(sub["models"].values()))
            for n, s in _shapes_of(sm).items():
                out[f"{mn}.{n}"] = s
        else:
            for n, s in _shapes_of(m).items():
                out[f"{mn}.{n}"] = s
                if mn == model_name:
                    out[n] = s
    return out


def _rank(v):
    r = 0
    while isinstance(v, list):
        r += 1
        v = v[0] if v else None
    return r


def cell_axes(doc, model_name, tests, doc_path):
    shapes = doc_shapes(doc, model_name, doc_path)
    for t in tests:
        for block in ("parameter_overrides", "initial_conditions"):
            for k, v in t.get(block, {}).items():
                s = shapes.get(k)
                if s is None or not isinstance(v, list):
                    continue
                for _ in range(len(s) - _rank(v)):
                    v = [v]
                t[block][k] = v
        for a in t.get("assertions", []):
            s = shapes.get(a["variable"])
            if not s:
                continue
            lead = [ax for ax in s if ax in CELL]
            if lead:
                a["coords"] = {**{ax: 1 for ax in lead}, **a.get("coords", {})}
    return tests
