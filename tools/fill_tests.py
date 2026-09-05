#!/usr/bin/env python3
"""Fill the numbers of a hand-written esm inline-test skeleton from WRF kernel dumps.

Phase 0.4 of PLAN.md.  The author writes the component and, next to it, a
sidecar ``<component>.esm.fill.json`` that says, per test, which dump supplies
the inputs, which dump(s) supply the reference outputs, and how every esm name
maps onto a dump variable (a NumPy expression over the column-sliced dump
variables, see ``tools/esm_dump.py``).  This script then

* writes the per-test **input-field library** ``tests/<id>_inputs.esm`` -- a
  template-library file whose ``match`` rules lower the component's
  ``input_<name>`` rewrite targets to the dump's column profiles (esm-spec 6.6
  ``parameter_overrides`` admit scalars only, so shaped inputs are injected
  through the test's ``expression_template_imports``, esm-spec 9.7.10); and
* splices the filled **test block** (``parameter_overrides``, the injection,
  and the ``assertions`` with their expected values) into the component's
  ``tests`` array, leaving every other byte of the file untouched.

It never writes an equation: the physics is authored by hand, this tool only
transcribes numbers.  Numbers are written with Python's shortest round-trip
``repr`` (full binary64 precision).

Sidecar schema (JSON)::

    {
      "dump_dir": "<directory holding the dumps>",
      "model": "<model name in the .esm>",
      "tests": [
        {
          "id": "...", "description": "...",
          "inputs": "<dump.json>",                       # kernel inputs (+ dt-free outputs)
          "reference": "<dump.json>" | ["<dt>", "<dt/2>", "<dt/4>"],
                                                          # outputs; a 3-list is Richardson-
                                                          # extrapolated to dt -> 0 (2nd order)
          "column": 0,                                    # horizontal index in the tile
          "time_span": {"start": 0, "end": 1},
          "tolerance": {"rel": 1e-6},                     # test-level default
          "fields_library": "tests/<id>_inputs.esm",      # relative to the .esm file
          "fields": {"input_theta": {"expr": "tx/pi2d", "axis": "lev"}, ...},
          "parameter_overrides": {"hfx": "hfx", ...},     # esm parameter -> expression
          "assertions": [
            {"variable": "hpbl", "expr": "hpbl", "tolerance": {"rel": 1e-6}},
            {"variable": "dthdt", "expr": "ttnp", "axis": "lev",
             "reduce": "Linf_error", "tolerance": {"abs": 1e-10}},
            {"variable": "dthdt", "expr": "ttnp", "coords": {"lev": [1, 2, 15]},
             "tolerance": {"rel": 1e-6}}
          ]
        }
      ]
    }

Expressions are evaluated with NumPy in a namespace holding every column-sliced
variable of the inputs dump and of the (extrapolated) reference, plus the
helpers ``cumsum0(a)`` (``[0, a1, a1+a2, ...]``), ``concat``, ``zeros`` and
``np``.  A ``reduce`` assertion (``Linf_error`` / ``L2_error``) compares the
whole column against the expression's array, written inline as an aggregate
over the axis of a ``const`` array (esm-spec 6.6.5 inline reference); a
``coords`` assertion with a list of 1-based levels expands to one scalar
assertion per level.

Usage:
    python3 tools/fill_tests.py <component.esm> [--sidecar <file>] [--only <test id>]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import esm_dump  # noqa: E402


# --------------------------------------------------------------------------- dumps
def _column(path: str, col: int) -> dict[str, Any]:
    return esm_dump.column(esm_dump.load(path), col)


def _reference(dump_dir: str, spec: Any, col: int) -> dict[str, Any]:
    """Reference outputs; a 3-list [dt, dt/2, dt/4] is Richardson-extrapolated."""
    if isinstance(spec, str):
        return _column(os.path.join(dump_dir, spec), col)
    cols = [_column(os.path.join(dump_dir, f), col) for f in spec]
    if len(cols) == 2:
        w = [-1.0, 2.0]
    elif len(cols) == 3:
        w = [1.0 / 3.0, -2.0, 8.0 / 3.0]
    else:
        raise ValueError("reference must be a file or a list of 2 or 3 files")
    out: dict[str, Any] = {}
    for name, a in cols[-1].items():
        if isinstance(a, np.ndarray) and np.issubdtype(a.dtype, np.floating):
            out[name] = sum(wi * np.asarray(c[name], dtype=float) for wi, c in zip(w, cols))
        else:
            out[name] = a
    return out


def _namespace(inputs: dict[str, Any], ref: dict[str, Any]) -> dict[str, Any]:
    ns: dict[str, Any] = {"np": np, "concat": lambda *a: np.concatenate([np.atleast_1d(np.asarray(x, dtype=float)) for x in a]),
                          "zeros": np.zeros, "cumsum0": lambda a: np.concatenate([[0.0], np.cumsum(np.asarray(a, dtype=float))])}
    ns.update(inputs)
    ns.update(ref)
    return ns


def _eval(expr: str, ns: dict[str, Any]) -> Any:
    return eval(expr, {"__builtins__": {}}, ns)  # noqa: S307 -- author-written sidecar


# --------------------------------------------------------------------------- esm AST
def _const_gather(values: np.ndarray, axis: str) -> dict[str, Any]:
    """aggregate(k from axis; index(const [...], k)) -- an inline column literal."""
    return {"op": "aggregate", "output_idx": ["k"], "args": [], "ranges": {"k": {"from": axis}},
            "expr": {"op": "index", "args": [{"op": "const", "args": [], "value": [float(x) for x in values]}, "k"]}}


def _fields_library(test: dict[str, Any], ns: dict[str, Any], esm_name: str) -> dict[str, Any]:
    templates = {}
    for op_name, spec in test["fields"].items():
        vals = np.asarray(_eval(spec["expr"], ns), dtype=float).reshape(-1)
        templates[op_name] = {
            "params": [],
            "match": {"op": op_name, "args": []},
            "body": _const_gather(vals, spec["axis"]),
        }
    return {
        "esm": "1.0.0",
        "metadata": {
            "name": f"{test['id']}_inputs",
            "description": (f"Input column fields of inline test `{test['id']}` of {esm_name}, extracted by "
                            f"tools/fill_tests.py from the WRF kernel dump {test['inputs']} (tile column "
                            f"{test.get('column', 0)}). Each match rule lowers one `input_<name>` rewrite "
                            f"target of the component to the dumped profile; the test injects this library "
                            f"through its expression_template_imports (esm-spec 9.7.10). Generated data -- do not edit."),
            "license": "MIT",
            "tags": ["test-inputs", "generated"],
        },
        "expression_templates": templates,
    }


def _assertions(test: dict[str, Any], ns: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for a in test["assertions"]:
        val = _eval(a["expr"], ns)
        base = {"variable": a["variable"], "time": float(a.get("time", 0.0))}
        if "reduce" in a:
            arr = np.asarray(val, dtype=float).reshape(-1)
            d = dict(base, reduce=a["reduce"], expected=0.0, reference=_const_gather(arr, a["axis"]))
            if "tolerance" in a:
                d["tolerance"] = a["tolerance"]
            out.append(d)
        elif "coords" in a:
            arr = np.asarray(val, dtype=float).reshape(-1)
            (axis, levels), = a["coords"].items()
            for lev in levels:
                d = dict(base, coords={axis: int(lev)}, expected=float(arr[int(lev) - 1]))
                if "tolerance" in a:
                    d["tolerance"] = a["tolerance"]
                out.append(d)
        else:
            d = dict(base, expected=float(np.asarray(val).reshape(-1)[0]))
            if "tolerance" in a:
                d["tolerance"] = a["tolerance"]
            out.append(d)
    return out


def _test_block(test: dict[str, Any], ns: dict[str, Any]) -> dict[str, Any]:
    block: dict[str, Any] = {"id": test["id"], "description": test["description"]}
    block["parameter_overrides"] = {k: float(np.asarray(_eval(e, ns)).reshape(-1)[0])
                                    for k, e in test.get("parameter_overrides", {}).items()}
    block["expression_template_imports"] = [{"ref": "./" + test["fields_library"].lstrip("./")}]
    block["time_span"] = test.get("time_span", {"start": 0.0, "end": 1.0})
    if "tolerance" in test:
        block["tolerance"] = test["tolerance"]
    block["assertions"] = _assertions(test, ns)
    return block


# --------------------------------------------------------------------------- splicing
def _find_tests_span(text: str, model: str) -> tuple[int, int]:
    """Character span of the JSON array value of ``"tests"`` inside model ``model``."""
    mkey = text.index(f'"{model}"')
    tkey = text.index('"tests"', mkey)
    start = text.index("[", tkey)
    depth, i, in_str, esc = 0, start, False, False
    while i < len(text):
        c = text[i]
        if in_str:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
        elif c == '"':
            in_str = True
        elif c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                return start, i + 1
        i += 1
    raise ValueError("unbalanced tests array")


def _indent_of(text: str, pos: int) -> str:
    line_start = text.rfind("\n", 0, pos) + 1
    return text[line_start:pos][: len(text[line_start:pos]) - len(text[line_start:pos].lstrip())]


def _compact_number_arrays(text: str) -> str:
    """Put every JSON array that holds only numbers on one line (the inline reference columns)."""
    import re
    return re.sub(r"\[\s*((?:-?\d[^\[\]{}\"]*?,\s*)*-?\d[^\[\]{}\"]*?)\s*\]",
                  lambda m: "[" + re.sub(r"\s+", " ", m.group(1)).strip() + "]", text)


def fill(esm_path: str, sidecar_path: str, only: str | None = None) -> None:
    side = json.load(open(sidecar_path))
    text = open(esm_path).read()
    doc = json.loads(text)
    model = side["model"]
    existing = {t["id"]: t for t in doc["models"][model].get("tests", [])}
    esm_dir = os.path.dirname(os.path.abspath(esm_path))
    for test in side["tests"]:
        if only and test["id"] != only:
            continue
        col = int(test.get("column", 0))
        inputs = _column(os.path.join(side["dump_dir"], test["inputs"]), col)
        ref = _reference(side["dump_dir"], test["reference"], col)
        ns = _namespace(inputs, ref)
        lib_path = os.path.join(esm_dir, test["fields_library"])
        os.makedirs(os.path.dirname(lib_path), exist_ok=True)
        with open(lib_path, "w") as f:
            f.write(_compact_number_arrays(json.dumps(_fields_library(test, ns, doc["metadata"]["name"]), indent=1)))
            f.write("\n")
        existing[test["id"]] = _test_block(test, ns)
        print(f"filled test {test['id']}: {len(existing[test['id']]['assertions'])} assertions, inputs -> {lib_path}")
    # splice the tests array back in sidecar order (tests the sidecar does not
    # know about keep their place first), preserving everything else byte-for-byte
    order = [t["id"] for t in doc["models"][model].get("tests", []) if t["id"] not in {s["id"] for s in side["tests"]}]
    order += [t["id"] for t in side["tests"] if t["id"] in existing]
    existing = {i: existing[i] for i in order}
    start, end = _find_tests_span(text, model)
    indent = _indent_of(text, text.rfind('"tests"', 0, start))
    body = _compact_number_arrays(json.dumps(list(existing.values()), indent=2))
    body = ("\n" + indent).join(body.split("\n"))
    open(esm_path, "w").write(text[:start] + body + text[end:])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("esm", help="component .esm file whose tests to fill")
    ap.add_argument("--sidecar", help="mapping file (default: <esm>.fill.json)")
    ap.add_argument("--only", help="fill only this test id")
    args = ap.parse_args(argv)
    fill(args.esm, args.sidecar or args.esm + ".fill.json", args.only)
    return 0


if __name__ == "__main__":
    sys.exit(main())
