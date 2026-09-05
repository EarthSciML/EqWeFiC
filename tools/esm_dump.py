#!/usr/bin/env python3
"""Read the JSON dumps written by the WRF fork's ``module_esm_dump``.

A dump file (``esm_dump_<scheme>_<call>.json``) holds every input and output of
one physics-kernel call for one WRF tile.  ``load`` returns the variables as
NumPy arrays reshaped to their Fortran extents (column-major order preserved),
and ``column`` slices the tile down to a single column (first horizontal index
by default) so the arrays line up with a single-column esm component: 1-D
``[lev]``/``[lev_nodes]`` fields and scalars.

This module reads dumps only; it never writes equations.  The companion
``fill_tests.py`` (Phase 0.4) uses it to fill ``parameter_overrides`` and
``assertions`` of a hand-written test skeleton.

Usage:
    python3 tools/esm_dump.py <dump.json>            # summary of every variable
    python3 tools/esm_dump.py <dump.json> tx exch_hx # print those columns
    python3 tools/esm_dump.py <dump.json> --flat in.txt  # kernel-driver input
"""
from __future__ import annotations

import json
import sys
from typing import Any

import numpy as np


def load(path: str) -> dict[str, Any]:
    """Load a dump; returns {"scheme", "call", "vars": {name: array|scalar}}."""
    with open(path) as f:
        doc = json.load(f)
    out: dict[str, Any] = {"scheme": doc["scheme"], "call": doc["call"], "vars": {}}
    for name, v in doc["vars"].items():
        kind, shape, data = v["kind"], v["shape"], v["data"]
        if kind == "str":
            out["vars"][name] = data[0]
            continue
        dtype = {"real64": np.float64, "real32": np.float64, "int": np.int64, "bool": bool}[kind]
        arr = np.asarray(data, dtype=dtype)
        if shape:
            arr = arr.reshape(shape, order="F")
        else:
            arr = arr.reshape(())[()]
        out["vars"][name] = arr
    return out


def column(dump: dict[str, Any], i: int = 0) -> dict[str, Any]:
    """Slice every tile array to horizontal index ``i`` (0-based).

    Rank-2 ``(its:ite, k)`` arrays become 1-D over k; rank-1 ``(its:ite)``
    arrays become scalars; rank-3 ``(its:ite, k, n)`` become 2-D; scalars and
    strings pass through.  Arrays whose first extent is not the tile width are
    passed through unchanged.
    """
    v = dump["vars"]
    its, ite = int(v.get("its", 1)), int(v.get("ite", 1))
    width = ite - its + 1
    out: dict[str, Any] = {}
    for name, a in v.items():
        if isinstance(a, np.ndarray) and a.ndim >= 1 and a.shape[0] == width:
            out[name] = a[i, ...] if a.ndim > 1 else a[i]
        else:
            out[name] = a
    return out


def to_flat(dump: dict[str, Any], path: str) -> None:
    """Write the dump as the flat text the fork's kernel drivers read.

    One record per variable::

        <name> <kind> <rank> <extent_1> ... <extent_rank>
        <values, whitespace separated, Fortran (column-major) order>

    Kinds: real (written with 17 significant digits), int, bool (T/F), str.
    """
    with open(path, "w") as f:
        for name, a in dump["vars"].items():
            if isinstance(a, str):
                f.write(f"{name} str 0\n{a}\n")
                continue
            arr = np.asarray(a)
            dims = " ".join(str(n) for n in arr.shape)
            flat = arr.reshape(-1, order="F")
            if arr.dtype == bool:
                f.write(f"{name} bool {arr.ndim} {dims}\n" + " ".join("T" if x else "F" for x in flat) + "\n")
            elif np.issubdtype(arr.dtype, np.integer):
                f.write(f"{name} int {arr.ndim} {dims}\n" + " ".join(str(int(x)) for x in flat) + "\n")
            else:
                f.write(f"{name} real {arr.ndim} {dims}\n" + " ".join(f"{float(x):.17g}" for x in flat) + "\n")


def _summary(dump: dict[str, Any]) -> None:
    print(f"scheme={dump['scheme']} call={dump['call']}")
    for name, a in dump["vars"].items():
        if isinstance(a, np.ndarray):
            print(f"  {name:22s} {str(a.shape):14s} min={a.min():.6g} max={a.max():.6g}")
        else:
            print(f"  {name:22s} scalar          {a}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    d = load(sys.argv[1])
    if len(sys.argv) == 4 and sys.argv[2] == "--flat":
        to_flat(d, sys.argv[3])
    elif len(sys.argv) == 2:
        _summary(d)
    else:
        col = column(d)
        for name in sys.argv[2:]:
            print(name, np.array2string(np.asarray(col[name]), precision=10, max_line_width=120))
