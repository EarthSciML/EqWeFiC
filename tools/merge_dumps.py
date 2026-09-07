#!/usr/bin/env python3
"""Merge several WRF dumps of one model step into a single column dump.

``tools/fill_tests.py`` gives a test's sidecar exactly two dumps (``inputs``
and ``reference``) to evaluate its expressions against.  A coupled document
that assembles a whole physics suite needs more than two: the driver-level
``subassembly_all``/``_pbl``/``_rad`` dumps plus the kernel-level dumps of the
individual schemes, all for the same WRF step.  This script slices each source
dump to one column (``tools/esm_dump.column``) and writes them out as ONE dump
whose variables are prefixed by a per-source tag, so a sidecar can name any of
them (``rad_pavel``, ``pbl_exch_h``, ``all_rthblten``, ...) in one namespace.

It moves numbers only; it never touches an ``.esm`` file.

Usage:
    python3 tools/merge_dumps.py OUT.json TAG=PATH [TAG=PATH ...] [--column N]

Example (the em_scm_xy suite at step 1410, kernel call 2820)::

    python3 tools/merge_dumps.py $D/subassembly/esm_dump_scm_suite_1410.json \\
        all=$D/subassembly/esm_dump_subassembly_all_1410.json \\
        pbl=$D/subassembly/esm_dump_subassembly_pbl_1410.json \\
        rad=$D/subassembly/esm_dump_subassembly_rad_1410.json \\
        sfc=$D/sfclayrev/scm_ref_2820.json \\
        sw=$D/sw/scm_ref_2820.json \\
        ysu=$D/ysu/scm_ref_2820.json \\
        rrtm=$D/rrtm/scm_ref2_2820.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import esm_dump  # noqa: E402


def _entry(a) -> dict:
    """One variable in the dump wire format (kind / shape / flat Fortran data)."""
    if isinstance(a, str):
        return {"kind": "str", "shape": [], "data": [a]}
    arr = np.asarray(a)
    if arr.dtype == bool:
        kind = "bool"
    elif np.issubdtype(arr.dtype, np.integer):
        kind = "int"
    else:
        kind = "real64"
    flat = arr.reshape(-1, order="F")
    conv = bool if kind == "bool" else (int if kind == "int" else float)
    return {"kind": kind, "shape": list(arr.shape), "data": [conv(x) for x in flat]}


def merge(out_path: str, sources: dict[str, str], column: int = 0) -> None:
    merged: dict[str, dict] = {}
    scheme, call = [], None
    for tag, path in sources.items():
        dump = esm_dump.load(path)
        scheme.append(f"{tag}={dump['scheme']}@{dump['call']}")
        if call is None:
            call = dump["call"]
        for name, val in esm_dump.column(dump, column).items():
            merged[f"{tag}_{name}"] = _entry(val)
    # keep the tile trivial so a later esm_dump.column() is a no-op
    merged["its"] = {"kind": "int", "shape": [], "data": [1]}
    merged["ite"] = {"kind": "int", "shape": [], "data": [1]}
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({"scheme": "merged(" + ", ".join(scheme) + ")",
                   "call": call, "vars": merged}, f)
    print(f"wrote {out_path}: {len(merged)} variables from {len(sources)} dumps")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("sources", nargs="+", metavar="TAG=PATH")
    ap.add_argument("--column", type=int, default=0, help="tile column to slice (0-based)")
    args = ap.parse_args(argv)
    srcs = {}
    for s in args.sources:
        tag, _, path = s.partition("=")
        if not path:
            ap.error(f"expected TAG=PATH, got {s!r}")
        srcs[tag] = path
    merge(args.out, srcs, args.column)
    return 0


if __name__ == "__main__":
    sys.exit(main())
