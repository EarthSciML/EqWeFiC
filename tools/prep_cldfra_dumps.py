#!/usr/bin/env python3
"""Build the single-column input dumps and real64 reference replays for CalCldfra.

cal_cldfra1 is a pointwise diagnostic, so its reference cases are just columns
pulled out of dumps that were written for other schemes: the RRTM longwave hooks
already dump every argument cal_cldfra1 reads (t3d, p3d, qv3d, qc3d, qi3d, qs3d)
alongside the CLDFRA the model itself computed.  For each (dump, column) this
script

* writes ``<name>.json`` -- the same dump narrowed to the six input fields of
  ONE column, so the component's test can use ``column: 0`` uniformly; and
* writes ``<name>.flat`` and runs the fork's real64 ``cldfra_driver`` on it,
  giving ``<name>_driver.json`` with ``cldfra3d_out`` and ``cldfra1_flag_out``.

The real64 replay is the reference, not the in-model real32 CLDFRA: the Randall
branch divides by (RHGRID*QVS_WEIGHT - QV)**GAMMA, a cancellation that costs the
single-precision model up to 5e-5 absolute near saturation (measured; evaluating
this component's own expressions in binary32 reproduces the model's CLDFRA to
6e-8, so the gap is WRF's arithmetic, not the transcription).

It never writes an equation, and it writes no .esm file: tools/fill_tests.py
turns these dumps into the numbers of the hand-written test skeleton.

Usage:  python3 tools/prep_cldfra_dumps.py [--driver <path>] [--out <dir>]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import esm_dump  # noqa: E402

DUMPS = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/dumps"

# (name, source dump, 0-based column, mp_physics)  -- every case is WSM6 with icloud = 1
CASES = [
    ("scm_clear_1", f"{DUMPS}/scm_ref_raw/esm_dump_rrtm_1.json", 0, 6),
    ("scm_cloud_120", f"{DUMPS}/scm_cloud_raw/esm_dump_rrtm_120.json", 0, 6),
    ("scm_cloud_360", f"{DUMPS}/scm_cloud_raw/esm_dump_rrtm_360.json", 0, 6),
    ("qss_warm_61", f"{DUMPS}/rrtm/qss_i20_61.json", 19, 6),
    ("qss_mixed_268", f"{DUMPS}/rrtm/qss_i20_268.json", 19, 6),
    ("qss_ice_396", f"{DUMPS}/rrtm/qss_i40_396.json", 39, 6),
]

FIELDS = ("t3d", "p3d", "qv3d", "qc3d", "qi3d", "qs3d")
FLAGS = ("f_qv", "f_qc", "f_qi", "f_qs")


def one_column(path: str, col: int, mp_physics: int) -> dict:
    """The dump narrowed to `col`, as a one-column tile (its = ite = 1)."""
    v = esm_dump.load(path)["vars"]
    out = {"its": np.int64(1), "ite": np.int64(1),
           "kts": np.int64(v["kts"]), "kte": np.int64(v["kte"]),
           "mp_physics": np.int64(mp_physics)}
    for f in FLAGS:
        out[f] = np.asarray(v[f], dtype=bool).reshape(())[()]
    for f in FIELDS:
        out[f] = np.asarray(v[f][col], dtype=np.float64).reshape(1, -1)
    return out


def write_dump(vars_: dict, scheme: str, call: int, path: str) -> None:
    kinds = {np.dtype(bool): "bool", np.dtype(np.int64): "int", np.dtype(np.float64): "real64"}
    doc = {"scheme": scheme, "call": call, "vars": {}}
    for name, a in vars_.items():
        arr = np.asarray(a)
        doc["vars"][name] = {"kind": kinds[arr.dtype], "shape": list(arr.shape),
                             "data": [bool(x) if arr.dtype == bool else
                                      int(x) if arr.dtype == np.int64 else float(x)
                                      for x in arr.reshape(-1, order="F")]}
    with open(path, "w") as f:
        json.dump(doc, f)


def main(argv=None) -> int:
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--driver", default=os.path.join(here, "..", "WRF", "kernels", "build", "cldfra_driver"))
    ap.add_argument("--out", default=f"{DUMPS}/cldfra")
    ap.add_argument("--flat-dir", default=None, help="where to keep the intermediate .flat files")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    flat_dir = args.flat_dir or args.out
    os.makedirs(flat_dir, exist_ok=True)
    for name, src, col, mp in CASES:
        v = one_column(src, col, mp)
        write_dump(v, "cldfra", 0, os.path.join(args.out, f"{name}.json"))
        flat = os.path.join(flat_dir, f"{name}.flat")
        esm_dump.to_flat({"scheme": "cldfra", "call": 0, "vars": v}, flat)
        subprocess.run([args.driver, flat, os.path.join(args.out, f"{name}_driver.json")], check=True)
        print(f"{name}: column {col} of {os.path.basename(src)} -> {name}.json + {name}_driver.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
