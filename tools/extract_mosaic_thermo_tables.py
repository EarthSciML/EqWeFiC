#!/usr/bin/env python3
"""Transcribe MOSAIC's thermodynamic parameter TABLES into the hand-authored
components/aerosol/mosaic/mosaic_thermo_parameters.esm.

The tables are read from the instrumented model's own arrays (fields th_a_zsr,
th_b_zsr, th_aw_min, th_b_mtem, th_kappa, th_mw_aer, th_dens_aer of a
dumps/chem170_mosaic_therm dump), because load_mosaic_parameters assigns them from
default-real literals: each entry is the binary32 rounding widened to double, which
is what the double-precision thermodynamics actually computes with.  Only the
`default` of the seven shaped parameters is written; every other field of the file
is left as authored.

    python3 tools/extract_mosaic_thermo_tables.py <dump.json>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402

PARAMS = os.path.join(os.path.dirname(__file__), "..", "components", "aerosol", "mosaic",
                      "mosaic_thermo_parameters.esm")
FIELDS = {  # parameter -> dump field
    "aw_min": "th_aw_min", "b_zsr": "th_b_zsr", "a_zsr": "th_a_zsr", "b_mtem": "th_b_mtem",
    "kappa": "th_kappa", "mw_aer": "th_mw_aer", "dens_aer": "th_dens_aer",
}


def nested(a: np.ndarray):
    """A Fortran-ordered array as nested lists indexed [i][j][k] = a(i,j,k)."""
    return a.tolist()


def main(dump: str) -> None:
    v = load(dump)["vars"]
    doc = json.load(open(PARAMS))
    var = doc["models"]["MosaicThermoParameters"]["variables"]
    for p, f in FIELDS.items():
        a = np.asarray(v[f], dtype=float)
        # the model assigned these from default-real literals, so each entry should
        # be a binary32 value; say so where one is not (it is still transcribed)
        n = int(np.sum(a != a.astype(np.float32).astype(float)))
        if n:
            print(f"note: {p} has {n} entries that are not binary32 values")
        var[p]["default"] = nested(a)
    with open(PARAMS, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print("tables written:", ", ".join(f"{p}{list(np.shape(v[f]))}" for p, f in FIELDS.items()))


if __name__ == "__main__":
    main(sys.argv[1])
