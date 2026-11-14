#!/usr/bin/env python3
"""Fill the inline tests of components/aerosol/mosaic/coagulation_assembly.esm from
the chem170 MOSAIC coagulation capture and its real64 coagsolv replay.

It transcribes numbers only; the component's equations are hand-written.  Each
case is one (model step, column level, composition quantity): the INPUTS are what
the instrumented model recorded at coagsolv's entry (dumps/chem170_mosaic_coag,
fork f311091e: coag_dpdry, coag_distrib) plus the real64 replay's kernel; the
EXPECTED values are the real64 replay's splitting tables and rates
(chem170_scm/coag64/r64_<step>.json, kernels/coagsolv_rate_driver).  The script
replaces the component's `tests` array wholesale, so the case list below IS the
test list.

    python3 tools/fill_mosaic_coag_tests.py <dump_dir> <replay_dir>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402

COMPONENT = os.path.join(os.path.dirname(__file__), "..", "components", "aerosol", "mosaic",
                         "coagulation_assembly.esm")

# coagsolv's composition columns, in mosaic_coag_1clm's order: the ten
# ncomp_plustracer_aer masses, then hysw, water and the dry volume
COMP = ["so4", "no3", "cl", "co3", "nh4", "na", "ca", "oin", "oc", "bc", "hysw", "water", "dryvol"]
UNIT = {"hysw": "mol/mol-air", "water": "mol/mol-air", "dryvol": "cm3/mol-air"}

# (model step, level 1-based, composition quantity, what the case is for)
CASES = [
    (1, 1, "so4", "first step, surface layer: the initial, nearly unimodal Aitken/accumulation distribution"),
    (1, 1, "water", None),
    (1, 30, "nh4", "first step, mid-column (level 30): five times the surface number, mostly in the two smallest bins"),
    (1, 30, "dryvol", None),
    (689, 1, "oin", "step 689, surface layer: after the morning nucleation burst, the smallest bins hold 5x the initial number"),
    (689, 1, "hysw", None),
    (1621, 15, "oc", "step 1621, level 15 (inside the boundary layer late on day 2)"),
    (1621, 15, "so4", None),
    (2134, 45, "water", "step 2134, level 45 (free troposphere)"),
    (2134, 45, "nh4", None),
    (2521, 59, "dryvol", "last dumped step, top model layer"),
    (2521, 59, "so4", None),
]


# The component is written on the cell axes (lev, bin, cq); each case is one cell,
# so every input carries a leading lev axis of size 1 (and q a trailing cq axis of
# size 1), and the splitting-table references broadcast over lev.
def const3(a: np.ndarray) -> dict:
    return {"op": "faq", "output_idx": ["l", "i", "j", "k"], "args": [],
            "ranges": {"l": {"from": "lev"}, "i": {"from": "bin"}, "j": {"from": "bin"}, "k": {"from": "bin"}},
            "expr": {"op": "index", "args": [{"op": "const", "args": [], "value": a.tolist()}, "i", "j", "k"]}}


def const2(a: np.ndarray) -> dict:
    return {"op": "faq", "output_idx": ["l", "i", "j"], "args": [],
            "ranges": {"l": {"from": "lev"}, "i": {"from": "bin"}, "j": {"from": "bin"}},
            "expr": {"op": "index", "args": [{"op": "const", "args": [], "value": a.tolist()}, "i", "j"]}}


def main(dump_dir: str, replay_dir: str) -> None:
    doc = json.load(open(COMPONENT))
    tests = []
    for step, lev, comp, purpose in CASES:
        d = load(os.path.join(dump_dir, f"esm_dump_mosaic_{step}.json"))["vars"]
        r = load(os.path.join(replay_dir, f"r64_{step}.json"))["vars"]
        k = lev - 1
        c = COMP.index(comp)
        nbin = int(d["rate_nbin"])
        dp = d["coag_dpdry"][:, k]
        num = d["coag_distrib"][:, 0, k]
        q = d["coag_distrib"][:, 1 + c, k]
        cbr = r["coag_cbr"][:, :, k]
        fink = r["coag_fink_flat"][:, :, nbin * k:nbin * (k + 1)]
        floss = r["coag_floss"][:, :, k]
        dn = r["rate_coag_num"][:, k]
        dq = r["rate_coag_vol"][:, c, k]
        unit = UNIT.get(comp, "mol/mol-air")
        desc = (f"Reference: the em_scm_xy chem_opt = 170 column (data/eqwefic/chem170_scm, reference namelist "
                f"unchanged), MOSAIC coagulation capture dumps/chem170_mosaic_coag (WRF fork earthsciml-instrumented-chem170 "
                f"f311091e; bit-identical to dumps/chem170_wrf on every shared field), model step {step}, column level {lev}. "
                f"Inputs are coagsolv's own entry state: dry diameters coag_dpdry and the distribution coag_distrib "
                f"(number in #/cm3; the composition quantity here is column {c + 1} of 13, {comp}, in {unit}). "
                f"The kernel cbr and every expected value come from the REAL64 replay of coagsolv on those inputs "
                f"(kernels/coagsolv_rate_driver, chem170_scm/coag64/r64_{step}.json), because WRF's real32 coagsolv forms "
                f"the small high-side split shares as 1 - finklow and keeps only a few digits of them (FORTRAN_BUGS N117); "
                f"the real32 replay of the same inputs agrees with the in-model rates to real32 round-off. "
                f"Tolerances measured over all 16 steps x 59 levels, not assumed: splitting tables 4.4e-16 absolute "
                f"(bound 1e-12), number rates 1.6e-10 relative per bin (bound 1e-9), composition rates 1.9e-10 of the "
                f"column's largest bin rate (bound 1e-9 relative with an absolute floor of 1e-9 of that maximum, because "
                f"production and loss cancel to 1.4e-8 relative in the smallest entries).")
        if purpose:
            desc = f"{purpose[0].upper()}{purpose[1:]}. " + desc
        asserts = []
        if purpose:  # the splitting tables do not depend on the composition quantity: pin them once per state
            asserts.append({"variable": "fink", "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
                            "tolerance": {"abs": 1e-12, "rel": 0.0}, "reference": const3(fink)})
            asserts.append({"variable": "floss", "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
                            "tolerance": {"abs": 1e-12, "rel": 0.0}, "reference": const2(floss)})
            for b in range(nbin):
                asserts.append({"variable": "dnum_dt", "coords": {"lev": 1, "bin": b + 1}, "time": 0.0, "expected": float(dn[b]),
                                "tolerance": {"rel": 1e-9, "abs": 0.0}})
        floor = 1e-9 * float(np.abs(dq).max())
        for b in range(nbin):
            asserts.append({"variable": "dq_dt", "coords": {"lev": 1, "bin": b + 1, "cq": 1}, "time": 0.0, "expected": float(dq[b]),
                            "tolerance": {"rel": 1e-9, "abs": floor}})
        tests.append({
            "id": f"coag_s{step}_lev{lev}_{comp}",
            "description": desc,
            "parameter_overrides": {"dp_dry": [dp.tolist()], "cbr": [cbr.tolist()], "num": [num.tolist()],
                                    "q": [[[x] for x in q.tolist()]]},
            "time_span": {"start": 0, "end": 1},
            "assertions": asserts,
        })
    doc["models"]["MosaicCoagulationAssembly"]["tests"] = tests
    with open(COMPONENT, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    print(f"{len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
