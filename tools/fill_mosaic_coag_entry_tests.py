#!/usr/bin/env python3
"""Fill the inline tests of components/aerosol/mosaic/coagulation_entry_state.esm from
the chem170 MOSAIC coagulation capture (dumps/chem170_mosaic_coag, fork f311091e).

It transcribes numbers only; the component's equations are hand-written.  Each case
is one (model step, column level, size bin).  INPUTS: the bin's rows of the column
state entering coagulation (rsub_newnuc), the carried dry diagnostics (aqmassdry_sub,
aqvoldry_sub, adrydens_sub), cairclm and the section geometry.  EXPECTED: what the
instrumented model recorded at coagsolv's entry (coag_distrib, coag_dpdry, coag_dpwet,
coag_denswet).  The carried diagnostics are dumped after aerchemistry, BEFORE
nucleation, so a case is admitted only where nucleation left every row of the bin
unchanged (checked below); there the dumped diagnostics are the ones coagulation saw.
The script replaces the component's `tests` array wholesale, so the case list below
IS the test list.

    python3 tools/fill_mosaic_coag_entry_tests.py <dump_dir>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402

COMPONENT = os.path.join(os.path.dirname(__file__), "..", "components", "aerosol", "mosaic",
                         "coagulation_entry_state.esm")
COMP = ["so4", "no3", "cl", "co3", "nh4", "na", "ca", "oin", "oc", "bc"]
REL = 2e-6   # measured worst 4.7e-7 over 7426 cells (real32 reference)

# (model step, level 1-based, bin 1-based, regime)
CASES = [
    (918, 1, 1, "surface layer, smallest section, a step with no nucleation anywhere in the column"),
    (918, 30, 4, "mid-column, section 4 (the accumulation-mode peak)"),
    (1080, 10, 8, "inside the boundary layer, the largest section"),
    (2400, 1, 2, "surface layer late on day 2, section 2"),
    (2521, 59, 1, "top model layer at the last dumped step, smallest section"),
    (1, 1, 8, "first step, surface layer, largest section (nucleation active in the column but not in this section)"),
    (689, 20, 5, "step 689, level 20, section 5"),
    (2134, 45, 3, "step 2134, level 45 (free troposphere), section 3"),
]


def main(dump_dir: str) -> None:
    doc = json.load(open(COMPONENT))
    tests = []
    for step, lev, b1, regime in CASES:
        d = load(os.path.join(dump_dir, f"esm_dump_mosaic_{step}.json"))["vars"]
        k, b = lev - 1, b1 - 1
        mp = d["massptr_aer"][:, b] - 1
        npt, wp, hp = d["numptr_aer"][b] - 1, d["waterptr_aer"][b] - 1, d["hyswptr_aer"][b] - 1
        rows = list(mp) + [npt, wp, hp]
        rs, rt = d["rsub_newnuc"], d["rsub_therm"]
        assert np.array_equal(rs[rows, k], rt[rows, k]), (step, lev, b1, "nucleation touched this bin")
        dist = d["coag_distrib"][b, :, k]
        po = {
            "r": rs[mp, k].tolist(), "mw": d["mw_aer"].tolist(), "dens": d["dens_aer"].tolist(),
            "n_in": float(rs[npt, k]), "water_in": float(rs[wp, k]), "hysw_in": float(rs[hp, k]),
            "mass_in": float(d["aqmassdry_sub"][b, k]), "vol_in": float(d["aqvoldry_sub"][b, k]),
            "rho_in": float(d["adrydens_sub"][b, k]), "cair": float(d["cairclm"][k]),
            "v_lo": float(d["volumlo_sect"][b]), "v_hi": float(d["volumhi_sect"][b]), "v_cen": float(d["volumcen_sect"][b]),
            "d_lo": float(d["dlo_sect"][b]), "d_hi": float(d["dhi_sect"][b]), "d_cen": float(d["dcen_sect"][b]),
        }
        exact = {"abs": 0.0, "rel": 0.0}
        tol = {"rel": REL, "abs": 0.0}
        a = [
            {"variable": "implausible", "time": 0.0, "expected": 0.0, "tolerance": exact},
            {"variable": "empty", "time": 0.0, "expected": 0.0, "tolerance": exact},
            {"variable": "num_conc", "time": 0.0, "expected": float(dist[0]), "tolerance": tol},
            {"variable": "dp_dry", "time": 0.0, "expected": float(d["coag_dpdry"][b, k]), "tolerance": tol},
            {"variable": "dp_wet", "time": 0.0, "expected": float(d["coag_dpwet"][b, k]), "tolerance": tol},
            {"variable": "dens_wet", "time": 0.0, "expected": float(d["coag_denswet"][b, k]), "tolerance": tol},
            {"variable": "hysw", "time": 0.0, "expected": float(dist[11]), "tolerance": tol if dist[11] else exact},
            {"variable": "water", "time": 0.0, "expected": float(dist[12]), "tolerance": tol if dist[12] else exact},
            {"variable": "vol_dry", "time": 0.0, "expected": float(dist[13]), "tolerance": tol},
        ]
        for l, name in enumerate(COMP):
            v = float(dist[1 + l])
            a.append({"variable": "q_comp", "coords": {"comp": l + 1}, "time": 0.0, "expected": v,
                      "tolerance": tol if v else exact})
        desc = (f"{regime[0].upper()}{regime[1:]}: model step {step}, level {lev}, section {b1} of 8. Reference: the em_scm_xy "
                f"chem_opt = 170 column (data/eqwefic/chem170_scm, reference namelist unchanged), MOSAIC coagulation capture "
                f"dumps/chem170_mosaic_coag (WRF fork earthsciml-instrumented-chem170 f311091e; bit-identical to dumps/chem170_wrf "
                f"on every shared field). Inputs are the bin's rows of the state entering coagulation (rsub_newnuc) and the carried "
                f"dry diagnostics, which were dumped before nucleation and are admitted only because nucleation left every row of "
                f"this bin unchanged; expected values are what the model recorded at coagsolv's entry (coag_distrib, coag_dpdry, "
                f"coag_dpwet, coag_denswet), GNU gfortran 13.3 REAL*4. Plain path: neither the implausible-diagnostics recompute "
                f"nor the empty-bin reset fires (pinned exactly), and the number is inside the section. Tolerance rel {REL:g}: "
                f"measured 4.7e-7 worst over the 7426 admissible cells of the 16 dumped steps, the real32 rounding of the "
                f"reference; components the bin does not carry are pinned at exactly 0.")
        tests.append({"id": f"entry_s{step}_lev{lev}_bin{b1}", "description": desc, "parameter_overrides": po,
                      "time_span": {"start": 0, "end": 1}, "assertions": a})
    doc["models"]["MosaicCoagulationEntryState"]["tests"] = tests
    with open(COMPONENT, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    print(f"{len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main(sys.argv[1])
