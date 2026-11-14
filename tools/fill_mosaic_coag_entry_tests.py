#!/usr/bin/env python3
"""Fill the inline tests of components/aerosol/mosaic/coagulation_entry_state.esm from
the chem170 MOSAIC coagulation capture (dumps/chem170_mosaic_coag, fork f311091e) and
its REAL64 replays (chem170_scm/coag_entry64, kernels/mosaic_subproc_driver).

It transcribes numbers only; the component's equations are hand-written.  Each case
is one (model step, column level, size bin).  INPUTS: the bin's rows of the column
state entering coagulation (rsub_newnuc), the carried dry diagnostics (aqmassdry_sub,
aqvoldry_sub, adrydens_sub), cairclm and the section geometry.  EXPECTED: what the
instrumented model recorded at coagsolv's entry (coag_distrib, coag_dpdry, coag_dpwet,
coag_denswet).  The carried diagnostics are dumped after aerchemistry, BEFORE
nucleation, so a case is admitted only where nucleation left every row of the bin
unchanged (checked below); there the dumped diagnostics are the ones coagulation saw.

Since 2026-11-12 the EXPECTED values come from the real64 replay of mosaic_coag_1clm on
the same inputs (base_<step>.json, rel 1e-9), not from the real32 in-model record
(rel 2e-6), and seven PERTURBED cases drive the branches no dumped step takes
(tools/perturb_mosaic_coag_entry.py writes their inputs; pert_<case>.json are the
replays).  The in-model record is still read and must agree with each unperturbed
case's replay to real32 rounding (checked below).
The script replaces the component's `tests` array wholesale, so the case list below
IS the test list.

    python3 tools/fill_mosaic_coag_entry_tests.py <dump_dir> <replay_dir>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402
import perturb_mosaic_coag_entry as pert  # noqa: E402

COMPONENT = os.path.join(os.path.dirname(__file__), "..", "components", "aerosol", "mosaic",
                         "coagulation_entry_state.esm")
COMP = ["so4", "no3", "cl", "co3", "nh4", "na", "ca", "oin", "oc", "bc"]
REL = 1e-9   # double-precision convention against the real64 replay
REL32 = 2e-6  # in-model real32 record vs the replay: measured 4.7e-7 worst

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


def case_test(v, rep, step, lev, b1, regime, flags, kind, edit_note=None):
    """One test: inputs from the (possibly edited) state ``v`` entering coagulation,
    expected from the real64 replay ``rep``."""
    k, b = lev - 1, b1 - 1
    mp = v["massptr_aer"][:, b] - 1
    npt, wp, hp = v["numptr_aer"][b] - 1, v["waterptr_aer"][b] - 1, v["hyswptr_aer"][b] - 1
    rows = list(mp) + [npt, wp, hp]
    rs = rep["rsub_newnuc_replay"]
    assert np.array_equal(rs[rows, k], v["rsub_therm"][rows, k]), (step, lev, b1, "nucleation touched this bin")
    dist = rep["coag_distrib_replay"][b, :, k]
    po = {
        "r": [[rs[mp, k].tolist()]], "mw": v["mw_aer"].tolist(), "dens": v["dens_aer"].tolist(),
        "n_in": float(rs[npt, k]), "water_in": float(rs[wp, k]), "hysw_in": float(rs[hp, k]),
        "mass_in": float(v["aqmassdry_sub"][b, k]), "vol_in": float(v["aqvoldry_sub"][b, k]),
        "rho_in": float(v["adrydens_sub"][b, k]), "cair": float(v["cairclm"][k]),
        "v_lo": float(v["volumlo_sect"][b]), "v_hi": float(v["volumhi_sect"][b]), "v_cen": float(v["volumcen_sect"][b]),
        "d_lo": float(v["dlo_sect"][b]), "d_hi": float(v["dhi_sect"][b]), "d_cen": float(v["dcen_sect"][b]),
    }
    exact = {"abs": 0.0, "rel": 0.0}
    tol = {"rel": REL, "abs": 0.0}
    implausible, empty = flags
    a = [
        {"variable": "implausible", "time": 0.0, "expected": implausible, "tolerance": exact},
        {"variable": "empty", "time": 0.0, "expected": empty, "tolerance": exact},
        {"variable": "num_conc", "time": 0.0, "expected": float(dist[0]), "tolerance": tol},
        {"variable": "dp_dry", "time": 0.0, "expected": float(rep["coag_dpdry_replay"][b, k]), "tolerance": tol},
        {"variable": "dp_wet", "time": 0.0, "expected": float(rep["coag_dpwet_replay"][b, k]), "tolerance": tol},
        {"variable": "dens_wet", "time": 0.0, "expected": float(rep["coag_denswet_replay"][b, k]), "tolerance": tol},
        {"variable": "hysw", "time": 0.0, "expected": float(dist[11]), "tolerance": tol if dist[11] else exact},
        {"variable": "water", "time": 0.0, "expected": float(dist[12]), "tolerance": tol if dist[12] else exact},
        {"variable": "vol_dry", "time": 0.0, "expected": float(dist[13]), "tolerance": tol},
    ]
    for x in a:  # the component is on the cell axes (lev, bin): each case is one cell, r is [1][1][comp]
        x["coords"] = {"lev": 1, "bin": 1}
    for l, name in enumerate(COMP):
        x = float(dist[1 + l])
        a.append({"variable": "q_comp", "coords": {"lev": 1, "bin": 1, "comp": l + 1}, "time": 0.0, "expected": x,
                  "tolerance": tol if x else exact})
    common = (f"Reference: the em_scm_xy chem_opt = 170 column (data/eqwefic/chem170_scm, reference namelist "
              f"unchanged), MOSAIC coagulation capture dumps/chem170_mosaic_coag (WRF fork earthsciml-instrumented-chem170 "
              f"f311091e). Expected values are what coagsolv is handed in a REAL64 replay of mosaic_coag_1clm "
              f"(kernels/mosaic_subproc_driver built with the coag entry hooks, fork 66a1b222, -fdefault-real-8; "
              f"chem170_scm/coag_entry64/{'base_%d' % step if kind == 'base' else 'pert_' + kind}.json) on these same inputs; "
              f"tolerance rel {REL:g}, the double-precision convention. Inputs are the bin's rows of the state entering "
              f"coagulation (nucleation leaves this section untouched, checked) and the carried dry diagnostics. ")
    if kind == "base":
        desc = (f"{regime[0].upper()}{regime[1:]}: model step {step}, level {lev}, section {b1} of 8. " + common +
                f"Plain path: neither the implausible-diagnostics recompute nor the empty-bin reset fires (pinned exactly), "
                f"and the number is inside the section. The real32 in-model record at coagsolv's entry agrees with this "
                f"replay to {REL32:g} (4.7e-7 worst over the 7426 admissible cells of the 16 dumped steps); components the "
                f"bin does not carry are pinned at exactly 0.")
        tid = f"entry_s{step}_lev{lev}_bin{b1}"
    else:
        desc = (f"BRANCH CASE, PERTURBED INPUT ({kind}): {regime}. Model step {step}, level {lev}, section {b1} of 8, with "
                f"that one section edited before the replay (EqWeFiC tools/perturb_mosaic_coag_entry.py) -- no dumped step "
                f"takes this branch, so the input is synthetic and the expected value is still the Fortran's. " + common +
                f"Flags pinned exactly: implausible = {implausible:g}, empty = {empty:g}.")
        tid = f"entry_branch_{kind}_s{step}_lev{lev}_bin{b1}"
    return {"id": tid, "description": desc, "parameter_overrides": po,
            "time_span": {"start": 0, "end": 1}, "assertions": a}


FLAGS = {"dens_high": (1.0, 0.0), "dens_low": (1.0, 0.0), "mass_mismatch": (1.0, 0.0),
         "negative_comp": (1.0, 0.0), "empty": (0.0, 1.0), "num_high": (0.0, 0.0), "num_low": (0.0, 0.0)}


def main(dump_dir: str, replay_dir: str) -> None:
    doc = json.load(open(COMPONENT))
    tests = []
    for step, lev, b1, regime in CASES:
        v = load(os.path.join(dump_dir, f"esm_dump_mosaic_{step}.json"))["vars"]
        rep = load(os.path.join(replay_dir, f"base_{step}.json"))["vars"]
        k, b = lev - 1, b1 - 1
        for a, m in [("coag_dpdry_replay", "coag_dpdry"), ("coag_dpwet_replay", "coag_dpwet"),
                     ("coag_denswet_replay", "coag_denswet")]:
            assert abs(rep[a][b, k] - v[m][b, k]) <= REL32 * abs(v[m][b, k]), (step, lev, b1, a)
        tests.append(case_test(v, rep, step, lev, b1, regime, (0.0, 0.0), "base"))
    for case in pert.CASES:
        cid, step, lev, b1, edit, regime = case
        v = pert.edited(dump_dir, case)["vars"]
        rep = load(os.path.join(replay_dir, f"pert_{cid}.json"))["vars"]
        tests.append(case_test(v, rep, step, lev, b1, regime, FLAGS[cid], cid))
    doc["models"]["MosaicCoagulationEntryState"]["tests"] = tests
    with open(COMPONENT, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    print(f"{len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
