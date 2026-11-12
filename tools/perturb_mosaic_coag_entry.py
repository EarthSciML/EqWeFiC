#!/usr/bin/env python3
"""Write perturbed mosaic_subproc_driver inputs that drive the branches of MOSAIC's
coagulation entry block (mosaic_coag_1clm, chem/module_mosaic_coag.F:155-293) that the
16 dumped chem170 steps never take.

Each case starts from one step of the MOSAIC coagulation capture
(dumps/chem170_mosaic_coag) and edits ONE section at ONE level -- never section 1,
which nucleation writes, so the state entering coagulation at the edited section is
exactly the edited input.  The real64 replay (kernels/mosaic_subproc_driver, built
with the coag entry hooks enabled) then records what coagsolv is handed.  This script
writes inputs only; the expected values come from the Fortran.

    python3 tools/perturb_mosaic_coag_entry.py <dump_dir> <out_dir>

Writes <out_dir>/<case>.flat and <out_dir>/cases.json (case -> step, level, section,
what was edited and why).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import esm_dump  # noqa: E402

COMP = ["so4", "no3", "cl", "co3", "nh4", "na", "ca", "oin", "oc", "bc"]

# (case id, step, level 1-based, section 1-based, edit, regime)
CASES = [
    ("dens_high", 918, 1, 4, "dens_high",
     "carried dry density 25 g/cm3, above the 20 g/cm3 plausibility bound: mass and volume recomputed from the composition"),
    ("dens_low", 1080, 10, 3, "dens_low",
     "carried dry density 0.05 g/cm3, below the 0.1 g/cm3 bound: recomputed from the composition"),
    ("mass_mismatch", 918, 30, 5, "mass_mismatch",
     "carried dry mass 1e-3 relative above volume x density (plausible density, inconsistent triple): recomputed"),
    ("negative_comp", 2400, 1, 6, "negative_comp",
     "implausible density AND a negative mixing ratio for the section's second-largest component by mass, which the recompute clips at zero"),
    ("empty", 918, 1, 3, "empty",
     "every component zeroed and the carried dry mass 1e-31 g/mol (below 1.01e-30), in a section that carries water and hysteresis water: empty-bin reset, both waters dropped, number put on the section centre"),
    ("num_high", 2134, 45, 3, "num_high",
     "number x 1e3, so the mean dry volume per particle is below the section's lower bound: number lowered to V/(f v_lo)"),
    ("num_low", 689, 20, 5, "num_low",
     "number x 1e-3, so the mean dry volume per particle is above the section's upper bound: number raised to V f/v_hi"),
]


def edited(dump_dir: str, case: tuple) -> dict:
    """The capture step of ``case`` with the case's edit applied (vars, in place)."""
    cid, step, lev, b1, edit, regime = case
    d = esm_dump.load(os.path.join(dump_dir, f"esm_dump_mosaic_{step}.json"))
    d["vars"].pop("coag_fink")          # rank 4; the flat reader stops at rank 3
    v = d["vars"]
    k, b = lev - 1, b1 - 1
    mp = v["massptr_aer"][:, b] - 1
    npt = v["numptr_aer"][b] - 1
    rs = v["rsub_therm"]
    mass, vol, rho = v["aqmassdry_sub"], v["aqvoldry_sub"], v["adrydens_sub"]
    if edit == "dens_high":
        rho[b, k] = 25.0
    elif edit == "dens_low":
        rho[b, k] = 0.05
    elif edit == "mass_mismatch":
        mass[b, k] = vol[b, k] * rho[b, k] * 1.001
    elif edit == "negative_comp":
        rho[b, k] = 25.0
        # the component with the SECOND-largest mass is made negative, so that
        # clipping it (or not) moves the recomputed dry mass well beyond round-off
        r, mw = rs[mp, k], v["mw_aer"]
        i_neg = sorted(range(len(mp)), key=lambda i: r[i] * mw[i])[-2]
        rs[mp[i_neg], k] = -abs(r[i_neg])
    elif edit == "empty":
        rs[mp, k] = 0.0
        mass[b, k] = 1e-31
        vol[b, k] = 1e-31 / rho[b, k]
    elif edit == "num_high":
        rs[npt, k] *= 1e3
    elif edit == "num_low":
        rs[npt, k] *= 1e-3
    else:
        raise ValueError(edit)
    return d


def main(dump_dir: str, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    meta = {}
    for case in CASES:
        cid, step, lev, b1, edit, regime = case
        esm_dump.to_flat(edited(dump_dir, case), os.path.join(out_dir, f"{cid}.flat"))
        meta[cid] = {"step": step, "lev": lev, "bin": b1, "edit": edit, "regime": regime}
    with open(os.path.join(out_dir, "cases.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"{len(meta)} perturbed inputs in {out_dir}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
