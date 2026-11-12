#!/usr/bin/env python3
"""Fill the inline tests of couplings/wetdep_ls_column.esm from the chem170
wet-scavenging capture (data/eqwefic/chem170_scm/run_wetscav.sh) and its real64
replays (kernels/wetdep_ls_driver).  Numbers only; the coupling is hand-written.

Each case is one model step: ONE column state (qc, rho, w, ze, rain, dt) and all 169
scavenged species -- sulf, then every aerosol species (nv > numgas) in registry order,
the order WetDepColumn.c states -- feed the single WetDepLS mount.  The script also
checks the claim the species axis rests on -- that wetdep_ls changes nothing off it --
and refuses a step where another gas changes.

    python3 tools/fill_wetdep_ls_column_tests.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402
from fill_wetdep_ls_tests import (PARAM_FIRST_SCALAR, column_inputs, species_assertions,  # noqa: E402
                                  species_names)

DATA = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic"
DUMPS = f"{DATA}/dumps/chem170_wetscav"
REGISTRY = f"{DATA}/wrf-chem170/Registry/registry.chem"
DOC = os.path.join(os.path.dirname(__file__), "..", "couplings", "wetdep_ls_column.esm")
STEPS = [1300, 1405, 1441, 1453, 1468, 1480]


def main():
    names = species_names(REGISTRY)
    doc = json.load(open(DOC))
    tests = []
    for step in STEPS:
        d = load(f"{DUMPS}/esm_dump_wetdep_ls_{step}.json")["vars"]
        r = load(f"{DUMPS}/replay/wetdep_ls_{step}.json")["vars"]
        numgas = int(d["numgas"])
        cin = np.asarray(d["chem_in"], dtype=float)
        rout = np.asarray(r["chem_out"], dtype=float)
        wout = np.asarray(d["chem_out"], dtype=float)
        axis = [names.index("sulf") + PARAM_FIRST_SCALAR] + list(range(numgas + 1, cin.shape[1] + 1))
        assert len(axis) == doc["metaparameters"]["NSPEC"]["default"] == 169
        changed = set(np.nonzero(np.any(rout != cin, axis=0) | np.any(wout != cin, axis=0))[0] + 1)
        off = sorted(changed - set(axis))
        if off:
            raise SystemExit(f"step {step}: species off the axis changed: "
                             f"{[names[nv - PARAM_FIRST_SCALAR] for nv in off]}")
        cols = [nv - 1 for nv in axis]
        c_in, c_out = cin[:, cols], rout[:, cols]
        over = {f"WetDepColumn.{k}": v for k, v in column_inputs(d).items()}
        over["WetDepColumn.c"] = c_in.tolist()
        dt = over["WetDepColumn.dt"]
        nscav = int(np.any(c_out != c_in, axis=0).sum())
        tests.append({
            "id": f"chem170_wetscav_step{step}",
            "description": (
                f"Reference: the chem_opt = 170 wet-scavenging capture (data/eqwefic/chem170_scm/run_wetscav.sh: the "
                f"reference column with wetscav_onoff = -1 and cu_physics = 0; fork ctessum-claude/WRF branch "
                f"earthsciml-instrumented-chem170, GNU gfortran real32), the `wetdep_ls` dump at model step {step} "
                f"(dt = {dt:g} s), column (its, jts); grid-scale rain {float(d['rain']):.6g} mm. All 169 scavenged "
                f"species on one species axis; {nscav} of them are scavenged at this step and the other "
                f"{169 - nscav} must come out unchanged. Expected values are the real64 replay of the dumped inputs "
                f"through chem/module_wetdep_ls.F compiled unchanged (kernels/wetdep_ls_driver): activation flags "
                f"exactly, removal fraction and tendency at 1e-9 of their largest magnitudes. {len(changed)} of the "
                f"233 chem species change in the whole dump, every one on the axis (checked by the fill script)."),
            "parameter_overrides": over, "time_span": {"start": 0, "end": 0},
            "assertions": species_assertions("WetDep.", c_in, c_out, dt)})
    doc["models"]["WetDepColumn"]["tests"] = tests
    with open(DOC, "w") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"{len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main()
