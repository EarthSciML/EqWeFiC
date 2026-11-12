#!/usr/bin/env python3
"""Fill the inline tests of couplings/wetdep_ls_column.esm from the chem170
wet-scavenging capture (data/eqwefic/chem170_scm/run_wetscav.sh) and its real64
replays (kernels/wetdep_ls_driver).  Numbers only; the coupling is hand-written.

Each case is one model step: ONE column state (qc, rho, w, ze, rain, dt) feeds every
mounted WetDepLS instance, and each instance's species comes from the same dump.  The
script also checks the claim the coupling rests on -- that wetdep_ls changes only
sulf and the aerosol species (nv > numgas) -- and refuses a step where a gas changes.

    python3 tools/fill_wetdep_ls_column_tests.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402
from fill_wetdep_ls_tests import PARAM_FIRST_SCALAR, const_ref, species_names  # noqa: E402

DATA = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic"
DUMPS = f"{DATA}/dumps/chem170_wetscav"
REGISTRY = f"{DATA}/wrf-chem170/Registry/registry.chem"
DOC = os.path.join(os.path.dirname(__file__), "..", "couplings", "wetdep_ls_column.esm")
STEPS = [1300, 1405, 1441, 1453, 1468, 1480]


def linf(var, values):
    scale = float(np.max(np.abs(values)))
    tol = {"abs": 1e-9 * scale, "rel": 0.0} if scale > 0 else {"abs": 0.0, "rel": 0.0}
    return {"variable": var, "time": 0.0, "reduce": "Linf_error", "expected": 0.0, "tolerance": tol,
            "reference": const_ref(values)}


def main():
    names = species_names(REGISTRY)
    doc = json.load(open(DOC))
    mounts = {k: None for k, v in doc["models"].items() if v.get("ref", "").endswith("/wetdep_ls.esm")}
    # the species of each mount is named in its alpha edge's description
    for c in doc["coupling"]:
        if c["to"].endswith(".alpha"):
            mounts[c["to"].split(".")[0]] = c["description"].split()[0]
    tests = []
    for step in STEPS:
        d = load(f"{DUMPS}/esm_dump_wetdep_ls_{step}.json")["vars"]
        r = load(f"{DUMPS}/replay/wetdep_ls_{step}.json")["vars"]
        numgas, dt = int(d["numgas"]), float(d["dt"])
        cin = np.asarray(d["chem_in"], dtype=float)
        rout = np.asarray(r["chem_out"], dtype=float)
        wout = np.asarray(d["chem_out"], dtype=float)
        changed = set(np.nonzero(np.any(rout != cin, axis=0) | np.any(wout != cin, axis=0))[0] + 1)
        bad = [names[nv - PARAM_FIRST_SCALAR] for nv in changed
               if nv <= numgas and names[nv - PARAM_FIRST_SCALAR] != "sulf"]
        if bad:
            raise SystemExit(f"step {step}: gases {bad} changed")
        dz = np.asarray(d["dz8w"], dtype=float)
        over = {"qc": [float(v) for v in np.asarray(d["qc"], dtype=float)],
                "rho": [float(v) for v in np.asarray(d["rho"], dtype=float)],
                "w": [float(v) for v in np.asarray(d["vvel"], dtype=float)],
                "ze": [float(v) for v in np.concatenate([[0.0], np.cumsum(dz)])],
                "rain": float(d["rain"]), "dt": dt}
        a, nscav = [], 0
        for key, sp in mounts.items():
            nv = names.index(sp) + PARAM_FIRST_SCALAR
            c_in, c_out = cin[:, nv - 1], rout[:, nv - 1]
            over[f"{key}.c"] = [float(v) for v in c_in]
            scav = bool(np.any(c_out != c_in))
            nscav += scav
            a.append({"variable": f"{key}.active", "time": 0.0, "expected": 1.0 if scav else 0.0,
                      "tolerance": {"abs": 0.0, "rel": 0.0}})
            a.append(linf(f"{key}.c_out", c_out))
            a.append(linf(f"{key}.dc_dt", (c_out - c_in) / dt))
        tests.append({
            "id": f"chem170_wetscav_step{step}",
            "description": (
                f"Reference: the chem_opt = 170 wet-scavenging capture (data/eqwefic/chem170_scm/run_wetscav.sh: the "
                f"reference column with wetscav_onoff = -1 and cu_physics = 0; fork ctessum-claude/WRF branch "
                f"earthsciml-instrumented-chem170, GNU gfortran real32), the `wetdep_ls` dump at model step {step} "
                f"(dt = {dt:g} s), column (its, jts); grid-scale rain {float(d['rain']):.6g} mm. ONE column state feeds "
                f"all {len(mounts)} mounted instances; only each instance's species differs, and its efficiency "
                f"comes from the shared parameter mount. Expected values are the real64 replay of the dumped "
                f"inputs through chem/module_wetdep_ls.F compiled unchanged (kernels/wetdep_ls_driver), at 1e-9 of "
                f"each column's maximum; {nscav} of the {len(mounts)} mounted species are scavenged at this step, "
                f"and {len(changed)} of the 233 chem species change in the whole dump, every one of them sulf or an "
                f"aerosol species (checked by the fill script)."),
            "parameter_overrides": over, "time_span": {"start": 0, "end": 0}, "assertions": a})
    doc["models"]["WetDepColumn"]["tests"] = tests
    with open(DOC, "w") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"{len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main()
