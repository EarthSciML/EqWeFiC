#!/usr/bin/env python3
"""Fill the inline tests of components/atmospheric_deposition/wrf_wetdep_ls/wetdep_ls.esm
from the `wetdep_ls` dumps of the chem170 wet-scavenging capture and their real64
kernel replays (wrf-chem170/kernels/wetdep_ls_driver).

It transcribes numbers only; the component's equations are hand-written.  Each
case is one (model step, species) pair: the inputs are the in-model column the
hook dumped, the expected values are the replay's scavenged column.  The script
replaces the component's `tests` array wholesale, so the case list below IS the
test list.

    python3 tools/fill_wetdep_ls_tests.py <dump_dir> <replay_dir>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402

COMPONENT = os.path.join(os.path.dirname(__file__), "..", "components", "atmospheric_deposition",
                         "wrf_wetdep_ls", "wetdep_ls.esm")
REGISTRY_LINE = "package   cbmz_mosaic_kpp chem_opt==170"
PARAM_FIRST_SCALAR = 2

# (model step, species name); the step is also the wetdep_ls call index
CASES = [
    (1300, "so4_a03"),      # dry step (rain 4e-22 mm, under the gate): nothing scavenged
    (1405, "sulf"),         # onset: 2.7e-9 mm of rain, just over the 1e-10 gate
    (1405, "num_a03"),
    (1429, "so4_a04"),
    (1441, "sulf"),
    (1441, "water_a05"),
    (1453, "so4_a03"),      # most layers scavenged in the event
    (1453, "oc_a07"),
    (1468, "sulf"),         # peak rain
    (1468, "nh4_a04"),
    (1468, "num_a06"),
    (1480, "so4_a02"),      # tail: rain still falling, cloud water left in few layers
]


def species_names(registry: str) -> list[str]:
    for line in open(registry):
        if line.startswith(REGISTRY_LINE):
            lst = line.split()[-1]
            assert lst.startswith("chem:")
            return lst[len("chem:"):].split(",")
    raise SystemExit("cbmz_mosaic_kpp package line not found")


def const_ref(values: np.ndarray) -> dict:
    return {"op": "faq", "output_idx": ["k"], "args": [], "ranges": {"k": {"from": "lev"}},
            "expr": {"op": "index", "args": [
                {"op": "const", "args": [], "value": [float(v) for v in values]}, "k"]}}


def linf(var: str, values: np.ndarray, rel: float) -> dict:
    scale = float(np.max(np.abs(values)))
    return {"variable": var, "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
            "tolerance": {"abs": rel * scale}, "reference": const_ref(values)}


def main(dump_dir: str, replay_dir: str, registry: str) -> None:
    names = species_names(registry)
    tests = []
    for step, sp in CASES:
        nv = names.index(sp) + PARAM_FIRST_SCALAR       # 1-based chem index
        d = load(os.path.join(dump_dir, f"esm_dump_wetdep_ls_{step}.json"))["vars"]
        r = load(os.path.join(replay_dir, f"wetdep_ls_{step}.json"))["vars"]
        numgas = int(d["numgas"])
        alpha_name = "alpha_sulf" if sp == "sulf" else "alpha_aerosol"
        assert sp == "sulf" or nv > numgas
        c_in = np.asarray(d["chem_in"])[:, nv - 1].astype(float)
        c_out = np.asarray(r["chem_out"])[:, nv - 1].astype(float)
        c_wrf = np.asarray(d["chem_out"])[:, nv - 1].astype(float)
        dz = np.asarray(d["dz8w"], dtype=float)
        ze = np.concatenate([[0.0], np.cumsum(dz)])
        dt = float(d["dt"])
        rain = float(d["rain"])
        dcdt = (c_out - c_in) / dt
        scav = np.abs(c_out - c_in) > 0
        wrf_dev = float(np.max(np.abs(c_wrf - c_out)) / max(np.max(np.abs(c_out)), 1e-300))
        desc = (f"Reference: the chem_opt = 170 wet-scavenging capture of data/eqwefic/chem170_scm/run_wetscav.sh "
                f"(the reference column with wetscav_onoff = -1 and cu_physics = 0, so WSM6 makes the rain; fork "
                f"ctessum-claude/WRF branch earthsciml-instrumented-chem170, GNU gfortran real32), dumped by the "
                f"`wetdep_ls` hook in chem_driver.F at model step {step} (dt = {dt:g} s), column (its, jts), species "
                f"`{sp}` (chem index {nv}, alpha = par.{alpha_name}). Grid-scale rain this step {rain:.6g} mm; "
                f"{int(scav.sum())} of {len(c_in)} layers scavenged. Expected values from the real64 replay of the "
                f"dumped inputs through chem/module_wetdep_ls.F compiled unchanged (kernels/wetdep_ls_driver), so "
                f"the tolerance is 1e-9 of the column maximum; the in-model real32 result differs from the replay "
                f"by {wrf_dev:.2g} of the column maximum.")
        a = [{"variable": "active", "time": 0.0, "expected": 1.0 if scav.any() else 0.0,
              "tolerance": {"abs": 0.0, "rel": 0.0}},
             linf("c_out", c_out, 1e-9)]
        if scav.any():
            a.append(linf("dc_dt", dcdt, 1e-9))
        else:
            a.append({"variable": "dc_dt", "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
                      "tolerance": {"abs": 0.0, "rel": 0.0}, "reference": const_ref(dcdt)})
        tests.append({
            "id": f"chem170_wetscav_step{step}_{sp}",
            "description": desc,
            "parameter_overrides": {
                "c": [float(v) for v in c_in],
                "qc": [float(v) for v in np.asarray(d["qc"], dtype=float)],
                "rho": [float(v) for v in np.asarray(d["rho"], dtype=float)],
                "w": [float(v) for v in np.asarray(d["vvel"], dtype=float)],
                "ze": [float(v) for v in ze],
                "rain": rain, "dt": dt, "alpha": 1.0 if sp == "sulf" else 0.5,
            },
            "time_span": {"start": 0, "end": 1},
            "assertions": a,
        })
    doc = json.load(open(COMPONENT))
    doc["models"]["WetDepLS"]["tests"] = tests
    with open(COMPONENT, "w") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"{len(tests)} tests written")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
         "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/wrf-chem170/Registry/registry.chem")
