#!/usr/bin/env python3
"""Fill the inline tests of components/atmospheric_deposition/wrf_wetdep_ls/wetdep_ls.esm
from the `wetdep_ls` dumps of the chem170 wet-scavenging capture and their real64
kernel replays (wrf-chem170/kernels/wetdep_ls_driver).

It transcribes numbers only; the component's equations are hand-written.  Each
case is one model step on a two-species axis, `sulf` (the one gas, NGAS = 1) and
one aerosol species: the inputs are the in-model columns the hook dumped, the
expected values are the replay's scavenged columns.  The script
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

# (model step, aerosol species); sulf is always the first entry of the species axis
CASES = [
    (1300, "so4_a03"),      # dry step (rain 4e-22 mm, under the gate): nothing scavenged
    (1405, "num_a03"),      # onset: 2.7e-9 mm of rain, just over the 1e-10 gate
    (1429, "so4_a04"),
    (1441, "water_a05"),
    (1453, "so4_a03"),      # most layers scavenged in the event
    (1453, "oc_a07"),
    (1468, "nh4_a04"),      # peak rain
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
    """A reference over the asserted field's own axes: (lev) or (lev, spc)."""
    values = np.asarray(values, dtype=float)
    out, ranges = (["k"], {"k": {"from": "lev"}}) if values.ndim == 1 else \
        (["k", "s"], {"k": {"from": "lev"}, "s": {"from": "spc"}})
    return {"op": "faq", "output_idx": out, "args": [], "ranges": ranges,
            "expr": {"op": "index", "args": [
                {"op": "const", "args": [], "value": values.tolist()}, *out]}}


def linf(var: str, values: np.ndarray, rel: float, ref: dict | None = None) -> dict:
    """Linf over the whole field at `rel` of its largest magnitude; exact where it is all zero."""
    scale = float(np.max(np.abs(values)))
    tol = {"abs": rel * scale, "rel": 0.0} if scale > 0 else {"abs": 0.0, "rel": 0.0}
    return {"variable": var, "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
            "tolerance": tol, "reference": ref if ref is not None else const_ref(values)}


def active_ref(flags: np.ndarray) -> dict:
    return {"op": "faq", "output_idx": ["s"], "args": [], "ranges": {"s": {"from": "spc"}},
            "expr": {"op": "index", "args": [
                {"op": "const", "args": [], "value": [float(v) for v in flags]}, "s"]}}


def removal_fraction(c_in: np.ndarray, c_out: np.ndarray) -> np.ndarray:
    """The fraction each layer lost, read from the replay where it changed; refuses a floored cell,
    where c_out = 1e-16 and the fraction is no longer 1 - c_out/c_in."""
    changed = c_out != c_in
    assert not np.any(changed & (c_out == 1e-16)), "floored cell: removal not recoverable"
    return np.where(changed, 1.0 - c_out / np.where(changed, c_in, 1.0), 0.0)


def column_inputs(d: dict) -> dict:
    dz = np.asarray(d["dz8w"], dtype=float)
    return {"qc": [float(v) for v in np.asarray(d["qc"], dtype=float)],
            "rho": [float(v) for v in np.asarray(d["rho"], dtype=float)],
            "w": [float(v) for v in np.asarray(d["vvel"], dtype=float)],
            "ze": [float(v) for v in np.concatenate([[0.0], np.cumsum(dz)])],
            "rain": float(d["rain"]), "dt": float(d["dt"])}


def species_assertions(prefix: str, c_in: np.ndarray, c_out: np.ndarray, dt: float) -> list[dict]:
    """active (exact), removal and the tendency (1e-9 of their largest magnitudes) on a
    (lev, spc) axis."""
    scav = np.any(c_out != c_in, axis=0)
    return [{"variable": f"{prefix}active", "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
             "tolerance": {"abs": 0.0, "rel": 0.0}, "reference": active_ref(scav.astype(float))},
            linf(f"{prefix}removal", removal_fraction(c_in, c_out), 1e-9),
            linf(f"{prefix}dc_dt", (c_out - c_in) / dt, 1e-9)]


def main(dump_dir: str, replay_dir: str, registry: str) -> None:
    names = species_names(registry)
    nv_sulf = names.index("sulf") + PARAM_FIRST_SCALAR       # 1-based chem index
    tests = []
    for step, sp in CASES:
        nv = names.index(sp) + PARAM_FIRST_SCALAR
        d = load(os.path.join(dump_dir, f"esm_dump_wetdep_ls_{step}.json"))["vars"]
        r = load(os.path.join(replay_dir, f"wetdep_ls_{step}.json"))["vars"]
        assert nv > int(d["numgas"]) >= nv_sulf
        cols = [nv_sulf - 1, nv - 1]
        c_in = np.asarray(d["chem_in"], dtype=float)[:, cols]
        c_out = np.asarray(r["chem_out"], dtype=float)[:, cols]
        c_wrf = np.asarray(d["chem_out"], dtype=float)[:, cols]
        over = column_inputs(d)
        dt, rain = over["dt"], over["rain"]
        scav = np.abs(c_out - c_in) > 0
        wrf_dev = float(np.max(np.max(np.abs(c_wrf - c_out), axis=0) /
                               np.maximum(np.max(np.abs(c_out), axis=0), 1e-300)))
        desc = (f"Reference: the chem_opt = 170 wet-scavenging capture of data/eqwefic/chem170_scm/run_wetscav.sh "
                f"(the reference column with wetscav_onoff = -1 and cu_physics = 0, so WSM6 makes the rain; fork "
                f"ctessum-claude/WRF branch earthsciml-instrumented-chem170, GNU gfortran real32), dumped by the "
                f"`wetdep_ls` hook in chem_driver.F at model step {step} (dt = {dt:g} s), column (its, jts), on the "
                f"species axis [sulf (chem index {nv_sulf}, alpha = par.alpha_sulf), {sp} (chem index {nv}, alpha = "
                f"par.alpha_aerosol)]. Grid-scale rain this step {rain:.6g} mm; {int(scav[:, 0].sum())} and "
                f"{int(scav[:, 1].sum())} of {len(c_in)} layers scavenged. Expected values from the real64 replay of "
                f"the dumped inputs through chem/module_wetdep_ls.F compiled unchanged (kernels/wetdep_ls_driver), so "
                f"the tolerance is 1e-9 of each field's largest magnitude; the in-model real32 result differs from the "
                f"replay by {wrf_dev:.2g} of the column maximum.")
        over["c"] = c_in.tolist()
        tests.append({
            "id": f"chem170_wetscav_step{step}_sulf_{sp}",
            "description": desc,
            "parameter_overrides": over,
            "time_span": {"start": 0, "end": 0},
            "assertions": species_assertions("", c_in, c_out, dt),
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
