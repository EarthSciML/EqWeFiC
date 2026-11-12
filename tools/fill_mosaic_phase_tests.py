#!/usr/bin/env python3
"""Fill the inline tests of components/aerosol/mosaic/aerosol_phase_state.esm.

Sources (numbers only; the component's equations are hand-written):
  - the chem170 aerosol_phase_state capture, dumps/chem170_mosaic_phase (WRF fork
    earthsciml-instrumented-chem170, esm_ps_* hooks): FIRST (_c1) and LAST (_c2) call
    per (bin, level) of each dumped step, selected below;
  - hand-specified bins run through WRF's own compiled aerosol_phase_state by
    kernels/mosaic_therm_driver (mode 2), chem170_scm/therm64/{cases_ps.json,ps.json}.
Both record real(kind=8) values (the thermodynamics is double precision), so the
reference is double precision: rel 1e-9.

What is asserted depends on the MESA branch the Fortran took (ps_sc path):
  path 2 (liquid): everything the component computes;
  paths 1, 4, 5 (solid): state, path, niter, ah2o_a, kelvin, water = 0, the partition;
  path 6 (mesa_ptc, not transcribed): only path, niter = 1 is NOT implied, so only path.

    python3 tools/fill_mosaic_phase_tests.py <dump_dir|-> <ps.json>
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402

COMP = os.path.join(os.path.dirname(__file__), "..", "components", "aerosol", "mosaic", "aerosol_phase_state.esm")
REL = 1e-9
SC = {n: i for i, n in enumerate("water_in hyst_in up_in num jhyst_in state phase jhyst path niter ah2o_a kelvin water "
                                 "water_up mass_dry vol_dry mass_wet vol_wet dpmv sigma growth mdrh rh_pc".split())}
IONS = {1: "SO4", 2: "NO3", 3: "Cl", 4: "NH4", 5: "MSA", 107: "CO3", 108: "Na", 109: "Ca", 110: "OIN", 111: "OC", 112: "BC"}
SRC_COL = ("Reference: the em_scm_xy chem_opt = 170 column (data/eqwefic/chem170_scm, reference namelist unchanged), "
           "aerosol_phase_state capture dumps/chem170_mosaic_phase (WRF fork earthsciml-instrumented-chem170, the esm_ps_* "
           "hooks); the thermodynamics runs in real(kind=8), a double-precision reference")
SRC_DRV = ("Reference: WRF's own compiled aerosol_phase_state (main/libwrflib.a of the chem170 build, the model's flags) run "
           "on a hand-specified bin by kernels/mosaic_therm_driver mode 2 (conform_electrolytes, then aerosol_phase_state), "
           "cases in chem170_scm/therm64/make_cases.py, recorded by the same esm_ps_* hooks; double precision")
# column cases: (step, level, bin, slot, why)
COL = [
    (918, 1, 1, "c1", "upper leg at the surface, mid-afternoon, one Kelvin pass"),
    (918, 21, 1, "c1", "upper leg, smallest section, two Kelvin passes"),
    (1621, 49, 5, "c2", "upper leg at 71% relative humidity, the last call of the step"),
    (1621, 56, 1, "c2", "upper leg in the upper troposphere, two Kelvin passes, the last call"),
    (1, 53, 1, "c1", "effloresced below crh in the upper troposphere"),
    (1621, 59, 8, "c2", "effloresced below crh at the model top, largest section, the last call"),
    (1, 48, 1, "c1", "lower leg below the mutual deliquescence point"),
    (1621, 56, 8, "c2", "lower leg below the mutual deliquescence point, largest section, the last call"),
    (918, 57, 2, "c1", "excess-sulfate bin handed to mesa_ptc, which converges on mass (all solid)"),
    (1, 11, 1, "c1", "handed to mesa_ptc, which converges on the fluxes"),
    (1, 50, 2, "c2", "handed to mesa_ptc, which does NOT converge in nmax_mesa = 80 iterations"),
]


def tol(x):
    return {"rel": REL, "abs": 0.0} if x != 0 else {"rel": 0.0, "abs": 0.0}


EXACT = {"rel": 0.0, "abs": 0.0}


def build(v, k, b, slot, tid, why, src):
    s = v[f"ps_sc_{slot}"][:, b, k]
    path = int(s[SC["path"]])
    path = path % 10 if path >= 6 else path   # 6 + 10*outcome + 100*iterations for mesa_ptc
    po = {"rh_pc": float(s[SC["rh_pc"]]), "t_k": float(v["ps_t_k"][k]), "num": float(s[SC["num"]]),
          "aer_t": v[f"ps_aert_in_{slot}"][:, b, k].tolist(), "elec": v[f"ps_elet_in_{slot}"][:22, b, k].tolist(),
          "water_hyst": float(s[SC["hyst_in"]])}
    assert float(v[f"ps_elet_in_{slot}"][22, b, k]) == 0.0
    a = [{"variable": "path", "time": 0.0, "expected": float(path), "tolerance": EXACT}]
    def put(var, x, coords=None):
        d = {"variable": var, "time": 0.0, "expected": float(x), "tolerance": tol(float(x))}
        if coords:
            d["coords"] = coords
        a.append(d)
    def put_exact(var, x):
        a.append({"variable": var, "time": 0.0, "expected": float(x), "tolerance": EXACT})
    put("water_up", s[SC["water_up"]])
    put("mass_dry", s[SC["mass_dry"]])
    if path != 6:
        put_exact("state", s[SC["state"]])
        put_exact("jhyst", s[SC["jhyst"]])
        put_exact("niter", s[SC["niter"]])
        put("ah2o_a", s[SC["ah2o_a"]])
        put("kelvin", s[SC["kelvin"]])
        put("water", s[SC["water"]])
        put("vol_dry", s[SC["vol_dry"]])
        if path == 5:
            put("mdrh", s[SC["mdrh"]])
        if path == 2:
            for n in ("mass_wet", "vol_wet", "dpmv"):
                put(n, s[SC[n]])
            put("sigma_soln", s[SC["sigma"]])
            for g in range(5):
                put("kel", v[f"ps_kel_{slot}"][g, b, k], {"gas": g + 1})
        for e in range(22):
            put("elec_l", v[f"ps_elel_{slot}"][e, b, k], {"el": e + 1})
            put("elec_s", v[f"ps_eles_{slot}"][e, b, k], {"el": e + 1})
        for ia in IONS:
            put("aer_l", v[f"ps_aerl_{slot}"][ia - 1, b, k], {"aer": ia})
            put("aer_s", v[f"ps_aers_{slot}"][ia - 1, b, k], {"aer": ia})
    names = ["", "efflorescence below crh, all solid", "upper hysteresis leg, fully liquid", "", "no salt present, all solid",
             "below the mutual deliquescence point, all solid", "handed to mesa_ptc (only the branch, water_up and dry mass asserted)"]
    desc = (f"{why[0].upper()}{why[1:]}: MESA path {path} ({names[path]}), {int(s[SC['niter']])} Kelvin pass(es), relative "
            f"humidity {s[SC['rh_pc']]:.4f}%. {src}; tolerance rel {REL:g}, integer flags exactly.")
    return {"id": tid, "description": desc, "parameter_overrides": po, "time_span": {"start": 0, "end": 1}, "assertions": a}


def main(dump_dir, drv):
    tests = []
    if dump_dir != "-":
        cache = {}
        for step, lev, b1, slot, why in COL:
            if step not in cache:
                cache[step] = load(os.path.join(dump_dir, f"esm_dump_mosaic_{step}.json"))["vars"]
            tests.append(build(cache[step], lev - 1, b1 - 1, slot, f"col_s{step}_lev{lev}_bin{b1}_{slot}",
                               f"{why}; model step {step}, level {lev}, section {b1} of 8, the {'first' if slot == 'c1' else 'last'} call of the step", SRC_COL))
    v = load(drv)["vars"]
    for i, c in enumerate(json.load(open(os.path.join(os.path.dirname(drv), "cases_ps.json")))):
        tests.append(build(v, i, 0, "c1", f"r_{c['id']}", f"hand-specified bin '{c['id']}': {c['intent']}", SRC_DRV))
    doc = json.load(open(COMP))
    doc["models"]["MosaicAerosolPhaseState"]["tests"] = tests
    with open(COMP, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    print(f"aerosol_phase_state.esm: {len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
