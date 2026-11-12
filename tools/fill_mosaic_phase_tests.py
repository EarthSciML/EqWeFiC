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

The component is evaluated at a given water activity ah2o_a (the unknown of the Kelvin equilibrium; the
Kelvin iteration itself is solver and is not written, CLAUDE.md).  Each call gives a test at the Fortran's
FINAL iterate (ah2o_a as it left it) and, for a call that took a second pass, a test at the FIRST iterate
aw0, where aw_next must be the second iterate the Fortran moved to.  What is asserted at the final iterate
depends on the MESA branch the Fortran took (ps_sc path):
  path 2 (liquid): everything the component computes, and the residual aw0/K - x (the update the Fortran's
                   1e-2 tolerance declined to apply);
  paths 1, 4, 5 (solid): state, path, water = 0, the partition, and kelvin = 1 / residual 0 when the first
                   pass was already solid (a solid reached on the SECOND pass carries the first iterate's K,
                   asserted at the first iterate instead);
  path 6 (mesa_ptc, not evaluated here): only path.

    python3 tools/fill_mosaic_phase_tests.py <dump_dir|-> <ps.json>
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402
from mosaic_cell_axes import cell_axes  # noqa: E402

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
    (1, 1, 2, "c1", "handed to mesa_ptc, which dissolves every salt (all liquid), one Kelvin pass"),
    (1, 1, 1, "c1", "handed to mesa_ptc, all liquid, two Kelvin passes"),
    (1, 1, 6, "c1", "handed to mesa_ptc, which leaves the bin mixed, one Kelvin pass"),
    (1, 43, 1, "c1", "handed to mesa_ptc, mixed, two Kelvin passes"),
    (60, 59, 2, "c1", "handed to mesa_ptc, which converges on mass (all solid) on the second Kelvin pass"),
]


def tol(x):
    return {"rel": REL, "abs": 0.0} if x != 0 else {"rel": 0.0, "abs": 0.0}


EXACT = {"rel": 0.0, "abs": 0.0}
PC_TO_FRAC = float(np.float32(0.01))   # pp.pc_to_frac, the binary32 0.01 of ah2o = rh_pc*0.01


def build(v, k, b, slot, tid, why, src):
    s = v[f"ps_sc_{slot}"][:, b, k]
    path = int(s[SC["path"]])
    path = path % 10 if path >= 6 else path   # 6 + 10*outcome + 100*iterations for mesa_ptc
    niter = int(s[SC["niter"]])
    rh = float(s[SC["rh_pc"]])
    aw0 = rh * PC_TO_FRAC
    x = float(s[SC["ah2o_a"]])
    base = {"rh_pc": rh, "t_k": float(v["ps_t_k"][k]), "num": float(s[SC["num"]]),
            "aer_t": v[f"ps_aert_in_{slot}"][:, b, k].tolist(), "elec": v[f"ps_elet_in_{slot}"][:22, b, k].tolist(),
            "water_hyst": float(s[SC["hyst_in"]])}
    if path == 6:
        # mesa_ptc's outcome, the consumer's inputs on this branch (the MESA equilibrium is solver, CLAUDE.md)
        base["mesa_state"] = float(s[SC["state"]])
        base["water_mesa"] = float(s[SC["water"]])
    assert float(v[f"ps_elet_in_{slot}"][22, b, k]) == 0.0
    a = [{"variable": "path", "time": 0.0, "expected": float(path), "tolerance": EXACT}]
    def put(var, val, coords=None, aa=a):
        d = {"variable": var, "time": 0.0, "expected": float(val), "tolerance": tol(float(val))}
        if coords:
            d["coords"] = coords
        aa.append(d)
    def put_exact(var, val, aa=a):
        aa.append({"variable": var, "time": 0.0, "expected": float(val), "tolerance": EXACT})
    put("water_up", s[SC["water_up"]])
    put("mass_dry", s[SC["mass_dry"]])
    mesa_solid = path == 6 and int(s[SC["state"]]) == 1
    solid_late = (path in (1, 4, 5) or mesa_solid) and niter > 1
    if path == 6:
        # given MESA's phase and water: the state, the water and the Kelvin factors formed from them
        put_exact("state", s[SC["state"]])
        put_exact("jhyst", s[SC["jhyst"]])
        put("water", s[SC["water"]])
        for g in range(5):
            put("kel", v[f"ps_kel_{slot}"][g, b, k], {"gas": g + 1})
        if not solid_late:
            put("kelvin", s[SC["kelvin"]])
            a.append({"variable": "kelvin_residual", "time": 0.0, "expected": aw0 / float(s[SC["kelvin"]]) - x,
                      "tolerance": {"rel": REL, "abs": 1e-14}})
        if not mesa_solid:
            for n in ("mass_wet", "vol_wet", "dpmv"):
                put(n, s[SC[n]])
            put("sigma_soln", s[SC["sigma"]])
    if path != 6:
        put_exact("state", s[SC["state"]])
        put_exact("jhyst", s[SC["jhyst"]])
        if not solid_late:
            put("kelvin", s[SC["kelvin"]])
            # the update the Fortran's tolerance declined to apply (0 exactly on a first-pass solid)
            # a difference of two O(1) numbers: an absolute floor of 1e-14 for the cancellation
            a.append({"variable": "kelvin_residual", "time": 0.0, "expected": aw0 / float(s[SC["kelvin"]]) - x,
                      "tolerance": {"rel": REL, "abs": 1e-14}})
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
             "below the mutual deliquescence point, all solid",
             "handed to mesa_ptc, whose outcome (phase and water) is given; asserted: the branch, the state, water and the Kelvin factors formed from them"]
    head = f"{why[0].upper()}{why[1:]}: MESA path {path} ({names[path]}), {niter} Kelvin pass(es), relative humidity {rh:.4f}%."
    tail = f" {src}; tolerance rel {REL:g}, integer flags exactly."
    po = dict(base, ah2o_a=x)
    out = [{"id": tid, "description": head + " Evaluated at the water activity the Fortran's Kelvin iteration stopped at"
            + (" (kelvin_residual is the update its tolerance declined to apply)." if path == 2 or (path == 6 and not solid_late) else
               "; on this solid branch K = 1 and the residual is exactly 0." if path != 6 and not solid_late else
               ("; kelvin is not asserted here: the Fortran leaves the first iterate's K in it" + (" (asserted at that iterate)." if path != 6 else ".")) if solid_late else ".")
            + tail, "parameter_overrides": po, "time_span": {"start": 0, "end": 1}, "assertions": a}]
    if niter == 2 and path != 6:
        a1 = []
        put("aw_next", x, aa=a1)
        if solid_late:
            put("kelvin", s[SC["kelvin"]], aa=a1)
        out.append({"id": tid + "_it1", "description": head + " Evaluated at the FIRST iterate aw0 = rh_pc*0.01 the "
                    "Fortran started from: aw_next = aw0/K(aw0) must be the second iterate it moved to (its final "
                    "ah2o_a)" + ("; and kelvin is the K the Fortran left in kelvin(ibin) when the second iterate "
                    "turned out solid" if solid_late else "") + "." + tail,
                    "parameter_overrides": dict(base, ah2o_a=aw0), "time_span": {"start": 0, "end": 1},
                    "assertions": a1})
    return out


# Converged states (2026-11-12): the Kelvin root x = aw0/K(x) found by WRF's own compiled
# aerosol_phase_state with its tolerance tightened from 1e-2 to 1e-13 (ESM_KEL_TOL, up to 1000 updates;
# the tolerance-study build chem170_scm/solver_tol, out_kel_tight).  The Fortran stops when the update it
# would apply is within 1e-13 of x, so at its final x |aw0/K(x) - x| <= 1e-13 x; KEL_TOL is twice that
# bound at x <= 1, the margin being the component's round-off against the Fortran's K.
KEL_TOL = 2e-13
# (solver_tol cases file, case): liquid upper-leg calls (path 2), chosen by relative humidity and by how many
# updates the model's own 1e-2 tolerance took (1 or 2) and the root needs (5-8); (88, 97) is the call
# farthest from its root (1.0e-2 in ah2o_a), (50, 57) the driest (36.5%)
KEL_SEL = [(50, 57), (83, 8), (57, 86), (88, 97), (93, 33)]


def converged_tests(sdir):
    idx = {(r["file"], r["case"]): r for r in json.load(open(os.path.join(sdir, "cases", "index.json")))}
    out, cache = [], {}
    for f, c in KEL_SEL:
        if f not in cache:
            cache[f] = (load(os.path.join(sdir, "out_kel_tight", f"all_{f:03d}.json"))["vars"],
                        load(os.path.join(sdir, "out_default", f"all_{f:03d}.json"))["vars"])
        v, v0 = cache[f]
        s, s0 = v["ps_sc_c1"][:, 0, c - 1], v0["ps_sc_c1"][:, 0, c - 1]
        assert int(s[SC["path"]]) == 2 and int(s0[SC["path"]]) == 2, (f, c)
        rh = float(s[SC["rh_pc"]]); aw0 = rh * PC_TO_FRAC; x = float(s[SC["ah2o_a"]]); K = float(s[SC["kelvin"]])
        assert abs(aw0 / K - x) <= 1e-13 * x, (f, c, aw0 / K - x)
        r = idx[(f, c)]
        base = {"rh_pc": rh, "t_k": float(v["ps_t_k"][c - 1]), "num": float(s[SC["num"]]),
                "aer_t": v["ps_aert_in_c1"][:, 0, c - 1].tolist(), "elec": v["ps_elet_in_c1"][:22, 0, c - 1].tolist(),
                "water_hyst": float(s[SC["hyst_in"]])}
        a = [{"variable": "kelvin_residual", "time": 0.0, "expected": 0.0, "tolerance": {"rel": 0.0, "abs": KEL_TOL}},
             {"variable": "path", "time": 0.0, "expected": 2.0, "tolerance": EXACT}]
        for n, val in (("aw_next", x), ("kelvin", K), ("water", s[SC["water"]])):
            a.append({"variable": n, "time": 0.0, "expected": float(val), "tolerance": tol(float(val))})
        dx = abs(float(s0[SC["ah2o_a"]]) - x) / x
        out.append({"id": f"kelvin_root_s{r['step']}_lev{r['k'] + 1}_bin{r['b'] + 1}_{'c' + str(r['slot'])}",
                    "description": (f"THE KELVIN ROOT, not a visited iterate: chem170 column step {r['step']} level {r['k'] + 1} "
                                    f"section {r['b'] + 1} (slot {r['slot']}), a liquid upper-leg bin at relative humidity "
                                    f"{rh:.4f}%, solved by WRF's own compiled aerosol_phase_state with its Kelvin tolerance "
                                    f"tightened from 1e-2 to 1e-13 (ESM_KEL_TOL; chem170_scm/solver_tol, out_kel_tight), "
                                    f"{int(s[SC['niter']])} updates (the model's own setting took {int(s0[SC['niter']])}, and its "
                                    f"ah2o_a is {dx:.2e} relative from this root). Evaluated at the root x: the residual "
                                    f"aw0/K(x) - x is zero within {KEL_TOL:g} absolute (twice the tightened stopping test "
                                    "|aw0/K(x) - x| <= 1e-13 x), and K, the water and aw_next = x at rel 1e-9. This is the "
                                    "algebraic equation the EarthSciIntegrator `solve` (newton) is given "
                                    "(couplings/mosaic_kelvin_equilibrium.esm)."),
                    "parameter_overrides": dict(base, ah2o_a=x), "time_span": {"start": 0, "end": 1}, "assertions": a})
    return out


def main(dump_dir, drv, tight=None):
    tests = []
    if dump_dir != "-":
        cache = {}
        for step, lev, b1, slot, why in COL:
            if step not in cache:
                cache[step] = load(os.path.join(dump_dir, f"esm_dump_mosaic_{step}.json"))["vars"]
            tests.extend(build(cache[step], lev - 1, b1 - 1, slot, f"col_s{step}_lev{lev}_bin{b1}_{slot}",
                               f"{why}; model step {step}, level {lev}, section {b1} of 8, the {'first' if slot == 'c1' else 'last'} call of the step", SRC_COL))
    v = load(drv)["vars"]
    for i, c in enumerate(json.load(open(os.path.join(os.path.dirname(drv), "cases_ps.json")))):
        tests.extend(build(v, i, 0, "c1", f"r_{c['id']}", f"hand-specified bin '{c['id']}': {c['intent']}", SRC_DRV))
    if tight is not None:
        tests.extend(converged_tests(tight))
    doc = json.load(open(COMP))
    # the mdrh_table_* tests belong to tools/fill_mosaic_mdrh_table_tests.py: keep them
    keep = [t for t in doc["models"]["MosaicAerosolPhaseState"]["tests"] if t["id"].startswith("mdrh_table_")]
    cell_axes(doc, "MosaicAerosolPhaseState", tests, COMP)
    doc["models"]["MosaicAerosolPhaseState"]["tests"] = tests + keep
    with open(COMP, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    print(f"aerosol_phase_state.esm: {len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions "
          f"(plus {len(keep)} mdrh_table tests kept)")


if __name__ == "__main__":
    main(*sys.argv[1:4])
