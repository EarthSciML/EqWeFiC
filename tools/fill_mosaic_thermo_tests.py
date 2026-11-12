#!/usr/bin/env python3
"""Fill the inline tests of three MOSAIC thermodynamics components from the chem170
compute_activities capture (dumps/chem170_mosaic_therm, WRF fork
earthsciml-instrumented-chem170; chem170_scm/run_mosaic_therm.sh):

  aerosol_water_zsr.esm      water_a        <- ah2o_a, electrolyte(:,jliquid), aer(:,jtotal)
  mtem_binary_log_gamma.esm  log_gamz       <- relhumclm (see below)
  mtem_activities.esm        mc, ma, gam, activity, gam_ratio
                                            <- water_a, aer(ions,jliquid),
                                               electrolyte(:,jliquid), log_gamz, Keq_ll

It transcribes numbers only; the components' equations are hand-written.  The
capture records compute_activities' inputs and outputs at the FIRST (_c1) and the
LAST (_c2) call per (bin, level) of each dumped step.  The thermodynamics runs in
real(kind=8), so the in-model values are the double-precision reference (rel 1e-9).

Two things the capture taught, kept here because they decide what is compared:
- log_gamz is computed in update_thermodynamic_constants from ah2o = relhumclm(k),
  but aerosol_phase_state later resets ah2o = rh_pc*0.01 with rh_pc = 100.*relhumclm
  formed in REAL*4, which differs by up to 8e-8 relative; th_ah2o records the later
  value, so the MTEM input is the dump's relhumclm (FORTRAN_BUGS N118).
- compute_activities assigns only some entries of gam and activity on each path and
  leaves the rest stale from an earlier call; only the entries it sets on that path
  are asserted (the component's description lists them).

    python3 tools/fill_mosaic_thermo_tests.py <dump_dir> [<chem170_scm/therm64/act.json>]

The optional second argument adds the hand-specified branch cases (kernels/mosaic_therm_driver,
chem170_scm/therm64/make_cases.py; their ids and intents from cases_act.json beside it).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402

HERE = os.path.join(os.path.dirname(__file__), "..", "components", "aerosol", "mosaic")
WATER = os.path.join(HERE, "aerosol_water_zsr.esm")
MTEM = os.path.join(HERE, "mtem_binary_log_gamma.esm")
ACT = os.path.join(HERE, "mtem_activities.esm")
PARAMS = os.path.join(HERE, "mosaic_thermo_parameters.esm")
REL = 1e-9
SRC = ("Reference: the em_scm_xy chem_opt = 170 column (data/eqwefic/chem170_scm, reference namelist unchanged), "
       "MOSAIC thermodynamics capture dumps/chem170_mosaic_therm (WRF fork earthsciml-instrumented-chem170, the "
       "compute_activities hooks), whose fields shared with dumps/chem170_mosaic_coag are bit-identical. The "
       "thermodynamics runs in real(kind=8), so the in-model values are a double-precision reference")
SRC_R = ("Reference: WRF's own compiled thermodynamics (main/libwrflib.a of the chem170 build, the model's flags) run "
         "on a hand-specified bin by kernels/mosaic_therm_driver (mode 1: conform_electrolytes, do_full_deliquescence, "
         "compute_activities with ah2o_a = the stated water activity), cases in chem170_scm/therm64/make_cases.py, "
         "recorded through the same compute_activities hooks as the column capture; a double-precision reference "
         "for a branch the chem170 column never reaches")
ELNAMES = ["(NH4)2SO4", "letovicite", "NH4HSO4", "NH4MSA", "NH4NO3", "NH4Cl", "Na2SO4", "Na3H(SO4)2", "NaHSO4",
           "NaMSA", "NaNO3", "NaCl", "Ca(NO3)2", "CaCl2", "Ca(MSA)2", "H2SO4", "MSA", "HNO3", "HCl", "H.HSO4",
           "CaSO4", "CaCO3"]

# (step, level, bin, call slot, regime)
CASES = [
    (1, 7, 7, "c1", "sulfate-POOR, charge-balanced (a_c = 0 exactly): B27 makes H+ = 2 sqrt(K_w), twice the neutral value"),
    (1, 13, 7, "c1", "sulfate-POOR with a round-off cation excess (a_c < 0)"),
    (1, 10, 7, "c2", "sulfate-POOR with a round-off anion excess (a_c > 0), the last call of the step"),
    (918, 47, 7, "c1", "sulfate-POOR in the upper troposphere at mid-afternoon"),
    (1, 1, 1, "c2", "sulfate-RICH, X_T < 1 (bisulfate-dominated), surface layer, smallest section, last call"),
    (1, 1, 1, "c1", "sulfate-RICH, 1 < X_T < 1.5, the same bin at its first call"),
    (1, 53, 1, "c2", "sulfate-RICH at low relative humidity (a_w < 0.3), upper troposphere"),
    (1, 14, 1, "c1", "sulfate-RICH at high relative humidity (a_w > 0.6)"),
    (689, 20, 5, "c1", "sulfate-RICH, accumulation-mode section 5, step 689"),
    (2134, 45, 3, "c2", "sulfate-RICH, free troposphere, step 2134, last call"),
    (2521, 30, 4, "c1", "sulfate-RICH at the last dumped step"),
]
# hand-specified cases whose water activity also gets a binary-MTEM matrix test
MTEM_REPLAY = {"dilute98": "very wet", "dilute995": "near saturation", "dry05": "very dry, below every a_min"}
# levels for the binary MTEM matrix: (step, level, why)
MTEM_LEVELS = [
    (918, 1, "surface layer, mid-afternoon (relative humidity 0.50)"),
    (1, 14, "a moist layer at the first step (relative humidity 0.61)"),
    (1, 53, "upper troposphere, dry (relative humidity 0.30)"),
    (2400, 30, "mid-column late on day 2"),
]


def params():
    v = json.load(open(PARAMS))["models"]["MosaicThermoParameters"]["variables"]
    return {k: np.array(v[k]["default"], dtype=float) for k in ("nu_poor", "nu_rich", "aw_min")}


def tol_for(x):
    return {"rel": REL, "abs": 0.0} if x != 0 else {"rel": 0.0, "abs": 0.0}


def main(dump_dir: str, replay: str | None = None) -> None:
    P = params()
    cache = {}

    def dump(step):
        if step not in cache:
            path = (replay if step == "R" else os.path.join(dump_dir, f"esm_dump_mosaic_{step}.json"))
            cache[step] = load(path)["vars"]
        return cache[step]

    cases = list(CASES)
    if replay:
        for i, c in enumerate(json.load(open(os.path.join(os.path.dirname(replay), "cases_act.json")))):
            cases.append(("R", i + 1, 1, "c1", c["intent"], c["id"]))
    wtests, atests, mtests = [], [], []
    for case in cases:
        step, lev, b1, slot, regime = case[:5]
        v = dump(step)
        k, b = lev - 1, b1 - 1
        assert v["th_ncall"][b, k] >= (1 if slot == "c1" else 2), (step, lev, b1, slot)
        el = v[f"th_elec_{slot}"][:, b, k]
        water = float(v[f"th_water_{slot}"][b, k])
        assert water > 0.0, case
        which = "first" if slot == "c1" else "last"
        if step == "R":
            where = f"hand-specified bin '{case[5]}', T = {v['th_t_k'][k]:.2f} K"
            src, tid = SRC_R, f"r_{case[5]}"
        else:
            where = f"model step {step}, level {lev}, section {b1} of 8, the {which} compute_activities call of the step"
            src, tid = SRC, f"s{step}_lev{lev}_bin{b1}_{slot}"
        # ---- ZSR water
        wtests.append({
            "id": f"zsr_{tid}",
            "description": (f"{regime[0].upper()}{regime[1:]}: {where}; bin water activity "
                            f"{v[f'th_ah2o_a_{slot}'][b, k]:.4f}. {src}; tolerance rel {REL:g}."),
            "parameter_overrides": {"a_w": float(v[f"th_ah2o_a_{slot}"][b, k]), "elec": el.tolist(),
                                    "aer_t": v[f"th_aert_{slot}"][:, b, k].tolist()},
            "time_span": {"start": 0, "end": 1},
            "assertions": [{"variable": "water", "time": 0.0, "expected": water, "tolerance": tol_for(water)}],
        })
        # ---- activities
        ia = v["th_iaer"] - 1
        ions = v[f"th_aerl_{slot}"][ia, b, k]
        xt = float(v[f"th_xt_{slot}"][b, k])
        poor = xt > 2 or xt < 0
        k1, _, kw = v["th_keq_ll"][:, k]
        po = {"water": water, "so4": float(ions[0]), "no3": float(ions[1]), "cl": float(ions[2]),
              "nh4": float(ions[3]), "msa": float(ions[4]), "na": float(ions[5]), "ca": float(ions[6]),
              "elec": el.tolist(), "log_gamz": v["th_log_gamz"][:, :, k].tolist(),
              "keq_hso4": float(k1), "k_w": float(kw)}
        G, A = v[f"th_gam_{slot}"][:, b, k], v[f"th_act_{slot}"][:, b, k]
        mc, ma = v[f"th_mc_{slot}"][:, b, k], v[f"th_ma_{slot}"][:, b, k]
        exact = {"rel": 0.0, "abs": 0.0}
        a = [{"variable": "poor", "time": 0.0, "expected": 1.0 if poor else 0.0, "tolerance": exact},
             {"variable": "x_t", "time": 0.0, "expected": xt, "tolerance": tol_for(xt)}]
        for name, val in zip(["mc_h", "mc_nh4", "mc_na", "mc_ca"], mc):
            a.append({"variable": name, "time": 0.0, "expected": float(val), "tolerance": tol_for(val)})
        for name, val in zip(["ma_hso4", "ma_so4", "ma_no3", "ma_cl", "ma_msa"], ma):
            a.append({"variable": name, "time": 0.0, "expected": float(val), "tolerance": tol_for(val)})
        if poor:
            xm = P["nu_poor"] * el
            mixed = [e for e in range(22) if P["nu_poor"][e] > 0 and xm[e] > 0 and e not in (17, 18)]
            setg = mixed + [17, 18] + [1, 2, 3, 7, 8, 9, 14]
            seta = setg
            ratio_set = xm[4] > 0
        else:
            setg = [e for e in range(22) if e not in (16, 20, 21)]
            seta = [e for e in range(22) if e not in (15, 16, 19, 20, 21)]
            ratio_set = True
        for e in sorted(setg):
            a.append({"variable": "gam", "coords": {"el": e + 1}, "time": 0.0, "expected": float(G[e]), "tolerance": tol_for(G[e])})
        for e in sorted(seta):
            a.append({"variable": "activity", "coords": {"el": e + 1}, "time": 0.0, "expected": float(A[e]), "tolerance": tol_for(A[e])})
        if ratio_set:
            gr = float(v[f"th_gam_ratio_{slot}"][b, k])
            a.append({"variable": "gam_ratio", "time": 0.0, "expected": gr, "tolerance": tol_for(gr)})
        skipped = sorted(set(range(22)) - set(setg))
        note = ("B27 is in force: the H+ asserted is the Fortran's 0.5 a_c + sqrt(a_c^2 + 4 K_w), not the "
                "charge-balance root. " if poor else "")
        atests.append({
            "id": f"act_{tid}",
            "description": (f"{regime[0].upper()}{regime[1:]}: {where}; X_T = {xt:.4g}. {note}{src}; tolerance rel "
                            f"{REL:g} (exact zeros pinned exactly). "
                            f"Not asserted, because compute_activities leaves them stale from an earlier call on this "
                            f"path: gam of {', '.join(ELNAMES[e] for e in skipped) or 'none'}"
                            + ("" if ratio_set else "; gam_ratio (NH4NO3 absent)") + "."),
            "parameter_overrides": po, "time_span": {"start": 0, "end": 1}, "assertions": a,
        })
    for step, lev, why in MTEM_LEVELS:
        v = dump(step)
        k = lev - 1
        rh = float(v["relhumclm"][k])
        L = v["th_log_gamz"][:, :, k]
        a = []
        for i in range(22):
            for j in range(22):
                x = float(L[i, j])
                a.append({"variable": "log_gamz", "coords": {"el": i + 1, "solv": j + 1}, "time": 0.0, "expected": x, "tolerance": tol_for(x)})
        mtests.append({
            "id": f"mtem_s{step}_lev{lev}",
            "description": (f"{why[0].upper()}{why[1:]}: model step {step}, level {lev}; relative humidity {rh:.6f} "
                            f"(relhumclm, the value update_thermodynamic_constants evaluates log_gamz at -- not the "
                            f"rh_pc*0.01 that aerosol_phase_state later stores in ah2o, N118). {SRC}; all 484 entries, "
                            f"the {int(np.count_nonzero(L))} pairs mtem_compute_log_gamz lists at rel {REL:g} (measured "
                            f"1.3e-10 worst over 944 level-steps: the polynomial cancels near its zeros) and the "
                            f"{int(np.sum(L == 0))} it leaves at zero exactly."),
            "parameter_overrides": {"a_w": rh}, "time_span": {"start": 0, "end": 1}, "assertions": a,
        })
    if replay:
        v = dump("R")
        for i, c in enumerate(json.load(open(os.path.join(os.path.dirname(replay), "cases_act.json")))):
            if c["id"] not in MTEM_REPLAY:
                continue
            L = v["th_log_gamz"][:, :, i]
            aw = float(np.float32(c["rh"]))   # the driver reads relhum as default REAL, as relhumclm is
            a = [{"variable": "log_gamz", "coords": {"el": p + 1, "solv": q + 1}, "time": 0.0,
                  "expected": float(L[p, q]), "tolerance": tol_for(float(L[p, q]))}
                 for p in range(22) for q in range(22)]
            mtests.append({
                "id": f"mtem_r_{c['id']}",
                "description": (f"Water activity {aw:.9g} ({MTEM_REPLAY[c['id']]}), outside the 0.30-0.77 the column "
                                f"reaches: the fit-range clamp of mtem_compute_log_gamz. {SRC_R}; all 484 entries at "
                                f"rel {REL:g}, zeros exactly."),
                "parameter_overrides": {"a_w": aw}, "time_span": {"start": 0, "end": 1}, "assertions": a,
            })
    for path, model, tests in ((WATER, "MosaicAerosolWaterZSR", wtests), (ACT, "MosaicMtemActivities", atests),
                               (MTEM, "MosaicMtemBinaryLogGamma", mtests)):
        doc = json.load(open(path))
        doc["models"][model]["tests"] = tests
        with open(path, "w") as f:
            json.dump(doc, f, indent=2)
            f.write("\n")
        print(f"{os.path.basename(path)}: {len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
