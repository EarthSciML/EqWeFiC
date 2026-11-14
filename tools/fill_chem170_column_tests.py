#!/usr/bin/env python3
"""Fill the inline tests of couplings/chem170_column.esm, the joint CBM-Z + MOSAIC column.

Numbers only; the document's equations and edges are hand-written. One test per model step s
of data/eqwefic/chem170_scm/run_full_mosaic.sh (the run_full.sh configuration -- emissions,
Gunther, wetdep_ls, Wesely -- with the MOSAIC hooks and the ASTEM tape on), at ONE gas +
aerosol state: the state ASTEM's first semi-volatile sub-step forms its fluxes at (see
data/eqwefic/chem170_full/replay_joint.py, which builds it and replays every process at it).
The three parts reuse the fill scripts of the columns this one joins:

  CBMZColumn   tools/fill_cbmz_column_tests.build_test at the state X (CBM-Z's gases, molec/cm3),
               chemistry from vdot_rconst_replay (RCONST recomputed at X)
  GasParticle  tools/fill_mosaic_column_tests.build on the in-model ASTEM tape and dump
  AerChem      tools/fill_mosaic_aerchem_column_tests.build with rsub_therm := the state's rsub

and the join adds the exchange: ASTEM's gas tendencies (GasParticle's Rates.dgas, expected
from the tape) and nucleation's H2SO4/NH3 uptake (the nucleation replay) handed back to
CBM-Z's H2SO4, HNO3, HCl and NH3 by factgas*conv/conv1a, and ASTEM's condensation added to the
aerosol components. ASTEM's rates have a reference only on the levels whose eight fluxes were all
formed at the state (fill_mosaic_column_tests' `clean`); on every other level the four
exchanging gases' column tendencies, their exchange terms and the summed aerosol tendency are
not asserted (the other 62 gases are, everywhere).

    python3 tools/fill_chem170_column_tests.py
"""
import collections
import copy
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from esm_dump import load  # noqa: E402
from esm_compact_json import dumps  # noqa: E402

DATA = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic"
DUMPS = f"{DATA}/dumps/chem170_full_mosaic"
REPLAY = f"{DATA}/chem170_full/replay_mosaic"
TAPE = f"{DATA}/chem170_full/mosaic/astem.tape.gz"
DOC = os.path.join(HERE, "..", "couplings", "chem170_column.esm")
sys.path.insert(0, f"{DATA}/chem170_scm/mosaic_astem")
from astem_tape import read_tape  # noqa: E402


def mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


fcz = mod("fcz", f"{HERE}/fill_cbmz_column_tests.py")
fmc = mod("fmc", f"{HERE}/fill_mosaic_column_tests.py")
fac = mod("fac", f"{HERE}/fill_mosaic_aerchem_column_tests.py")

STEPS = [(364, "local noon of the first day: the strongest photochemistry and emissions; nucleation at the levels the replay finds"),
         (790, "an evening shower: WSM6 rain, so wetdep_ls removes sulf at every rainy layer")]
REL = 1e-9
FACTGAS, CM3_PER_M3, NMOL_PER_MOL, FACTDENS = 1.0e6, 1.0e6, 1.0e9, 28966.0
EXCH = {"H2SO4": 1, "HNO3": 2, "HCl": 3, "NH3": 4}      # ASTEM gas axis
COMP_OF_GAS = {1: 1, 2: 2, 3: 3, 4: 5}                  # gas -> rsub component (so4 no3 cl nh4)


def renamed(a, old, new):
    a = copy.deepcopy(a)
    if a["variable"].startswith(old + "."):
        a["variable"] = new + a["variable"][len(old):]
    return a


def build(step, what, calls_by, cz_states):
    st = json.load(open(f"{REPLAY}/state_{step}.json"))
    X = np.asarray(st["X"], dtype=float)
    conv = np.asarray(st["conv"], dtype=float)
    nlev = X.shape[1]
    # --- CBM-Z part at X
    tc, expected, dep_term = fcz.build_test(step, what, False, cz_states, dumps=DUMPS, replay=REPLAY,
                                             X=X, vdot_key="vdot_rconst_replay")
    # --- MOSAIC gas-particle part: the in-model tape at the state
    v = load(f"{DUMPS}/esm_dump_mosaic_{step}.json")["vars"]
    po_m, asr_m, census = fmc.build(step - 1, calls_by[step - 1], v, skip_unconverged_state=True)
    # --- MOSAIC aerosol chemistry at the state's rsub
    m = dict(v)
    m["rsub_therm"] = np.load(f"{REPLAY}/rsub_S_{step}.npy")
    dd = load(f"{DUMPS}/esm_dump_mosaic_drydep_{step - 2}.json")["vars"]
    vx = load(f"{DUMPS}/esm_dump_chem_vertmx_{step - 2}.json")["vars"]
    assert int(dd["ktau"]) == int(vx["ktau"]) == int(v["ktau"]) == step
    nuc = load(f"{REPLAY}/nuc_{step}.json")["vars"]
    ta = fac.build(step, m=m, dd=dd, vx=vx, nuc=nuc, coag=load(f"{REPLAY}/coag_{step}.json")["vars"],
                   dep=load(f"{REPLAY}/dd_{step}.json")["vars"], desc_in="")

    # --- the join
    rho = np.asarray(v["rho_phy"], dtype=float)
    conv1a = rho / FACTDENS * CM3_PER_M3 * NMOL_PER_MOL
    gas = np.asarray(st["tape_gas"], dtype=float)                      # [lev, 5]
    rs = m["rsub_therm"]
    rnames = [str(v["name_%04d" % i]) for i in range(1, int(v["ltot2"]) + 1)]
    q_h2so4, q_nh3 = rs[rnames.index("h2so4")], rs[rnames.index("nh3")]
    rnuc = np.asarray(nuc["rate_newnuc_replay"], dtype=float)
    dt = float(v["dtchem"])
    dgas, daer = {}, {}
    for a in asr_m:
        c = a.get("coords", {})
        if a["variable"] == "Rates.dgas":
            dgas[(c["lev"], c["gas"])] = (a["expected"], a.get("tolerance", {}).get("abs", 0.0))
        elif a["variable"] == "Rates.daer":
            daer[(c["lev"], c["bin"], c["gas"])] = (a["expected"], a.get("tolerance", {}).get("abs", 0.0))
    clean = sorted({k for (k, g) in dgas})

    asr = []

    def put(var, val, rel, absb, **coords):
        a = {"variable": var, "time": 0.0, "expected": float(val), "tolerance": {"rel": rel, "abs": float(absb)}}
        if coords:
            a["coords"] = coords
        asr.append(a)

    for k in range(1, nlev + 1):
        for g in range(1, 6):
            put("Join.gas_nmol", gas[k - 1, g - 1], REL, 0.0, lev=k, gas=g)
        put("Join.q_h2so4", q_h2so4[k - 1], REL, 0.0, lev=k)
        put("Join.q_nh3", q_nh3[k - 1], REL, 0.0, lev=k)
    exch = {}
    for n, g in EXCH.items():
        exch[n] = {}
        for k in clean:
            dg, fl = dgas[(k, g)]
            nu = (rnuc[0, k - 1] if n == "H2SO4" else rnuc[1, k - 1] if n == "NH3" else 0.0) / dt
            val = (dg / conv1a[k - 1] + nu) * FACTGAS * conv[k - 1]
            scale = (fl / conv1a[k - 1] + 1e-9 * abs(nu)) * FACTGAS * conv[k - 1]
            exch[n][k] = (val, scale)
            put(f"Join.aerosol_{n}", val, REL, scale, lev=k)
    # aerosol: the aerosol-chemistry column's dr_dt plus ASTEM's condensation, on the clean levels
    for a in ta["assertions"]:
        if a["variable"] != "Column.dr_dt":
            continue
        c = a["coords"]
        k, b, comp = c["lev"], c["bin"], c["comp"]
        if k not in clean:
            continue
        add, addfl = 0.0, 0.0
        for g, cc in COMP_OF_GAS.items():
            if cc == comp:
                x, fl = daer[(k, b, g)]
                add, addfl = x / conv1a[k - 1], fl / conv1a[k - 1]
        put("Join.dr_dt", a["expected"] + add, REL, a["tolerance"]["abs"] + addfl, lev=k, bin=b, comp=comp)

    # CBM-Z assertions: the four exchanging gases carry the exchange, and only where it is known
    # a radical in steady state (O3P, OH) has a net tendency that is the roundoff-level difference of
    # its gross production and loss; at the state the column and the replay agree to 1-3 ulps of the
    # gross loss L X (replay_joint.gross_loss), so every species' bound has a floor of 1e-14 L X
    gross = json.load(open(f"{REPLAY}/gross_{step}.json"))
    cz_asr = []
    for a in tc["assertions"]:
        n, lev = a["variable"][len("CBMZColumn.d"):].split("_dt[")
        lev = int(lev.rstrip("]"))
        if n in gross:
            a = copy.deepcopy(a)
            a["tolerance"]["abs"] = max(a["tolerance"]["abs"], 1e-14 * gross[n][lev - 1])
        if n in EXCH:
            if lev not in clean:
                continue
            val, scale = exch[n][lev]
            a = copy.deepcopy(a)
            a["expected"] = a["expected"] + val
            a["tolerance"] = {"rel": REL, "abs": a["tolerance"]["abs"] + scale}
        cz_asr.append(a)

    po = dict(tc["parameter_overrides"])
    for key, val in po_m.items():
        if key != "Column.gas":                      # the join supplies ASTEM's gases
            po["GasParticle" + key[len("Column"):]] = val
    for key, val in ta["parameter_overrides"].items():
        if key not in ("Column.q_h2so4", "Column.q_nh3"):  # the join supplies nucleation's gases
            po["AerChem" + key[len("Column"):]] = val
    all_asr = (cz_asr + [renamed(a, "Column", "GasParticle") for a in asr_m]
               + [renamed(a, "Column", "AerChem") for a in ta["assertions"]] + asr)
    nnuc = int(np.sum(rnuc[4] != 0))
    desc = (f"The joint chem_opt = 170 column at model step {step} ({fcz.lst(step)} LST), {what}. Reference: "
            "data/eqwefic/chem170_scm/run_full_mosaic.sh (the run_full.sh configuration -- emiss_opt = 4, "
            "emiss_inpt_opt = 101, bio_emiss_opt = 1, wetscav_onoff = -1, cu_physics = 0 -- with the MOSAIC hooks and "
            "the ASTEM tape; fork ctessum-claude/WRF branch earthsciml-instrumented-chem170, binary chem170_full/bin). "
            "THE STATE is the one ASTEM's first semi-volatile sub-step forms its fluxes at: per section MOSAIC's last "
            "aerosol_phase_state call of the step (ps_*_c2), per level the condensing gases of the tape's first "
            "sub-step (SVT), every other gas KPP's output of the step; CBM-Z's four exchanging gases are put in molec/cm3 "
            "by inverting the document's own down chain, so ASTEM is handed the tape's gases. Inputs the Fortran's "
            "solvers decided are taken as the Fortran left them, as in couplings/mosaic_column.esm: ah2o_a (Kelvin), "
            "the MESA solid on path 6, dp_wet. Expected values, all at that one state, real64: CBM-Z's Vdot with RCONST "
            "recomputed at the state (kernels/cbmz_driver), conv times the emissions_driver increment over dt, vertmx "
            "with K = 0 (Romberg dt -> 0) and wetdep_ls (sulf) replayed at the state; ASTEM's phase state, "
            "thermodynamics and rates from the in-model tape (fill_mosaic_column_tests.py's selection); nucleation, "
            "coagulation and aerosol dry deposition replayed at the state's rsub (chem170_full/replay_joint.py: "
            f"nucleation fires at {nnuc} of {nlev} levels); and the join -- ASTEM's dgas and nucleation's uptake handed "
            "back to H2SO4, HNO3, HCl and NH3, ASTEM's condensation added to the components -- on the "
            f"{len(clean)} levels whose eight ASTEM fluxes were formed at the state (levels {clean}); there the four "
            "exchanging gases' column tendencies are asserted with the exchange in them, elsewhere not at all. "
            f"Census of the gas-particle part: {dict(sorted((str(a), b) for a, b in census.items()))}. Tolerances: "
            "the three parts' own (rel 1e-10 and the deposition bound for CBM-Z, rel 1e-9 and the parts' floors for "
            "MOSAIC), rel 1e-9 plus the parts' floors on the join. Filled by tools/fill_chem170_column_tests.py.")
    return {"id": f"chem170_joint_step{step}", "description": desc, "parameter_overrides": po,
            "initial_conditions": tc["initial_conditions"], "time_span": {"start": 0, "end": 0},
            "tolerance": fmc.TEST_TOL, "assertions": all_asr}


def main():
    doc = json.load(open(DOC))
    cz = json.load(open(fcz.DOC))["models"]["CBMZColumn"]
    cz_states = [n for n in cz["variables"] if cz["variables"][n]["type"] == "unknown"
                 and cz["variables"][n].get("units") == "molec/cm^3" and n not in ("n_air",)]
    calls_by = collections.defaultdict(dict)
    for c in read_tape(TAPE):
        calls_by[c["step"]][c["lev"]] = c
    tests = [build(s, w, calls_by, cz_states) for s, w in STEPS]
    doc["models"]["Join"]["tests"] = tests
    open(DOC, "w").write(dumps(doc) + "\n")
    print(f"{len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main()
