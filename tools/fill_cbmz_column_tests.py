#!/usr/bin/env python3
"""Fill the inline tests of couplings/cbmz_column.esm from the chem170_full capture.

Numbers only; the coupling's equations and edges are hand-written.  Each test is one
model step S of data/eqwefic/chem170_scm/run_full.sh (emissions, Gunther, wetdep_ls and
Wesely all on), at ONE column state: CBM-Z's state at KPP entry.  Inputs are that state
(the column species' initial values, molec/cm3) and the dumped column/surface fields each
WRF process reads at step S -- the cbmz_kpp, emissions and wetdep_ls hooks at their call S,
the wesely and chem_vertmx hooks at their call S-2 (they are not called at the first two
steps, so their call N is ktau N+2; checked on the dumped ktau).  Expected values are the
SUM of every process's instantaneous tendency at that state, each from a real64 kernel
replay (data/eqwefic/chem170_full/replay_column.py):

    d<X>/dt = Vdot_X (cbmz_mosaic_Fun at var_in)
              + conv * [ emissions_driver increment / dt            (emitted species)
                         + vertmx with K = 0, Romberg dt -> 0       (every chem gas)
                         + wetdep_ls increment / dt at the state ]  (sulf)

with conv = 1e-6 dens2con_a rho_phy, the KPP interface's own.  Every species and every
layer is asserted.  `--mixing` fills the vertmx variant instead (the vertmix branch): the
replay with the dumped ekmfull, mixing and deposition together, with WRF's discard of the
top layer applied (dry_dep_driver.F writes vertmx's result back for kts..kte-1 only).

    python3 tools/fill_cbmz_column_tests.py [--mixing] [--doc PATH]
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from esm_dump import load  # noqa: E402
from fill_emissions_tests import pe  # noqa: E402

DATA = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic"
DUMPS = f"{DATA}/dumps/chem170_full"
REPLAY = f"{DATA}/chem170_full/replay"
sys.path.insert(0, f"{DATA}/chem170_full")
from replay_column import KPP2CHEM, DENS2CON_A, kpp_indices, chem_names  # noqa: E402

DOC = os.path.join(HERE, "..", "couplings", "cbmz_column.esm")
STEPS = [(364, "local noon of the first day, the strongest photochemistry and emissions of the run"),
         (790, "an evening shower: WSM6 rain, so wetdep_ls removes sulf at every rainy layer"),
         (810, "the same shower twenty minutes on, the rain easing"),
         (1080, "local midnight: every photolysis frequency zero and Gunther dark"),
         (1621, "the second morning, an aged state after a day of emissions and chemistry")]
PHOTO = {"j_o31d": 1, "j_o33p": 2, "j_no2": 3, "j_no3o": 5, "j_hno2": 6, "j_hno3": 7, "j_hno4": 8,
         "j_h2o2": 9, "j_ch2or": 10, "j_ch2om": 11}
GAS_FIELDS = ["so2", "no", "no2", "co", "eth", "hc3", "hc5", "hc8", "xyl", "ol2", "olt", "oli", "tol",
              "csl", "hcho", "ald", "ket", "ora2", "nh3", "ch3oh", "c2h5oh"]
AER_FIELDS = ["so4i", "so4j", "so4c", "no3i", "no3j", "no3c", "orgi", "orgj", "orgc", "eci", "ecj",
              "ecc", "pm25i", "pm25j", "pm_10"]
KPP_ONLY = ["O1D", "O3P"]
RTOL = 1e-10
# The deposition term is the one piece whose reference is single precision: the replayed
# vertmx reads the real32 ddvel WRF formed, while the column forms the velocity through the
# Wesely components in double precision -- the same comparison those components' own tests
# make at rel 1e-5. So each assertion's absolute bound is rel 1e-5 of its deposition term
# (and never below 1e-12 of the species' column maximum).
DEP_RTOL = 1e-5


def fl(a):
    return [float(x) for x in np.asarray(a, dtype=float).reshape(-1)]


def lst(step):
    h = (6.0 + step / 60.0) % 24.0      # the run starts at 12 UTC = 06 LST (96.7 W, solar)
    return f"{int(h):02d}:{int(round((h % 1) * 60)) % 60:02d}"


def main(mixing, doc_path):
    doc = json.load(open(doc_path))
    model = doc["models"]["CBMZColumn"]
    names = chem_names()
    ind = kpp_indices()
    states = [n for n in model["variables"] if model["variables"][n]["type"] == "unknown"
              and model["variables"][n].get("units") == "molec/cm^3" and n not in ("n_air",)]
    tests, worst = [], {}
    for step, what in STEPS:
        k = load(f"{DUMPS}/esm_dump_cbmz_kpp_{step}.json")["vars"]
        e = load(f"{DUMPS}/esm_dump_emissions_{step}.json")["vars"]
        wl = load(f"{DUMPS}/esm_dump_wetdep_ls_{step}.json")["vars"]
        we = load(f"{DUMPS}/esm_dump_wesely_{step - 2}.json")["vars"]
        vx = load(f"{DUMPS}/esm_dump_chem_vertmx_{step - 2}.json")["vars"]
        assert int(e["ktau"]) == step and int(wl["ktau"]) == step, step
        assert int(we["ktau"]) == step and int(vx["ktau"]) == step, step
        assert np.array_equal(np.asarray(k["temp"]), np.asarray(e["t_phy"], dtype=float)), step
        rho = np.asarray(k["rho_phy"], dtype=float)
        # one column geometry and one density serve every process at a given model step
        assert np.array_equal(np.asarray(e["z_at_w"]), np.asarray(vx["zzfull"])), step
        assert np.array_equal(np.asarray(vx["dryrho"], dtype=float), rho), step
        conv = 1e-6 * DENS2CON_A * rho
        var = np.asarray(k["var_in"], dtype=float)
        vdot = np.asarray(load(f"{REPLAY}/kp_{step}.out.json")["vars"]["vdot"], dtype=float)
        emr = load(f"{DUMPS}/replay/emissions_{step}.json")["vars"]
        dt = float(e["dt"])
        emis = (np.asarray(emr["d_anthro_gas"], dtype=float) + np.asarray(emr["d_bio"], dtype=float)) / dt
        wlo = load(f"{REPLAY}/wl_{step}.json")["vars"]
        st = json.load(open(f"{REPLAY}/state_{step}.json"))
        wet = (np.asarray(wlo["chem_out"], dtype=float)[:, names.index("sulf") + 1]
               - np.asarray(st["S"]["sulf"])) / float(wl["dt"])
        rom = json.load(open(f"{REPLAY}/rom_{step}_{'k' if mixing else 'k0'}.json"))
        dd = dict(zip(rom["species"], np.asarray(rom["tend"])))
        nlev = var.shape[1]

        expected, dep_term = {}, {}
        for n in states:
            c = KPP2CHEM.get(n)
            tot = np.zeros(nlev) if n in ("HCl", "NH3") else vdot[ind[n] - 1, :].copy()
            dep_term[n] = np.zeros(nlev)
            if c is not None:
                dep_term[n] = conv * dd[c]
                ppmv = emis[:, names.index(c) + 1] + dd[c]
                if mixing:
                    ppmv[-1] -= dd[c][-1]     # WRF discards vertmx's top layer
                if c == "sulf":
                    ppmv = ppmv + wet
                tot = tot + conv * ppmv
            expected[n] = tot

        ea = np.asarray(e["emis_ant"], dtype=float)
        po = {"CBMZColumn.T": fl(k["temp"]), "CBMZColumn.M_air": fl(k["c_m"]),
              "CBMZColumn.H2O": fl(np.asarray(k["fix"])[0]), "CBMZColumn.M_fix": fl(np.asarray(k["fix"])[1]),
              "CBMZColumn.rho": fl(rho), "CBMZColumn.p": fl(e["p_phy"]), "CBMZColumn.alt": fl(e["alt"]),
              "CBMZColumn.ze": fl(e["z_at_w"]),
              "CBMZColumn.qc": fl(wl["qc"]), "CBMZColumn.w": fl(wl["vvel"]),
              "CBMZColumn.rain": float(wl["rain"]), "CBMZColumn.dt": float(wl["dt"]),
              "CBMZColumn.gsw": float(e["gsw"]), "CBMZColumn.lu": float(e["ivgtyp"]),
              "CBMZColumn.julday": float(e["julday"]), "CBMZColumn.emiss_inpt_opt": float(e["emiss_inpt_opt"]),
              "CBMZColumn.volumcen": fl(e["volumcen_sect"]),
              "CBMZColumn.T_sfc": float(we["ta"]), "CBMZColumn.rh": float(we["rhchem"]),
              "CBMZColumn.wetflag": float(bool(we["wetflag"])), "CBMZColumn.rainflag": float(bool(we["rainflag"])),
              "CBMZColumn.iland": float(we["iland"]), "CBMZColumn.iseason": float(we["iseason"]),
              "CBMZColumn.ustar": float(we["ustar"]), "CBMZColumn.rmol": float(we["rmol"]),
              "CBMZColumn.z0": float(we["zntt"])}
        if mixing:
            po["CBMZColumn.exch_h"] = fl(vx["ekmfull"])
        for j, i in PHOTO.items():
            po[f"CBMZColumn.{j}"] = fl(np.asarray(k["jv"])[i - 1])
        for f in GAS_FIELDS:
            po[f"Anthro.E_{f}"] = float(ea[0, pe(f"e_{f}") - 1])
        for f in AER_FIELDS:
            po[f"Aer.E_{f}"] = float(ea[0, pe(f"e_{f}") - 1])
        ic = {}
        for n in states:
            ic[f"CBMZColumn.{n}"] = fl(var[ind[n] - 1])
        # the flag the column forms from its own state must be the one wesely_driver formed
        hi = var[ind["NH3"] - 1, 0] > 2 * var[ind["SO2"] - 1, 0]
        assert hi == bool(we["highnh3"]), (step, hi, we["highnh3"])

        a = []
        for n in states:
            ex = expected[n]
            floor = 1e-12 * float(np.max(np.abs(ex)))
            for lev in range(nlev):
                a.append({"variable": f"CBMZColumn.d{n}_dt[{lev + 1}]", "time": 0.0,
                          "expected": float(ex[lev]),
                          "tolerance": {"rel": RTOL,
                                        "abs": max(floor, DEP_RTOL * abs(float(dep_term[n][lev])))}})
        rain = float(wl["rain"])
        tests.append({
            "id": f"chem170_full_step{step}" + ("_mixing" if mixing else ""),
            "description": (
                f"Reference: the chem_opt = 170 full-column capture (data/eqwefic/chem170_scm/run_full.sh: the reference "
                f"namelist with emiss_opt = 4, emiss_inpt_opt = 101, bio_emiss_opt = 1, wetscav_onoff = -1 and "
                f"cu_physics = 0; fork ctessum-claude/WRF branch earthsciml-instrumented-chem170) at model step "
                f"{step} (dt = {dt:g} s, {lst(step)} LST), column (its, jts): {what}. Grid-scale rain over the step "
                f"{rain:.3g} kg/m2; j(NO2) at the surface {float(np.asarray(k['jv'])[2, 0]):.3g} /s; Wesely "
                f"iland {int(we['iland'])}, iseason {int(we['iseason'])}, rh {float(we['rhchem']):.1f} %"
                f"{', rain flag' if we['rainflag'] else ''}{', dew flag' if we['wetflag'] else ''}"
                f"{', high NH3' if we['highnh3'] else ''}. The column state is CBM-Z's at KPP entry (var_in), and "
                f"every expected value is the SUM, at that one state, of each process's instantaneous tendency from "
                f"a real64 replay of the Fortran: Vdot of cbmz_mosaic_Fun (kernels/cbmz_driver); conv times the "
                f"emissions_driver increment over dt (kernels/emissions_driver; independent of the chemical state); "
                + ("conv times vertmx with the dumped ekmfull and ddvel, mixing and deposition together, extrapolated "
                   "to dt -> 0 (Romberg on dt = 0.2 ... 0.025 s), with vertmx's top layer discarded as "
                   "dry_dep_driver discards it; "
                   if mixing else
                   "conv times vertmx with K = 0 -- the deposition flux through the ground alone -- with the "
                   "ddvel dry_dep_driver hands it, extrapolated to dt -> 0 (Romberg on dt = 20 ... 2.5 s); ")
                + f"and, for sulf, conv times wetdep_ls's increment over dt at this state "
                f"(kernels/wetdep_ls_driver with the dumped call's gases replaced by it). conv = 1e-6 dens2con_a "
                f"rho_phy is the KPP interface's own. Every species and every layer is asserted, rel {RTOL:g}, with an absolute bound of {DEP_RTOL:g} of the "
                f"layer's deposition term: the replayed vertmx reads WRF's real32 velocity, the column forms it through "
                f"the Wesely components in double precision (their own tests compare at rel 1e-5)."),
            "parameter_overrides": po, "initial_conditions": ic,
            "time_span": {"start": 0, "end": 0}, "assertions": a})
    model["tests"] = tests
    with open(doc_path, "w") as f:
        f.write(dumps(doc) + "\n")
    print(f"{len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


from esm_compact_json import dumps  # noqa: E402

if __name__ == "__main__":
    args = sys.argv[1:]
    mixing = "--mixing" in args
    doc = DOC
    if "--doc" in args:
        doc = args[args.index("--doc") + 1]
    main(mixing, doc)
