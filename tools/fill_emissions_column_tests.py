#!/usr/bin/env python3
"""Fill the inline tests of couplings/emissions_column.esm from the chem170 emissions
capture (data/eqwefic/chem170_scm/run_emis.sh) and its real64 kernel replays.

Numbers only; the coupling's equations and edges are hand-written.  Inputs are the
column state the `emissions` hook dumped (layer fields, land use, day, gsw, the
dumped emis_ant of layer 1 into the mounted anthropogenic schemes); expected values
are the WHOLE emissions_driver increment of the replay, (chem_out - chem_in)/dt, over
every layer -- not the per-stage increments the component tests use.  The `tests`
array is replaced wholesale, so CASES is the test list.

    python3 tools/fill_emissions_column_tests.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402
from fill_emissions_tests import DATA, BUILD, EMIS, p, pe  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.join(HERE, "..", "couplings", "emissions_column.esm")

# dawn, solar noon, late afternoon, local midnight and the last dumped step (RADM2
# variety); dawn, noon and midnight for the CBM-Z variety
CASES = [("default", 1), ("default", 361), ("default", 721), ("default", 1081), ("default", 2879),
         ("pnnlcm", 1), ("pnnlcm", 361), ("pnnlcm", 1081)]
GAS_FIELDS = ["so2", "no", "no2", "co", "eth", "hc3", "hc5", "hc8", "xyl", "ol2", "olt", "oli", "tol", "csl",
              "hcho", "ald", "ket", "ora2", "nh3", "ch3oh", "c2h5oh"]
AER_FIELDS = ["so4i", "so4j", "so4c", "no3i", "no3j", "no3c", "orgi", "orgj", "orgc", "eci", "ecj", "ecc",
              "pm25i", "pm25j", "pm_10"]
AER = {"so4": "lptr_so4_aer", "no3": "lptr_no3_aer", "oc": "lptr_oc_aer", "bc": "lptr_bc_aer",
       "oin": "lptr_oin_aer", "num": "numptr_aer"}


def const(values):
    return {"op": "const", "args": [], "value": values}


def linf(var, arr, axes):
    a = np.asarray(arr, dtype=float)
    idx = ["k", "b"][: a.ndim]
    ref = {"op": "faq", "output_idx": idx, "args": [], "ranges": {i: {"from": ax} for i, ax in zip(idx, axes)},
           "expr": {"op": "index", "args": [const(a.tolist())] + idx}}
    scale = float(np.max(np.abs(a)))
    return {"variable": var, "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
            "tolerance": {"abs": 1e-9 * scale, "rel": 0.0}, "reference": ref}


def main():
    doc = json.load(open(DOC))
    model = doc["models"]["EmissionsColumn"]
    gas = [v[1:-3] for v in model["variables"] if v.startswith("d") and v.endswith("_dt") and v[1:-3] not in AER]
    tests = []
    for mode, step in CASES:
        d = load(f"{DATA}/dumps/chem170_emis_{mode}/esm_dump_emissions_{step}.json")["vars"]
        r = load(f"{DATA}/dumps/chem170_emis_{mode}/replay/emissions_{step}.json")["vars"]
        dt = float(d["dt"])
        ci = np.asarray(d["chem_in"], dtype=float)
        inc = (np.asarray(r["chem_out"], dtype=float) - ci) / dt  # (lev, species)
        winc = (np.asarray(d["chem_out"], dtype=float) - ci) / dt
        ea = np.asarray(d["emis_ant"], dtype=float)
        col = lambda n: [float(v) for v in np.asarray(d[n], dtype=float)]
        over = {"T": col("t_phy"), "p": col("p_phy"), "rho": col("rho"), "alt": col("alt"), "ze": col("z_at_w"),
                "gsw": float(d["gsw"]), "lu": int(d["ivgtyp"]), "julday": int(d["julday"]),
                "emiss_inpt_opt": int(d["emiss_inpt_opt"]), "volumcen": col("volumcen_sect")}
        for f in GAS_FIELDS:
            over[f"Anthro.E_{f}"] = float(ea[0, pe(f"e_{f}") - 1])
        for f in AER_FIELDS:
            over[f"Aer.E_{f}"] = float(ea[0, pe(f"e_{f}") - 1])
        a, covered = [], set()
        for s in gas:
            a.append(linf(f"d{s}_dt", inc[:, p(s) - 1], ["lev"]))
            covered.add(p(s) - 1)
        for c, ptr in AER.items():
            cols = np.asarray(d[ptr], dtype=int) - 1
            a.append(linf(f"d{c}_dt", inc[:, cols], ["lev", "bin"]))
            covered |= set(int(i) for i in cols)
        touched = set(np.nonzero(np.any(inc != 0, axis=0))[0]) | set(np.nonzero(np.any(winc != 0, axis=0))[0])
        stray = touched - covered
        if stray:
            raise SystemExit(f"{mode} {step}: species {sorted(stray)} change but are not exposed")
        nz = np.abs(inc) > 0
        dev = float(np.max(np.abs(winc[nz] - inc[nz]) / np.abs(inc[nz]))) if nz.any() else 0.0
        tests.append({
            "id": f"chem170_emis_{mode}_step{step}",
            "description": (
                f"Reference: the chem_opt = 170 emissions capture (data/eqwefic/chem170_scm/run_emis.sh, MODE={mode}: "
                f"emiss_opt = 4, emiss_inpt_opt = {int(d['emiss_inpt_opt'])}, bio_emiss_opt = 1; {BUILD}) at model "
                f"step {step} (dt = {dt:g} s), column (its, jts); USGS category {over['lu']}, day {over['julday']}, "
                f"gsw {over['gsw']:.6g} W/m2. Inputs are the dumped column state and the dumped emission fields of "
                f"layer 1. Expected values are the WHOLE emissions_driver increment over the step divided by dt, "
                f"(chem_out - chem_in)/dt at every layer, from the real64 replay of bio_emissions, "
                f"cbmz_addemiss_anthro, add_biogenics and mosaic_addemiss compiled unchanged and called in "
                f"emissions_driver's order (kernels/emissions_driver), at 1e-9 of each column's maximum; so e_bio "
                f"reaches add_biogenics through the coupling, not as an input. {int(np.count_nonzero(np.any(inc != 0, axis=0)))} "
                f"of the 233 chem species change, all in layer 1, and every one of them is asserted here (the fill "
                f"script refuses a case where a species changes that the document does not expose). The in-model "
                f"real32 increments differ from the replay by at most {dev:.2g} relative."),
            "parameter_overrides": over, "time_span": {"start": 0, "end": 0}, "assertions": a})
    model["tests"] = tests
    with open(DOC, "w") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"{len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


if __name__ == "__main__":
    main()
