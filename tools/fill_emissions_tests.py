#!/usr/bin/env python3
"""Fill the inline tests of components/emissions/wrf_chem/*.esm from the chem170
emissions capture (data/eqwefic/chem170_scm/run_emis.sh) and its real64 kernel
replays (wrf-chem170/kernels/emissions_driver, via chem170_scm/emis_flat.py and
chem170_scm/emis_synth.py).

It transcribes numbers only; the components' equations are hand-written.  The
inputs of each case are the in-model state the `emissions` hook dumped (or, for
the synthetic Gunther cases, that state with the biogenic inputs rewritten);
the expected values are the replay's per-stage increments divided by dt:
`d_anthro_gas` (cbmz_addemiss_anthro), `d_bio` (add_biogenics), `d_aer`
(mosaic_addemiss) and the Gunther `e_bio`.  Each component's `tests` array is
replaced wholesale, so the case lists below ARE the test lists.

    python3 tools/fill_emissions_tests.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = os.path.join(HERE, "..", "components", "emissions", "wrf_chem")
DATA = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic"
REGISTRY = f"{DATA}/wrf-chem170/Registry/registry.chem"
PFS = 2  # PARAM_FIRST_SCALAR

# (input variety, model step), dt = 60 s from 1997-06-20 12Z (about 05:30 local solar
# time at -96.7): the first step, the first hourly reread (13Z), local midnight (06Z)
# and the last dumped step, for each variety
ANTHRO = [("default", 1), ("default", 61), ("default", 1081), ("default", 2879),
          ("pnnlcm", 1), ("pnnlcm", 61), ("pnnlcm", 1081), ("pnnlcm", 2879)]
# in-model Gunther steps (cropland, June): dawn, solar noon, late afternoon, midnight,
# the second solar noon
BIO_STEPS = [1, 361, 721, 1081, 1801]
BIO_SYNTH = ["lu01_urban", "lu05_cropgrass", "lu06_cropwood", "lu07_grass", "lu08_shrub",
             "lu11_decid", "lu11_decid_night", "lu12_decneedle", "lu14_evneedle", "lu15_mixed",
             "lu15_mixed_winter", "lu02_crop_winter", "lu16_water", "lu24_ice"]

BUILD = ("GNU gfortran real32 in model; fork ctessum-claude/WRF branch earthsciml-instrumented-chem170, "
         "the `emissions` hook around emissions_driver in chem/chem_driver.F")


def package(cond: str, array: str) -> list[str]:
    for line in open(REGISTRY):
        if line.startswith("package") and f" {cond} " in f" {line} ":
            for part in line.split()[-1].split(";"):
                name, _, members = part.partition(":")
                if name == array:
                    return members.split(",")
    raise SystemExit(f"no package {cond}")


CHEM = package("chem_opt==170", "chem")
EMIS = package("emiss_opt==4", "emis_ant")


def p(name: str) -> int:
    return CHEM.index(name) + PFS


def pe(name: str) -> int:
    """WRF's p_e_<name>; 1 (the dummy slot WRF leaves every absent pointer at) when the
    emiss_opt = 4 package does not carry the field -- E_PM_10 is such a field."""
    return EMIS.index(name) + PFS if name in EMIS else 1


def ref(values, axis: str) -> dict:
    return {"op": "faq", "output_idx": ["k"], "args": [], "ranges": {"k": {"from": axis}},
            "expr": {"op": "index", "args": [{"op": "const", "args": [], "value": [float(v) for v in values]}, "k"]}}


def scalar(var: str, value: float) -> dict:
    if value == 0.0:
        return {"variable": var, "time": 0.0, "expected": 0.0, "tolerance": {"abs": 0.0, "rel": 0.0}}
    return {"variable": var, "time": 0.0, "expected": float(value), "tolerance": {"rel": 1e-9}}


def vector(var: str, values, axis: str) -> dict:
    scale = float(np.max(np.abs(values)))
    return {"variable": var, "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
            "tolerance": {"abs": 1e-9 * scale, "rel": 0.0}, "reference": ref(values, axis)}


def case(mode: str, step: int):
    d = load(f"{DATA}/dumps/chem170_emis_{mode}/esm_dump_emissions_{step}.json")["vars"]
    r = load(f"{DATA}/dumps/chem170_emis_{mode}/replay/emissions_{step}.json")["vars"]
    return d, r


def write(fname: str, model: str, tests: list) -> None:
    path = os.path.join(COMP, fname)
    doc = json.load(open(path))
    doc["models"][model]["tests"] = tests
    with open(path, "w") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"{fname}: {len(tests)} tests, {sum(len(t['assertions']) for t in tests)} assertions")


def wrf_dev(d, r, cols) -> float:
    """Largest in-model (real32) departure from the replay over the given chem columns, relative."""
    ci = np.asarray(d["chem_in"], dtype=float)[0]
    dw = np.asarray(d["chem_out"], dtype=float)[0] - ci
    dr = np.asarray(r["chem_out"], dtype=float)[0] - ci
    dev = [abs(dw[c] - dr[c]) / abs(dr[c]) for c in cols if dr[c] != 0.0]
    return max(dev) if dev else 0.0


def anthro_gas() -> list:
    direct = ["csl", "so2", "no", "ald", "hcho", "ora2", "nh3", "eth", "co", "ol2", "olt", "oli", "tol", "xyl", "ket"]
    fields = ["so2", "no", "no2", "co", "eth", "hc3", "hc5", "hc8", "xyl", "ol2", "olt", "oli", "tol", "csl",
              "hcho", "ald", "ket", "ora2", "nh3", "ch3oh", "c2h5oh"]
    tests = []
    for mode, step in ANTHRO:
        d, r = case(mode, step)
        dt = float(d["dt"])
        ea = np.asarray(d["emis_ant"], dtype=float)  # (kemit, num_emis_ant)
        inc = np.asarray(r["d_anthro_gas"], dtype=float)[0]
        over = {"emiss_inpt_opt": int(d["emiss_inpt_opt"]), "rho": float(np.asarray(d["rho"])[0]),
                "dz": float(np.asarray(d["dz8w"])[0])}
        for f in fields:
            over[f"E_{f}"] = float(ea[0, pe(f"e_{f}") - 1])
        a = [scalar(f"d_{s}", inc[p(s) - 1] / dt) for s in direct + ["par", "no2", "ch3oh", "c2h5oh"]]
        a.append(scalar("radm2_variety", 1.0 if int(d["emiss_inpt_opt"]) in (1, 102) else 0.0))
        cols = [p(s) - 1 for s in direct + ["par", "no2", "ch3oh", "c2h5oh"]]
        tests.append({
            "id": f"chem170_emis_{mode}_step{step}",
            "description": (f"Reference: the chem_opt = 170 emissions capture (data/eqwefic/chem170_scm/run_emis.sh, "
                            f"MODE={mode}: emiss_opt = 4, emiss_inpt_opt = {int(d['emiss_inpt_opt'])}, synthetic hourly "
                            f"wrfchemi files from make_wrfchemi.py; {BUILD}) at model step {step} (dt = {dt:g} s), "
                            f"surface layer of column (its, jts). Inputs are the dumped emission fields of layer 1 and "
                            f"the layer's density and depth; expected values are the real64 replay of "
                            f"cbmz_addemiss_anthro compiled unchanged (kernels/emissions_driver), its increment over "
                            f"the step divided by dt, at 1e-9 relative. The in-model real32 increments differ from the "
                            f"replay by at most {wrf_dev(d, r, cols):.2g} relative. Isoprene carries no tendency "
                            f"although E_ISO is non-zero (FORTRAN_BUGS.md N116)."),
            "parameter_overrides": over, "time_span": {"start": 0, "end": 1}, "assertions": a})
    return tests


def anthro_aer() -> list:
    comps = {"so4": "lptr_so4_aer", "no3": "lptr_no3_aer", "oc": "lptr_oc_aer", "bc": "lptr_bc_aer",
             "oin": "lptr_oin_aer", "num": "numptr_aer"}
    fields = ["so4i", "so4j", "so4c", "no3i", "no3j", "no3c", "orgi", "orgj", "orgc", "eci", "ecj", "ecc",
              "pm25i", "pm25j", "pm_10"]
    tests = []
    for mode, step in ANTHRO:
        d, r = case(mode, step)
        dt = float(d["dt"])
        ea = np.asarray(d["emis_ant"], dtype=float)
        inc = np.asarray(r["d_aer"], dtype=float)[0]
        over = {"emiss_inpt_opt": int(d["emiss_inpt_opt"]), "alt": float(np.asarray(d["alt"])[0]),
                "dz": float(np.asarray(d["dz8w"])[0]),
                "volumcen": [float(v) for v in np.asarray(d["volumcen_sect"], dtype=float)]}
        for f in fields:
            over[f"E_{f}"] = float(ea[0, pe(f"e_{f}") - 1])
        a, cols = [], []
        for c, ptr in comps.items():
            idx = np.asarray(d[ptr], dtype=int)
            cols += [int(i) - 1 for i in idx]
            a.append(vector(f"d_{c}", inc[idx - 1] / dt, "bin"))
        tests.append({
            "id": f"chem170_emis_{mode}_step{step}",
            "description": (f"Reference: the chem_opt = 170 emissions capture (data/eqwefic/chem170_scm/run_emis.sh, "
                            f"MODE={mode}: emiss_opt = 4, emiss_inpt_opt = {int(d['emiss_inpt_opt'])}; {BUILD}) at "
                            f"model step {step} (dt = {dt:g} s), surface layer. Inputs are the dumped aerosol emission "
                            f"fields, the layer's specific volume and depth and MOSAIC's run-time section centre "
                            f"volumes; expected values are the real64 replay of mosaic_addemiss (kernels/emissions_driver), "
                            f"each section's increment over the step divided by dt, at 1e-9 of each vector's maximum. "
                            f"The in-model real32 increments differ from the replay by at most "
                            f"{wrf_dev(d, r, cols):.2g} relative (cancellation in the large number mixing ratios)."),
            "parameter_overrides": over, "time_span": {"start": 0, "end": 1}, "assertions": a})
    return tests


def bio_cases():
    """(id, description stem, dump vars, replay vars) for the in-model and synthetic Gunther cases."""
    out = []
    for step in BIO_STEPS:
        d, r = case("default", step)
        out.append((f"chem170_emis_default_step{step}",
                    f"the chem_opt = 170 emissions capture (MODE=default; {BUILD}) at model step {step}", d, r, None))
    base = load(f"{DATA}/dumps/chem170_emis_default/esm_dump_emissions_1.json")["vars"]
    for tag in BIO_SYNTH:
        r = load(f"{DATA}/dumps/chem170_emis_synth/replay/{tag}.json")["vars"]
        flat = open(f"{DATA}/dumps/chem170_emis_synth/flat/{tag}.flat").read().split("\n")
        rec = {}
        for i in range(0, len(flat) - 1):
            parts = flat[i].split()
            if len(parts) >= 3 and parts[1] in ("int", "real", "bool") and parts[0] in ("ivgtyp", "julday", "gsw", "t_phy"):
                rec[parts[0]] = flat[i + 1].split()
        d = dict(base)
        d["ivgtyp"] = int(rec["ivgtyp"][0]); d["julday"] = int(rec["julday"][0]); d["gsw"] = float(rec["gsw"][0])
        t = np.asarray(base["t_phy"], dtype=float).copy(); t[0] = float(rec["t_phy"][0]); d["t_phy"] = t
        out.append((f"synthetic_{tag}",
                    f"a synthetic case: the step-1 dump of the chem_opt = 170 emissions capture (MODE=default) with "
                    f"its biogenic inputs rewritten (data/eqwefic/chem170_scm/emis_synth.py, case {tag})", d, r, tag))
    return out


def gunther() -> list:
    tests = []
    for tid, stem, d, r, _ in bio_cases():
        eb = np.asarray(r["e_bio"], dtype=float)
        over = {"lu": int(d["ivgtyp"]), "julday": int(d["julday"]), "gsw": float(d["gsw"]),
                "T": float(np.asarray(d["t_phy"])[0]), "p": float(np.asarray(d["p_phy"])[0])}
        a = [vector("e_bio", eb, "radm2")] if np.any(eb != 0) else [
            {"variable": "e_bio", "time": 0.0, "reduce": "Linf_error", "expected": 0.0,
             "tolerance": {"abs": 0.0, "rel": 0.0}, "reference": ref(eb, "radm2")}]
        tests.append({
            "id": tid,
            "description": (f"Reference: {stem}; USGS category {over['lu']}, day {over['julday']}, gsw "
                            f"{over['gsw']:.6g} W/m2, T {over['T']:.6g} K. Expected e_bio from the real64 replay of "
                            f"bio_emissions (chem/module_bioemi_simple.F compiled unchanged, kernels/emissions_driver), "
                            f"at 1e-9 of its maximum."),
            "parameter_overrides": over, "time_span": {"start": 0, "end": 1}, "assertions": a})
    return tests


def gunther_to_cbmz() -> list:
    species = ["iso", "oli", "xyl", "olt", "ald", "hcho", "eth", "ora2", "co", "no", "ol2", "ket", "par"]
    tests = []
    for tid, stem, d, r, _ in bio_cases():
        dt = float(d["dt"])
        inc = np.asarray(r["d_bio"], dtype=float)[0]
        over = {"e_bio": [float(v) for v in np.asarray(r["e_bio"], dtype=float)],
                "dz": float(np.asarray(d["dz8w"])[0])}
        a = [scalar(f"d_{s}", inc[p(s) - 1] / dt) for s in species]
        tests.append({
            "id": tid,
            "description": (f"Reference: {stem}. Inputs are the replay's e_bio and the surface layer's depth; expected "
                            f"values are the real64 replay of add_biogenics (GUNTHER1 branch, compiled unchanged, "
                            f"kernels/emissions_driver), its increment over the step divided by dt, at 1e-9 relative. "
                            f"The zero tendencies of ISO, OLI, XYL and ETH are FORTRAN_BUGS.md B25."),
            "parameter_overrides": over, "time_span": {"start": 0, "end": 1}, "assertions": a})
    return tests


if __name__ == "__main__":
    write("cbmz_anthropogenic.esm", "CBMZAnthroEmissions", anthro_gas())
    write("mosaic_anthropogenic.esm", "MosaicAnthroEmissions", anthro_aer())
    write("gunther_biogenic.esm", "GuntherBiogenicEmissions", gunther())
    write("gunther_to_cbmz.esm", "GuntherToCBMZ", gunther_to_cbmz())
