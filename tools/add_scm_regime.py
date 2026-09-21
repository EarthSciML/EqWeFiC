#!/usr/bin/env python3
"""Transcribe one more WRF step of the em_scm_xy reference run into an inline
test of couplings/eqweather_scm.esm.

The document's test SHAPE is hand-authored: which variables are asserted, with
which reduction, over which axis, and why, is fixed by the existing step-1410
test and is not generated here.  What this script does is what CLAUDE.md allows
a script to do -- read the numbers for another step out of the Fortran dumps and
put them in the already-authored places.  Every expression below is the one the
existing fill sidecars (couplings/scm_physics_column*.esm.fill.json,
couplings/scm_dynamics_column*.esm.fill.json) already use for that name, and
`--verify` re-derives the two EXISTING tests from the same table and reports any
disagreement, so a wrong expression cannot slip in unnoticed.

The PBL tendencies are NOT taken from the driver-level rthblten/rublten/...:
those are WRF's implicit dt = 60 s step, not a derivative.  They are the
dt -> 0 Richardson limit of the YSU kernel replays, exactly as
tools/reref_pbl_dt0.py computes them for the existing tests.

Usage:
    python3 tools/add_scm_regime.py --verify
    python3 tools/add_scm_regime.py --step 900 --id eqweather_scm_step900_night
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "tools"))
from reref_pbl_dt0 import pbl_limit  # noqa: E402

DUMPS = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/dumps"
DOC = os.path.join(REPO, "couplings/eqweather_scm.esm")
MODEL = "EqWeatherScm"
NL = 59

STEPS = {60: 120, 300: 600, 900: 1800, 1410: 2820}   # WRF step -> physics call


# --------------------------------------------------------------------------- namespace
def _load(path, tag=None):
    out = {}
    for k, x in json.load(open(path))["vars"].items():
        if not (isinstance(x, dict) and "data" in x):
            continue
        a = np.array(x["data"], dtype=float)
        # the tile is two columns wide with i fastest; take column 1
        if a.size == 2 * NL:
            a = a.reshape(NL, 2)[:, 0]
        elif a.size == 2 * (NL + 1):
            a = a.reshape(NL + 1, 2)[:, 0]
        elif a.size == 2:
            a = a[:1]
        out[(tag + "_" if tag else "") + k] = a
    return out


def namespace(step):
    call = STEPS[step]
    n = {}
    n.update(_load(f"{DUMPS}/dyn_column/esm_dump_dyn_suite_{step}.json"))
    n.update(_load(f"{DUMPS}/subassembly/esm_dump_scm_suite_{step}.json"))
    n.update(_load(f"{DUMPS}/subassembly/esm_dump_subassembly_mp_{step}.json", "mp"))
    n.update(_load(f"{DUMPS}/solthm_raw/esm_dump_theta_m_conv_{step}.json", "thm"))
    n.update({"pbl_" + k: v for k, v in pbl_limit(call).items()})   # dt -> 0 YSU limit
    return n


def s(x):
    """A dump entry used as a scalar."""
    return float(np.ravel(x)[0])


# --------------------------------------------------------------------------- values
def parameters(N):
    return {
        "psfc0": s(N["all_psfc"]), "xland0": s(N["all_xland"]),
        "mavail0": s(N["pbl_mavail"]), "pblh_prev0": s(N["pbl_pblh"]),
        "dx0": s(N["sfc_dx"]), "lakemask0": 0.0, "water_depth0": 0.0,
        "znt_prev0": s(N["sfc_znt"]), "ust_prev0": s(N["sfc_ust"]),
        "mol_prev0": s(N["sfc_mol"]), "qsfc_prev0": s(N["sfc_qsfc"]),
        "hfx_prev0": s(N["sfc_hfx"]), "qfx_prev0": s(N["sfc_qfx"]),
        "ustm_prev0": s(N["sfc_ustm"]), "chs_prev0": s(N["sfc_chs"]),
        "chs2_prev0": s(N["sfc_chs2"]), "cqs2_prev0": s(N["sfc_cqs2"]),
        "isfflx0": s(N["pbl_isfflx"]), "shalwater_z00": 0.0, "scm_force_flux0": 0.0,
        "isftcflx0": 0.0, "iz0tlnd0": 0.0,
        "topdown0": 1.0 if s(N["pbl_ysu_topdown_pblmix"]) else 0.0,
        "co2vmr0": s(N["rrtm_co2vmr"]), "n2ovmr0": s(N["rrtm_n2ovmr"]),
        "ch4vmr0": s(N["rrtm_ch4vmr"]), "emiss0": s(N["all_emiss"]),
        "csza0": s(N["sw_coszen"]), "albedo0": s(N["sw_albedo"]),
        "solcon0": s(N["sw_solcon"]), "obscur0": s(N["sw_obscur"]),
        "cssca0": s(N["sw_cssca"]), "icloud0": s(N["sw_icloud"]),
        "tmn0": s(N["pbl_tmn"]), "thc0": 0.03999999910593033,
        "snowc0": s(N["all_snowc"]), "radiation0": 1.0, "ifsnow0": s(N["pbl_ifsnow"]),
        "T_deep0": float(N["pbl_tslb"][4]), "z_deep0": float(N["pbl_zs"][4]),
        "mu_p0": s(N["dyn_mu"]), "mu_b0": s(N["dyn_mub"]), "p_top0": s(N["grid_p_top"]),
        "f0": s(N["dyn_f"]), "e0": s(N["dyn_e"]), "sina0": s(N["dyn_sina"]),
        "cosa0": s(N["dyn_cosa"]), "dampcoef0": s(N["dyn_dampcoef"]),
        "zdamp0": s(N["dyn_zdamp"]),
    }


def initial_conditions(N):
    return {
        "u": N["dyn_u"], "v": N["dyn_v"], "w": N["dyn_w"], "phi_p": N["dyn_ph"],
        "theta_m": N["dyn_t"] + 300.0, "qv": N["rad_qv"], "qc": N["rad_qc"],
        "qi": N["rad_qi"], "qr": N["rad_qr"], "qs": N["rad_qs"], "qg": N["rad_qg"],
        "T_s": N["pbl_tslb"][:4], "T_sk": N["pbl_tsk"],
    }


def references(N):
    """Every asserted name -> its reference value (array or scalar)."""
    g = N.get
    mut = s(N["thm_mut"])
    dt = s(N["pbl_dt"])

    du_dyn = (g("dyn_dru_coriolis") + g("dyn_dru_curvature") + g("dyn_dru_rayleigh")) / g("dyn_muu")
    dv_dyn = (g("dyn_drv_coriolis") + g("dyn_drv_curvature") + g("dyn_drv_rayleigh")) / g("dyn_muv")
    dw_nonpg = (g("dyn_drw_coriolis") + g("dyn_drw_curvature") + g("dyn_drw_rayleigh")) / g("dyn_mut")
    dw = dw_nonpg + g("dyn_drw_pg_buoy") / g("dyn_mut")
    dth_dyn = g("dyn_dt_rayleigh") / g("dyn_mut")

    rthraten = g("all_rthraten")
    pbl_dthdt = g("pbl_pbl_dthdt")
    rth_sum = rthraten + pbl_dthdt

    # conv_t_tendf_to_moist, applied as a DELTA on WRF's own output so that the
    # dumped real32 baseline is preserved and only the PBL re-reference moves.
    qv = g("rad_qv")
    rvovrd = s(N["thm_rvovrd"])
    A = 1.0 + rvovrd * qv
    B = rvovrd * ((g("dyn_t") + 300.0) / A)
    dthm_phys = (g("thm_t_tendf_out") / mut
                 + A * (rth_sum - g("thm_t_tendf_in") / mut)
                 + B * (g("pbl_pbl_dqvdt") - g("thm_qv_tend") / mut))

    mp = lambda n: (g("mp_" + n + "_out") - g("mp_" + n)) / s(N["mp_dt"])
    zero = np.zeros(NL)

    return {
        "du_dt_dyn": du_dyn, "dv_dt_dyn": dv_dyn, "dw_dt_nonpg": dw_nonpg,
        "dw_dt": dw, "dtheta_m_dt_dyn": dth_dyn, "p_p": g("dyn_p"),
        "dphi_dt": g("dyn_w"),                       # reference is wrf.g * this
        "sfc_hfx": s(N["sfc_hfx_out"]), "sfc_qfx": s(N["sfc_qfx_out"]),
        "sfc_lh": s(N["sfc_lh"]), "sfc_qsfc": s(N["sfc_qsfc_out"]),
        "sfc_wspd": s(N["sfc_wspd_out"]), "sfc_fm": s(N["sfc_fm_out"]),
        "sfc_fh": s(N["sfc_fh_out"]), "sfc_gz1oz0": s(N["sfc_gz1oz0_out"]),
        "sfc_ust": s(N["sfc_ust_out"]), "sfc_znt": s(N["sfc_znt_out"]),
        "sfc_psim": s(N["sfc_psim_out"]), "sfc_psih": s(N["sfc_psih_out"]),
        "sfc_br": s(N["sfc_br_out"]), "sfc_u10": s(N["sfc_u10"]),
        "sfc_v10": s(N["sfc_v10"]), "sfc_th2": s(N["sfc_th2"]),
        "sfc_t2": s(N["sfc_t2"]), "sfc_q2": s(N["sfc_q2"]),
        "sfc_flhc": s(N["sfc_flhc_out"]), "sfc_flqc": s(N["sfc_flqc_out"]),
        "sfc_regime": s(N["sfc_regime_out"]), "sfc_rmol": s(N["sfc_rmol_out"]),
        "sfc_mol": s(N["sfc_mol_out"]),
        "pbl_hpbl": s(N["all_pblh"]), "pbl_kpbl": s(N["pbl_kpbl_out"]),
        "pbl_K_h": np.append(g("pbl_exch_h"), 0.0),
        "pbl_K_m": np.append(g("pbl_exch_m"), 0.0),
        "pbl_dthdt": pbl_dthdt, "pbl_dudt": g("pbl_pbl_dudt"),
        "pbl_dvdt": g("pbl_pbl_dvdt"), "pbl_dqvdt": g("pbl_pbl_dqvdt"),
        "pbl_dqcdt": g("pbl_pbl_dqcdt"), "pbl_dqidt": g("pbl_pbl_dqidt"),
        "F_up": g("rad_totuflux"), "F_dn": g("rad_totdflux"),
        "dthdt_lw": g("rad_rthratenlw"), "dthdt_sw": g("rad_rthratensw"),
        "rthraten": rthraten, "rth_phys_sum": rth_sum,
        "gsw": s(N["all_gsw"]), "glw": s(N["all_glw"]), "olr": s(N["all_olr"]),
        "cldfra": g("rad_cldfra"),
        "lsm_hfx": s(N["pbl_hfx_out"]), "lsm_qfx": s(N["pbl_qfx_out"]),
        "lsm_lh": s(N["pbl_lh_out"]), "lsm_qsfc": s(N["pbl_qsfc_out"]),
        "lsm_capg": s(N["pbl_capg_out"]), "lsm_land": 1.0,
        "dTsk_dt_wrf_step": (s(N["pbl_tsk_out"]) - s(N["pbl_tsk"])) / dt,
        "lsm_dTs_dt": (g("pbl_tslb_out")[:4] - g("pbl_tslb")[:4]) / dt,
        "dtheta_m_dt_phys": dthm_phys,
        "dtheta_m_conv_delta": dthm_phys - rth_sum,
        "dtheta_m_dt": dth_dyn + dthm_phys,
        "du_dt": du_dyn + g("pbl_pbl_dudt"), "dv_dt": dv_dyn + g("pbl_pbl_dvdt"),
        "mp_dqv_dt": mp("qv"), "mp_dqc_dt": mp("qc"), "mp_dqi_dt": mp("qi"),
        "mp_dqr_dt": mp("qr"), "mp_dqs_dt": mp("qs"), "mp_dqg_dt": mp("qg"),
        "mp_dtheta_dt": mp("th"),
        "mp_melt_dqi_dt": zero, "mp_rate_dqi_dt": zero, "mp_sed_dqi_dt": zero,
    }


# --------------------------------------------------------------------------- test I/O
def const_node(assertion):
    r = assertion.get("reference")
    if isinstance(r, dict) and r.get("op") == "faq":
        e = r["expr"]
        while isinstance(e, dict) and e.get("op") in ("*", "index"):
            cands = [a for a in e["args"] if isinstance(a, dict)]
            hit = [c for c in cands if c.get("op") == "const"]
            if hit:
                return hit[0]
            e = cands[0] if cands else None
    return None


def assertion_value(a):
    c = const_node(a)
    if c is not None:
        return np.array(c["value"], dtype=float)
    if "coords" in a or a.get("reduce") is None:
        return float(a["expected"])
    return None


def set_assertion(a, val):
    c = const_node(a)
    if c is not None:
        c["value"] = [float(x) for x in np.ravel(val)]
    else:
        a["expected"] = float(np.ravel(val)[0])


def coord_value(a, val):
    """The scalar a `coords` assertion should carry, 1-based indexing."""
    i = list(a["coords"].values())[0]
    return float(np.ravel(val)[int(i) - 1])


def verify(doc):
    bad = 0
    for ti, step in ((0, 1410), (1, 60)):
        t = doc["models"][MODEL]["tests"][ti]
        N = namespace(step)
        P, IC, R = parameters(N), initial_conditions(N), references(N)
        for k, v in t["parameter_overrides"].items():
            if k not in P:
                print(f"  test{ti} param {k}: NO EXPRESSION"); bad += 1
            elif float(v) != P[k]:
                print(f"  test{ti} param {k}: {v} != {P[k]}"); bad += 1
        for k, v in t.get("initial_conditions", {}).items():
            a = np.array(v, dtype=float).ravel()
            if k not in IC:
                print(f"  test{ti} ic {k}: NO EXPRESSION"); bad += 1
            elif not np.array_equal(a, np.ravel(IC[k])):
                print(f"  test{ti} ic {k}: max diff {np.max(np.abs(a - np.ravel(IC[k]))):.3e}"); bad += 1
        for a in t["assertions"]:
            v = a["variable"]
            got = assertion_value(a)
            if v not in R:
                print(f"  test{ti} assert {v}: NO EXPRESSION"); bad += 1
                continue
            want = coord_value(a, R[v]) if "coords" in a else R[v]
            want = np.ravel(np.asarray(want, dtype=float))
            got = np.ravel(np.asarray(got, dtype=float))
            if got.shape != want.shape:
                print(f"  test{ti} assert {v}: shape {got.shape} vs {want.shape}"); bad += 1
            else:
                d = np.max(np.abs(got - want))
                amp = max(np.max(np.abs(want)), 1e-300)
                # 1e-7 relative, not exact: the four rth_phys_sum-derived
                # references were transcribed from WRF's real32 all_rth_phys_sum
                # while this table forms rthraten + rthblten itself, which
                # differs in the last real32 bits (<= 2e-11 on a 4e-4 signal).
                if d > 3e-7 * amp:
                    print(f"  test{ti} assert {v}: max diff {d:.3e} (amp {amp:.3e})"); bad += 1
    print("VERIFY:", "clean" if bad == 0 else f"{bad} disagreements")
    return bad


MP_NAMES = ["mp_dqv_dt", "mp_dqc_dt", "mp_dqi_dt", "mp_dqr_dt", "mp_dqs_dt",
            "mp_dqg_dt", "mp_dtheta_dt"]


def build(doc, step, test_id, description, template=0, with_mp=True):
    import copy
    t = copy.deepcopy(doc["models"][MODEL]["tests"][template])
    if with_mp:
        have = {a["variable"] for a in t["assertions"]}
        donor = {a["variable"]: a for a in doc["models"][MODEL]["tests"][1]["assertions"]}
        for n in MP_NAMES:
            if n not in have:
                t["assertions"].append(copy.deepcopy(donor[n]))
    N = namespace(step)
    P, IC, R = parameters(N), initial_conditions(N), references(N)
    t["id"] = test_id
    t["description"] = description
    t["parameter_overrides"] = P
    t["initial_conditions"] = {k: (float(np.ravel(v)[0]) if np.size(v) == 1
                                   else [float(x) for x in np.ravel(v)])
                               for k, v in IC.items()}
    kept = []
    for a in t["assertions"]:
        v = a["variable"]
        if v not in R:
            raise SystemExit(f"no expression for {v}")
        set_assertion(a, coord_value(a, R[v]) if "coords" in a else R[v])
        kept.append(a)
    t["assertions"] = kept
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--step", type=int)
    ap.add_argument("--id")
    ap.add_argument("--template", type=int, default=0)
    ap.add_argument("--description", default="")
    args = ap.parse_args()
    src = open(DOC).read()
    doc = json.loads(src)
    if args.verify:
        sys.exit(1 if verify(doc) else 0)
    t = build(doc, args.step, args.id, args.description, args.template)
    doc["models"][MODEL]["tests"].append(t)
    open(DOC, "w").write(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
    print(f"appended {args.id} ({len(t['assertions'])} assertions)")


if __name__ == "__main__":
    main()
