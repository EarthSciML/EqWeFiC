#!/usr/bin/env python3
"""Re-reference the coupled documents' YSU PBL assertions to the dt -> 0 limit.

WRF's driver-level rthblten/rublten/rvblten/rqvblten/rqcblten/rqiblten are the
IMPLICIT tridiagonal step's (c_new - c_old)/dt at dt = 60 s.  They are not
instantaneous derivatives, and TENDENCY_ERROR_BUDGET.md measured that this --
not any ESM defect -- is the whole of the 2-31 % coupled PBL residual.  The
component tests in components/atmospheric_dynamics/ysu/ysu.esm already assert
against the dt -> 0 limit; this script moves the COUPLED documents onto the
same reference.

The limit is the second-order Richardson extrapolation

    T(0) = (8 T(dt/4) - 6 T(dt/2) + T(dt)) / 3,   dt = 0.04 s

of real64 replays of the SAME dumped column through kernels/ysu_driver of the
ctessum-claude/WRF fork (data/eqwefic/dumps/ysu/scmref_<call>_dt<dt>.json,
produced from scm_ref_<call>.flat by `ESM_DT=<dt> ./bin/ysu_driver`).  The
kernel is built with NEED_B4B_DURING_CCPP_TESTING = 1, under which ttnp is a
TEMPERATURE tendency, so the theta reference is ttnp / pi2d.

Assertions whose reference is a SUM containing a PBL term are recomputed from
the new PBL arrays; the moist-theta conversion is applied as a DELTA on WRF's
own conv_t_tendf_to_moist output so that the dumped baseline is preserved.

This script only transcribes numbers into hand-authored tests blocks; it writes
no equation and adds no assertion.  The .esm files round-trip through
json.dumps(indent=1), so nothing else in them moves.

Usage:  python3 tools/reref_pbl_dt0.py [--check]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
YSU = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/dumps/ysu"
DYN = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/dumps/dyn_column"
DTS = ("0.04", "0.02", "0.01")
NLEV = 59

# (file, model, test index, YSU kernel call, WRF step)
TARGETS = [
    ("couplings/eqweather_scm.esm", "EqWeatherScm", 0, 2820, 1410),
    ("couplings/eqweather_scm.esm", "EqWeatherScm", 1, 120, 60),
    ("couplings/scm_physics_column.esm", "ScmPhysicsColumn", 0, 2820, 1410),
    ("couplings/scm_physics_column_night.esm", "ScmPhysicsColumn", 0, 1800, 900),
    ("couplings/geometry_surface_pbl_column.esm", "GeometrySurfacePBLColumn", 0, 2820, 1410),
    ("couplings/geometry_surface_pbl_column_night.esm", "GeometrySurfacePBLColumn", 0, 1800, 900),
]

# esm assertion name -> kernel output name
PBL_MAP = {
    "pbl_dudt": "utnp",
    "pbl_dvdt": "vtnp",
    "pbl_dthdt": "ttnp",        # divided by pi2d below
    "pbl_dqvdt": "qvtnp",
    "pbl_dqcdt": "qctnp",
    "pbl_dqidt": "qitnp",
}


def read_flat(path):
    lines = open(path).read().split("\n")
    out = {}
    for i in range(0, len(lines) - 1, 2):
        head = lines[i].split()
        if not head:
            break
        out[head[0]] = lines[i + 1]
    return out


def flat_col(flat, name):
    a = np.array(flat[name].split(), dtype=float)
    if a.size == 2 * NLEV:          # (its:ite, kts:kte), i fastest -> column 1
        a = a.reshape(NLEV, 2)[:, 0]
    return a


def replay(call, dt):
    v = json.load(open(f"{YSU}/scmref_{call}_dt{dt}.json"))["vars"]
    out = {}
    for k, val in v.items():
        if isinstance(val, dict) and "data" in val:
            a = np.array(val["data"], dtype=float)
            if a.size == 2 * NLEV:
                a = a.reshape(NLEV, 2)[:, 0]
            out[k] = a
    return out


def richardson(call, dts=DTS):
    a, b, c = (replay(call, d) for d in dts)
    return {k: (8 * c[k] - 6 * b[k] + a[k]) / 3 for k in PBL_MAP.values()}


def pbl_limit(call):
    """dt -> 0 PBL tendency column for every esm assertion name."""
    r = richardson(call)
    pi2d = flat_col(read_flat(f"{YSU}/scm_ref_{call}.flat"), "pi2d")
    out = {}
    for name, kern in PBL_MAP.items():
        out[name] = r[kern] / pi2d if kern == "ttnp" else r[kern]
    return out


def moist_weights(step):
    """A = 1 + (Rv/Rd) q_v and B = (Rv/Rd) theta for conv_t_tendf_to_moist."""
    v = json.load(open(f"{DYN}/esm_dump_dyn_suite_{step}.json"))["vars"]
    g = lambda k: np.array(v[k]["data"], dtype=float)
    rvovrd = float(v["dyn_rvovrd"]["data"][0])
    theta_m = g("dyn_t") + float(v["dyn_t0"]["data"][0])
    qv = g("rad_qv")
    A = 1.0 + rvovrd * qv
    return A, rvovrd * (theta_m / A)


def const_of(assertion):
    r = assertion.get("reference")
    if not (isinstance(r, dict) and r.get("op") == "faq"):
        return None
    e = r.get("expr")
    if not (isinstance(e, dict) and e.get("op") == "index"):
        return None
    c = e["args"][0]
    if isinstance(c, dict) and c.get("op") == "const" and isinstance(c.get("value"), list):
        return c
    return None


def dump_style(src, doc):
    """The (indent, ensure_ascii) pair under which this file round-trips byte for byte.

    The couplings were not all written with the same json.dumps settings, and
    re-emitting one in the other's style would rewrite every non-ASCII
    description line for nothing.
    """
    for ensure in (False, True):
        if json.dumps(doc, indent=1, ensure_ascii=ensure) + "\n" == src:
            return ensure
    raise SystemExit(f"{path}: does not round-trip through json.dumps(indent=1); "
                     "refusing to rewrite it")


def process(path, model, ti, call, step, check):
    src = open(os.path.join(REPO, path)).read()
    doc = json.loads(src)
    ensure = dump_style(src, doc)
    test = doc["models"][model]["tests"][ti]
    nodes = {}
    for a in test["assertions"]:
        c = const_of(a)
        if c is not None and len(c["value"]) == NLEV:
            nodes.setdefault(a["variable"], []).append(c)
    old = {k: np.array(v[0]["value"], dtype=float) for k, v in nodes.items()}

    new = dict(old)
    lim = pbl_limit(call)
    for name in PBL_MAP:
        if name in nodes:
            new[name] = lim[name]

    # sums that contain a PBL term
    if "rth_phys_sum" in nodes:
        new["rth_phys_sum"] = old["rth_phys_sum"] + (new["pbl_dthdt"] - old["pbl_dthdt"])
    if "du_dt" in nodes:
        new["du_dt"] = old["du_dt"] + (new["pbl_dudt"] - old["pbl_dudt"])
    if "dv_dt" in nodes:
        new["dv_dt"] = old["dv_dt"] + (new["pbl_dvdt"] - old["pbl_dvdt"])
    if "dtheta_m_dt_phys" in nodes:
        A, B = moist_weights(step)
        new["dtheta_m_dt_phys"] = (old["dtheta_m_dt_phys"]
                                   + A * (new["rth_phys_sum"] - old["rth_phys_sum"])
                                   + B * (new["pbl_dqvdt"] - old["pbl_dqvdt"]))
    if "dtheta_m_dt" in nodes:
        new["dtheta_m_dt"] = old["dtheta_m_dt"] + (new["dtheta_m_dt_phys"] - old["dtheta_m_dt_phys"])
    if "dtheta_m_conv_delta" in nodes:
        new["dtheta_m_conv_delta"] = new["dtheta_m_dt_phys"] - new["rth_phys_sum"]

    print(f"{path} [{test.get('id')}] call {call}")
    for name in sorted(new):
        if np.array_equal(new[name], old[name]):
            continue
        amp_o, amp_n = np.max(np.abs(old[name])), np.max(np.abs(new[name]))
        shift = np.max(np.abs(new[name] - old[name]))
        print(f"   {name:22s} shift {shift:.4e}  amp {amp_o:.4e} -> {amp_n:.4e}")
        if not check:
            for node in nodes[name]:
                node["value"] = [float(x) for x in new[name]]

    if not check:
        open(os.path.join(REPO, path), "w").write(json.dumps(doc, indent=1, ensure_ascii=ensure) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="report only, write nothing")
    args = ap.parse_args()
    for path, model, ti, call, step in TARGETS:
        process(path, model, ti, call, step, args.check)


if __name__ == "__main__":
    main()
