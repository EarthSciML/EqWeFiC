#!/usr/bin/env python3
"""Compare MOSAIC's nucleation stage before and after the B20 fix.

Usage:
    python3 tools/newnuc_b20_compare.py <stock_dump_dir> <fixed_dump_dir> [step ...]

Both directories hold `esm_dump_mosaic_<step>.json` from the same em_scm_xy
chem_opt = 170 case, one built from stock WRF 4.8 and one with
`qh2so4_avail = qh2so4_cur - qh2so4_crit` (FORTRAN_BUGS B20).

For each step it reports, over the levels where nucleation fires:
  * the H2SO4 vapour left after the call, as a fraction of what entered and of
    the critical concentration the scheme relaxes to,
  * the (mole NH4)/(mole SO4) of the new particles,
  * the new particle number, and
  * how many levels the guard let through.

Only step 1 is a like-for-like model comparison: the two runs' trajectories
separate afterwards, because the aerosol feeds Dudhia's scattering (N86).  For
the later steps use the kernel replay, which runs both schemes on the same
input column.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np


def load(path):
    d = json.load(open(path))["vars"]

    def A(k):
        v = d[k]
        if not v["shape"]:
            return v["data"][0]
        return np.array(v["data"]).reshape(v["shape"], order="F")

    return d, A


def stage(path):
    """Per-level nucleation summary of one dump."""
    d, A = load(path)
    names = [d["name_%04d" % i]["data"][0] for i in range(1, A("ltot2") + 1)]
    kh = names.index("h2so4")
    dt = A("dtchem")
    before = A("rsub_therm")[kh]          # vapour entering the nucleation stage
    after = A("rsub_newnuc")[kh]          # and leaving it
    fd = A("rate_newnuc_finitedt")        # increments / dtnuc
    act = np.abs(fd[4]) > 0
    dq_so4 = fd[2] * dt
    dq_nh4 = fd[3] * dt
    with np.errstate(divide="ignore", invalid="ignore"):
        nh4_per_so4 = np.where(dq_so4 > 0, dq_nh4 / dq_so4, np.nan)
    return dict(act=act, before=before, after=after, num=fd[4] * dt,
                nh4_per_so4=nh4_per_so4, dt=dt,
                crit=A("newnuc_composition")[2] if "newnuc_composition" in d else None)


def report(stock_dir, fixed_dir, steps):
    print("nucleation stage, stock WRF vs the B20 fix")
    print("step | levels firing        | vapour left / vapour in   | NH4/SO4 of new particles | number made (per mol-air)")
    print("     | stock  fixed         | stock      fixed          | stock     fixed          | stock        fixed")
    for n in steps:
        ps, pf = (os.path.join(d, f"esm_dump_mosaic_{n}.json") for d in (stock_dir, fixed_dir))
        if not (os.path.exists(ps) and os.path.exists(pf)):
            continue
        s, f = stage(ps), stage(pf)
        if not (s["act"].any() or f["act"].any()):
            print("%5d |   0      0           | (nucleation fires at no level in either run)" % n)
            continue
        def frac(x):
            m = x["act"]
            return np.median(x["after"][m] / x["before"][m]) if m.any() else float("nan")
        def med(x, key):
            m = x["act"]
            return np.nanmedian(x[key][m]) if m.any() else float("nan")
        print("%5d | %4d   %4d          | %9.4g  %9.4g      | %8.4g  %8.4g       | %10.4g  %10.4g"
              % (n, s["act"].sum(), f["act"].sum(), frac(s), frac(f),
                 med(s, "nh4_per_so4"), med(f, "nh4_per_so4"),
                 med(s, "num"), med(f, "num")))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    steps = [int(x) for x in sys.argv[3:]] or [1, 60, 230, 364, 601, 639, 689, 918,
                                               1080, 1621, 1730, 1801, 2101, 2134, 2400, 2521]
    report(sys.argv[1], sys.argv[2], steps)
