#!/usr/bin/env python3
"""Compare MOSAIC's instantaneous sub-process rates with the finite-step increments.

The WRF fork's `mosaic` dump (branch earthsciml-instrumented-chem170) carries, for
one column, both the state at every stage boundary and the rate each sub-process
forms internally (module_esm_mosaic_rates.F).  This script reports, per process,
the rate, the increment/dt recovered from the state, and their ratio, so a
stage-2 tolerance can be set on the quantity that is actually resolved.

    python3 tools/mosaic_rate_check.py <esm_dump_mosaic_<call>.json> [...]

Units.  rsub holds mol/mol-air for gases and aerosol mass, particles/mol-air for
number.  MOSAIC's own gas arrays are nmol/m^3, so the condensation rates are
converted with conv1a = cairclm * 1e15 (cairclm is mol/cm^3), the same factor
map_mosaic_species uses.  Coagulation works on num_distrib = rsub*cairclm
(#/cm^3), so its number rates are divided by cairclm to compare with rsub.
"""
from __future__ import annotations

import json
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


def report(path):
    d, A = load(path)
    names = [d["name_%04d" % i]["data"][0] for i in range(1, A("ltot2") + 1)]
    dt = A("dtchem")
    cair = A("cairclm")            # mol/cm^3
    conv1a = cair * 1.0e15         # mol/mol-air -> nmol/m^3
    nbin = A("rate_nbin")
    print(f"== {path}   ktau {A('ktau')}  dtchem {dt} s  nbin {nbin}")

    # --- nucleation: dt -> 0 rate vs the increment the model applies
    rate = A("rate_newnuc")
    fd = A("rate_newnuc_finitedt")
    rows = ["h2so4", "nh3", "so4_a", "nh4_a", "num_a"]
    act = np.abs(fd[:5]).max(axis=0) > 0
    if act.any():
        k = int(np.argmax(np.abs(fd[4])))
        print("   nucleation (level %d of %d active):" % (k + 1, act.sum()))
        for i, nm in enumerate(rows):
            r, f = rate[i, k], fd[i, k]
            ratio = r / f if f else float("nan")
            print("      %-6s dt->0 rate %12.5g   increment/dt %12.5g   ratio %9.3g"
                  % (nm, r, f, ratio))
    else:
        print("   nucleation: inactive at every level of this column")

    # --- coagulation: rate vs the increment recovered from rsub
    num = A("rate_coag_num")                      # (#/cm^3)/s
    prev, post = A("rsub_newnuc"), A("rsub_coag")
    numptr = A("numptr_aer").astype(int)
    print("   coagulation, number, at the level of the largest increment:")
    inc_all = np.array([(post[numptr[b] - 1] - prev[numptr[b] - 1]) / dt for b in range(nbin)])
    k = int(np.unravel_index(np.abs(inc_all).argmax(), inc_all.shape)[1])
    for b in range(nbin):
        r = num[b, k] / cair[k]                   # (#/cm^3)/s -> (#/mol-air)/s
        f = inc_all[b, k]
        ratio = r / f if f else float("nan")
        print("      bin %d  rate %12.5g   increment/dt %12.5g   ratio %9.3g"
              % (b + 1, r, f, ratio))

    # --- gas-particle transfer: the non-volatile rate against the H2SO4 loss
    nv = A("rate_cond_nonvol")                    # (ngas, nbin, nk), nmol/m^3/s
    ih2so4 = A("rate_igas_h2so4")
    tot = nv[ih2so4 - 1].sum(axis=0) / conv1a     # -> mol/mol-air/s, gas loss
    kh2so4 = names.index("h2so4")
    inc = (A("rsub_therm")[kh2so4] - A("rsub_in")[kh2so4]) / dt
    k = int(np.argmax(np.abs(tot)))
    print("   gas-particle transfer, H2SO4 (level %d):" % (k + 1))
    print("      condensation rate %12.5g   d(h2so4)/dt from the state %12.5g   ratio %9.3g"
          % (-tot[k], inc[k], (-tot[k] / inc[k]) if inc[k] else float("nan")))
    sv = A("rate_cond_semivol")
    for gas_name, idx, row in (("HNO3", A("rate_igas_hno3"), "hno3"),
                               ("HCl", A("rate_igas_hcl"), "hcl"),
                               ("NH3", A("rate_igas_nh3"), "nh3")):
        tot = sv[idx - 1].sum(axis=0) / conv1a
        kk = names.index(row)
        inc = (A("rsub_therm")[kk] - A("rsub_in")[kk]) / dt
        k = int(np.argmax(np.abs(tot)))
        ratio = (-tot[k] / inc[k]) if inc[k] else float("nan")
        print("      %-4s rate %12.5g   d/dt from the state %12.5g   ratio %9.3g  (level %d)"
              % (gas_name, -tot[k], inc[k], ratio, k + 1))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        report(p)
