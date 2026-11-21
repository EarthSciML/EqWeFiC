#!/usr/bin/env python3
"""Fill the reference profiles of couplings/eqweather_scm_trajectory.esm.

Rewrites ONLY the numbers in the hand-authored test `eqweather_scm_traj_step1410_60s`:
each assertion's `reference` const (or scalar `expected`) at t = 60 s. Which
reference each variable is compared against is fixed here and stated in the
test description:

- u, v, theta_m, q_v: WRF's dt -> 0 limit, the first-order extrapolation
  (15 X(dt=6) - 6 X(dt=15)) / 9 of two restarts of the reference run from its
  18:29 state (data/eqwefic/scm_ref_dtsens/dt6, dt15).
- w, phi', q_i, T_s, T_sk: the reference run itself at its own dt = 60 s
  (data/eqwefic/scm_ref_perstep, wrfout frame 1410), because for these the
  ESM-vs-WRF difference is inside WRF's own dt spread and the extrapolation
  does not converge (w, phi': acoustic mode; q_i: trace; T_s, T_sk: see the
  description).
- csza at t = 0, 30, 45, 54 s: WRF's COSZEN at 18:30 in the dt = 60 / 30 / 15 / 6
  runs, i.e. calc_coszen at the START of each run's last step (xtime 1409 +
  0, 30, 45, 54 s): four samples of the sun WRF moves, which pin the clock
  that drives the document's mounted solar geometry.

All are real32 WRF builds (../WRF/main/wrf.exe, em_scm_xy). The initial state
(frame 1409) is the one couplings/eqweather_scm.esm's step-1410 test uses,
bitwise. Run:  python3 tools/fill_eqweather_scm_trajectory.py
"""
import json
import sys

import numpy as np
import xarray as xr

DOC = "couplings/eqweather_scm_trajectory.esm"
DATA = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/"
PERSTEP = DATA + "scm_ref_perstep/wrfout_d01_1999-10-22_19:00:00"
DTSENS = DATA + "scm_ref_dtsens/dt{}/wrfout_d01_1999-10-23_18:30:00"
FRAME = 1410          # wrfout frame of t = 60 s (frame 1409 is t = 0)
WRF = {"u": "U", "v": "V", "w": "W", "phi_p": "PH", "theta_m": "THM",
       "qv": "QVAPOR", "qi": "QICE", "T_s": "TSLB", "T_sk": "TSK"}
OFFSET = {"theta_m": 300.0}  # THM is theta_m - 300
DT0 = {"u", "v", "theta_m", "qv"}
CSZA_DT = {0.0: 60, 30.0: 30, 45.0: 15, 54.0: 6}  # last step of run dt starts at 60 - dt


def column(ds, name, frame):
    a = np.asarray(ds[name].values[frame], dtype=np.float64)
    return a[..., 1, 1]


def main():
    ps = xr.open_dataset(PERSTEP)
    r6 = xr.open_dataset(DTSENS.format(6))
    r15 = xr.open_dataset(DTSENS.format(15))
    doc = json.load(open(DOC))
    test = doc["models"]["EqWeatherScm"]["tests"][0]
    for a in test["assertions"]:
        v = a["variable"]
        if v == "csza":
            r = xr.open_dataset(DTSENS.format(CSZA_DT[a["time"]]))
            a["expected"] = float(column(r, "COSZEN", 0))
            continue
        if v not in WRF:
            continue
        w = WRF[v]
        if v in DT0:
            x = (15.0 * column(r6, w, 0) - 6.0 * column(r15, w, 0)) / 9.0
        else:
            x = column(ps, w, FRAME)
        x = x + OFFSET.get(v, 0.0)
        if v == "T_s":
            x = x[:4]  # the 5th WRF layer is the fixed reservoir T_deep
        if "reference" in a:
            a["reference"]["expr"]["args"][0]["value"] = [float(y) for y in np.atleast_1d(x)]
        else:
            a["expected"] = float(x)
    open(DOC, "w").write(json.dumps(doc, indent=1, ensure_ascii=True) + "\n")
    sys.stdout.write("filled %d assertions\n" % len(test["assertions"]))


if __name__ == "__main__":
    main()
