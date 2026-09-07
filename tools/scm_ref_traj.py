#!/usr/bin/env python3
"""Single-column reference trajectory from the WRF SCM run.

The wrfout domain is 2x2 mass points (3x3 staggered) with periodic BCs and a
horizontally uniform initial state, so every mass point is bitwise identical.
We extract the mass point (j=1, i=1) 0-based.

Regenerate (needs netCDF4):
    python3 tools/scm_ref_traj.py extract
Load without netCDF4:
    from tools.scm_ref_traj import load
    tr = load(); tr.at('Tk', hours=12)[0]   # lowest model level Tk at t=12h
"""
import os
import sys

import numpy as np

SRC = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/scm_ref"
OUT = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/scm_ref_traj"
WRFOUT = os.path.join(SRC, "wrfout_d01_1999-10-22_19:00:00")
WRFINPUT = os.path.join(SRC, "wrfinput_d01")
J, I = 1, 1  # 0-based mass point
RCP = 2.0 / 7.0  # R_d/c_p used for the Exner function
G = 9.81  # gravity used for z = (PH+PHB)/g

# name -> (kind,) ; kind: 'sfc' scalar/time, 'lev' mass levels, 'levs' w levels, 'soil'
MOIST = ["QVAPOR", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP"]
SFC = ["TSK", "HFX", "LH", "QFX", "GRDFLX", "SWDOWN", "GLW", "OLR", "RAINNC",
       "RAINC", "SNOWNC", "T2", "Q2", "TH2", "PSFC", "U10", "V10", "PBLH",
       "UST", "ZNT", "COSZEN", "ALBEDO", "EMISS", "SNOW", "CANWAT", "SFROFF"]
SOIL = ["TSLB", "SMOIS", "SH2O"]
STATIC = ["ZNU", "ZNW", "ZS", "DZS", "DNW", "P_TOP", "XLAT", "XLONG", "HGT",
          "TMN", "VEGFRA", "XLAND", "IVGTYP", "ISLTYP", "LU_INDEX", "F", "E"]


def _column(ds, name):
    """Return the (time, ...) column at mass point (J, I), destaggering U/V."""
    v = ds.variables[name]
    a = np.asarray(v[:], dtype=np.float64)
    dims = v.dimensions
    if "west_east_stag" in dims:  # U: average the two x-faces of the cell
        return 0.5 * (a[..., J, I] + a[..., J, I + 1])
    if "south_north_stag" in dims:  # V: average the two y-faces
        return 0.5 * (a[..., J, I] + a[..., J + 1, I])
    if "west_east" in dims:
        return a[..., J, I]
    return a  # (Time,) or (Time, lev) already


def _derive(d):
    """Add derived fields to a dict of raw columns (in place) and return it."""
    d["theta"] = d["T"] + 300.0
    d["p"] = d["P"] + d["PB"]
    d["Tk"] = d["theta"] * (d["p"] / 1.0e5) ** RCP
    d["z_stag"] = (d["PH"] + d["PHB"]) / G
    d["z"] = 0.5 * (d["z_stag"][:, :-1] + d["z_stag"][:, 1:])
    d["W_m"] = 0.5 * (d["W"][:, :-1] + d["W"][:, 1:])
    return d


def _read(path):
    import netCDF4
    ds = netCDF4.Dataset(path)
    have = set(ds.variables)
    d = {}
    for n in (["T", "P", "PB", "U", "V", "W", "PH", "PHB", "CLDFRA"]
              + MOIST + SFC + SOIL + STATIC):
        if n in have:
            d[n] = np.atleast_1d(_column(ds, n))
    raw = ds.variables["Times"][:]
    d["Times"] = np.array(["".join(np.char.decode(np.asarray(r, dtype="S1")))
                           for r in raw])
    d["time_s"] = np.asarray(ds.variables["XTIME"][:], dtype=np.float64) * 60.0 \
        if "XTIME" in have else np.zeros(len(d["Times"]))
    ds.close()
    return _derive(d)


def extract(wrfout=WRFOUT, wrfinput=WRFINPUT, outdir=OUT):
    os.makedirs(outdir, exist_ok=True)
    out = _read(wrfout)
    ini = {"init_" + k: v for k, v in _read(wrfinput).items()}
    np.savez_compressed(os.path.join(outdir, "scm_ref_traj.npz"), **out, **ini)

    def csv(name, header, rows):
        with open(os.path.join(outdir, name), "w") as f:
            f.write(",".join(header) + "\n")
            for r in rows:
                f.write(",".join("%s" % x if isinstance(x, str) else "%.9g" % x
                                 for x in r) + "\n")

    nt = len(out["Times"])
    csv("times.csv", ["index", "iso", "time_s", "time_h"],
        [(i, out["Times"][i], out["time_s"][i], out["time_s"][i] / 3600.0)
         for i in range(nt)])
    sfc = [k for k in SFC if k in out]
    csv("surface.csv", ["time_s"] + sfc,
        [[out["time_s"][t]] + [out[k][t] for k in sfc] for t in range(nt)])
    lev = ["z", "p", "P", "PB", "theta", "T", "Tk", "U", "V", "W_m", "CLDFRA"] \
        + [q for q in MOIST if q in out]
    nz = out["theta"].shape[1]
    csv("profiles_mass.csv", ["time_s", "k"] + lev,
        [[out["time_s"][t], k] + [out[v][t, k] for v in lev]
         for t in range(nt) for k in range(nz)])
    csv("profiles_stag.csv", ["time_s", "k", "z_stag", "PH", "PHB", "W"],
        [[out["time_s"][t], k, out["z_stag"][t, k], out["PH"][t, k],
          out["PHB"][t, k], out["W"][t, k]]
         for t in range(nt) for k in range(nz + 1)])
    soil = [k for k in SOIL if k in out]
    csv("soil.csv", ["time_s", "layer", "depth_m"] + soil,
        [[out["time_s"][t], k, out["ZS"][t, k]] + [out[v][t, k] for v in soil]
         for t in range(nt) for k in range(out["ZS"].shape[1])])
    csv("static.csv", ["name", "k", "value"],
        [[n, k, float(np.ravel(out[n][0])[k])]
         for n in STATIC if n in out
         for k in range(np.ravel(np.atleast_1d(out[n][0])).size)])
    ini_lev = [v for v in lev if "init_" + v in ini]
    csv("init_profile.csv", ["k"] + ini_lev,
        [[k] + [ini["init_" + v][0, k] for v in ini_lev] for k in range(nz)])
    return outdir


class Traj:
    """Tiny accessor over the .npz extract (no netCDF needed)."""

    def __init__(self, z):
        self.z = z
        self.time_s = z["time_s"]

    def __getitem__(self, name):
        return self.z[name]

    @property
    def names(self):
        return sorted(self.z.files)

    def index(self, hours=None, seconds=None):
        t = seconds if seconds is not None else hours * 3600.0
        return int(np.argmin(np.abs(self.time_s - t)))

    def at(self, name, hours=None, seconds=None):
        """Field `name` at the history time nearest the requested time."""
        return self.z[name][self.index(hours, seconds)]


def load(outdir=OUT):
    return Traj(np.load(os.path.join(outdir, "scm_ref_traj.npz")))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "extract":
        print("wrote", extract())
    else:
        tr = load()
        print(len(tr.time_s), "frames;", len(tr.names), "fields")
