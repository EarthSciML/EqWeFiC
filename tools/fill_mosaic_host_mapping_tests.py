"""Fill components/aerosol/mosaic/mosaic_host_mapping.esm's tests from dumps/chem170_wrf.

Per step: the forward map (imap = 0), host columns chem_in / t_phy / p_phy / rho_phy /
moist_in(p_qv) against rsub_in / cairclm / ptotclm / relhumclm, and the reverse map
(imap > 0), MOSAIC's final state rsub_coag against chem_out. The dumps are real32, so the
expected values are WRF's real32 results and the tolerance is rel 1e-6.

    python3 tools/fill_mosaic_host_mapping_tests.py [steps...]
"""
import json
import os
import sys
from collections import OrderedDict as OD

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import esm_dump  # noqa: E402

W = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/dumps/chem170_wrf/"
DOC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "components", "aerosol", "mosaic", "mosaic_host_mapping.esm")
STEPS = [int(a) for a in sys.argv[1:]] or [1, 918, 1730]
GASES = ["so2", "h2so4", "no2", "no", "o3", "hno3", "h2o2", "hcho", "nh3", "n2o5", "no3", "pan", "hono", "oh", "ho2", "hcl"]
QV = 1          # moist_in column of p_qv
AER_LEVELS = [1, 2, 20, 40, 57, 59]
TOL = OD([("rel", 1e-6), ("abs", 0.0)])
# relhumclm comes out of a real32 exp of an argument near -22: one ulp of the argument is
# 2e-6 of the result, so it is bound at rel 1e-5 (measured worst between 5e-6 and 1e-5).
TOL_RH = OD([("rel", 1e-5), ("abs", 0.0)])


def build(s):
    m = esm_dump.load(W + f"esm_dump_mosaic_{s}.json")["vars"]
    names = [str(m["name_%04d" % i]) for i in range(1, int(m["ltot2"]) + 1)]
    gi = [names.index(g) for g in GASES]          # 0-based rsub row == chem column
    ci = np.asarray(m["chem_in"], dtype=float)    # [lev, species]
    co = np.asarray(m["chem_out"], dtype=float)
    ri = np.asarray(m["rsub_in"], dtype=float)    # [row, lev]
    rc = np.asarray(m["rsub_coag"], dtype=float)
    mptr = np.asarray(m["massptr_aer"], dtype=int) - 1   # [comp, bin]
    nptr = np.asarray(m["numptr_aer"], dtype=int) - 1
    wptr = np.asarray(m["waterptr_aer"], dtype=int) - 1
    hptr = np.asarray(m["hyswptr_aer"], dtype=int) - 1
    nlev, nbin, ncomp = ci.shape[0], len(nptr), mptr.shape[0]
    aer = lambda X, lev_first: np.stack([np.stack([X[:, mptr[c, b]] if lev_first else X[mptr[c, b], :] for c in range(ncomp)], axis=1) for b in range(nbin)], axis=1)  # noqa: E731
    per_bin = lambda X, ptr, lev_first: np.stack([X[:, ptr[b]] if lev_first else X[ptr[b], :] for b in range(nbin)], axis=1)  # noqa: E731
    po = OD()
    po["t_phy"] = np.asarray(m["t_phy"], dtype=float).tolist()
    po["p_phy"] = np.asarray(m["p_phy"], dtype=float).tolist()
    po["rho_phy"] = np.asarray(m["rho_phy"], dtype=float).tolist()
    po["qv"] = np.asarray(m["moist_in"], dtype=float)[:, QV].tolist()
    po["chem_gas"] = ci[:, gi].tolist()
    po["chem_aer"] = aer(ci, True).tolist()
    po["chem_num"] = per_bin(ci, nptr, True).tolist()
    po["chem_water"] = per_bin(ci, wptr, True).tolist()
    po["chem_hysw"] = per_bin(ci, hptr, True).tolist()
    po["mw"] = np.asarray(m["mw_aer"], dtype=float).tolist()
    po["gas_out"] = rc[gi, :].T.tolist()
    po["r_out"] = aer(rc, False).tolist()
    po["n_out"] = per_bin(rc, nptr, False).tolist()
    po["water_out"] = per_bin(rc, wptr, False).tolist()
    po["hysw_out"] = per_bin(rc, hptr, False).tolist()
    asr = []

    def put(var, coords, val):
        a = OD([("variable", var), ("coords", OD(coords)), ("time", 0.0), ("expected", float(val)),
                ("tolerance", TOL_RH if var == "relhum" else TOL)])
        asr.append(a)

    exp_lev = {"cair": np.asarray(m["cairclm"]), "ptot": np.asarray(m["ptotclm"]), "relhum": np.asarray(m["relhumclm"]),
               "temp": ri[int(m["ktemp"]) - 1], "h2o": ri[int(m["kh2o"]) - 1]}
    for v, arr in exp_lev.items():
        for k in range(nlev):
            put(v, [("lev", k + 1)], arr[k])
    for k in range(nlev):
        for g, row in enumerate(gi):
            put("gas", [("lev", k + 1), ("gas", g + 1)], ri[row, k])
            put("chem_gas_out", [("lev", k + 1), ("gas", g + 1)], co[k, row])
    for k in [x - 1 for x in AER_LEVELS]:
        for b in range(nbin):
            for c in range(ncomp):
                put("r", [("lev", k + 1), ("bin", b + 1), ("comp", c + 1)], ri[mptr[c, b], k])
                put("chem_aer_out", [("lev", k + 1), ("bin", b + 1), ("comp", c + 1)], co[k, mptr[c, b]])
            for v, ptr in [("n", nptr), ("water", wptr), ("hysw", hptr)]:
                put(v, [("lev", k + 1), ("bin", b + 1)], ri[ptr[b], k])
            for v, ptr in [("chem_num_out", nptr), ("chem_water_out", wptr), ("chem_hysw_out", hptr)]:
                put(v, [("lev", k + 1), ("bin", b + 1)], co[k, ptr[b]])
    tcold = int(np.sum(np.asarray(m["t_phy"]) < 233.0))
    desc = (f"The chem_opt = 170 column at model step {s} (dumps/chem170_wrf esm_dump_mosaic_{s}, WRF fork earthsciml-instrumented-chem170, "
            f"real32): forward, the host columns entering mosaic_aerchem_driver against MOSAIC's state rsub_in, cairclm, ptotclm and "
            f"relhumclm (every level; components, number and waters at levels {AER_LEVELS}); reverse, the state after nucleation and "
            f"coagulation rsub_coag against chem_out. {tcold} levels are colder than 233 K, where MOSAIC's temperature is the floor. "
            f"The {len(GASES)} gases: {', '.join(GASES)} (h2so4 is chem(p_sulf)). WRF maps in real32: rel 1e-6. "
            f"Filled by tools/fill_mosaic_host_mapping_tests.py.")
    return OD([("id", f"host_map_s{s}"), ("description", desc), ("parameter_overrides", po),
               ("time_span", OD([("start", 0), ("end", 1)])), ("assertions", asr)])


def main():
    doc = json.load(open(DOC), object_pairs_hook=OD)
    doc["models"]["MosaicHostMapping"]["tests"] = [build(s) for s in STEPS]
    # the test binds the column sizes of chem_opt = 170
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from esm_compact_json import dumps
    open(DOC, "w").write(dumps(doc) + "\n")


if __name__ == "__main__":
    main()
