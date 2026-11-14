"""Fill couplings/mosaic_aerchemistry_column.esm's tests from the chem170 dumps and the
real64 replays at ONE state per step (data/eqwefic/chem170_scm/aerchem64, run_aerchem64.sh).

The state is MOSAIC's nucleation-entry state rsub_therm of model step s (dumps/chem170_wrf
esm_dump_mosaic_<s>), with the carried dry diagnostics, the section geometry and the
component properties of the same dump. Dry deposition is evaluated in the same state's
lowest layer with u*, the aerodynamic resistance (esm_dump_mosaic_drydep_<s>) and dz8w(1)
(esm_dump_chem_vertmx_<s>) of the dry-deposition call numbered s (model step s+2: those
hooks count from the third step), a constructed state that the replay dd_<s>.json is the
reference for. References, all real64 replays of WRF's own routines at that state:
  nuc_<s>.json   rate_newnuc_replay: the nucleation increments (mosaic_newnuc_1clm)
  coag_<s>.json  rate_coag_*_replay: coagsolv's instantaneous rates with nucleation
                 skipped, i.e. at rsub_therm (ESM_SKIP_NEWNUC=1)
  dd_<s>.json    vdep_aer: mosaic_drydep_1clm on level 1 of rsub_therm
The summed tendencies are formed here from those references exactly as the coupling sums
them; each is asserted at rel 1e-9 with an absolute floor of 1e-9 of the largest term
that enters it (coagulation's production and loss cancel, so its net is bounded against
the column's largest coagulation rate of that quantity at that level).

    python3 tools/fill_mosaic_aerchem_column_tests.py [steps...]
"""
import json
import os
import sys
from collections import OrderedDict as OD

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import esm_dump  # noqa: E402
from esm_compact_json import dumps as compact_dumps  # noqa: E402

DATA = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic"
W = DATA + "/dumps/chem170_wrf/"
R = DATA + "/chem170_scm/aerchem64/"
DOC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "couplings", "mosaic_aerchemistry_column.esm")
STEPS = [int(a) for a in sys.argv[1:]] or [1, 364, 918, 1730]
REL = 1e-9
FLOOR = 1e-9
NCOMP = 10


def A(v):
    return np.asarray(v, dtype=float)


def tol(expected, scale):
    return OD([("rel", REL), ("abs", float(FLOOR * scale))])


def build(s, m=None, dd=None, vx=None, nuc=None, coag=None, dep=None, desc_in=None):
    """One test at model step s. By default the chem170_wrf dumps and the aerchem64 replays; an
    importer passes its own dumps (m with rsub_therm = its state) and replays, and its description."""
    if m is None:
        m = esm_dump.load(W + f"esm_dump_mosaic_{s}.json")["vars"]
        dd = esm_dump.load(W + f"esm_dump_mosaic_drydep_{s}.json")["vars"]
        vx = esm_dump.load(W + f"esm_dump_chem_vertmx_{s}.json")["vars"]
        assert int(dd["ktau"]) == int(m["ktau"]) + 2 == int(vx["ktau"])
        nuc = esm_dump.load(R + f"nuc_{s}.json")["vars"]
        coag = esm_dump.load(R + f"coag_{s}.json")["vars"]
        dep = esm_dump.load(R + f"dd_{s}.json")["vars"]
    rs = A(m["rsub_therm"])  # [235, 59], rows are 1-based rsub indices
    nlev = rs.shape[1]
    mptr = np.asarray(m["massptr_aer"], dtype=int)  # [10, 8]
    nptr = np.asarray(m["numptr_aer"], dtype=int)
    wptr = np.asarray(m["waterptr_aer"], dtype=int)
    hptr = np.asarray(m["hyswptr_aer"], dtype=int)
    nbin = len(nptr)
    row = lambda p: rs[p - 1, :]  # noqa: E731
    temp = row(int(m["ktemp"]))
    n = np.stack([row(nptr[b]) for b in range(nbin)], axis=1)          # [lev, bin]
    water = np.stack([row(wptr[b]) for b in range(nbin)], axis=1)
    hysw = np.stack([row(hptr[b]) for b in range(nbin)], axis=1)
    r = np.stack([np.stack([row(mptr[c, b]) for c in range(NCOMP)], axis=1) for b in range(nbin)], axis=1)  # [lev, bin, comp]
    cair = A(m["cairclm"])
    po = OD()
    po["Column.temp"] = temp.tolist()
    po["Column.ptot"] = A(m["ptotclm"]).tolist()
    po["Column.cair"] = cair.tolist()
    po["Column.rh"] = A(m["relhumclm"]).tolist()
    po["Column.q_h2so4"] = row(int(m["kh2so4"])).tolist()
    po["Column.q_nh3"] = row(int(m["knh3"])).tolist()
    po["Column.n"] = n.tolist()
    po["Column.r"] = r.tolist()
    po["Column.water"] = water.tolist()
    po["Column.hysw"] = hysw.tolist()
    po["Column.mass_in"] = A(m["aqmassdry_sub"]).T.tolist()
    po["Column.vol_in"] = A(m["aqvoldry_sub"]).T.tolist()
    po["Column.rho_in"] = A(m["adrydens_sub"]).T.tolist()
    po["Column.mw"] = A(m["mw_aer"]).tolist()
    po["Column.dens"] = A(m["dens_aer"]).tolist()
    for k, f in [("v_lo", "volumlo_sect"), ("v_hi", "volumhi_sect"), ("v_cen", "volumcen_sect"),
                 ("d_lo", "dlo_sect"), ("d_hi", "dhi_sect"), ("d_cen", "dcen_sect")]:
        po["Column." + k] = A(m[f]).tolist()
    po["Column.ustar"] = float(dd["ust"])
    po["Column.aer_res"] = float(dd["aer_res"])
    po["Column.dz1"] = float(A(vx["dz8w"])[0])
    po["Column.dt_nuc"] = float(m["dtchem"])

    # references
    rnuc = A(nuc["rate_newnuc_replay"])        # [6, lev]: h2so4, nh3, so4a, nh4a, num increments, isize
    assert np.all((rnuc[5] == 0) | (rnuc[5] == 1)), "nucleation outside section 1"
    cnum = A(coag["rate_coag_num_replay"]).T   # [lev, bin]
    cvol = np.transpose(A(coag["rate_coag_vol_replay"]), (2, 0, 1))  # [lev, bin, 13]
    vdep = A(dep["vdep_aer"])
    dt = float(m["dtchem"])
    dz1 = po["Column.dz1"]
    sfc = np.zeros(nlev); sfc[0] = 1.0 / dz1

    asr = []

    def put(var, coords, val, scale):
        a = OD([("variable", var)])
        if coords:
            a["coords"] = OD(coords)
        a["time"] = 0.0
        a["expected"] = float(val)
        a["tolerance"] = tol(val, scale)
        asr.append(a)

    for b in range(nbin):
        put("Column.vdep", [("bin", b + 1)], vdep[b], 0.0)
    for k, i in [("nuc_h2so4", 0), ("nuc_nh3", 1), ("nuc_so4a", 2), ("nuc_nh4a", 3), ("nuc_num", 4)]:
        mx = np.max(np.abs(rnuc[i]))
        for l in range(nlev):
            put("Column." + k, [("lev", l + 1)], rnuc[i, l], mx * 0 + 0.0)
    for l in range(nlev):
        mx = np.max(np.abs(cnum[l]))
        for b in range(nbin):
            put("Column.dnum_coag", [("lev", l + 1), ("bin", b + 1)], cnum[l, b], mx)
    # summed tendencies
    for l in range(nlev):
        mxn = np.max(np.abs(cnum[l] / cair[l]))
        for b in range(nbin):
            tn = cnum[l, b] / cair[l]
            tnu = (rnuc[4, l] / dt) if b == 0 else 0.0
            td = sfc[l] * vdep[b] * n[l, b]
            put("Column.dn_dt", [("lev", l + 1), ("bin", b + 1)], tn + tnu - td, mxn + abs(tnu) + abs(td))
        for c in range(NCOMP + 2):
            mxq = np.max(np.abs(cvol[l, :, c]))
            for b in range(nbin):
                tq = cvol[l, b, c]
                if c < NCOMP:
                    st = r[l, b, c]
                    tnu = 0.0
                    if b == 0 and c == 0:
                        tnu = rnuc[2, l] / dt
                    if b == 0 and c == 4:
                        tnu = rnuc[3, l] / dt
                    var, co = "Column.dr_dt", [("lev", l + 1), ("bin", b + 1), ("comp", c + 1)]
                else:
                    st = hysw[l, b] if c == NCOMP else water[l, b]
                    tnu = 0.0
                    var = "Column.dhysw_dt" if c == NCOMP else "Column.dwater_dt"
                    co = [("lev", l + 1), ("bin", b + 1)]
                td = sfc[l] * vdep[b] * st
                put(var, co, tq + tnu - td, mxq + abs(tnu) + abs(td))
        put("Column.dh2so4_dt", [("lev", l + 1)], rnuc[0, l] / dt, 0.0)
        put("Column.dnh3_dt", [("lev", l + 1)], rnuc[1, l] / dt, 0.0)
    nnuc = int(np.sum(rnuc[4] != 0))
    desc = (f"The whole chem_opt = 170 column (59 levels x 8 sections), model step {s}, at MOSAIC's nucleation-entry state "
            f"rsub_therm (after gas-particle transfer and water uptake) of dumps/chem170_wrf esm_dump_mosaic_{s}; nucleation "
            f"fires at {nnuc} of {nlev} levels. Dry deposition reads u* = {po['Column.ustar']:.4g} m/s, aerodynamic resistance "
            f"{po['Column.aer_res']:.5g} s/m and dz1 = {dz1:.5g} m from the dry-deposition and vertmx calls numbered {s} "
            f"(model step {s + 2}): a constructed state, whose reference is a replay at exactly these inputs. Expected values: "
            f"real64 replays of WRF's own routines at this one state (data/eqwefic/chem170_scm/aerchem64: nuc_{s}, coag_{s} with "
            f"nucleation skipped so coagulation sees rsub_therm, dd_{s}), summed as the coupling sums them; rel 1e-9 with an "
            f"absolute floor of 1e-9 of the largest term entering each assertion. Filled by tools/fill_mosaic_aerchem_column_tests.py.")
    if desc_in is not None:
        desc = desc_in
    return OD([("id", f"aerchem_s{s}"), ("description", desc), ("parameter_overrides", po),
               ("time_span", OD([("start", 0), ("end", 0)])), ("assertions", asr)])


def main():
    doc = json.load(open(DOC), object_pairs_hook=OD)
    doc["models"]["Column"]["tests"] = [build(s) for s in STEPS]
    open(DOC, "w").write(compact_dumps(doc) + "\n")


if __name__ == "__main__":
    main()
