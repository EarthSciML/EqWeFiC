#!/usr/bin/env python3
"""Append chem_opt = 170 (CBM-Z) tests to the two Wesely components.

The Wesely components were filled from the RADM2 run (chem_opt = 101); their per-gas
tables now carry a CBM-Z row (wesely_data.esm, mech = 2).  This writes, for each chosen
step of the chem170 reference run (dumps/chem170_wrf, wesely hook in wesely_driver):

  surface_resistance.esm   rc (= rcx) and srfres of every CBM-Z gas slot 2..65, in
                           parts of NSPEC = 10 (the last part overlaps the one before)
  deposition_velocity.esm  vd (= ddvel0d) of every slot except 3 (sulf, which takes the
                           particle path through deppart; dry_dep_driver then gives it
                           HNO3's velocity anyway), in parts of NSPEC = 6, plus aer_res

Only the test blocks are written; the equations are hand-authored.  Re-running replaces
the tests it wrote before (ids starting wesely_cbmz_).

    python3 tools/fill_wesely_cbmz_tests.py [steps...]      # default 364 1080 2101
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import esm_dump  # noqa: E402

ROOT = os.path.dirname(HERE)
P = os.path.join(ROOT, 'components/atmospheric_deposition/wrf_wesely')
D = '/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/dumps/chem170_wrf'
REF = ("Reference: the chem_opt = 170 (CBMZ_MOSAIC_KPP) em_scm_xy reference run "
       "(data/eqwefic/chem170_scm/run.sh, dumps/chem170_wrf; GNU gfortran REAL*4; fork "
       "ctessum-claude/WRF branch earthsciml-instrumented-chem170, wesely hook in "
       "wesely_driver) at model step {step}, column (its, jts): iland {iland}, iseason "
       "{iseason}, ta {ta:.2f} K, rad {rad:.1f} W/m2, rh {rh:.1f} %, {wet}. Mechanism row 2 "
       "of wesely_data (the CBM-Z gas slots). Single-precision reference, rel 1e-5.")


def parts(slots, n):
    out = [slots[i:i + n] for i in range(0, len(slots), n)]
    if len(out[-1]) < n:
        out[-1] = slots[-n:]
    return out


def a(var, s, val):
    return {"variable": var, "time": 0.0, "coords": {"spc": s}, "expected": float(val),
            "tolerance": {"rel": 1e-05, "abs": 1e-30}}


def main(steps):
    sr_path, dv_path = f'{P}/surface_resistance.esm', f'{P}/deposition_velocity.esm'
    sr, dv = json.load(open(sr_path)), json.load(open(dv_path))
    srm = sr['models']['WRFWeselySurfaceResistance']
    dvm = dv['models']['WRFWeselyDepositionVelocity']
    srm['tests'] = [t for t in srm['tests'] if not t['id'].startswith('wesely_cbmz_')]
    dvm['tests'] = [t for t in dvm['tests'] if not t['id'].startswith('wesely_cbmz_')]
    for step in steps:
        v = esm_dump.load(f'{D}/esm_dump_wesely_{step}.json')['vars']
        numgas = int(v['numgas'])
        wet = ', '.join(n for n in ('rainflag', 'wetflag', 'highnh3') if v[n]) or 'dry'
        ref = REF.format(step=step, iland=int(v['iland']), iseason=int(v['iseason']),
                         ta=float(v['ta']), rad=float(v['rad']), rh=float(v['rhchem']),
                         wet=wet)
        for i, sl in enumerate(parts(list(range(2, numgas + 1)), 10), 1):
            srm['tests'].append({
                "id": f"wesely_cbmz_rc_step{step}_part{i}",
                "description": ref + f" Gas slots {', '.join(map(str, sl))} on the species axis. "
                               "srfres is zero for the LAST TWO slots (64 lim, 65 isopo2), "
                               "not HO and HO2 (FORTRAN_BUGS N124).",
                "parameter_overrides": {
                    "iland": float(v['iland']), "iseason": float(v['iseason']),
                    "T": float(v['ta']), "rad": float(v['rad']), "rh": float(v['rhchem']),
                    "wetflag": float(v['wetflag']), "rainflag": float(v['rainflag']),
                    "highnh3": float(v['highnh3']), "mech": 2.0, "n_gas": float(numgas),
                    "igas": [float(s) for s in sl]},
                "time_span": {"start": 0, "end": 1},
                "tolerance": {"rel": 1e-05},
                "assertions": [a('rc', j + 1, v['rcx'][s - 1]) for j, s in enumerate(sl)]
                + [a('srfres', j + 1, v['srfres'][s - 1]) for j, s in enumerate(sl)]})
        dslots = [s for s in range(2, numgas + 1) if s != 3]
        for i, sl in enumerate(parts(dslots, 6), 1):
            t = {"id": f"wesely_cbmz_vd_step{step}_part{i}",
                 "description": ref + f" ustar {float(v['ustar']):.4f} m/s, rmol "
                 f"{float(v['rmol']):.5f} 1/m, z0 {float(v['zntt']):.3f} m, first layer "
                 f"{float(v['dz8w1']):.2f} m. Gas slots {', '.join(map(str, sl))}; Rc from "
                 "the same dump, so this isolates the Ra/Rb/cellvg chain. Expected value is "
                 "ddvel0d, before dry_dep_driver's CBM-Z reassignments. Slot 3 (sulf) "
                 "is excluded: it takes deppart's particle path.",
                 "parameter_overrides": {
                     "ustar": float(v['ustar']), "rmol_in": float(v['rmol']),
                     "z0": float(v['zntt']), "dz1": float(v['dz8w1']), "mech": 2.0,
                     "igas": [float(s) for s in sl],
                     "srfres": [float(v['srfres'][s - 1]) for s in sl]},
                 "time_span": {"start": 0, "end": 1},
                 "tolerance": {"rel": 1e-05},
                 "assertions": [a('vd', j + 1, v['ddvel0d'][s - 1]) for j, s in enumerate(sl)]}
            if i == 1:
                t["assertions"].append({"variable": "aer_res", "time": 0.0,
                                        "expected": float(v['aer_res_def']),
                                        "tolerance": {"rel": 1e-05}})
            dvm['tests'].append(t)
    for path, doc in ((sr_path, sr), (dv_path, dv)):
        with open(path, 'w') as f:
            json.dump(doc, f, indent=2)
            f.write('\n')


if __name__ == '__main__':
    main([int(s) for s in sys.argv[1:]] or [364, 1080, 2101])
