#!/usr/bin/env python3
"""Fill the inline tests of components/gaschem/wrf_chem_interface/chem_met_state.esm.

Reference: data/eqwefic/chem170_scm/run_full_met.sh -- the chem170 full-chemistry column (the
configuration of dumps/chem170_full_mosaic, whose shared dumps it reproduces byte for byte) with the
meteorology hooks and a one-minute history.  The wrfout frame S is the state after step S, the state
chem_driver sees at the end of step S (p + pb reproduces KPP's p_phy exactly).

Inputs:  theta_m = THM + t0, p = P + PB, qv, qr from frame S; alt from the emissions dump of call S
         (chem_driver's grid%alt -- NOT 1/rho of the frame's geopotential, FORTRAN_BUGS N129).
Expected: cbmz_kpp call S (temp, rho_phy, c_m, fix) and wesely call S-2 (= ktau S: rhchem and the
         flags) where that call was dumped.

  python3 tools/fill_chem_met_state_tests.py   (rewrites only the file's `tests`)
"""
import json, sys
import numpy as np
from scipy.io import netcdf_file

D = '/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic'
DUMPS = f'{D}/dumps/chem170_full_met'
WRFOUT = f'{D}/chem170_full/run_met/wrfout_d01_1997-06-20_12:00:00'
ESM = 'components/gaschem/wrf_chem_interface/chem_met_state.esm'
STEPS = [362, 364, 1024, 1026]
REL = 1e-6
T0 = 300.0
EXACT = {'rel': 0.0, 'abs': 0.0}


def frame(v, name, S):
    a = np.asarray(v[name][S], dtype=float)
    return a[..., 0, 0]          # tile column (its, jts), the column every chemistry hook dumps


def dump(s, call):
    try:
        return json.load(open(f'{DUMPS}/esm_dump_{s}_{call}.json'))['vars']
    except FileNotFoundError:
        return None


def arr(x):
    return np.asarray(x['data'], dtype=float).reshape(x['shape'][::-1]).T if x['shape'] else float(x['data'][0])


def fl(a):
    return [float(x) for x in np.asarray(a).ravel()]


def main():
    v = netcdf_file(WRFOUT, 'r', mmap=False).variables
    tests = []
    for S in STEPS:
        k, e, we = dump('cbmz_kpp', S), dump('emissions', S), dump('wesely', S - 2)
        times = v['Times'][S].tobytes().decode()
        po = {'theta_m': fl(frame(v, 'THM', S) + T0), 'p': fl(frame(v, 'P', S) + frame(v, 'PB', S)),
              'alt': fl(arr(e['alt'])), 'qv': fl(frame(v, 'QVAPOR', S)), 'qr': fl(frame(v, 'QRAIN', S)),
              'raincv': 0.0}
        asr = []
        fix = arr(k['fix'])                     # [2, lev]: FIX(indf_H2O), FIX(indf_M)
        for var, ref in (('T', arr(k['temp'])), ('rho', arr(k['rho_phy'])), ('M_air', arr(k['c_m'])),
                         ('H2O', fix[0])):
            for lev, x in enumerate(ref, 1):
                asr.append({'variable': var, 'coords': {'lev': lev}, 'time': 0.0, 'expected': float(x)})
        sfc = ''
        if we is not None:
            asr.append({'variable': 'rhchem', 'time': 0.0, 'expected': float(arr(we['rhchem']))})
            for f in ('wetflag', 'rainflag'):
                asr.append({'variable': f, 'time': 0.0, 'expected': float(bool(we[f]['data'][0])),
                            'tolerance': EXACT})
            sfc = (f"; the wesely dump of call {S - 2} (ktau {S}): rhchem {arr(we['rhchem']):.4f} %, "
                   f"rain flag {int(bool(we['rainflag']['data'][0]))}, dew flag {int(bool(we['wetflag']['data'][0]))}")
        else:
            sfc = f'; wesely_driver was not dumped at ktau {S}, so its outputs are not asserted here'
        tests.append({
            'id': f'chem170_full_met_step{S}',
            'description': (f"WRF-Chem's view of the chem170 column at model step {S} ({times} UTC): inputs from "
                            f"wrfout frame {S} of data/eqwefic/chem170_scm/run_full_met.sh (THM + t0, P + PB, QVAPOR, "
                            f"QRAIN; the state after step {S}) and grid%alt from the emissions dump of call {S}; "
                            f"expected values the cbmz_kpp dump of call {S} (KPP's TEMP, rho_phy, C_M, FIX(indf_H2O))"
                            f"{sfc}. WRF forms them in real32, this component in real64: rel {REL:g}. T is theta_m "
                            f"times the Exner function (FORTRAN_BUGS B34): with chem_theta_is_moist = 0 every T "
                            f"assertion fails. Filled by tools/fill_chem_met_state_tests.py."),
            'parameter_overrides': po,
            'time_span': {'start': 0, 'end': 1},
            'tolerance': {'rel': REL, 'abs': 0.0},
            'assertions': asr})
    doc = json.load(open(ESM))
    doc['models']['WRFChemMetState']['tests'] = tests
    with open(ESM, 'w') as f:
        json.dump(doc, f, indent=2)
        f.write('\n')
    print(sum(len(t['assertions']) for t in tests), 'assertions in', len(tests), 'tests')


if __name__ == '__main__':
    main()
