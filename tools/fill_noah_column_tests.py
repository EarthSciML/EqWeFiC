#!/usr/bin/env python3
"""Fill the NoahColumn tests of couplings/noah_column.esm (numbers only; no equations).

One test per chemistry step S of the full chem170 column (362, 364, 1024, 1026), all at
t = 0 with an EMPTY time_span (the document has ODE states and every assertion is at the
initial state). Inputs: the noah dump of step S and the sfclayrev dump of call 2S-1 (row 1
of ktau S) from data/eqwefic/dumps/chem170_full_met_sfc. References: the real64 replays
written by tools/noah_column_replay.py (run that first):
  sfc64x_S.json        sfclayrev, exact-decimal constants   -> the Sfc assertions
  noah64_S_1e-6.json   SFLX at ESM_DT = 1e-6 s              -> the Noah fluxes and diagnostics
  noah64_S_0.01.json   SFLX at ESM_DT = 0.01 s              -> the soil and canopy tendencies
    python3 tools/fill_noah_column_tests.py [--measure]
--measure writes every assertion at rel 1e-14 / abs 0 so a run reports the raw residuals.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import esm_dump  # noqa: E402

D = '/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic'
DUMPS = f'{D}/dumps/chem170_full_met_sfc'
R = f'{D}/noah_column'
DOC = os.path.join(os.path.dirname(__file__), '..', 'couplings', 'noah_column.esm')
STEPS = (362, 364, 1024, 1026)
MEASURE = '--measure' in sys.argv

# Sfc output -> sfclayrev_driver name, tolerance
SFC = {'ust': 'ust_out', 'ustm': 'ustm_out', 'rmol': 'rmol_out', 'zol': 'zol_out', 'mol': 'mol_out',
       'br': 'br_out', 'wspd': 'wspd_out', 'psim': 'psim_out', 'psih': 'psih_out', 'fm': 'fm_out',
       'fh': 'fh_out', 'gz1oz0': 'gz1oz0_out', 'cpm': 'cpm_out', 'chs': 'chs_out', 'chs2': 'chs2_out',
       'cqs2': 'cqs2_out', 'qgh': 'qgh_out', 'flhc': 'flhc_out', 'flqc': 'flqc_out', 'hfx': 'hfx_out',
       'qfx': 'qfx_out', 'lh': 'lh', 'qsfc': 'qsfc_out', 'znt': 'znt_out', 'u10': 'u10', 'v10': 'v10',
       'th2': 'th2', 't2': 't2', 'q2': 'q2'}
# Noah mount.variable -> noah_driver name (ESM_DT = 1e-6 s)
NOAH = {'Heat.t1': 't1_out', 'Heat.sheat': 'sheat', 'Heat.ssoil': 'ssoil', 'Heat.ssoil_pre': 'ssoil_pre',
        'Heat.yy': 'yy', 'Heat.zz1': 'zz1', 'Pen.rch': 'rch', 'Pen.rr': 'rr', 'Pen.epsca': 'epsca',
        'Pen.etp_energy': 'etp', 'Can.rcs': 'rcs', 'Can.rct': 'rct', 'Can.rcq': 'rcq', 'Can.rcsoil': 'rcsoil',
        'Can.rc': 'rc', 'Can.pc': 'pc', 'Et.edir': 'edir', 'Et.ec': 'ec', 'Et.ett': 'ett', 'Et.eta': 'eta',
        'Et.eta_kinematic': 'eta_kinematic', 'Et.beta': 'beta', 'Drv.fdown': 'fdown', 'Drv.q1': 'q1',
        'Drv.albedo': 'albedo_out', 'Drv.emissi': 'emissi_out', 'Drv.z0': 'z0_out'}
# Tolerances (measured with --measure on esm 00808ba51; data/eqwefic/noah_column/run_measure_sfclit.log).
# Against the LITERAL-EXACT surface-layer replay (tools/noah_column_replay.py) every Sfc output agrees to
# <= 9e-13 relative and every Noah quantity to <= 3.1e-10 (sheat), so everything is at rel 1e-9.
# (Until 2026-11-20 the reference kernel kept the CCPP source's default-real literals .33/.32/0.001 in
# binary32, which put br/zol/rmol 2.2e-5 and mol/hfx 1.9e-4 off by night and held these at rel 5e-5;
# that run is data/eqwefic/noah_column/run_measure.log.)
CHS_FREE = {'Heat.ssoil_pre', 'Can.rcs', 'Can.rct', 'Can.rcq', 'Can.rcsoil', 'Can.rc', 'Drv.albedo', 'Drv.emissi', 'Drv.z0'}
# tendencies (ESM_DT = 0.01 s): document variable -> (state in, state out)
TEND = {'dstc_dt': ('stc', 'stc_out'), 'dsh2o_dt': ('sh2o', 'sh2o_out'), 'dsmc_dt': ('smc', 'smc_out')}


def v(d, k):
    return d[k]


def f(x):
    x = float(np.asarray(x).reshape(-1)[0])
    return 0.0 if np.isnan(x) else x


def const_col(vals, axis='soil'):
    return {'op': 'faq', 'output_idx': ['k'], 'args': [], 'ranges': {'k': {'from': axis}},
            'expr': {'op': 'index', 'args': [{'op': 'const', 'args': [], 'value': [float(x) for x in vals]}, 'k']}}


def tol(t):
    return {'rel': 1e-14, 'abs': 0.0} if MEASURE else t


def main() -> None:
    doc = json.load(open(DOC))
    tests = []
    for S in STEPS:
        n = esm_dump.column(esm_dump.load(f'{DUMPS}/esm_dump_noah_{S}.json'))
        s = esm_dump.column(esm_dump.load(f'{DUMPS}/esm_dump_sfclayrev_{2 * S - 1}.json'))
        sx = esm_dump.column(esm_dump.load(f'{R}/sfc64x_{S}.json'))
        s32 = esm_dump.column(esm_dump.load(f'{R}/sfc64_{S}.json'))
        nf = esm_dump.column(esm_dump.load(f'{R}/noah64_{S}_1e-6.json'))
        nt = esm_dump.column(esm_dump.load(f'{R}/noah64_{S}_0.01.json'))
        # the same column, read from the two dumps
        assert f(n['t3d1']) == f(s['t1d']) and f(n['tsk']) == f(s['tsk']) and f(n['chs']) == f(s['chs_out'])
        assert f(n['p8w1']) == f(s['psfcpa']) and f(n['dz8w1']) == f(s['dz8w1d'])
        assert int(n['nroot']) == 3 and np.allclose(n['sldpth'], [0.1, 0.3, 0.6, 1.0])
        po = {
            'shdfac': f(n['shdfac']), 'soiltyp': float(n['soiltyp']), 'vegtyp': float(n['vegtyp']),
            'slopetyp': float(n['slopetyp']), 'shdmin': f(n['shdmin']), 'shdmax': f(n['shdmax']),
            'tbot': f(n['tbot']), 'z_deep': 8.0, 'z_root': float(np.sum(np.asarray(n['sldpth'], float)[:3])), 'dt_canopy': 60.0,
            'tsk_in': f(n['tsk']), 'swdown': f(n['swdown']), 'glw': f(n['glw']),
            'albedo_in': f(n['albedo_out']), 'emiss_in': f(n['emissi']), 'rainbl': f(n['rainbl']),
            'dt_rain': f(n['dt']), 'sneqv': f(n['sneqv']), 'ffrozp': f(n['ffrozp']),
            'u1': f(s['ux']), 'v1': f(s['vx']), 't_air': f(n['t3d1']), 'qv_air': f(n['qv3d1']),
            'p_air': f(s['p1d']), 'p_sfc': f(n['p8w1']), 'p_w2': 2.0 * f(n['sfcprs']) - f(n['p8w1']),
            'dz1': f(s['dz8w1d']), 'mavail': f(s['mavail']), 'xland': f(s['xland']),
            'lakemask': f(s['lakemask']), 'water_depth': f(s['water_depth']), 'dx': f(s['dx']),
            'isfflx': float(bool(s['isfflx'])), 'shalwater_z0': float(bool(s['shalwater_z0'])),
            'scm_force_flux': float(bool(s['scm_force_flux'])), 'isftcflx': float(s['isftcflx']),
            'iz0tlnd': float(s['iz0tlnd']), 'pblh_prev': f(s['pblh']), 'znt_prev': f(s['znt']),
            'ust_prev': f(s['ust']), 'mol_prev': f(s['mol']), 'qsfc_prev': f(s['qsfc']),
            'hfx_prev': f(s['hfx']), 'qfx_prev': f(s['qfx']), 'ustm_prev': f(s['ustm']),
            'chs_prev': f(s['chs']), 'chs2_prev': f(s['chs2']), 'cqs2_prev': f(s['cqs2'])}
        ic = {'stc': [float(x) for x in n['stc']], 'smc': [float(x) for x in n['smc']],
              'sh2o': [float(x) for x in n['sh2o']], 'cmc': f(n['cmc'])}
        A = []
        for name, key in SFC.items():
            t = {'rel': 1e-9}
            A.append({'variable': f'Sfc.{name}', 'time': 0.0, 'expected': f(sx[key]), 'tolerance': tol(t)})
        A.append({'variable': 'Sfc.regime', 'time': 0.0, 'expected': f(sx['regime_out']), 'tolerance': {'abs': 0.0, 'rel': 0.0}})
        for name, key in NOAH.items():
            t = {'rel': 1e-9} if name in CHS_FREE else {'rel': 1e-9, 'abs': 1e-12}
            A.append({'variable': name, 'time': 0.0, 'expected': f(nf[key]), 'tolerance': tol(t)})
        A.append({'variable': 'Et.et', 'time': 0.0, 'reduce': 'Linf_error', 'expected': 0.0,
                  'reference': const_col(nf['et']), 'tolerance': tol({'abs': 1e-10, 'rel': 0.0})})
        h = 0.01
        tends = {}
        for var, (a, b) in TEND.items():
            tends[var] = (np.asarray(nt[b], float) - np.asarray(n[a], float)) / h
        tends['dcmc_dt'] = np.array([(f(nt['cmc_out']) - f(n['cmc'])) / h])
        TT = {'dstc_dt': 1e-9, 'dsh2o_dt': 3e-11, 'dsmc_dt': 3e-11}
        for var in ('dstc_dt', 'dsh2o_dt', 'dsmc_dt'):
            A.append({'variable': var, 'time': 0.0, 'reduce': 'Linf_error', 'expected': 0.0,
                      'reference': const_col(tends[var]), 'tolerance': tol({'abs': TT[var], 'rel': 0.0})})
        A.append({'variable': 'dcmc_dt', 'time': 0.0, 'expected': float(tends['dcmc_dt'][0]), 'tolerance': tol({'abs': 3e-14, 'rel': 0.0})})
        for name in ('Drv.snow_branch_active', 'Sm.supersaturation_active'):
            A.append({'variable': name, 'time': 0.0, 'expected': 0.0, 'tolerance': {'abs': 0.0, 'rel': 0.0}})
        lag = f(n['tsk']) - f(nf['t1_out'])
        c32 = abs(f(s32['chs_out']) - f(sx['chs_out'])) / f(sx['chs_out'])
        desc = (f"Step {S} of the FULL chem_opt = 170 column (em_scm_xy Kansas 1997-06-20, cu_physics = 0; "
                f"data/eqwefic/chem170_scm/run_full_met_sfc.sh, the same run as dumps/chem170_full_met: every noah dump "
                f"byte-identical), the step whose end state the chemistry sees. Inputs: the noah dump of step {S} and the "
                f"sfclayrev dump of call {2 * S - 1} (row j = 1 of ktau {S}; the surface layer is dumped per j-row call), "
                f"real32 values widened; the soil column (stc, smc, sh2o, cmc) is the initial condition. References "
                f"(tools/noah_column_replay.py): the surface layer is the real64 kernels/sfclayrev_driver with WRF's "
                f"thirteen real32 constant arguments replaced by their exact decimal values (with WRF's real32 copies "
                f"chs moves by {c32:.1e} relative here), built literal-exact (kernels/build_sfclit: the CCPP source's "
                f"default-real literals .33, .32 and 0.001 promoted to real64, as sfclayrev.esm writes them); Noah is the real64 kernels/noah_driver fed the SAME forcing this "
                f"document forms -- the driver preamble rebuilt in real64 from the raw inputs and the exact-constant "
                f"chs and qgh -- at ESM_DT = 1e-6 s for the fluxes and diagnostics (instantaneous: the soil state and the "
                f"canopy are the initial ones) and at 0.01 s for the four tendencies (x_out - x)/dt (the minimum of the "
                f"V: 1e-6 s loses the tendency to cancellation). The incoming skin temperature tsk_in is WRF's lagged "
                f"TSK; it differs from the skin temperature this state diagnoses (Heat.t1) by {lag:+.4f} K, which is the "
                f"size of the algebraic loop WRF breaks by lagging. Snow-free and unfrozen: the snow branch and the "
                f"supersaturation projection are asserted inactive at zero tolerance.")
        tests.append({'id': f'chem170_full_step{S}', 'description': desc, 'parameter_overrides': po,
                      'initial_conditions': ic, 'time_span': {'start': 0, 'end': 0}, 'assertions': A})
    doc['models']['NoahColumn']['tests'] = tests
    with open(DOC, 'w') as fh:
        fh.write(json.dumps(doc, indent=1, ensure_ascii=False) + '\n')
    print(f'filled {len(tests)} tests, {sum(len(t["assertions"]) for t in tests)} assertions')


if __name__ == '__main__':
    main()
