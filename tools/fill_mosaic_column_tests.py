#!/usr/bin/env python3
"""Fill the inline tests of couplings/mosaic_column.esm (MOSAIC's phase state, liquid
thermodynamics and gas-particle rates over a whole column, lev x bin).

    tools/fill_mosaic_column_tests.py <data/eqwefic> [step ...]

Sources, all from ONE capture (the chem170 ASTEM capture, data/eqwefic/notes/
chem170_mosaic_astem_capture.md), real(kind=8) values:
  * dumps/chem170_mosaic_astem/esm_dump_mosaic_<step+1>.json: the LAST aerosol_phase_state
    call of each bin of the step (ps_*_c2: the recompute at the end of ASTEM_non_volatiles,
    which is the state ASTEM's first semi-volatile sub-step forms its fluxes at -- checked
    bit-identical against the tape's SVB records on all 7552 bins: state, water, totals,
    Kelvin factors), relhumclm;
  * chem170_scm/mosaic_astem/astem.tape.gz: AIN (the cell), BIN (num, dp_wet, kg), the first
    sub-step's SVT gas, SVB state and FXO records (each bin's flux, activities, ions).
Inputs: the cell (T, p, RH, the gases after the non-volatile stage), and per
bin the totals, number, dp_wet, water_a_hyst, and the two equilibrium states the solvers
reached -- ah2o_a (Kelvin) and the solid electrolytes (MESA) -- as the Fortran left them.
Only numbers are written; the document's equations and structure are hand-authored.

What is asserted, per bin: the branch and jaerosolstate; on every bin whose MESA call (if
any) converged, the water and the five Kelvin factors; on the analytic branches the solid
electrolytes; on a bin ASTEM formed a liquid flux for at its own state, the liquid ions,
the H+/NH4+/NO3-/Cl- molalities, gam(HNO3), gam(HCl), gam_ratio and the NH4NO3/NH4Cl
activities ASTEM reads; kg of HNO3/HCl/NH3; and on every level whose eight fluxes were all
formed at one state (tools/fill_mosaic_astem_rates_tests.py `clean`), each bin's rates and
the gas tendencies and balances. Never: a quantity the Fortran left from a state it had
moved off (a mesa_ptc call that hit its 80-pass cap updates the solid after its last water;
compute_activities' stale gam_ratio in the sulfate-poor regime without NH4NO3), or ASTEM's
solver decisions (see the astem rates fill script).
"""
import collections
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / 'couplings' / 'mosaic_column.esm'
sys.path.insert(0, str(ROOT / 'tools'))
from esm_dump import load  # noqa: E402
from esm_compact_json import dumps  # noqa: E402

_spec = importlib.util.spec_from_file_location('far', ROOT / 'tools' / 'fill_mosaic_astem_rates_tests.py')
far = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(far)
fat = far.fat

NLEV, NBIN, NEL = 59, 8, 22
REL = 1e-9
SC = {n: i for i, n in enumerate("water_in hyst_in up_in num jhyst_in state phase jhyst path niter ah2o_a kelvin water "
                                 "water_up mass_dry vol_dry mass_wet vol_wet dpmv sigma growth mdrh rh_pc".split())}
EXACT = {'rel': 0.0, 'abs': 0.0}
ION = ['so4', 'no3', 'cl', 'nh4', 'msa', 'co3', 'na', 'ca', 'oin']
JNH4NO3 = 4                                   # 0-based electrolyte


TEST_TOL = {'rel': REL, 'abs': 0.0}   # the test-level tolerance every assertion inherits unless it says otherwise


def tol(x, floor=0.0):
    """rel 1e-9 (the test's level, so left out), with an absolute floor where given; an exact
    zero is bounded by the inherited abs 0 (esm-spec 6.6.3: rel cannot loosen a zero)."""
    return {'rel': REL, 'abs': floor} if floor else None


def A(var, val, t, **coords):
    a = {'variable': var, 'time': 0.0, 'expected': float(val)}
    if t is not None and t != TEST_TOL:
        a['tolerance'] = t
    a['coords'] = coords
    return a


def build(step, call_of_lev, v, skip_unconverged_state=False):
    """One test: the whole column at tape step `step` (model step step + 1). With
    skip_unconverged_state, jaerosolstate is not asserted on a bin whose mesa_ptc call hit its 80-pass
    cap (outcome 3): there it is the solver's stopping point, which CLAUDE.md says not to assert."""
    po = collections.defaultdict(list)
    asr = []
    census = collections.Counter()
    for k in range(1, NLEV + 1):
        c = call_of_lev[k]
        ain, sb = c['ain'], c['subs'][0]
        po['Column.t_k'].append(ain['t_k'])
        po['Column.p_atm'].append(ain['p_atm'])
        po['Column.rh_pc'].append(ain['rh_pc'])
        po['Column.relhum'].append(float(v['relhumclm'][k - 1]))
        po['Column.gas'].append(list(sb['gas']))
        rows = collections.defaultdict(list)
        for b in range(1, NBIN + 1):
            s = v['ps_sc_c2'][:, b - 1, k - 1]
            S, B = sb['bins'][b], c['bins'][b]
            rows['num'].append(B['num_a'])
            rows['dp_wet'].append(B['dp_wet'])
            rows['aer_t'].append(v['ps_aert_in_c2'][:, b - 1, k - 1].tolist())
            rows['elec_t'].append(v['ps_elet_in_c2'][:NEL, b - 1, k - 1].tolist())
            rows['water_hyst'].append(float(s[SC['hyst_in']]))
            rows['ah2o_a'].append(float(s[SC['ah2o_a']]))
            rows['elec_s_mesa'].append(list(S['E_s']))
        for n, r in rows.items():
            po[f'Column.{n}'].append(r)

        # ---- per bin: phase state, water, Kelvin, partition, activities, kg
        for b in range(1, NBIN + 1):
            s = v['ps_sc_c2'][:, b - 1, k - 1]
            S, B = sb['bins'][b], c['bins'][b]
            path = int(s[SC['path']])
            outcome = (path // 10) % 10 if path >= 6 else 0
            path = 6 if path >= 6 else path
            census[('path', path, outcome)] += 1
            co = dict(lev=k, bin=b)
            asr.append(A('Ps.path', path, EXACT, **co))
            if skip_unconverged_state and path == 6 and outcome == 3:
                census['state_skipped_unconverged'] += 1
            else:
                asr.append(A('Ps.state', S['jaerosolstate'], EXACT, **co))
            asr.append(A('Column.mesa', 1.0 if path == 6 else 0.0, EXACT, **co))
            fx = sb['fxo'].get(b)
            # mesa_ptc's water and Kelvin factors belong to the liquid of its last mesa_flux_salt
            # call; where that call then snapped a salt into the liquid (or, unconverged, the
            # last step moved the solid) the state left behind has another water, which
            # ASTEM_flux_wet recomputes (FXO) while keeping the stale Kelvin factors
            stale = outcome == 3 or (path == 6 and fx is not None and int(fx['route']) in (2, 4)
                                     and fx['water'] != S['water'])
            if stale:
                census['kel_skipped_' + ('unconverged' if outcome == 3 else 'mesa_snap')] += 1
            if not stale:
                asr.append(A('Ps.water', S['water'], tol(S['water']), **co))
                for g in range(5):
                    kv = float(v['ps_kel_c2'][g, b - 1, k - 1])
                    asr.append(A('Ps.kel', kv, tol(kv), lev=k, bin=b, gas=g + 1))
            if path != 6:
                fl = 1e-12 * max(max(abs(x) for x in S['E_t']), 1e-300)
                for e in range(NEL):
                    asr.append(A('Column.elec_s', S['E_s'][e], tol(S['E_s'][e], fl if S['E_s'][e] == 0 else 0.0), lev=k, bin=b, el=e + 1))
            kg = B['kg']
            for gname, i in (('hno3', 1), ('hcl', 2), ('nh3', 3)):
                asr.append(A(f'Rates.kg_{gname}', kg[i], tol(kg[i]), **co))
            own = (fx is not None and fx['ions_s'] == S['ions_s'] and fx['ions_l'] == S['ions_l']
                   and fx['E_s'] == S['E_s'])
            if fx is not None and int(fx['route']) in (2, 4) and own:
                census['liquid_thermo'] += 1
                lmax = max(max(abs(x) for x in fx['ions_l']), 1e-300)
                for i, nm in enumerate(ION):
                    if nm in ('co3', 'oin'):
                        continue
                    x = fx['ions_l'][i]
                    asr.append(A('Column.aer_l', x, tol(x, 1e-12 * lmax if x == 0 else 0.0), lev=k, bin=b, ion=i + 1))
                asr.append(A('Zsr.water', fx['water'], tol(fx['water']), **co))
                mols = [('mc_nh4', fx['mc'][fat.JC_NH4]), ('ma_no3', fx['ma'][fat.JA_NO3]), ('ma_cl', fx['ma'][fat.JA_CL])]
                if fat.wet_case(fx, fx['gas']) != 'case3':   # case 3 replaces mc(jc_h) by its own root
                    mols.insert(0, ('mc_h', fx['mc'][fat.JC_H]))
                else:
                    census['mc_h_case3'] += 1
                for nm, x in mols:
                    asr.append(A(f'Act.{nm}', x, tol(x), **co))
                for e in (18, 19):
                    asr.append(A('Act.gam', fx['gam'][e - 1], tol(fx['gam'][e - 1]), lev=k, bin=b, el=e))
                for e in (5, 6):
                    asr.append(A('Act.activity', fx['activity'][e - 1], tol(fx['activity'][e - 1]), lev=k, bin=b, el=e))
                poor = fat.xt(fx['ions_l']) > 2 or fat.xt(fx['ions_l']) < 0
                if not poor or fx['E_l'][JNH4NO3] > 0:
                    asr.append(A('Act.gam_ratio', fx['gam_ratio'], tol(fx['gam_ratio']), **co))
                else:
                    census['gam_ratio_stale'] += 1

        # ---- the rates, where all eight fluxes were formed at one state
        if not far.clean(sb):
            census['level_edited_between_bins'] += 1
            continue
        bad = None
        for b in sorted(sb['fxo']):
            _, why = far.bin_branch(c, sb, b)
            if why:
                bad = why
        unconv = [b for b in range(1, NBIN + 1) if int(v['ps_sc_c2'][SC['path'], b - 1, k - 1]) >= 6
                  and (int(v['ps_sc_c2'][SC['path'], b - 1, k - 1]) // 10 % 10 == 3
                       or (int(sb['fxo'][b]['route']) in (2, 4)
                           and sb['fxo'][b]['water'] != sb['bins'][b]['water']))]
        if bad or unconv:
            census['level_excluded:' + (bad or 'mesa_ptc water stale (unconverged or snapped)')] += 1
            continue
        census['levels_with_rates'] += 1
        g = sb['gas']
        sums = {n: [0.0, 0.0] for n in far.GASES}
        for b in range(1, NBIN + 1):
            fx, kg = sb['fxo'][b], c['bins'][b]['kg']
            for gname, i in (('h2so4', 0), ('msa', 4)):
                x = g[i] * kg[i]
                asr.append(A('Rates.daer', x, {'rel': REL, 'abs': 1e-12 * abs(x)}, lev=k, bin=b, gas=i + 1))
                sums[gname][0] += x
                sums[gname][1] += 1e-12 * abs(x)
            for gname, (x, fl) in far.semi_expected(c, b, fx).items():
                asr.append(A('Rates.daer', x, {'rel': REL, 'abs': fl}, lev=k, bin=b, gas=far.GAX[gname]))
                sums[gname][0] += x
                sums[gname][1] += fl
        for gname in far.GASES:
            x, fl = sums[gname]
            asr.append(A('Rates.dgas', -x, {'rel': REL, 'abs': fl}, lev=k, gas=far.GAX[gname]))
            asr.append(A('Rates.balance', 0.0, {'rel': 0.0, 'abs': max(fl, 1e-300) * 1e-3}, lev=k, gas=far.GAX[gname]))
    return dict(po), asr, census


def main():
    d = Path(sys.argv[1])
    steps = [int(x) for x in sys.argv[2:]] or [0, 917]
    calls = far.load(d / 'chem170_scm' / 'mosaic_astem')[0]
    by = collections.defaultdict(dict)
    for c in calls:
        by[c['step']][c['lev']] = c
    tests = []
    for st in steps:
        v = load(d / 'dumps' / 'chem170_mosaic_astem' / f'esm_dump_mosaic_{st + 1}.json')['vars']
        po, asr, census = build(st, by[st], v)
        print(f'step {st}: {len(asr)} assertions;', dict(sorted(census.items(), key=str)))
        desc = (f"The whole chem_opt = 170 column (59 levels x 8 sections), model step {st + 1}, at the state ASTEM's first "
                "semi-volatile sub-step forms its fluxes at (after ASTEM_non_volatiles and its phase-state recompute). "
                "Inputs: the cell (T, p, RH, the five gases after the non-volatile stage) from the ASTEM tape (AIN, SVT) and "
                "the dump (relhumclm), the cell's equilibrium constants -- Keq_ll(1) and Kp0 included -- computed from "
                "them; per bin the totals, number, dp_wet and water_a_hyst, and the two equilibrium states the "
                "Fortran's solvers stopped at, ah2o_a (Kelvin) and the solid electrolytes (MESA, read on path 6 only), "
                "from the ASTEM capture's last aerosol_phase_state call (dumps/chem170_mosaic_astem, ps_*_c2, bit-identical "
                "to the tape's SVB records). Expected, real64: "
                f"census {dict(sorted(((str(k), n) for k, n in census.items())))}. "
                "Per bin the branch and jaerosolstate (exactly), the water and the Kelvin factors (except where mesa_ptc hit "
                "its 80-pass cap: it updates the solid after its last water, so the Fortran's water is not the state's), "
                "the solid electrolytes on the analytic branches, kg; on a bin ASTEM formed a liquid flux for at its own "
                "state the liquid ions, molalities, gam(HNO3), gam(HCl), gam_ratio and the NH4NO3/NH4Cl activities; and on "
                "a level whose eight fluxes were all formed at one state, each bin's rates and the gas tendencies and "
                "balances (as couplings/mosaic_astem_rates.esm, now with the phase state and activities COMPUTED). "
                "rel 1e-9, flags exactly.")
        tests.append({'id': f'column_s{st + 1}', 'description': desc, 'parameter_overrides': po,
                      'time_span': {'start': 0, 'end': 0}, 'tolerance': TEST_TOL, 'assertions': asr})
    doc = json.loads(DOC.read_text())
    doc['models']['Column']['tests'] = tests
    DOC.write_text(dumps(doc) + '\n')
    print(f'{DOC.name}: {len(tests)} tests, {sum(len(t["assertions"]) for t in tests)} assertions')


if __name__ == '__main__':
    main()
