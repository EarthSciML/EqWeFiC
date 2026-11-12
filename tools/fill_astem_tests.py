#!/usr/bin/env python3
"""Fill the inline tests of MOSAIC's ASTEM physics components from the ASTEM
tape (data/eqwefic/chem170_scm/mosaic_astem/astem.tape.gz, written by the
esm_at_* hooks of the WRF fork branch earthsciml-instrumented-chem170 inside
the chem_opt = 170 column run; see data/eqwefic/notes/chem170_mosaic_astem_capture.md).

    tools/fill_astem_tests.py <astem_dir>

Writes the `tests` of
  components/aerosol/mosaic/mosaic_gas_equilibrium_constants.esm
  components/aerosol/mosaic/astem_liquid_surface_equilibrium.esm
  components/aerosol/mosaic/gas_particle_mass_transfer.esm  (real64 tests added
      beside the existing real32 ones, ids prefixed 'astem64_')
Only numbers are written: equations and structure are hand-authored.

Every reference is a real(kind=8) value ASTEM itself held, so tolerances are
rel 1e-9.  The surface equilibria are asserted at EVERY flux call the tape holds
(every bin, every ASTEM sub-step), not only the first: each is an instantaneous
evaluation at a state the Fortran visited, which is what the physics is.  The
sub-step at which ASTEM stops is never asserted (EqWeFiC CLAUDE.md, solvers).
"""
import json, sys, os, math, collections
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MOS = ROOT / 'components' / 'aerosol' / 'mosaic'
TOL = {'rel': 1e-9}

IH = dict(hno3=1, hcl=2, nh3=3)          # 0-based positions in the 5-gas vectors
JNH4NO3, JNH4CL, JHNO3, JHCL, JCACO3 = 4, 5, 17, 18, 21   # 0-based electrolyte
JC_H, JC_NH4 = 0, 1
JA_NO3, JA_CL = 2, 3
ION = dict(so4=0, no3=1, cl=2, nh4=3, msa=4, co3=5, na=6, ca=7, oin=8)


def a(var, val, tol=None, **coords):
    d = {'variable': var, 'time': 0.0, 'expected': float(val)}
    if coords:
        d['coords'] = coords
    if tol:
        d['tolerance'] = tol
    return d


HAND_NAMES = {}


def load(dirpath, hand=False):
    sys.path.insert(0, str(dirpath))
    from astem_tape import read_tape
    if hand:
        for line in open(Path(dirpath) / 'cases_hand_names.txt'):
            i, name, why = line.rstrip('\n').split('\t')
            HAND_NAMES[int(i)] = (name, why)
        return read_tape(Path(dirpath) / 'hand.tape')
    return read_tape(Path(dirpath) / 'astem.tape.gz')


def where(c, ibin, sb):
    if c['step'] < 0:
        name, why = HAND_NAMES[c['lev']]
        return (f"hand_{name}", f"Hand-specified bin '{name}' ({why}; chem170_scm/mosaic_astem/"
                "make_astem_hand.py), run through WRF's own compiled ASTEM_semi_volatiles by "
                "kernels/mosaic_therm_driver mode 6 at T = 290 K, p = 0.9 atm, 1000 particles/cm3")
    return (f"s{c['step']}_k{c['lev']}_b{ibin}_sub{sb['isteps']}",
            f"ASTEM tape (real64), chem_opt = 170 column, model step {c['step'] + 1}, level {c['lev']}, "
            f"bin {ibin}, ASTEM sub-step {sb['isteps']}")


def write_tests(fname, model, tests, keep=None):
    p = MOS / fname
    doc = json.loads(p.read_text())
    m = doc['models'][model]
    old = m.get('tests', [])
    m['tests'] = ([t for t in old if keep(t)] if keep else []) + tests
    p.write_text(json.dumps(doc, indent=2) + '\n')
    print(f'{fname}: {len(tests)} tests, '
          f'{sum(len(t["assertions"]) for t in m["tests"])} assertions')


def relhumclm(dirpath, step, lev, _cache={}):
    """The cell's relative humidity update_thermodynamic_constants ran at (relhumclm, N118),
    from the dump of the same capture (dumps/chem170_mosaic_astem, beside the tape's run dir)."""
    if step < 0:
        raise ValueError('hand calls carry no dump; use cell_relhum')
    if step not in _cache:
        sys.path.insert(0, str(ROOT / 'tools'))
        from esm_dump import load as dload
        dump = Path(dirpath).parents[1] / 'dumps' / 'chem170_mosaic_astem' / f'esm_dump_mosaic_{step + 1}.json'
        _cache[step] = dload(str(dump))['vars']['relhumclm']
    return float(_cache[step][lev - 1])


def cell_relhum(dirpath, c):
    """relhumclm of a column call (from the dump), or of a hand call: the driver holds it as a
    default REAL, so it is the binary32 of the case's RH, recovered from the taped rh_pc =
    100*relhumclm."""
    if c['step'] >= 0:
        return relhumclm(dirpath, c['step'], c['lev'])
    import numpy as np
    return float(np.float32(np.float32(c['ain']['rh_pc']) / np.float32(100.0)))


def constants_tests(calls, pick, dirpath):
    """One test per (step, level) picked: t_k and the cell's relative humidity in,
    keq_gl/keq_ll/kp/Kp0 out."""
    tests = []
    for c in pick:
        ain = c['ain']
        rh = cell_relhum(dirpath, c)
        if c['step'] >= 0:
            tid = f"astem64_s{c['step']}_k{c['lev']}"
            where_ = (f"ASTEM tape, chem_opt = 170 column, model step {c['step'] + 1}, level {c['lev']}: "
                      f"T = {ain['t_k']:.2f} K, relative humidity {rh:.4f} (relhumclm from dumps/chem170_mosaic_astem).")
        else:
            name = HAND_NAMES[c['lev']][0]
            tid = f"astem64_hand_{name}"
            where_ = (f"Hand-specified case '{name}' (chem170_scm/mosaic_astem/hand.tape, WRF's own compiled "
                      f"update_thermodynamic_constants via kernels/mosaic_therm_driver mode 6): T = {ain['t_k']:.2f} K, "
                      f"relative humidity {rh:.8f} (the driver's binary32 relhumclm).")
        tests.append({
            'id': tid,
            'description': (where_ + " The real64 constants ASTEM held after "
                            "update_thermodynamic_constants (AIN record). rel 1e-9."),
            'parameter_overrides': {'t_k': ain['t_k'], 'a_w': rh},
            'time_span': {'start': 0, 'end': 1},
            'tolerance': TOL,
            'assertions': [a('keq_nh3_gl', ain['keq_gl'][1]), a('keq_hno3_gl', ain['keq_gl'][2]),
                           a('keq_hcl_gl', ain['keq_gl'][3]), a('keq_hso4', ain['keq_ll'][0]),
                           a('keq_nh3_ll', ain['keq_ll'][1]),
                           a('k_w', ain['keq_ll'][2]), a('kp_nh4no3', ain['kp_nh4no3']),
                           a('kp_nh4cl', ain['kp_nh4cl']), a('kp_nh4no3_0', ain['kp_nh4no3_0']),
                           a('kp_nh4cl_0', ain['kp_nh4cl_0'])]})
    return tests


MW = dict(h2so4=98.0, hno3=63.0, hcl=36.5, nh3=17.0, msa=96.0)


def mtc_tests(calls, pick):
    """Real64 kg from the tape, the inputs aerosolmtc read (dp_wet, num_a, T, p)."""
    from_par = json.loads((MOS / 'mosaic_mass_transfer_parameters.esm').read_text())
    pv = from_par['models']['MosaicMassTransferParameters']['variables']
    gas_par = {g: (pv[f'mw_{g}']['default'], pv[f'vm_{g}']['default'])
               for g in ['h2so4', 'msa', 'hno3', 'hcl', 'nh3']}
    gidx = dict(h2so4=0, hno3=1, hcl=2, nh3=3, msa=4)
    tests = []
    for c, ibin, g in pick:
        b = c['bins'][ibin]
        mw, vm = gas_par[g]
        tests.append({
            'id': f"astem64_s{c['step']}_k{c['lev']}_b{ibin}_{g}",
            'description': (f"ASTEM tape (real64): chem_opt = 170 column, model step {c['step'] + 1}, "
                            f"level {c['lev']}, bin {ibin}, gas {g.upper()}; T = {c['ain']['t_k']:.2f} K, "
                            f"wet diameter {b['dp_wet']:.3e} cm, number {b['num_a']:.4g} /cm3. kg is the "
                            "real(kind=8) value aerosolmtc formed, read with its own real64 inputs, so the "
                            "bound is rel 1e-9 (the older tests above carry the real32 dump's 5e-7)."),
            'parameter_overrides': {'t_k': c['ain']['t_k'], 'p_atm': c['ain']['p_atm'],
                                    'dp_wet': b['dp_wet'], 'num_conc': b['num_a'],
                                    'mw_gas': mw, 'v_molar': vm, 'accom': pv['accom']['default'],
                                    'gas_conc': c['ain']['gas'][gidx[g]]},
            'time_span': {'start': 0, 'end': 1},
            'tolerance': TOL,
            'assertions': [a('kg', b['kg'][gidx[g]])]})
    return tests


def xt(ions):
    so4, msa = ions[ION['so4']], ions[ION['msa']]
    den = so4 + 0.5 * msa
    if den > 0:
        return (ions[ION['nh4']] + ions[ION['na']] + 2.0 * ions[ION['ca']]) / den
    return -1.0


XT_RICH = 1.999899983406067
NH3_NONE = 1.0000000195414814e-25
SERIES_X = 9.999999974752427e-07


def wet_case(fx, gas):
    """The branch ASTEM_flux_wet took, decided from the taped state (route 2/4)."""
    if fx['E_s'][JCACO3] > 0:
        return 'case1'
    x = xt(fx['ions_l'])
    if 0.0 <= x < XT_RICH:
        return 'case2'
    if gas[IH['nh3']] + fx['ions_l'][ION['nh4']] < NH3_NONE:
        return 'case4'
    return 'case3'


def wet_tests(calls, max_per=12):
    """Liquid-bin surface equilibria, every flux call of every sub-step.  Census
    first, then an even spread over steps/levels/bins within each branch."""
    groups = collections.defaultdict(list)
    for c in calls:
        ain = c['ain']
        for sb in c['subs']:
            for ibin, fx in sb['fxo'].items():
                if fx['route'] not in (2, 4):
                    continue
                g = fx['gas']
                case = wet_case(fx, g)
                kg = c['bins'][ibin]['kg']
                sub = case
                if case == 'case3':
                    kel = fx['kel']
                    knh4no3 = kel[3] * kel[1] * fx['activity'][JNH4NO3] * ain['kp_nh4no3']
                    knh4cl = kel[3] * kel[2] * fx['activity'][JNH4CL] * ain['kp_nh4cl']
                    qa = kg[3]
                    qb = -kg[3] * g[3] + kg[1] * g[1] + kg[2] * g[2]
                    qc = -(kg[1] * knh4no3 + kg[2] * knh4cl)
                    qx = 1e6 if qb == 0 else 4.0 * (qa / qb) * (qc / qb)
                    s_no3, s_cl = fx['sfc'][1], fx['sfc'][2]
                    # equilibrate_acids is reached when, after the tiny absorptions,
                    # the solution holds neither nitrate nor chloride: decide it from
                    # the state BEFORE the flux call (SVB ions, FXI gas), because the
                    # projection itself writes nitrate/chloride into the FXO record
                    pre, g0 = sb['bins'][ibin]['ions_l'], sb['fxi'][ibin]['gas']
                    no3_after = pre[ION['no3']] > 0 or (g0[3] * g0[1] > 0 and pre[ION['nh4']] * pre[ION['no3']] == 0)
                    cl_after = pre[ION['cl']] > 0 or (g0[3] * g0[2] > 0 and pre[ION['nh4']] * pre[ION['cl']] == 0)
                    if not (no3_after or cl_after):
                        hsrc = 'equilibrate_acids'
                    elif s_no3 > 0 and fx['ma'][JA_NO3] > 0:
                        hsrc = 'h_no3'
                    elif s_cl > 0 and fx['ma'][JA_CL] > 0:
                        hsrc = 'h_cl'
                    else:
                        hsrc = 'equilibrate_acids'
                    sub = f"case3_{'series' if abs(qx) < SERIES_X else 'closed'}_{hsrc}"
                elif case == 'case4':
                    num = fx['kel'][1] * ain['keq_gl'][3] * fx['ma'][JA_NO3] * fx['gam'][JHNO3] ** 2
                    den = fx['kel'][2] * ain['keq_gl'][2] * fx['ma'][JA_CL] * fx['gam'][JHCL] ** 2
                    sub = 'case4_' + ('off' if num == 0 or den == 0 else 'exchange')
                elif case == 'case1':
                    sub = 'case1'
                elif case == 'case2':
                    sub = ('case2' + ('_no3' if fx['ma'][JA_NO3] > 0 else '')
                           + ('_cl' if fx['ma'][JA_CL] > 0 else '')
                           + ('_nh4' if fx['mc'][JC_NH4] > 0 else ''))
                groups[sub].append((c, sb, ibin, fx))
    print('wet census:', {k: len(v) for k, v in sorted(groups.items())})
    tests = []
    for sub, rows in sorted(groups.items()):
        if sub.endswith('equilibrate_acids'):
            continue
        n = len(rows)
        idx = sorted({round(i * (n - 1) / max(max_per - 1, 1)) for i in range(min(max_per, n))})
        for i in idx:
            tests.append(wet_test(sub, *rows[i]))
    return tests, {k: len(v) for k, v in groups.items()}


def wet_test(sub, c, sb, ibin, fx):
    ain, b = c['ain'], c['bins'][ibin]
    g, kg = fx['gas'], b['kg']
    il = fx['ions_l']
    over = {
        'gas_hno3': g[1], 'gas_hcl': g[2], 'gas_nh3': g[3],
        'kg_hno3': kg[1], 'kg_hcl': kg[2], 'kg_nh3': kg[3],
        'kel_hno3': fx['kel'][1], 'kel_hcl': fx['kel'][2], 'kel_nh3': fx['kel'][3],
        'so4_l': il[ION['so4']], 'msa_l': il[ION['msa']], 'nh4_l': il[ION['nh4']],
        'na_l': il[ION['na']], 'ca_l': il[ION['ca']], 'caco3_s': fx['E_s'][JCACO3],
        'mc_h_act': fx['mc'][JC_H], 'mc_nh4': fx['mc'][JC_NH4],
        'ma_no3': fx['ma'][JA_NO3], 'ma_cl': fx['ma'][JA_CL],
        'gam_hno3': fx['gam'][JHNO3], 'gam_hcl': fx['gam'][JHCL], 'gam_ratio': fx['gam_ratio'],
        'act_nh4no3': fx['activity'][JNH4NO3], 'act_nh4cl': fx['activity'][JNH4CL],
        'keq_nh3_gl': ain['keq_gl'][1], 'keq_hno3_gl': ain['keq_gl'][2], 'keq_hcl_gl': ain['keq_gl'][3],
        'keq_nh3_ll': ain['keq_ll'][1], 'k_w': ain['keq_ll'][2],
        'kp_nh4no3': ain['kp_nh4no3'], 'kp_nh4cl': ain['kp_nh4cl']}
    asr = []
    flag = sub[:5]
    asr.append(a(flag, 1.0, {'abs': 0.0, 'rel': 0.0}))
    if sub == 'case1':
        for gname, i in (('hno3', 1), ('hcl', 2)):
            if g[i] > 9.999999747378752e-06:
                asr += [a(f'sfc_{gname}', fx['sfc'][i], {'abs': 0.0, 'rel': 0.0}),
                        a(f'df_{gname}', fx['df_s'][i], TOL),
                        a(f'cond_{gname}', fx['flux_s'][i], TOL)]
        asr.append(a('mc_h', fx['mc'][JC_H], TOL))
    else:
        exch = (1, 2, 3) if sub[:5] in ('case2', 'case3') else ((1, 2) if sub == 'case4_exchange' else ())
        for i in exch:
            gname = ['h2so4', 'hno3', 'hcl', 'nh3'][i]
            floor = 1e-12 * max(abs(g[i]), abs(fx['sfc'][i]), 1e-30)
            asr += [a(f'sfc_{gname}', fx['sfc'][i], {'rel': 1e-9, 'abs': floor}),
                    a(f'df_{gname}', fx['df_l'][i], {'rel': 1e-9, 'abs': floor}),
                    a(f'cond_{gname}', kg[i] * fx['df_l'][i], {'rel': 1e-9, 'abs': kg[i] * floor})]
        asr.append(a('mc_h', fx['mc'][JC_H], TOL))
    return {
        'id': f"astem64_{sub}_{where(c, ibin, sb)[0]}",
        'description': (f"Branch {sub}. {where(c, ibin, sb)[1]} (route {fx['route']}: "
                        f"{'all liquid' if fx['route'] == 2 else 'mixed, no solid NH4NO3/NH4Cl'}). Inputs are the "
                        "state the bin's ASTEM_flux_wet saw after its tiny-absorption adjustments and "
                        "compute_activities (FXO record), kg and the constants from the call's AIN/BIN records; "
                        "expected are the surface concentrations sfc_a, driving forces df_gas_l (df_gas_s for "
                        "CaCO3) and H+ it formed. cond_x is kg x df_gas, the product of two taped values. rel 1e-9, "
                        "absolute floor 1e-12 of the larger of gas and surface value."),
        'parameter_overrides': over,
        'time_span': {'start': 0, 'end': 1},
        'assertions': asr}


def spread(rows, n):
    if len(rows) <= n:
        return rows
    return [rows[round(i * (len(rows) - 1) / (n - 1))] for i in range(n)]



RTOL_EQB = 0.009999999776482582
JNACL, JCACL2 = 11, 13


def dry_rows(calls):
    groups = collections.defaultdict(list)
    skipped = collections.Counter()
    for c in calls:
        ain = c['ain']
        for sb in c['subs']:
            for ibin, fx in sb['fxo'].items():
                if fx['route'] != 1:
                    continue
                g, E = fx['gas'], fx['E_s']
                x = xt(fx['ions_s'])
                k0no3 = min(ain['kp_nh4no3_0'], ain['keq_sg'][0])
                k0cl = min(ain['kp_nh4cl_0'], ain['keq_sg'][1])
                if E[JCACO3] > 0:
                    sub = 'case1'
                elif 0.0 <= x < XT_RICH:
                    sub = 'case2' + ('_uptake' if g[3] > 9.999999747378752e-06 else '_none')
                elif E[JNACL] + E[JCACL2] > 0 and g[1] > 0:
                    prod = max(g[3] * g[2] - k0cl, 0.0) + E[JNH4CL]
                    sub = 'case3a'
                    if prod > 0:
                        phi = (g[3] * g[2] - ain['keq_sg'][1]) / max(g[3] * g[2], ain['keq_sg'][1])
                        if abs(phi) < RTOL_EQB:
                            skipped['case3b_deadband'] += 1
                            continue
                        sub = 'case3a_3b' if fx['int_s'][3] else 'case3a_3b_absent'
                else:
                    p1 = max(g[3] * g[1] - k0no3, 0.0) + E[JNH4NO3]
                    p2 = max(g[3] * g[2] - k0cl, 0.0) + E[JNH4CL]
                    if p1 <= 0 and p2 <= 0:
                        sub = 'none'
                    else:
                        ph1 = (g[3] * g[1] - ain['keq_sg'][0]) / max(g[3] * g[1], ain['keq_sg'][0])
                        ph2 = (g[3] * g[2] - ain['keq_sg'][1]) / max(g[3] * g[2], ain['keq_sg'][1])
                        if abs(ph1) < RTOL_EQB or abs(ph2) < RTOL_EQB:
                            skipped['case4_deadband'] += 1
                            continue
                        sub = 'case4_' + ('a' if fx['int_s'][1] else '') + ('b' if fx['int_s'][2] else '')
                        if sub == 'case4_':
                            sub = 'case4_absent'
                groups[sub].append((c, sb, ibin, fx))
    return groups, skipped


def dry_test(sub, c, sb, ibin, fx):
    ain, b = c['ain'], c['bins'][ibin]
    g, kg, E, i_s = fx['gas'], b['kg'], fx['E_s'], fx['ions_s']
    over = {
        'gas_hno3': g[1], 'gas_hcl': g[2], 'gas_nh3': g[3],
        'kg_hno3': kg[1], 'kg_hcl': kg[2], 'kg_nh3': kg[3],
        'so4_s': i_s[ION['so4']], 'msa_s': i_s[ION['msa']], 'nh4_s': i_s[ION['nh4']],
        'na_s': i_s[ION['na']], 'ca_s': i_s[ION['ca']], 'cl_s': i_s[ION['cl']],
        'caco3_s': E[JCACO3], 'nacl_s': E[JNACL], 'cacl2_s': E[JCACL2],
        'nh4no3_s': E[JNH4NO3], 'nh4cl_s': E[JNH4CL], 'salt_s': sum(E[:15]),
        'keq_sg_nh4no3': ain['keq_sg'][0], 'keq_sg_nh4cl': ain['keq_sg'][1],
        'kp_nh4no3_0': ain['kp_nh4no3_0'], 'kp_nh4cl_0': ain['kp_nh4cl_0']}
    flag = 'case' + sub[4]
    big = max(abs(v) for v in fx['flux_s'][1:4]) or 1.0
    tol = {'rel': 1e-9, 'abs': 1e-12 * big}
    asr = [a(flag, 1.0, {'abs': 0.0, 'rel': 0.0}),
           a('cond_hno3', fx['flux_s'][1], tol), a('cond_hcl', fx['flux_s'][2], tol),
           a('cond_nh3', fx['flux_s'][3], tol)]
    if sub in ('case4_a', 'case4_b', 'case3a_3b'):
        asr.append(a({'case4_a': 's4a_nh3', 'case4_b': 's4b_nh3', 'case3a_3b': 's3_nh3'}[sub],
                     fx['sfc'][3], {'rel': 1e-9, 'abs': 1e-12 * max(abs(fx['sfc'][3]), 1e-30)}))
    w = where(c, ibin, sb)
    return {
        'id': f"astem64_dry_{sub}_{w[0]}",
        'description': (f"Branch {sub}. {w[1]} (route 1: all solid). Inputs are the state the bin's "
                        "ASTEM_flux_dry saw (FXO record), kg and the constants from the call's AIN/BIN "
                        "records; expected are the solid fluxes flux_s it formed (and the surface NH3 its "
                        "quadratic gave, where that survives case 4ab's reconciliation). rel 1e-9, absolute "
                        "floor 1e-12 of the largest of the three fluxes."),
        'parameter_overrides': over,
        'time_span': {'start': 0, 'end': 1},
        'assertions': asr}


def dry_tests(calls, max_per=12):
    groups, skipped = dry_rows(calls)
    print('dry census:', {k: len(v) for k, v in sorted(groups.items())}, 'deadband (not asserted):', dict(skipped))
    tests = []
    for sub, rows in sorted(groups.items()):
        if sub in ('none',):
            continue
        for r in spread(rows, max_per):
            tests.append(dry_test(sub, *r))
    return tests

def main():
    d = Path(sys.argv[1])
    calls = load(d)
    # constants: one test per dumped step at three levels (surface, mid, top of the
    # aerosol column), so 240-310 K is spanned
    pick = [c for c in calls if c['lev'] in (1, 20, 40)]
    # plus the driest and the most humid cell of the capture, for Kp_nh4no3_0/Kp_nh4cl_0
    by_rh = sorted(calls, key=lambda c: relhumclm(d, c['step'], c['lev']))
    pick += [c for c in (by_rh[0], by_rh[-1]) if c not in pick]
    # plus one hand case per relative humidity the hand cases use (0.2-0.85: Kp0 above and
    # below keq_sg)
    hand = load(d, hand=True)
    seen = set()
    for c in hand:
        rh = cell_relhum(d, c)
        if rh not in seen:
            seen.add(rh)
            pick.append(c)
    t = constants_tests(calls, pick, d)
    for x, c in zip(t, pick):   # also the solid-gas constants
        x['assertions'] += [a('keq_sg_nh4no3', c['ain']['keq_sg'][0]),
                            a('keq_sg_nh4cl', c['ain']['keq_sg'][1])]
    write_tests('mosaic_gas_equilibrium_constants.esm', 'MosaicGasEquilibriumConstants', t)
    # mass transfer: real64 kg, spread over steps, levels, bins and the five gases
    rows = [(c, ib, g) for c in calls for ib in sorted(c['bins'])
            for g in ['h2so4', 'hno3', 'hcl', 'nh3', 'msa']
            if c['bins'][ib]['kg'][0] > 0]
    rows = spread(rows, 40)
    write_tests('gas_particle_mass_transfer.esm', 'MosaicGasParticleMassTransfer',
                mtc_tests(calls, rows), keep=lambda t: not t['id'].startswith('astem64_'))
    wt, cen = wet_tests(calls + hand)
    write_tests('astem_liquid_surface_equilibrium.esm', 'MosaicAstemLiquidSurfaceEquilibrium', wt)
    write_tests('astem_solid_surface_equilibrium.esm', 'MosaicAstemSolidSurfaceEquilibrium',
                dry_tests(calls + hand))


if __name__ == '__main__':
    main()
