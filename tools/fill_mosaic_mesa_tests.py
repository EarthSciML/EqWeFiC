#!/usr/bin/env python3
"""Fill the inline tests of components/aerosol/mosaic/{ions_to_electrolytes,mesa_ptc_step}.esm
and couplings/mosaic_mesa_pass.esm from the MESA tape.

Source (numbers only; the components' equations are hand-written): the MESA tape that
WRF fork earthsciml-instrumented-chem170 module_esm_mesa_tape writes from inside
mesa_ptc -- at the top of every pass, inside mesa_flux_salt once ions_to_electrolytes /
compute_activities / the first estimate are done, after the convergence test, after the
step, and at the routine's own exit -- when kernels/mosaic_therm_driver is given a tape
path.  Two sets, both real(kind=8) and both run through WRF's own compiled objects:
  - chem170_scm/mesa64/tape/all_NNN.tape.gz: every mesa_ptc hand-off of the chem170
    aerosol_phase_state capture (dumps/chem170_mosaic_phase), replayed in the driver
    (mode 5) and verified bit for bit against the column (verify_replay.py);
  - chem170_scm/mesa64/hand.tape: hand-specified bins (make_mesa_hand.py) for the
    branches the column never reaches.
Samples are chosen to cover branches (the census keys below), first occurrence in a
fixed file order, so the selection is deterministic.

  fill_mosaic_mesa_tests.py <data/eqwefic/chem170_scm/mesa64>
"""
import gzip, json, os, sys, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = os.path.join(HERE, '..', 'components', 'aerosol', 'mosaic')
NS, NE, NI = 15, 22, 9
# tape ion order so4 no3 cl nh4 msa co3 na ca oin  ->  component `sp` order so4 no3 cl nh4 msa na ca
SP = [0, 1, 2, 3, 4, 6, 7]
MW_E = [132.0, 247.0, 115.0, 113.0, 80.0, 53.5, 142.0, 262.0, 120.0, 118.0, 85.0, 58.5, 164.0, 111.0, 230.0]
DISSOLVABLE = [1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0]
SALT = 'nh4so4 lvcite nh4hso4 nh4msa nh4no3 nh4cl na2so4 na3hso4 nahso4 namsa nano3 nacl cano3 cacl2 camsa2'.split()
LAYOUT = {
    'ENTRY': [('t_k', 1), ('rh_pc', 1), ('ah2o_a', 1), ('mass_dry_a', 1), ('vol_dry_a', 1), ('mass_dry_salt', 1),
              ('keq_sl', NS), ('keq_ll', 3), ('jsalt_present', NS), ('epercent_t', NE), ('elec_t', NE), ('ions_t', NI)],
    'AERT': [('aer_t', None)], 'GAMZ': [('log_gamz', None)],
    'TOP': [('elec_s', NE), ('elec_l', NE), ('ions_s', NI), ('ions_l', NI), ('phi_old', NS), ('alpha', NS),
            ('phi', NS), ('water', 1)],
    'THERM': [('ions_l', NI), ('ions_s', NI), ('elec_l', NE), ('elec_s', NE), ('water', 1), ('activity', NE),
              ('eleliquid', NE)],
    'FLUX': [('iconv_mass', 1), ('iconv_flux', 1), ('idissolved', 1), ('sat', NS), ('phi', NS), ('flux', NS),
             ('frac_s', NS), ('frac_l', NS), ('elec_s', NE), ('ions_l', NI), ('ions_s', NI), ('eleliquid', NE),
             ('water', 1), ('activity', NE)],
    'STEP': [('phi_bar', NS), ('alpha', NS), ('phi_old', NS), ('hsalt', NS), ('flux', NS), ('phi', NS),
             ('elec_s', NE), ('ions_s', NI), ('ions_l', NI)],
    'EXIT': [('outcome', 1), ('state', 1), ('jhyst', 1), ('growth', 1), ('mass_wet', 1), ('vol_wet', 1),
             ('water', 1), ('mass_dry_a', 1), ('elec_s', NE), ('elec_l', NE), ('ions_s', NI), ('ions_l', NI),
             ('epercent_s', NE)],
}


def parse(tag, vals):
    out, i = {}, 0
    for name, n in LAYOUT[tag]:
        if n is None:
            out[name] = np.array(vals[i:]); i = len(vals)
        elif n == 1:
            out[name] = vals[i]; i += 1
        else:
            out[name] = np.array(vals[i:i + n]); i += n
    assert i == len(vals), (tag, i, len(vals))
    return out


def read(path):
    calls = {}
    with (gzip.open if path.endswith('.gz') else open)(path, 'rt') as f:
        for ln in f:
            p = ln.split()
            tag, case, call, it = p[0], int(p[1]), int(p[2]), int(p[3])
            r = parse(tag, [float(x) for x in p[4:]])
            c = calls.setdefault((case, call), {'its': {}})
            if tag in ('ENTRY', 'AERT', 'GAMZ', 'EXIT'):
                c[tag.lower()] = r
            else:
                c['its'].setdefault(it, {})[tag.lower()] = r
    return calls


def xt(a):  # a in tape ion order
    d = a[0] + 0.5 * a[4]
    return (a[3] + a[6] + 2.0 * a[7]) / d if d > 0 else -1.0


def keys_step(c, it, r):
    """Branch keys of one pass, read off the tape (selection only)."""
    k = set()
    th, fl, tp = r.get('therm'), r['flux'], r['top']
    if th is None:
        return {'no_water'}
    k.add('outcome_%d' % (1 if fl['iconv_mass'] else 2 if fl['iconv_flux'] else 0))
    if fl['iconv_mass'] == 0 and fl['iconv_flux'] == 1:
        k.add('flux_conv_dissolved_%d' % fl['idissolved'])
    k.add('it1' if it == 1 else 'it2' if it == 2 else 'itlate' if it >= 40 else 'itmid')
    if c['entry']['mass_dry_salt'] < np.float32(1e-30):
        k.add('no_salt')
    if th['ions_l'][7] > 0:
        k.add('calcium')
        for s in range(NS):
            if s in (12, 13) or (fl['sat'][s] < 1 and fl['elec_s'][s] == 0):
                continue
            if fl['frac_l'][s] < np.float32(0.01):
                k.add('ca_blocks_small_liquid')
            elif c['entry']['jsalt_present'][s] == 0:
                k.add('ca_blocks_absent')
    et = c['aert']['aer_t']
    if (2 * et[0] + et[1] + et[2] + et[4] + 2 * et[106]) - (2 * et[108] + et[107] + et[3]) < -1e-6 * np.max(np.abs(et)):
        k.add('charge_h_negative')
    mds = c['entry']['mass_dry_salt']
    if mds >= np.float32(1e-30):
        fr = sum(fl['elec_s'][s] * MW_E[s] * float(np.float32(1e-15)) for s in range(NS)) / mds
        if 0.97 <= fr < np.float32(0.98):
            k.add('mass_frac_just_below')
        if np.float32(0.98) <= fr < 0.99:
            k.add('mass_frac_just_above')
    if fl['iconv_mass'] == 0:
        pm = float(np.max(np.abs(fl['phi'])))
        if np.float32(0.01) < pm <= 0.02:
            k.add('phi_just_above_rtol')
        if 0.005 <= pm <= np.float32(0.01):
            k.add('phi_just_within_rtol')
        if fl['iconv_flux'] == 1 and np.sum(np.abs(fl['flux'])) == 0:
            if fl['elec_s'][21] + fl['elec_s'][20] + et[109] > 0:
                k.add('crustal_solids_%s' % ('liquid' if fl['idissolved'] == 1 else 'mixed'))
    for s in range(NS):
        if DISSOLVABLE[s] and fl['sat'][s] < 1 and np.float32(0.01) <= th['elec_s'][s] / max(np.sum(th['elec_s'][:NS]), 1e-300) < 0.02:
            k.add('small_salt_band')
    x = xt(fl['ions_l'])
    if np.float32(1.9999) <= x < 2:
        k.add('est_split_band')
    if x >= 2 or x < 0:
        k.add('est_poor')
        a = fl['ions_l']
        if (2 * a[0] + a[1] + a[2] + a[4]) - (2 * a[7] + a[3] + a[6]) > 1e-6 * np.max(np.abs(a)):
            k.add('est_poor_acid')
    elif fl['ions_l'][3] == 0:
        k.add('est_rich_no_nh4_' + ('gt1' if x > 1 else 'le1'))
    else:
        so4, na, nh4 = fl['ions_l'][0], fl['ions_l'][6], fl['ions_l'][3]
        xna = 1 + 0.5 * na / so4
        dum = 2 * so4 - na
        xn = 2 * nh4 / dum if dum > 0 else 0.0
        k.add('est_rich_nh4_' + (('xna_' + ('53' if xn >= np.float32(5 / 3) else '15' if xn >= 1.5 else '1'))
                                 if x >= xna else 'gt1' if x > 1 else 'le1'))
    for s in range(NS):
        if th['elec_s'][s] != fl['elec_s'][s]:
            k.add('dissolve_' + SALT[s])
        if fl['flux'][s] > 0:
            k.add('flux_pos_' + SALT[s])
        if fl['flux'][s] < 0:
            k.add('flux_neg_' + SALT[s])
        if 'step' in r:
            sp = r['step']
            if fl['flux'][s] > 0 and sp['flux'][s] == 0:
                k.add('taup_zero_pos')
            if fl['flux'][s] < 0:
                tpp = -fl['eleliquid'][s] / fl['flux'][s]; td = -fl['elec_s'][s] / fl['flux'][s]
                k.add('neg_' + ('taup0' if tpp == 0 else 'taud' if td < tpp else 'taup'))
            if sp['phi_bar'][s] < 0:
                k.add('phibar_neg')
            if sp['phi_bar'][s] == -10:
                k.add('phibar_floor')
            if sp['phi_bar'][s] > 0:
                k.add('phibar_pos')
            if sp['alpha'][s] == 0.9:
                k.add('alpha_09')
            if sp['alpha'][s] == 0.5:
                k.add('alpha_05')
    return k


def keys_i2e(r):
    t, th = r['top'], r['therm']
    a = np.maximum(t['ions_l'], 0.0)
    k = set()
    if th['elec_s'][20] != t['elec_s'][20]:
        k.add('caso4')
    x = xt(th['ions_l'])
    if np.float32(1.9999) <= x < 2:
        k.add('split_band')
    if x >= np.float32(1.9999) or x < 0:
        k.add('poor' + ('_ca' if th['ions_l'][7] > 0 else '') + ('_so4free' if th['ions_l'][0] == 0 else ''))
    else:
        L = th['elec_l']
        k.add('rich' + ('_acid' if L[15] > 0 else '') + ('_bisulf' if L[2] + L[8] > 0 else '') +
              ('_lvc' if L[1] + L[7] > 0 else '') + ('_sulf' if L[0] + L[6] > 0 else '') +
              ('_freemsa' if L[16] > 0 else '') + ('_msa' if L[3] + L[9] > 0 else '') +
              ('_hno3hcl' if L[17] + L[18] > 0 else ''))
    return k


def tol(rel=1e-9, ab=0.0):
    return {'rel': rel, 'abs': ab}


def put(A, var, x, coords=None, rel=1e-9, ab=0.0):
    d = {'variable': var}
    if coords:
        d['coords'] = coords
    d.update({'time': 0.0, 'expected': float(x), 'tolerance': tol(rel, ab)})
    A.append(d)


def put_exact(A, var, x, coords=None):
    put(A, var, x, coords, 0.0, 0.0)


def sp(v):
    return [float(v[i]) for i in SP]


def step_test(tid, why, src, c, it, r):
    e, t, th, fl = c['entry'], r['top'], r['therm'], r['flux']
    scale = float(np.max(np.abs(e['ions_t'])))  # largest total species amount of the bin
    fa = 1e-12 * scale                            # round-off floor of a total-minus-solid difference
    po = {'it': float(it), 't_k': float(e['t_k']), 'aer_t': [float(x) for x in c['aert']['aer_t']],
          'elec_t': [float(x) for x in e['elec_t']], 'elec_s_in': [float(x) for x in t['elec_s']],
          'phi_old_in': [float(x) for x in t['phi_old']], 'alpha_in': [float(x) for x in t['alpha']],
          'phi_in': [float(x) for x in t['phi']], 'a_l': sp(th['ions_l']),
          'caso4': float(th['elec_s'][20] - t['elec_s'][20]), 'elec_l': [float(x) for x in th['elec_l']],
          'water': float(th['water']), 'activity': [float(x) for x in th['activity']]}
    A = []
    for s in range(NS):
        put(A, 'keq_sl', e['keq_sl'][s], {'salt': s + 1})
        put_exact(A, 'present', e['jsalt_present'][s], {'salt': s + 1})
    put(A, 'mass_dry', e['mass_dry_a']); put(A, 'mass_dry_salt', e['mass_dry_salt'])
    for s in range(NS):
        for v, x in (('sat', fl['sat']), ('phi', fl['phi']), ('flux', fl['flux']), ('frac_s', fl['frac_s']),
                     ('frac_l', fl['frac_l'])):
            put(A, v, x[s], {'salt': s + 1})
    for i in range(NE):
        put(A, 'elec_s1', fl['elec_s'][i], {'el': i + 1}, ab=fa)
        put(A, 'eleliquid', fl['eleliquid'][i], {'el': i + 1}, ab=fa)
    for k, x in enumerate(sp(fl['ions_l'])):
        put(A, 'a_l1', x, {'sp': k + 1}, ab=fa)
    put_exact(A, 'conv_mass', fl['iconv_mass'])
    if fl['iconv_mass'] == 0:
        put_exact(A, 'conv_flux', fl['iconv_flux'])
        put_exact(A, 'dissolved', fl['idissolved'])
    oc = 1 if fl['iconv_mass'] else 2 if fl['iconv_flux'] else 0
    put_exact(A, 'outcome', oc)
    if 'step' in r:
        stp = r['step']
        for s in range(NS):
            for v, x in (('phi_bar', stp['phi_bar']), ('alpha', stp['alpha']), ('hsalt', stp['hsalt']),
                         ('flux_out', stp['flux']), ('phi_out', stp['phi'])):
                put(A, v, x[s], {'salt': s + 1})
        for i in range(NE):
            put(A, 'elec_s_new', stp['elec_s'][i], {'el': i + 1}, ab=fa)
        for k, x in enumerate(sp(stp['ions_s'])):
            put(A, 'a_s_new', x, {'sp': k + 1}, ab=fa)
        for k, x in enumerate(sp(stp['ions_l'])):
            put(A, 'a_l_new', x, {'sp': k + 1}, ab=fa)
    elif oc == 2:
        for i in range(NE):
            put(A, 'elec_s_exit', c['exit']['elec_s'][i], {'el': i + 1}, ab=fa)
    desc = (f"{why} Source: {src}, pass {it}. Inputs: the pass's starting state (TOP record) and its "
            "ions_to_electrolytes / water / activity results (THERM record); expected: the FLUX record "
            "(after mesa_convergence_criterion), the STEP record (after the step, when the pass does not "
            "exit) or the EXIT record (solid electrolytes on a flux-converged exit), and keq_sl / "
            "jsalt_present / dry masses from the ENTRY record. Real64 MESA tape (WRF's own compiled "
            "mesa_ptc via kernels/mosaic_therm_driver, WRF fork earthsciml-instrumented-chem170); rel 1e-9, "
            "amounts with an absolute floor of 1e-12 of the bin's largest total species amount (the "
            "round-off of a total-minus-solid difference), flags exactly.")
    return {'id': tid, 'description': desc, 'parameter_overrides': po, 'time_span': {'start': 0, 'end': 1},
            'assertions': A}


def i2e_test(tid, why, src, it, r):
    t, th = r['top'], r['therm']
    scale = float(max(np.max(np.abs(t['ions_l'])), np.max(np.abs(th['ions_l'])), 1e-300))
    fa = 1e-12 * scale
    A = []
    for k, x in enumerate(sp(th['ions_l'])):
        put(A, 'a', x, {'sp': k + 1}, ab=fa)
    put(A, 'caso4', th['elec_s'][20] - t['elec_s'][20], ab=fa)
    for i in range(NE):
        put(A, 'elec', th['elec_l'][i], {'el': i + 1}, ab=fa)
    desc = (f"{why} Source: {src}, pass {it}: the liquid at the top of the pass (TOP record) in, the "
            "liquid and electrolytes mesa_flux_salt's ions_to_electrolytes left (THERM record) out. Real64 "
            "MESA tape (WRF's own compiled ions_to_electrolytes via kernels/mosaic_therm_driver); rel 1e-9 "
            "with an absolute floor of 1e-12 of the largest liquid amount.")
    return {'id': tid, 'description': desc, 'parameter_overrides': {'a_in': sp(t['ions_l'])},
            'time_span': {'start': 0, 'end': 1}, 'assertions': A}


def pass_test(tid, why, src, c, it, r):
    e, t, th, fl = c['entry'], r['top'], r['therm'], r['flux']
    scale = float(np.max(np.abs(e['ions_t'])))
    fa = 1e-12 * scale
    gz = c['gamz']['log_gamz'].reshape((NE, NE), order='F')
    po = {'MesaPass.it': float(it), 'MesaPass.t_k': float(e['t_k']), 'MesaPass.a_w': float(e['ah2o_a']),
          'MesaPass.aer_t': [float(x) for x in c['aert']['aer_t']], 'MesaPass.elec_t': [float(x) for x in e['elec_t']],
          'MesaPass.elec_s_in': [float(x) for x in t['elec_s']], 'MesaPass.a_l_in': sp(t['ions_l']),
          'MesaPass.phi_old_in': [float(x) for x in t['phi_old']], 'MesaPass.alpha_in': [float(x) for x in t['alpha']],
          'MesaPass.phi_in': [float(x) for x in t['phi']], 'MesaPass.log_gamz': gz.tolist(),
          'MesaPass.keq_hso4': float(e['keq_ll'][0]), 'MesaPass.k_w': float(e['keq_ll'][2])}
    A = []
    put(A, 'Zsr.water', th['water'])
    for s in range(NS):
        if s != 7:  # Na3H(SO4)2: the tape holds it after mesa_flux_salt zeroes it; Act returns the MTEM value
            put(A, 'Act.activity', th['activity'][s], {'el': s + 1})
        for v, x in (('Step.sat', fl['sat']), ('Step.phi', fl['phi']), ('Step.flux', fl['flux'])):
            put(A, v, x[s], {'salt': s + 1})
    for i in range(NE):
        put(A, 'Step.eleliquid', fl['eleliquid'][i], {'el': i + 1}, ab=fa)
        put(A, 'Step.elec_s1', fl['elec_s'][i], {'el': i + 1}, ab=fa)
    oc = 1 if fl['iconv_mass'] else 2 if fl['iconv_flux'] else 0
    put_exact(A, 'Step.outcome', oc)
    if 'step' in r:
        for i in range(NE):
            put(A, 'Step.elec_s_new', r['step']['elec_s'][i], {'el': i + 1}, ab=fa)
        for k, x in enumerate(sp(r['step']['ions_l'])):
            put(A, 'Step.a_l_new', x, {'sp': k + 1}, ab=fa)
    desc = (f"{why} Source: {src}, pass {it}. The whole pass from its starting state (TOP record), the bin's "
            "totals and water activity (ENTRY), and its MTEM inputs (GAMZ, keq_ll): ZSR water and activities as "
            "mesa_flux_salt saw them (THERM), the fluxes and decisions (FLUX) and the next state (STEP). Real64 "
            "MESA tape; rel 1e-9, amounts with an absolute floor of 1e-12 of the bin's largest total amount, "
            "flags exactly.")
    return {'id': tid, 'description': desc, 'parameter_overrides': po, 'time_span': {'start': 0, 'end': 1},
            'assertions': A}


def main(d):
    hand_ids = json.load(open(os.path.join(d, 'cases_hand.json')))
    idx = json.load(open(os.path.join(d, 'replay', 'all_index.json')))
    byfile = {}
    for r in idx:
        byfile[(r['file'], r['case'])] = r
    sources = [('hand', os.path.join(d, 'hand.tape'))] + \
              [('col%03d' % n, f) for n, f in enumerate(sorted(glob.glob(os.path.join(d, 'tape', 'all_*.tape.gz'))))]
    step_sel, i2e_sel, pass_sel, seen_s, seen_i = [], [], [], set(), set()
    late_nc = 0
    for tag, f in sources:
        for (case, call), c in sorted(read(f).items()):
            if tag == 'hand':
                hid = hand_ids[case - 1]
                src = f"hand-specified bin '{hid['id']}' ({hid['intent']}; chem170_scm/mesa64/make_mesa_hand.py), mesa_ptc call {call}"
                short = f"hand_{hid['id']}_c{call}"
            else:
                m = byfile[(int(tag[3:]), case)]
                src = (f"chem170 column hand-off, step {m['step']} level {m['k'] + 1} bin {m['b'] + 1} slot "
                       f"{m['slot']} (aerosol_phase_state path {m['path']}), mesa_ptc call {call}")
                short = f"col_s{m['step']}_lev{m['k'] + 1}_bin{m['b'] + 1}_c{m['slot']}_m{call}"
            exit_oc = int(c['exit']['outcome'])
            for it, r in sorted(c['its'].items()):
                if 'therm' not in r:
                    continue
                ks = keys_step(c, it, r)
                new = ks - seen_s
                # also the last pass of the first two unconverged calls (the 80-pass cut-off)
                force = exit_oc == 3 and it == 80 and late_nc < 2
                if new or force:
                    if force:
                        late_nc += 1
                    seen_s |= ks
                    why = ('Covers: ' + ', '.join(sorted(new)) + '.') if new else 'The 80th and last pass of a call MESA does not converge (it then calls the bin mixed).'
                    step_sel.append(step_test(f"{short}_it{it}", why, src, c, it, r))
                    pass_sel.append(pass_test(f"{short}_it{it}", why, src, c, it, r))
                ki = keys_i2e(r) - seen_i
                if ki:
                    seen_i |= ki
                    i2e_sel.append(i2e_test(f"{short}_it{it}", 'Covers: ' + ', '.join(sorted(ki)) + '.', src, it, r))
    for name, model, tests in ((os.path.join(COMP, 'mesa_ptc_step.esm'), 'MosaicMesaPtcStep', step_sel),
                               (os.path.join(COMP, 'ions_to_electrolytes.esm'), 'MosaicIonsToElectrolytes', i2e_sel),
                               (os.path.join(HERE, '..', 'couplings', 'mosaic_mesa_pass.esm'), 'MesaPass', pass_sel)):
        p = name
        doc = json.load(open(p))
        doc['models'][model]['tests'] = tests
        with open(p, 'w') as f:
            json.dump(doc, f, indent=2); f.write('\n')
        print(name, len(tests), 'tests', sum(len(t['assertions']) for t in tests), 'assertions')
    print('step keys covered:', sorted(seen_s))
    print('i2e keys covered:', sorted(seen_i))


if __name__ == '__main__':
    main(sys.argv[1])
