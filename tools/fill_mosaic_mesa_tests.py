#!/usr/bin/env python3
"""Fill the inline tests of components/aerosol/mosaic/{ions_to_electrolytes,mesa_salt_equilibrium}.esm
and couplings/mosaic_mesa_salt_rates.esm from the MESA tape.

The MESA tests evaluate the solid-liquid equilibrium's RATE and RESIDUAL instantaneously at states
mesa_ptc actually visited (the pass's starting solid); they do not assert mesa_ptc's step, its
step-size heuristics or its stopping point, which are solver, not physics (CLAUDE.md).

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

Plus, when <mesa64>/tight/tight.tape.gz exists, tests at CONVERGED states (mesa_ptc run to rtol 1e-10 by the
tolerance-study build, mesa64/make_tight_cases.py): the rate zero and g in complementarity with E_s.
"""
import gzip, json, os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mosaic_cell_axes import cell_axes  # noqa: E402

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


def snapped(r):
    """Salts mesa_dissolve_small_salt returned to the liquid on this pass (a solver snap)."""
    th, fl = r['therm'], r['flux']
    return [s for s in range(NS) if th['elec_s'][s] != fl['elec_s'][s]]


def keys_step(c, it, r):
    """Physics-branch keys of the state a pass starts from, read off the tape (selection only)."""
    k = set()
    th, fl = r.get('therm'), r['flux']
    if th is None:
        return {'no_water'}
    k.add('all_solid_%d' % fl['iconv_mass'])
    if fl['iconv_mass'] == 0 and fl['iconv_flux'] == 1:
        k.add('all_liquid_%d' % fl['idissolved'])
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
    if fl['iconv_mass'] == 0 and fl['iconv_flux'] == 1 and np.sum(np.abs(fl['flux'])) == 0:
        if fl['elec_s'][21] + fl['elec_s'][20] + et[109] > 0:
            k.add('crustal_solids_%s' % ('liquid' if fl['idissolved'] == 1 else 'mixed'))
    if th['elec_s'][20] != r['top']['elec_s'][20]:
        k.add('caso4_to_solid')
    for s in snapped(r):
        k.add('snap_' + SALT[s])
    for s in range(NS):
        if fl['flux'][s] > 0:
            k.add('flux_pos_' + SALT[s])
        if fl['flux'][s] < 0:
            k.add('flux_neg_' + SALT[s])
        if fl['sat'][s] < 1 and fl['elec_s'][s] == 0 and s not in (12, 13):
            k.add('undersat_no_solid')
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


def rate_expected(c, r):
    """sat, phi, flux the equilibrium's rate gives at the pass's starting state.  Equal to the FLUX record
    except for a salt mesa_dissolve_small_salt snapped on that pass, whose taped flux_sl/phi_salt are 0 by
    the snap: there the rate is formed from the taped S under the same switch-offs (solver snap removed)."""
    fl, th = r['flux'], r['therm']
    sat, phi, flux = fl['sat'].copy(), fl['phi'].copy(), fl['flux'].copy()
    ca = th['ions_l'][7]
    for s in snapped(r):
        S = fl['sat'][s]
        off = ca > 0 and (fl['frac_l'][s] < np.float32(0.01) or c['entry']['jsalt_present'][s] == 0)
        flux[s] = 0.0 if off else S - 1.0
        phi[s] = 0.0 if off else (S - 1.0) / max(S, 1.0)
    return sat, phi, flux


def g_expected(c, r):
    """g, the equilibrium residual paired with E_s: S - 1 where MOSAIC lets the salt exchange (and the liquid
    holds water), 0 where its own switch-offs hold -- the taped flux without the complementarity switch-off."""
    fl, th = r['flux'], r['therm']
    ca = th['ions_l'][7]
    g = np.zeros(NS)
    if th['water'] > 0:
        for s in range(NS):
            off = s in (12, 13) or (ca > 0 and (fl['frac_l'][s] < np.float32(0.01) or c['entry']['jsalt_present'][s] == 0))
            g[s] = 0.0 if off else fl['sat'][s] - 1.0
    return g


def step_test(tid, why, src, c, it, r):
    e, t, th, fl = c['entry'], r['top'], r['therm'], r['flux']
    scale = float(np.max(np.abs(e['ions_t'])))  # largest total species amount of the bin
    fa = 1e-12 * scale                            # round-off floor of a total-minus-solid difference
    po = {'t_k': float(e['t_k']), 'aer_t': [float(x) for x in c['aert']['aer_t']],
          'elec_t': [float(x) for x in e['elec_t']], 'elec_s': [float(x) for x in t['elec_s']],
          'a_l': sp(th['ions_l']), 'caso4': float(th['elec_s'][20] - t['elec_s'][20]),
          'elec_l': [float(x) for x in th['elec_l']], 'water': float(th['water']),
          'activity': [float(x) for x in th['activity']]}
    A = []
    for s in range(NS):
        put(A, 'keq_sl', e['keq_sl'][s], {'salt': s + 1})
        put_exact(A, 'present', e['jsalt_present'][s], {'salt': s + 1})
    put(A, 'mass_dry', e['mass_dry_a']); put(A, 'mass_dry_salt', e['mass_dry_salt'])
    for k, x in enumerate(sp(t['ions_s'])):
        put(A, 'a_s', x, {'sp': k + 1}, ab=fa)
    for k, x in enumerate(sp(t['ions_l'])):
        put(A, 'a_l_state', x, {'sp': k + 1}, ab=fa)
    sat, phi, flux = rate_expected(c, r)
    g = g_expected(c, r)
    for s in range(NS):
        for v, x in (('sat', sat), ('phi', phi), ('flux', flux), ('g', g), ('frac_l', fl['frac_l'])):
            put(A, v, x[s], {'salt': s + 1})
    snap = snapped(r)
    if not snap:
        put_exact(A, 'all_solid', fl['iconv_mass'])
        if fl['iconv_mass'] == 0:
            put_exact(A, 'all_liquid', fl['idissolved'])
    note = (" mesa_dissolve_small_salt snapped " + ', '.join(SALT[s] for s in snap) + " on this pass: its rate is "
            "asserted as the taped S - 1 (the snap is solver), and the solid-mass flags are not asserted (the "
            "Fortran forms them after the snap).") if snap else ""
    desc = (f"{why} Source: {src}, pass {it}: the state mesa_ptc started this pass from.{note} Inputs: the pass's "
            "starting solid (TOP record) and the ions_to_electrolytes / water / activity results at it (THERM "
            "record); expected: S, phi, the rate and the liquid fractions (FLUX record), the mass balance against "
            "the TOP record's own solid and liquid ions, keq_sl / jsalt_present / dry masses (ENTRY record) and "
            "the phase flags (mesa_convergence_criterion). Real64 MESA tape (WRF's own compiled mesa_ptc via "
            "kernels/mosaic_therm_driver, WRF fork earthsciml-instrumented-chem170); rel 1e-9, amounts with an "
            "absolute floor of 1e-12 of the bin's largest total species amount (the round-off of a "
            "total-minus-solid difference), flags exactly.")
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
    po = {'MesaRates.t_k': float(e['t_k']), 'MesaRates.a_w': float(e['ah2o_a']),
          'MesaRates.aer_t': [float(x) for x in c['aert']['aer_t']], 'MesaRates.elec_t': [float(x) for x in e['elec_t']],
          'MesaRates.elec_s': [float(x) for x in t['elec_s']], 'MesaRates.log_gamz': gz.tolist(),
          'MesaRates.keq_hso4': float(e['keq_ll'][0]), 'MesaRates.k_w': float(e['keq_ll'][2])}
    A = []
    for k, x in enumerate(sp(t['ions_l'])):
        put(A, 'Eq.a_l_state', x, {'sp': k + 1}, ab=fa)
    put(A, 'Zsr.water', th['water'])
    sat, phi, flux = rate_expected(c, r)
    for s in range(NS):
        if s != 7:  # Na3H(SO4)2: the tape holds it after mesa_flux_salt zeroes it; Act returns the MTEM value
            put(A, 'Act.activity', th['activity'][s], {'el': s + 1})
        for v, x in (('Eq.sat', sat), ('Eq.phi', phi), ('Eq.flux', flux), ('Eq.g', g_expected(c, r))):
            put(A, v, x[s], {'salt': s + 1})
    if not snapped(r):
        put_exact(A, 'Eq.all_solid', fl['iconv_mass'])
    desc = (f"{why} Source: {src}, pass {it}. The rate and residual of the solid-liquid equilibrium at the state "
            "mesa_ptc started this pass from, assembled from that state alone: its solid (TOP record), the bin's "
            "totals and water activity (ENTRY), and its MTEM inputs (GAMZ, keq_ll); expected: the liquid by mass "
            "balance (TOP), the ZSR water and activities mesa_flux_salt saw (THERM), and S, phi and the rate "
            "(FLUX; a salt mesa_dissolve_small_salt snapped is asserted at the taped S - 1). Real64 MESA tape; "
            "rel 1e-9, amounts with an absolute floor of 1e-12 of the bin's largest total amount, flags exactly.")
    return {'id': tid, 'description': desc, 'parameter_overrides': po, 'time_span': {'start': 0, 'end': 1},
            'assertions': A}


# Converged states (2026-11-12): mesa_ptc run to rtol 1e-10 (ESM_RTOL_MESA, tolerance-study build
# chem170_scm/solver_tol) on hand-offs that converge there, taped (mesa64/tight, make_tight_cases.py).
# Its stopping test is |phi_s| <= 1e-10 on the last pass, so |S_s - 1| <= 1e-10 max(S_s, 1) for every
# salt that may exchange; G_TOL is twice that, the margin being the component's round-off against the
# Fortran's (the S values agree to ~1e-15).
G_TOL = 2e-10
TIGHT_WHY = {2: "a hand-off mesa_ptc converges at its own 1% test", 3: "a hand-off mesa_ptc does NOT converge in its 80 passes"}


def converged_asserts(A, c, r, pre):
    """At a converged state: the rate is zero for every salt; g_s is zero where solid is present and is the
    (non-positive) S_s - 1 where none is -- the complementarity E_s >= 0, g_s <= 0, E_s g_s = 0."""
    fl = r['flux']
    g = g_expected(c, r)
    zero = {'rel': 0.0, 'abs': G_TOL}
    for s in range(NS):
        A.append({'variable': pre + 'flux', 'coords': {'salt': s + 1}, 'time': 0.0, 'expected': 0.0, 'tolerance': zero})
        if fl['elec_s'][s] > 0:
            A.append({'variable': pre + 'g', 'coords': {'salt': s + 1}, 'time': 0.0, 'expected': 0.0, 'tolerance': zero})
        else:
            assert g[s] <= G_TOL, (s, g[s])
            put(A, pre + 'g', g[s], {'salt': s + 1}, ab=G_TOL)
    put_exact(A, pre + 'all_solid', fl['iconv_mass'])


def converged_desc(src, n_it):
    return (f"THE EQUILIBRIUM, not a visited pass: {src}, run with mesa_ptc's tolerance tightened from 1e-2 to "
            f"1e-10 (rtol_mesa; WRF's own compiled code, chem170_scm/solver_tol build), which converged in {n_it} "
            "passes; inputs are its final state (TOP record of the last pass) and the bin's ENTRY / GAMZ records. "
            "Asserted: the rate is zero for every salt, and g (the residual paired with E_s) is zero where solid is "
            "present and equal to the non-positive S - 1 where none is -- the complementarity E_s >= 0, g_s <= 0, "
            f"E_s g_s = 0 the bounded EarthSciIntegrator solve (D95) finds. Tolerance: absolute {G_TOL:g}, twice the "
            "tightened stopping test (|S - 1| <= 1e-10 max(S, 1)); S - 1 itself at rel 1e-9.")


def converged_step_test(tid, src, c, it, r):
    e, t, th = c['entry'], r['top'], r['therm']
    po = {'t_k': float(e['t_k']), 'aer_t': [float(x) for x in c['aert']['aer_t']],
          'elec_t': [float(x) for x in e['elec_t']], 'elec_s': [float(x) for x in t['elec_s']],
          'a_l': sp(th['ions_l']), 'caso4': float(th['elec_s'][20] - t['elec_s'][20]),
          'elec_l': [float(x) for x in th['elec_l']], 'water': float(th['water']),
          'activity': [float(x) for x in th['activity']]}
    A = []
    converged_asserts(A, c, r, '')
    return {'id': tid, 'description': converged_desc(src, it), 'parameter_overrides': po,
            'time_span': {'start': 0, 'end': 1}, 'assertions': A}


def converged_pass_test(tid, src, c, it, r):
    e, t, th = c['entry'], r['top'], r['therm']
    gz = c['gamz']['log_gamz'].reshape((NE, NE), order='F')
    po = {'MesaRates.t_k': float(e['t_k']), 'MesaRates.a_w': float(e['ah2o_a']),
          'MesaRates.aer_t': [float(x) for x in c['aert']['aer_t']], 'MesaRates.elec_t': [float(x) for x in e['elec_t']],
          'MesaRates.elec_s': [float(x) for x in t['elec_s']], 'MesaRates.log_gamz': gz.tolist(),
          'MesaRates.keq_hso4': float(e['keq_ll'][0]), 'MesaRates.k_w': float(e['keq_ll'][2])}
    A = []
    put(A, 'Zsr.water', th['water'])
    converged_asserts(A, c, r, 'Eq.')
    return {'id': tid, 'description': converged_desc(src, it) + " Assembled from the state alone (ions_to_electrolytes, "
            "ZSR water and MTEM activities at the converged solid).", 'parameter_overrides': po,
            'time_span': {'start': 0, 'end': 1}, 'assertions': A}


def tight_tests(d):
    tdir = os.path.join(d, 'tight')
    if not os.path.exists(os.path.join(tdir, 'tight.tape.gz')):
        return [], []
    meta = json.load(open(os.path.join(tdir, 'cases.json')))
    sidx = {(r['file'], r['case']): r for r in json.load(open(os.path.join(d, '..', 'solver_tol', 'cases', 'index.json')))}
    st, pt = [], []
    for (case, call), c in sorted(read(os.path.join(tdir, 'tight.tape.gz')).items()):
        m = meta[case - 1]
        assert int(c['exit']['outcome']) == 2, (case, call)
        it = max(c['its']); r = c['its'][it]
        assert not snapped(r) and np.max(np.abs(r['flux']['phi'])) <= 1e-10
        x = sidx[(m['file'], m['case'])]
        src = (f"{TIGHT_WHY[m['outcome_default']]} ({m['passes_default']} passes at the model's settings), chem170 "
               f"column step {x['step']} level {x['k'] + 1} bin {x['b'] + 1} slot {x['slot']}")
        tid = f"converged_s{x['step']}_lev{x['k'] + 1}_bin{x['b'] + 1}_c{x['slot']}"
        st.append(converged_step_test(tid, src, c, it, r))
        pt.append(converged_pass_test(tid, src, c, it, r))
    return st, pt


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
                    why = ('Covers: ' + ', '.join(sorted(new)) + '.') if new else 'The state of the 80th and last pass of a call mesa_ptc does not converge (it then calls the bin mixed): the residual there is what its 80-pass cap leaves.'
                    step_sel.append(step_test(f"{short}_it{it}", why, src, c, it, r))
                    pass_sel.append(pass_test(f"{short}_it{it}", why, src, c, it, r))
                ki = keys_i2e(r) - seen_i
                if ki:
                    seen_i |= ki
                    i2e_sel.append(i2e_test(f"{short}_it{it}", 'Covers: ' + ', '.join(sorted(ki)) + '.', src, it, r))
    ts, tp = tight_tests(d)
    step_sel += ts
    pass_sel += tp
    for name, model, tests in ((os.path.join(COMP, 'mesa_salt_equilibrium.esm'), 'MosaicMesaSaltEquilibrium', step_sel),
                               (os.path.join(COMP, 'ions_to_electrolytes.esm'), 'MosaicIonsToElectrolytes', i2e_sel),
                               (os.path.join(HERE, '..', 'couplings', 'mosaic_mesa_salt_rates.esm'), 'MesaRates', pass_sel)):
        p = name
        doc = json.load(open(p))
        doc['models'][model]['tests'] = cell_axes(doc, model, tests, p)
        with open(p, 'w') as f:
            json.dump(doc, f, indent=2); f.write('\n')
        print(name, len(tests), 'tests', sum(len(t['assertions']) for t in tests), 'assertions')
    print('step keys covered:', sorted(seen_s))
    print('i2e keys covered:', sorted(seen_i))


if __name__ == '__main__':
    main(sys.argv[1])
