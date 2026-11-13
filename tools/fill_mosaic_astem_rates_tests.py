#!/usr/bin/env python3
"""Fill the inline tests of couplings/mosaic_astem_rates.esm (MOSAIC's gas-particle
condensation rates over the whole size distribution) from the ASTEM tape.

    tools/fill_mosaic_astem_rates_tests.py <astem_dir>

<astem_dir> is data/eqwefic/chem170_scm/mosaic_astem (astem.tape.gz, hand.tape,
cases_hand_names.txt, mix.tape, cases_mix_names.txt, astem_tape.py; see
data/eqwefic/notes/chem170_mosaic_astem_capture.md).
Only numbers are written; the document's equations and structure are hand-authored.

A test is one state ASTEM formed its fluxes at, all bins at once:
  * column: the FIRST sub-step of an ASTEM call in which every bin saw the same gas and no
    bin's state was edited before its flux was formed (FXI gas = FXO gas = SVT gas, FXO
    ions/electrolytes = SVB); in the other sub-steps a solid bin's flux routine edits HCl
    and NH3 before the next bin forms its flux, so no single state reproduces them;
  * hand: each hand-specified bin of kernels/mosaic_therm_driver mode 6 in bin 1, bins 2..
    NBIN empty (state 0, number 0): hand.tape, and mix.tape (make_astem_mix_hand.py: mixed
    bins that reach ASTEM_flux_mix, and all-solid bins whose Kp0 is below Keq_sg).
Expected (all real64 values ASTEM held):
  * kg of HNO3, HCl, NH3 (BIN record);
  * each bin's semi-volatile rate: flux_s for a solid bin; flux_s + kg df_gas_l where ASTEM
    integrates the liquid (its liquid update is kg (gas - Heff aer_l), Heff aer_l = surface);
  * gas kg for H2SO4 and MSA in every bin holding aerosol;
  * the gas tendencies, minus the sum over bins of those rates; the balances, zero.
Never in a test (ASTEM's solver decisions, CLAUDE.md): a solid bin inside ASTEM_flux_dry's
or ASTEM_flux_mix's 1% deadband, a liquid (or the liquid half of a mixed bin) ASTEM skipped
although a driving force is nonzero (the signed equilibrium test, FORTRAN_BUGS.md B31), and
equilibrate_acids.
"""
import collections
import importlib.util
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / 'couplings' / 'mosaic_astem_rates.esm'
MODEL = 'AstemRates'
NBIN = 8
NEL = 22
REL = 1e-9
EMPTY_NUM = 1000.0   # 1/cm3 in the empty bins of the hand tests

_spec = importlib.util.spec_from_file_location('fat', ROOT / 'tools' / 'fill_astem_tests.py')
fat = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fat)
_spec = importlib.util.spec_from_file_location('fmx', ROOT / 'tools' / 'fill_astem_mix_tests.py')
fmx = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fmx)

GASES = ['h2so4', 'hno3', 'hcl', 'nh3', 'msa']
SEMI = [('hno3', 1, 'no3'), ('hcl', 2, 'cl'), ('nh3', 3, 'nh4')]   # gas, 0-based gas index, ion
JNH4NO3, JNH4CL = 4, 5                                               # 0-based electrolytes


def load(d):
    sys.path.insert(0, str(d))
    from astem_tape import read_tape

    def names_of(fname):
        out = {}
        for line in open(Path(d) / fname):
            i, name, why = line.rstrip('\n').split('\t')
            out[int(i)] = (name, why)
        return out
    hand = [(read_tape(Path(d) / 'hand.tape'), names_of('cases_hand_names.txt'), 'make_astem_hand.py'),
            (read_tape(Path(d) / 'mix.tape'), names_of('cases_mix_names.txt'), 'make_astem_mix_hand.py')]
    return read_tape(Path(d) / 'astem.tape.gz'), hand


def clean(sb):
    """Every bin formed its flux at the sub-step's own state (no edit between bins)."""
    g = sb['gas']
    for b, fx in sb['fxo'].items():
        pre = sb['bins'][b]
        if sb['fxi'][b]['gas'] != g or fx['gas'] != g:
            return False
        if fx['ions_s'] != pre['ions_s'] or fx['ions_l'] != pre['ions_l'] or fx['E_s'] != pre['E_s']:
            return False
    return True


def route_of(state, E_s):
    """ASTEM_semi_volatiles' dispatch, as the coupling writes it."""
    if state == 1:
        return 1
    if state == 2:
        return 2
    if state == 3:
        return 3 if (E_s[JNH4NO3] > 0 or E_s[JNH4CL] > 0) else 4
    return 0


def bin_branch(c, sb, b):
    """(branch label, excluded reason or None) for one bin's flux call."""
    fx = sb['fxo'][b]
    route = int(fx['route'])
    if route != route_of(sb['bins'][b]['jaerosolstate'], fx['E_s']):
        return None, 'route mismatch'
    if route == 3:
        case = fmx.mix_case(fx)
        if fmx.deadband(case, fx, c['ain']):
            return None, 'mixed solid 1% deadband'
        for i in (1, 2, 3):
            if fx['int_l'][i] == 0 and fx['df_l'][i] != 0:
                return None, 'liquid half skipped by the signed equilibrium test (B31)'
        return 'mix_' + case, None
    if route == 1:
        groups, skipped = fat.dry_rows([{'ain': c['ain'], 'subs': [{'fxo': {b: fx}}]}])
        if skipped:
            return None, 'solid 1% deadband'
        return 'dry_' + next(iter(groups)), None
    for i in (1, 2, 3):
        if fx['int_l'][i] == 0 and fx['df_l'][i] != 0:
            return None, 'liquid skipped by the signed equilibrium test (B31)'
    return 'wet_' + fat.wet_case(fx, fx['gas']), None


def semi_expected(c, b, fx):
    """Per-gas semi-volatile rate ASTEM applies to bin b, and its absolute floor."""
    kg = c['bins'][b]['kg']
    out = {}
    for gname, i, _ in SEMI:
        if fx['route'] == 1:
            v = fx['flux_s'][i]
            floor = 1e-12 * max(abs(x) for x in fx['flux_s'][1:4])
        else:
            v = fx['flux_s'][i] + (kg[i] * fx['df_l'][i] if fx['int_l'][i] else 0.0)
            floor = 1e-12 * kg[i] * max(abs(fx['gas'][i]), abs(fx['sfc'][i]))
        out[gname] = (v, floor)
    return out


GAX = {'h2so4': 1, 'hno3': 2, 'hcl': 3, 'nh3': 4, 'msa': 5, 'so4': 1, 'no3': 2, 'cl': 3, 'nh4': 4}


def a(var, val, tol, b=None):
    """One assertion on the (single-cell, NLEV = 1) assembly; the per-ion rates and the gas
    tendencies live on the gas axis (daer[lev, bin, gas], dgas/balance[lev, gas])."""
    coords = {'lev': 1}
    if b is not None:
        coords['bin'] = b
    for p in ('daer_', 'dgas_', 'balance_'):
        if var.startswith(p):
            coords['gas'] = GAX[var[len(p):]]
            var = p[:-1]
    d = {'variable': f'{MODEL}.{var}', 'time': 0.0, 'expected': float(val), 'tolerance': tol}
    d['coords'] = coords
    return d


EXACT = {'abs': 0.0, 'rel': 0.0}


def build_test(tid, desc, c, sb, bins):
    """bins: {b: fx} for the bins that hold aerosol; the others are empty."""
    ain = c['ain']
    g = next(iter(bins.values()))['gas']   # the gas the fluxes were formed at (= SVT gas in a clean column sub-step)
    assert all(fx['gas'] == g for fx in bins.values())
    over = {'t_k': ain['t_k'], 'p_atm': ain['p_atm'], 'relhum': fat.cell_relhum(ASTEM_DIR, c)}
    for k, name in enumerate(GASES):
        over[f'gas_{name}'] = g[k]
    arr = collections.defaultdict(list)
    for b in range(1, NBIN + 1):
        if b in bins:
            fx, B = bins[b], c['bins'][b]
            arr['state'].append(sb['bins'][b]['jaerosolstate'])
            arr['num'].append(B['num_a'])
            arr['dp_wet'].append(B['dp_wet'])
            arr['kel'].append(list(fx['kel']))
            arr['aer_s'].append(list(fx['ions_s']))
            arr['aer_l'].append(list(fx['ions_l']))
            arr['elec_s'].append(list(fx['E_s']))
            arr['mc_h'].append(fx['mc'][fat.JC_H])
            arr['mc_nh4'].append(fx['mc'][fat.JC_NH4])
            arr['ma_no3'].append(fx['ma'][fat.JA_NO3])
            arr['ma_cl'].append(fx['ma'][fat.JA_CL])
            arr['gam'].append(list(fx['gam']))
            arr['gam_ratio'].append(fx['gam_ratio'])
            arr['activity'].append(list(fx['activity']))
        else:
            # an empty bin (state 0) still given a number, so its kg is nonzero and the
            # has_aer / route gates are what keep its rates at zero (aerosolmtc skips such a
            # bin, so ASTEM's own kg there is the 0 initialize_mosaic_variables left)
            arr['state'].append(0.0)
            arr['num'].append(EMPTY_NUM)
            arr['dp_wet'].append(1e-5)
            arr['kel'].append([1.0] * 5)
            arr['aer_s'].append([0.0] * 9)
            arr['aer_l'].append([0.0] * 9)
            arr['elec_s'].append([0.0] * NEL)
            for k in ('mc_h', 'mc_nh4', 'ma_no3', 'ma_cl'):
                arr[k].append(0.0)
            arr['gam'].append([1.0] * NEL)
            arr['gam_ratio'].append(1.0)
            arr['activity'].append([0.0] * NEL)
    over.update({k: [v] for k, v in arr.items()})   # the single cell, lev = 1
    over['gas'] = [[over.pop(f'gas_{name}') for name in GASES]]

    asr, sums = [], {n: [0.0, 0.0] for n in GASES}
    for b in range(1, NBIN + 1):
        if b not in bins:
            for v in ('has_aer', 'route_solid', 'route_liquid', 'route_mix',
                      'daer_so4', 'daer_no3', 'daer_cl', 'daer_nh4', 'daer_msa'):
                asr.append(a(v, 0.0, EXACT, b))
            continue
        fx, kg = bins[b], c['bins'][b]['kg']
        route = int(fx['route'])
        asr.append(a('has_aer', 1.0, EXACT, b))
        asr.append(a('route_solid', 1.0 if route == 1 else 0.0, EXACT, b))
        asr.append(a('route_liquid', 1.0 if route in (2, 4) else 0.0, EXACT, b))
        asr.append(a('route_mix', 1.0 if route == 3 else 0.0, EXACT, b))
        for gname, i, _ in SEMI:
            asr.append(a(f'kg_{gname}', kg[i], {'rel': REL, 'abs': 0.0}, b))
        for gname, i, ion in (('h2so4', 0, 'so4'), ('msa', 4, 'msa')):
            v = g[i] * kg[i]
            fl = 1e-12 * abs(v)
            asr.append(a(f'daer_{ion}', v, {'rel': REL, 'abs': fl}, b))
            sums[gname][0] += v
            sums[gname][1] += fl
        for gname, (v, fl) in semi_expected(c, b, fx).items():
            ion = dict((s[0], s[2]) for s in SEMI)[gname]
            asr.append(a(f'daer_{ion}', v, {'rel': REL, 'abs': fl}, b))
            sums[gname][0] += v
            sums[gname][1] += fl
    for gname in GASES:
        s, fl = sums[gname]
        asr.append(a(f'dgas_{gname}', -s, {'rel': REL, 'abs': fl}))
        asr.append(a(f'balance_{gname}', 0.0, {'rel': 0.0, 'abs': max(fl, 1e-300) * 1e-3}))
    over = {f'{MODEL}.{k}': v for k, v in over.items()}
    return {'id': tid, 'description': desc, 'parameter_overrides': over,
            'time_span': {'start': 0, 'end': 0}, 'assertions': asr}


ASTEM_DIR = None
KP0 = '--kp0' in sys.argv


def main():
    global ASTEM_DIR
    d = Path([x for x in sys.argv[1:] if not x.startswith('--')][0])
    ASTEM_DIR = d
    calls, hands = load(d)

    # ---- column: census of usable first sub-steps
    usable, reasons = [], collections.Counter()
    for c in calls:
        sb = c['subs'][0]
        if not clean(sb):
            reasons['state edited between bins'] += 1
            continue
        br, bad = {}, None
        for b in sorted(sb['fxo']):
            lab, why = bin_branch(c, sb, b)
            if why:
                bad = why
                break
            br[b] = lab
        if bad:
            reasons[bad] += 1
            continue
        usable.append((c, sb, br))
    print(f'column first sub-steps: {len(calls)}; usable {len(usable)}; excluded {dict(reasons)}')
    seen = collections.Counter(l for _, _, br in usable for l in br.values())
    print('branches over usable bins:', dict(sorted(seen.items())))

    # greedy cover of (route pattern, branch) then an even spread over the run
    picked, covered = [], set()
    def sig(u):
        c, sb, br = u
        pat = tuple(int(sb['fxo'][b]['route']) for b in sorted(sb['fxo']))
        return {('pat', pat)} | {('br', l) for l in br.values()}
    rest = list(usable)
    while True:
        best = max(rest, key=lambda u: len(sig(u) - covered), default=None)
        if best is None or not (sig(best) - covered):
            break
        picked.append(best)
        covered |= sig(best)
        rest.remove(best)
    n_extra = 4
    for k in range(n_extra):
        u = rest[round(k * (len(rest) - 1) / max(n_extra - 1, 1))]
        if u not in picked:
            picked.append(u)
    picked.sort(key=lambda u: (u[0]['step'], u[0]['lev']))

    tests = []
    for c, sb, br in picked:
        pat = ''.join(str(int(sb['fxo'][b]['route'])) for b in sorted(sb['fxo']))
        desc = (f"ASTEM tape (real64), chem_opt = 170 column, model step {c['step'] + 1}, level {c['lev']}, "
                f"first ASTEM sub-step, all {NBIN} bins at one state (ASTEM routes by bin: {pat}; 1 solid, "
                f"2 liquid, 3 mixed with solid NH4NO3/NH4Cl, 4 mixed handled as liquid). Bin branches: "
                + ', '.join(f'{b}:{br[b]}' for b in sorted(br)) + '. '
                "Inputs: the FXO records (each bin's state where its flux was formed), kg's inputs from AIN/BIN; "
                "expected: kg (BIN), each bin's semi-volatile rate (flux_s solid; flux_s + kg df_gas_l liquid), "
                "gas kg for H2SO4 and MSA, the gas tendencies as minus the sum over bins of those taped rates, "
                "the balances zero. rel 1e-9; absolute floor 1e-12 of the scale of each rate (kg times the larger "
                "of gas and surface value for a liquid bin, the largest of the three solid fluxes for a solid one); "
                "route flags exactly.")
        tests.append(build_test(f"astem_s{c['step']}_k{c['lev']}", desc, c, sb, sb['fxo']))

    # ---- hand bins
    hreasons = collections.Counter()
    for hand, names, maker in hands:
        for c in hand:
            sb = c['subs'][0]
            name, why = names[c['lev']]
            if name == 'wet_poor_noacid':
                hreasons['equilibrate_acids'] += 1
                continue
            if name.startswith('kp0_') and not KP0:
                # the be1fafe6d binary evaluates Const's Kp_nh4no3_0 / Kp_nh4cl_0 to NaN once Const
                # is lifted over lev (an index into a subsystem array inside a lifted system;
                # data/eqwefic/phase3_probes/const_lift/), and Sol's min(NaN, Keq_sg) then silently
                # picks Keq_sg: these bins, whose rate is set by Kp0, cannot pass on that binary.
                # EarthSciAST main computes Kp0 correctly. Pass --kp0 to include them.
                hreasons['Kp0 bin (lifted-Const NaN on be1fafe6d)'] += 1
                continue
            lab, bad = bin_branch(c, sb, 1)
            if bad:
                hreasons[bad] += 1
                print(f'  hand {name}: excluded ({bad})')
                continue
            desc = (f"Hand-specified bin '{name}' ({why}; chem170_scm/mosaic_astem/{maker}), run through "
                    f"WRF's own compiled ASTEM by kernels/mosaic_therm_driver mode 6 at T = 290 K, p = 0.9 atm "
                    f"(ASTEM route {int(sb['fxo'][1]['route'])}, branch {lab}), placed in bin 1 of a distribution "
                    f"whose bins 2..{NBIN} are empty (state 0, but a number of {EMPTY_NUM:g} /cm3 so that their kg is nonzero: every rate there must still be exactly 0, as ASTEM's dispatch and non-volatile loop skip a no_aerosol bin). "
                    "Inputs: the bin's FXO record (its state where its flux was formed, after ASTEM's state edits), "
                    "kg's inputs from AIN/BIN; expected as in the column tests. rel 1e-9.")
            tests.append(build_test(f'astem_hand_{name}', desc, c, sb, {1: sb['fxo'][1]}))
    print('hand excluded:', dict(hreasons))

    doc = json.loads(DOC.read_text())
    doc['models']['Tests']['tests'] = tests
    DOC.write_text(json.dumps(doc, indent=2) + '\n')
    print(f'{DOC.name}: {len(tests)} tests, {sum(len(t["assertions"]) for t in tests)} assertions')


if __name__ == '__main__':
    main()
