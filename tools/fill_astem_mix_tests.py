#!/usr/bin/env python3
"""Fill the inline tests of MOSAIC's mixed-bin surface equilibrium (ASTEM_flux_mix) and the
solid-surface tests that need hand bins the chem170 column never produces, from the hand
ASTEM tape mix.tape (kernels/mosaic_therm_driver mode 6 on the bins of
chem170_scm/mosaic_astem/make_astem_mix_hand.py; see
data/eqwefic/notes/chem170_mosaic_astem_capture.md).

    tools/fill_astem_mix_tests.py <astem_dir>

Writes
  components/aerosol/mosaic/astem_mixed_surface_equilibrium.esm   all tests
  components/aerosol/mosaic/astem_solid_surface_equilibrium.esm   the tests with ids
      'astemmix_*' (tools/fill_astem_tests.py keeps them): each route-3 bin's SOLID half at
      mixed = 1 (ASTEM_flux_mix's Keq_sg reset), and the all-solid bins whose Kp0 lies
      below Keq_sg (ASTEM_flux_dry case 4 at Kp0)
Only numbers are written: equations and structure are hand-authored.

Every reference is a real(kind=8) value ASTEM held at the bin's FXO record (its state after
the flux routine's tiny absorptions / degassings, the state its fluxes were formed at), so
tolerances are rel 1e-9.  The rates asserted are the continuous ones, flux_s + kg df_gas_l:
where ASTEM's signed liquid equilibrium test (FORTRAN_BUGS B31) skips an evaporating
liquid the .esm rate still evaporates it, as in astem_liquid_surface_equilibrium.esm.  A
bin inside a solid salt's 1% deadband is not asserted (ASTEM's solver decision).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MOS = ROOT / 'components' / 'aerosol' / 'mosaic'
REL = 1e-9
EXACT = {'abs': 0.0, 'rel': 0.0}
RTOL_EQB = 0.009999999776482582
ION = dict(so4=0, no3=1, cl=2, nh4=3, msa=4, co3=5, na=6, ca=7, oin=8)
JNH4NO3, JNH4CL, JNACL, JCACL2, JCACO3 = 4, 5, 11, 13, 21
XT_RICH = 1.999899983406067
PREFIX = 'astemmix_'


def a(var, val, tol):
    return {'variable': var, 'time': 0.0, 'expected': float(val), 'tolerance': tol}


def tol(scale):
    return {'rel': REL, 'abs': 1e-12 * max(abs(scale), 1e-300)}


def xt(ions):
    den = ions[ION['so4']] + 0.5 * ions[ION['msa']]
    return (ions[ION['nh4']] + ions[ION['na']] + 2.0 * ions[ION['ca']]) / den if den > 0 else -1.0


def phi(prod, k):
    return (prod - k) / max(prod, k)


def mix_case(fx):
    E, g = fx['E_s'], fx['gas']
    if E[JCACO3] > 0 or 0.0 <= xt(fx['ions_l']) < XT_RICH:
        return 'wet'
    if E[JNACL] + E[JCACL2] > 0 and g[1] > 0:
        return 'case3'
    if E[JNH4NO3] > 0 and E[JNH4CL] > 0:
        return 'case4'
    if E[JNH4NO3] > 0:
        return 'case4a'
    if E[JNH4CL] > 0:
        return 'case4b'
    return 'none'


def deadband(case, fx, ain):
    """True where a solid salt's relative driving force is inside ASTEM's 1% band: the
    Fortran switches the salt off there (solver decision), the .esm rate does not."""
    g, ksg = fx['gas'], ain['keq_sg']
    p_no3 = phi(g[3] * g[1], ksg[0])
    p_cl = phi(g[3] * g[2], ksg[1])
    if case == 'case3':
        return fx['E_s'][JNH4CL] > 0 and abs(p_cl) < RTOL_EQB
    if case == 'case4':
        return abs(p_no3) < RTOL_EQB or abs(p_cl) < RTOL_EQB
    if case == 'case4a':
        return abs(p_no3) <= RTOL_EQB
    if case == 'case4b':
        return abs(p_cl) <= RTOL_EQB
    return False


def keqs(fx, ain):
    kel, act = fx['kel'], fx['activity']
    return (kel[3] * kel[1] * act[JNH4NO3] * ain['kp_nh4no3'],
            kel[3] * kel[2] * act[JNH4CL] * ain['kp_nh4cl'])


def solid_over(fx, ain, kg, mixed):
    g, E, i_s = fx['gas'], fx['E_s'], fx['ions_s']
    return {
        'gas_hno3': g[1], 'gas_hcl': g[2], 'gas_nh3': g[3],
        'kg_hno3': kg[1], 'kg_hcl': kg[2], 'kg_nh3': kg[3],
        'so4_s': i_s[ION['so4']], 'msa_s': i_s[ION['msa']], 'nh4_s': i_s[ION['nh4']],
        'na_s': i_s[ION['na']], 'ca_s': i_s[ION['ca']], 'cl_s': i_s[ION['cl']],
        'caco3_s': E[JCACO3], 'nacl_s': E[JNACL], 'cacl2_s': E[JCACL2],
        'nh4no3_s': E[JNH4NO3], 'nh4cl_s': E[JNH4CL], 'salt_s': sum(E[:15]),
        'keq_sg_nh4no3': ain['keq_sg'][0], 'keq_sg_nh4cl': ain['keq_sg'][1],
        'kp_nh4no3_0': ain['kp_nh4no3_0'], 'kp_nh4cl_0': ain['kp_nh4cl_0'],
        'mixed': float(mixed)}


WHERE = ("Hand-specified bin '{name}' ({why}; chem170_scm/mosaic_astem/make_astem_mix_hand.py), run "
         "through WRF's own compiled ASTEM_semi_volatiles by kernels/mosaic_therm_driver mode 6 at "
         "T = 290 K, p = 0.9 atm, 1000 particles/cm3 (mix.tape)")


def solid_test(name, why, c, fx, mixed, case):
    """The solid half: astem_solid_surface_equilibrium.esm at the bin's FXO state."""
    ain, kg = c['ain'], c['bins'][1]['kg']
    big = max(abs(v) for v in fx['flux_s'][1:4]) or 1.0
    t = {'rel': REL, 'abs': 1e-12 * big}
    flag = 'case3' if case == 'case3' else 'case4'
    asr = [a(flag, 1.0, EXACT)]
    if case == 'case4a':
        # ASTEM_flux_mix_case4a forms the solid NH4NO3 surface alone (ASTEM_flux_dry_case4a): the
        # solid component's f4a, whatever its own case 4 makes of NH4Cl, which is liquid here
        asr.append(a('f4a', fx['flux_s'][1], t))
    elif case == 'case4b':
        asr.append(a('f4b', fx['flux_s'][2], t))
    else:
        asr += [a('cond_hno3', fx['flux_s'][1], t), a('cond_hcl', fx['flux_s'][2], t),
                a('cond_nh3', fx['flux_s'][3], t)]
    k0 = ((ain['keq_sg'][0], ain['keq_sg'][1]) if mixed else
          (min(ain['kp_nh4no3_0'], ain['keq_sg'][0]), min(ain['kp_nh4cl_0'], ain['keq_sg'][1])))
    asr += [a('keq0_nh4no3', k0[0], {'rel': 0.0, 'abs': 0.0}), a('keq0_nh4cl', k0[1], {'rel': 0.0, 'abs': 0.0})]
    # the surface NH3 its quadratic gave, where nothing later overwrote sfc_a(NH3): a single
    # active salt with no liquid surface after it (case 4ab and 3a+3b recompute or keep it)
    s_nh3 = fx['sfc'][3]
    act_no3, act_cl = fx['int_s'][1] > 0, fx['int_s'][2] > 0
    liquid_after = mixed and (fx['int_l'][1] or fx['int_l'][2] or fx['df_l'][1] or fx['df_l'][2])
    if case == 'case3' and fx['int_s'][3] and not liquid_after:
        asr.append(a('s3_nh3', s_nh3, tol(s_nh3)))
    elif case.startswith('case4') and not liquid_after:
        if act_no3 and not act_cl:
            asr.append(a('s4a_nh3', s_nh3, tol(s_nh3)))
        elif act_cl and not act_no3:
            asr.append(a('s4b_nh3', s_nh3, tol(s_nh3)))
    over = solid_over(fx, ain, kg, mixed)
    what = ('the SOLID half of ASTEM_flux_mix (route 3, mixed = 1: the solid products at Keq_sg)' if mixed else
            'ASTEM_flux_dry case 4 on an all-solid bin whose Kp0 is below Keq_sg (keq0 = Kp0)')
    return {
        'id': f'{PREFIX}dry_{name}',
        'description': (WHERE.format(name=name, why=why) + f'. Tests {what}. Inputs are the state the bin\'s '
                        'flux routine saw (FXO record, after its tiny absorptions/degassings), kg and the constants '
                        'from the AIN/BIN records; expected are the solid fluxes flux_s it formed (and the surface NH3 '
                        'its quadratic gave where no later surface overwrote it). rel 1e-9, absolute floor 1e-12 of the '
                        'largest of the three fluxes; keq0 exactly.'),
        'parameter_overrides': over, 'time_span': {'start': 0, 'end': 1}, 'assertions': asr}


def mix_test(name, why, c, fx, case):
    ain, kg = c['ain'], c['bins'][1]['kg']
    g, E, il, fs = fx['gas'], fx['E_s'], fx['ions_l'], fx['flux_s']
    k_no3, k_cl = keqs(fx, ain)
    # the solid half enters as coupling targets; here it is the solid flux ASTEM formed
    # (the solid component's own values are tested in its astemmix_dry_* tests)
    sol = {'acid_hno3': 0.0, 'f3b': 0.0, 'f4a': 0.0, 'f4b': 0.0, 'scale4': 1.0, 'nh3_need': 0.0}
    if case == 'case3':
        sol.update(acid_hno3=fs[1], f3b=fs[3])
    elif case == 'case4':
        sol.update(f4a=fs[1], f4b=fs[2], nh3_need=fs[3])
    elif case == 'case4a':
        sol.update(f4a=fs[1])
    elif case == 'case4b':
        sol.update(f4b=fs[2])
    over = {
        'gas_hno3': g[1], 'gas_hcl': g[2], 'gas_nh3': g[3],
        'kg_hno3': kg[1], 'kg_hcl': kg[2], 'kg_nh3': kg[3],
        'so4_l': il[ION['so4']], 'msa_l': il[ION['msa']], 'nh4_l': il[ION['nh4']],
        'na_l': il[ION['na']], 'ca_l': il[ION['ca']],
        'caco3_s': E[JCACO3], 'nacl_s': E[JNACL], 'cacl2_s': E[JCACL2],
        'nh4no3_s': E[JNH4NO3], 'nh4cl_s': E[JNH4CL],
        'keq_nh4no3': k_no3, 'keq_nh4cl': k_cl,
        'cond_l_hno3': 0.0, 'cond_l_hcl': 0.0, 'cond_l_nh3': 0.0}
    over.update(sol)
    asr = [a(f, 1.0 if f == case else 0.0, EXACT) for f in ('wet', 'case3', 'case4', 'case4a', 'case4b')]
    liq = {'case4a': ('on_l_nh4cl', 'scl_nh3', 'scl_hcl', 'lcl_hcl', 'lcl_nh3', 2),
           'case4b': ('on_l_nh4no3', 'sno3_nh3', 'sno3_hno3', 'lno3_hno3', 'lno3_nh3', 1)}.get(case)
    if liq:
        on_v, s_n, s_a, r_a, r_n, ia = liq
        on = 1.0 if (g[3] * g[ia] != 0 and (k_cl if ia == 2 else k_no3) != 0) else 0.0
        asr.append(a(on_v, on, EXACT))
        if on:
            # the liquid surface is formed last, so sfc_a(NH3) and sfc_a(acid) are its own
            for v, val, sc in ((s_n, fx['sfc'][3], fx['sfc'][3]), (s_a, fx['sfc'][ia], fx['sfc'][ia]),
                               (r_a, kg[ia] * fx['df_l'][ia], kg[ia] * max(abs(g[ia]), abs(fx['sfc'][ia]))),
                               (r_n, kg[3] * fx['df_l'][3], kg[3] * max(abs(g[3]), abs(fx['sfc'][3])))):
                asr.append(a(v, val, tol(sc)))
        else:
            asr += [a(r_a, 0.0, EXACT), a(r_n, 0.0, EXACT)]
    big = max(abs(fs[i] + kg[i] * fx['df_l'][i]) for i in (1, 2, 3)) or 1.0
    for gname, i in (('hno3', 1), ('hcl', 2), ('nh3', 3)):
        asr.append(a(f'cond_{gname}', fs[i] + kg[i] * fx['df_l'][i], {'rel': REL, 'abs': 1e-12 * big}))
    skipped = [n for n, i in (('HNO3', 1), ('HCl', 2), ('NH3', 3)) if fx['df_l'][i] != 0 and not fx['int_l'][i]]
    note = ('' if not skipped else
            f" ASTEM's signed liquid equilibrium test skipped the evaporating liquid ({', '.join(skipped)}: "
            "FORTRAN_BUGS B31); the rate asserted is the continuous one, which evaporates it.")
    return {
        'id': f'{PREFIX}{case}_{name}',
        'description': (WHERE.format(name=name, why=why) + f' (route 3: ASTEM_flux_mix, mixed case {case[4:]}). '
                        'Inputs are the state ASTEM_flux_mix formed its fluxes at (FXO record, after its tiny '
                        'absorptions/degassings), kg (BIN record), the solution products kel kel a kp (AIN/FXO); the '
                        "solid half enters as the solid fluxes ASTEM formed (the solid component's astemmix_dry_* tests "
                        'check those). Expected: the case flags, the liquid single-salt surface (sfc_a, kg df_gas_l) '
                        'and the bin rates flux_s + kg df_gas_l. rel 1e-9, absolute floor 1e-12 of the largest rate.'
                        + note),
        'parameter_overrides': over, 'time_span': {'start': 0, 'end': 1}, 'assertions': asr}


def write_tests(fname, model, tests, keep):
    p = MOS / fname
    doc = json.loads(p.read_text())
    m = doc['models'][model]
    m['tests'] = [t for t in m.get('tests', []) if keep(t)] + tests
    p.write_text(json.dumps(doc, indent=2) + '\n')
    print(f'{fname}: {len(m["tests"])} tests, {sum(len(t["assertions"]) for t in m["tests"])} assertions '
          f'({len(tests)} written here)')


def main():
    d = Path(sys.argv[1])
    sys.path.insert(0, str(d))
    from astem_tape import read_tape
    names = {}
    for line in open(d / 'cases_mix_names.txt'):
        i, name, why = line.rstrip('\n').split('\t')
        names[int(i)] = (name, why)
    calls = read_tape(d / 'mix.tape')
    mixt, solt, skipped = [], [], []
    for c in calls:
        name, why = names[c['lev']]
        fx = c['subs'][0]['fxo'][1]
        route = int(fx['route'])
        if route == 3:
            case = mix_case(fx)
            if deadband(case, fx, c['ain']):
                skipped.append((name, 'solid 1% deadband'))
                continue
            mixt.append(mix_test(name, why, c, fx, case))
            if case in ('case3', 'case4', 'case4a', 'case4b'):
                solt.append(solid_test(name, why, c, fx, 1, case))
        elif route == 1 and name.startswith('kp0_'):
            solt.append(solid_test(name, why, c, fx, 0, 'case4'))
        else:
            skipped.append((name, f'route {route}'))
    print('skipped:', skipped)
    write_tests('astem_mixed_surface_equilibrium.esm', 'MosaicAstemMixedSurfaceEquilibrium', mixt,
                keep=lambda t: False)
    write_tests('astem_solid_surface_equilibrium.esm', 'MosaicAstemSolidSurfaceEquilibrium', solt,
                keep=lambda t: not t['id'].startswith(PREFIX))


if __name__ == '__main__':
    main()
