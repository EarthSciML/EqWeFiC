#!/usr/bin/env python3
"""Fill the `mdrh_table_*` tests of components/aerosol/mosaic/aerosol_phase_state.esm: all 63
mutual-deliquescence curves mdrh_t(1:63) at a given temperature.

Sources (numbers only):
  - kernels/mosaic_therm_driver run once per temperature of a sweep
    (chem170_scm/therm64/mdrh/run_mdrh_sweep.sh; the driver dumps the table of its last
    case as mdrh_t_last);
  - the chem170 aerosol_phase_state capture (dumps/chem170_mosaic_phase, ps_mdrh_t): the
    table as update_thermodynamic_constants left it at the step's last MOSAIC call, which is
    the top level (59) of the column (checked: that level's t_k reproduces it bit for bit;
    the fill refuses a step where it does not).
Other tests of the component are left untouched.

  fill_mosaic_mdrh_table_tests.py <data/eqwefic>
"""
import glob, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = os.path.join(HERE, '..', 'components', 'aerosol', 'mosaic', 'aerosol_phase_state.esm')
PP = os.path.join(HERE, '..', 'components', 'aerosol', 'mosaic', 'mosaic_phase_parameters.esm')


def table(t, pp):
    D = np.array(pp['d_mdrh']['default'])
    return np.array([pp['mdrh_ca']['default'] if (j in (7, 8) or 34 <= j <= 51)
                     else D[j - 1, 0] + t * (D[j - 1, 1] + t * (D[j - 1, 2] + t * D[j - 1, 3])) + 1.0
                     for j in range(1, 64)])


def test(tid, t, m, src):
    A = [{'variable': 'mdrh_table', 'coords': {'mdrhcurve': j + 1}, 'time': 0.0, 'expected': float(m[j]),
          'tolerance': {'rel': 1e-9, 'abs': 0.0}} for j in range(63)]
    return {'id': tid, 'description': (
        f"All 63 mutual deliquescence curves at t_k = {t!r} K: mdrh_t(1:63) as update_thermodynamic_constants "
        f"fills them (drh_mutual). {src} Real64 (mdrh_t is real(kind=8), the coefficients binary32); rel 1e-9. "
        "Only the table is asserted; the bin is left at its defaults."),
        'parameter_overrides': {'t_k': float(t)}, 'time_span': {'start': 0, 'end': 1}, 'assertions': A}


def main(data):
    pp = json.load(open(PP))['models']['MosaicPhaseParameters']['variables']
    tests, seen = [], set()
    for f in sorted(glob.glob(os.path.join(data, 'chem170_scm', 'therm64', 'mdrh', 'mdrh_*.json'))):
        v = json.load(open(f))['vars']
        t = float(os.path.basename(f)[5:-5])
        m = np.array(v['mdrh_t_last']['data'])[:63]
        assert np.array_equal(table(t, pp), m) or np.max(np.abs(table(t, pp) - m) / np.abs(m)) < 1e-12, f
        tests.append(test(f'mdrh_table_driver_{t:g}K', t, m,
                          f"Source: kernels/mosaic_therm_driver at {t:g} K (chem170_scm/therm64/mdrh, WRF fork "
                          "earthsciml-instrumented-chem170)."))
        seen.add(t)
    for f in sorted(glob.glob(os.path.join(data, 'dumps', 'chem170_mosaic_phase', '*.json')),
                    key=lambda p: int(p.split('_')[-1][:-5])):
        v = json.load(open(f))['vars']
        step = int(f.split('_')[-1][:-5])
        t = float(v['ps_t_k']['data'][58])
        m = np.array(v['ps_mdrh_t']['data'])[:63]
        if not np.array_equal(table(t, pp), m):
            sys.exit(f'step {step}: level 59 t_k does not reproduce ps_mdrh_t')
        if t in seen:
            continue
        seen.add(t)
        tests.append(test(f'mdrh_table_col_s{step}', t, m,
                          f"Source: the chem170 aerosol_phase_state capture, step {step}, the table as left by "
                          "the column's last MOSAIC call (level 59)."))
    doc = json.load(open(COMP))
    mdl = doc['models']['MosaicAerosolPhaseState']
    mdl['tests'] = [t for t in mdl['tests'] if not t['id'].startswith('mdrh_table_')] + tests
    with open(COMP, 'w') as fh:
        json.dump(doc, fh, indent=2); fh.write('\n')
    print(len(tests), 'table tests,', 63 * len(tests), 'assertions')


if __name__ == '__main__':
    main(sys.argv[1])
