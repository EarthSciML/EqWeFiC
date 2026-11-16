#!/usr/bin/env python3
"""Fill the inline tests of components/aerosol/mosaic/section_moving.esm.

Inputs per test (one capture step, the whole column): the state move_sections
starts from -- the PRE-growth section state rsub_in (number and the ten
components + hysteresis water + water), the pre-growth dry mass and density
(ms_grow rows 1-2) -- and the step's dry-volume growth rate
dvdt = (drymass_aft/drydens_aft - drymass_pre/drydens_pre)/dtchem.
Expected values: the dt -> 0 limit of WRF's own move_sections, from the real64
replays (wrf-chem170/kernels/movesect_driver) at growth fractions d1 and d2,
Richardson-extrapolated:  R = (d2 R(d1) - d1 R(d2))/(d2 - d1),
R(d) = (rsub_exit - rsub_entry)/(d dtchem).

usage: fill_mosaic_movesect_tests.py <dump_dir> <replay_dir> <component.esm> step [step ...]
The numbers are extracted mechanically; the component and the test structure
are hand-written (CLAUDE.md).
"""
import json
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import esm_dump  # noqa: E402
from esm_compact_json import dumps  # noqa: E402

D1, D2 = '1e-2', '2e-2'
REL = 1.0e-5          # relative bound
FLOOR = 1.0e-6        # absolute floor, as a fraction of the largest |rate| of that quantity in the test


def limit(rep_dir, step, dt):
    """Richardson limit of the frozen-composition replays (ESM_MS_FREEZE=1)."""
    r1 = esm_dump.load(f'{rep_dir}/f_{step}_{D1}.json')['vars']
    r2 = esm_dump.load(f'{rep_dir}/f_{step}_{D2}.json')['vars']
    d1, d2 = float(D1), float(D2)
    R1 = (r1['ms_rsub_aft_replay'] - r1['ms_rsub_bef_replay'])/(d1*dt)
    R2 = (r2['ms_rsub_aft_replay'] - r2['ms_rsub_bef_replay'])/(d2*dt)
    return (d2*R1 - d1*R2)/(d2 - d1)


DESC = (
    "Model step {step}, the whole 59-level column x 8 sections. Reference: the em_scm_xy chem_opt = 170 column "
    "(data/eqwefic/chem170_scm, reference namelist unchanged), move_sections capture dumps/chem170_mosaic_movesect "
    "(WRF fork earthsciml-instrumented-chem170 711bb56a; bit-identical to dumps/chem170_mosaic_coag on all 5216 shared "
    "fields). Inputs are move_sections' own entry state: the PRE-growth section number and the twelve volume-carried "
    "rows (rsub_in: the ten components, the hysteresis water, the water), the pre-growth dry mass and density "
    "(drymass_pregrow, drydens_pregrow) and the step's dry-volume growth rate dvdt = (drymass_aftgrow/drydens_aftgrow - "
    "drymass_pregrow/drydens_pregrow)/dtchem. Expected values: the dt -> 0 limit of WRF's own move_sections, from the "
    "REAL64 replay (wrf-chem170/kernels/movesect_driver, chem170_scm/movesect64) with the step's growth scaled to the "
    "fractions 1e-2 and 2e-2 and the rows held at their pre-growth values (ESM_MS_FREEZE=1), Richardson-extrapolated. "
    "Real64 because WRF's real32 call forms each transferred fraction as 1 minus a number close to 1 (FORTRAN_BUGS N125). "
    "Tolerance measured over all 16 capture steps (98,176 entries), not assumed: every entry within rel 1e-5 or 1e-6 of "
    "the largest rate of its quantity in the step (worst 0.13 of that bound); the bound is set by the replay's own "
    "round-off and O(delta^2) truncation, not by the operator."
)

HAND = (
    "PERTURBED entry state, built by tools/perturb_mosaic_movesect.py from capture step 918 (chem170_scm/movesect64/"
    "{step}.json), to reach the branches the chem170 column never takes -- every column section holds its mean volume "
    "in the lower or middle third, is active, and is clear of the 1 % edge band. Levels 1-10, sections 2-4 unless "
    "stated: (1) mean at 0.85 of the width, the upper-third triangle, growing; (2) mean at 0.5, the full linear "
    "reconstruction, growing; (3) mean 0.5 % below the upper edge and (4) 0.5 % above the lower edge, inside MOSAIC's "
    "1 % band (no transfer); (5) pre-growth dry density -1, MOSAIC's no_aerosol (no transfer); (6) dry mass 1e-21 "
    "g/mol, below smallmassaa (no transfer; its water set to 0, since MOSAIC's initial conform zeroes it, a projection "
    "outside this operator); (7) as 1 and (8) as 2 but shrinking by 1e-3 of the dry volume; (9) the top section "
    "growing and (10) the bottom section shrinking, both at mean 0.5 (closed ends). WRF's own move_sections, replayed "
    "in real64 on these inputs exactly as for the capture steps, is the reference; same extraction and tolerance."
)

def main():
    dump_dir, rep_dir, comp, steps = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    doc = json.load(open(comp))
    model = doc['models']['MosaicSectionMoving']
    tests = []
    for step in steps:
        src = (f'{rep_dir}/{step}.json' if step.startswith('hand')
               else f'{dump_dir}/esm_dump_mosaic_{step}.json')
        v = esm_dump.load(src)['vars']
        dt = float(v['dtchem'])
        nb, nc = int(v['nsize_aer']), int(v['ncomp_plustracer_aer'])
        mp = np.asarray(v['massptr_aer']).reshape(nc, nb, order='F') - 1
        npt = np.asarray(v['numptr_aer']) - 1
        rows = np.vstack([mp, (np.asarray(v['hyswptr_aer']) - 1)[None, :],
                          (np.asarray(v['waterptr_aer']) - 1)[None, :]])      # [trc, bin]
        rin, grow = v['rsub_in'], v['ms_grow']
        nk = rin.shape[1]
        R = limit(rep_dir, step, dt)
        num = rin[npt, :].T                                   # [lev, bin]
        q = np.transpose(rin[rows, :], (2, 1, 0))             # [lev, bin, trc]
        mpre, dpre, maft, daft = (grow[i].T for i in range(4))
        with np.errstate(divide='ignore', invalid='ignore'):
            dvdt = np.where((dpre > 0) & (daft > 0), (maft/daft - mpre/dpre)/dt, 0.0)
        dnum = R[npt, :].T
        dq = np.transpose(R[rows, :], (2, 1, 0))
        asserts = []
        fl = FLOOR*np.abs(dnum).max()
        for k in range(nk):
            for j in range(nb):
                asserts.append({'variable': f'dnum_dt[{k+1},{j+1}]', 'time': 0.0, 'expected': float(dnum[k, j]),
                                'tolerance': {'rel': REL, 'abs': float(fl)}})
        for s_ in range(q.shape[2]):
            fl = FLOOR*np.abs(dq[:, :, s_]).max()
            for k in range(nk):
                for j in range(nb):
                    asserts.append({'variable': f'dq_dt[{k+1},{j+1},{s_+1}]', 'time': 0.0,
                                    'expected': float(dq[k, j, s_]), 'tolerance': {'rel': REL, 'abs': float(fl)}})
        tests.append({
            'id': f'movesect_s{step}',
            'description': (HAND if step.startswith('hand') else DESC).format(step=step),
            'time_span': {'start': 0.0, 'end': 0.0},
            'parameter_overrides': {
                'num': num.tolist(), 'q': q.tolist(),
                'drymass': mpre.tolist(), 'drydens': dpre.tolist(), 'dvdt': dvdt.tolist()},
            'assertions': asserts})
    model['tests'] = tests
    with open(comp, 'w') as f:
        f.write(dumps(doc) + '\n')
    print(sum(len(t['assertions']) for t in tests), 'assertions')


if __name__ == '__main__':
    main()
