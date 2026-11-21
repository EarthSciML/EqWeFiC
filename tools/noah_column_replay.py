#!/usr/bin/env python3
"""Real64 references for couplings/noah_column.esm (sfclayrev -> Noah on the chem170 column).

For each chemistry step S of the full-met capture (dumps/chem170_full_met_sfc, the same run as
dumps/chem170_full_met: every noah dump byte-identical, the sfclayrev/ysu dumps differing only in
uninitialised scratch inputs) this

1. replays sfclayrev call 2S-1 (row j = 1 of ktau S: the surface layer is dumped per j-row call,
   two rows on em_scm_xy) through the real64 kernels/sfclayrev_driver, twice: with the dumped
   arguments as they are (sfc64_S.json), and with the thirteen physical constants WRF passes the
   kernel (cp, g, rovcp, r, xlv, svp1-3, svpt0, ep1, ep2, karman, p1000mb) set to their exact
   decimal values instead of WRF's real32 copies (sfc64x_S.json). The second is the reference the
   coupled document is tested against: sfclayrev.esm documents that the real32 copies (svpt0 =
   273.1499939 K above all) move the kernel by up to 4e-6 relative, which would otherwise reach
   every Noah quantity through chs and bound the whole column at that level. The two replays'
   difference, propagated through Noah, is reported rather than asserted. The second replay runs
   the LITERAL-EXACT build of the kernel (kernels/build_sfclit: sf_sfclayrev.F90 compiled with
   -fdefault-real-8 -fdefault-double-8, everything else as `make`). The CCPP kernel declares its
   variables kind_phys but writes its literals in default real, so in the ordinary real64 build
   `.33`, `.32` and `0.001` are still binary32 (0.330000013, 0.319999993, 0.00100000005). Those
   moved br by 2.0e-5 and mol by 1.9e-4 relative at the night steps (thx - thgb cancellation, via
   vconv**.33 in wspd), the offset first read as the zolri secant (2026-11-20). sfclayrev.esm writes
   the decimals, which is the scheme as published;
2. rebuilds the Noah SFLX arguments that the Noah driver (module_sf_noahdrv.F) forms from the
   atmosphere and the surface layer, IN REAL64, from the same raw inputs the coupled document
   reads (t3d1, qv3d1, p8w(1), sfcprs -> p8w(2), glw, swdown, the entering albedo and emissivity,
   rainbl) and the real64 surface layer's chs and qgh, so that the SFLX reference is the one the
   coupled document's own inputs imply and the real32/real64 surface-layer difference (up to
   4e-6, sfclayrev.esm) does not leak into the Noah references; the rebuilt values are compared
   with WRF's real32 ones and the worst relative difference is reported;
3. replays SFLX through the real64 kernels/noah_driver at several ESM_DT.

This script transcribes numbers only; it writes no equations.
    python3 tools/noah_column_replay.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import esm_dump  # noqa: E402

D = '/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic'
DUMPS = f'{D}/dumps/chem170_full_met_sfc'
OUT = f'{D}/noah_column'
K = f'{D}/wrf-chem170/kernels/build'
KL = f'{D}/wrf-chem170/kernels/build_sfclit'   # literal-exact sfclayrev (see the docstring)
STEPS = (362, 364, 1024, 1026)
DTS = ('1e-6', '1e-4', '0.01', '1', '60')

# module_sf_noahdrv.F: A2 = 17.67, A3 = 273.15, A4 = 29.65, A23M4 = A2 (A3 - A4); CAPA = R_D/CP.
A2, A3, A4 = 17.67, 273.15, 29.65
R_D, CP = 287.0, 7.0 * 287.0 / 2.0
R_V = 461.6
# lib/wrf_constants.esm's values of the constants WRF hands sf_sfclayrev_run (as real32 copies)
EXACT = {'cp': CP, 'g': 9.81, 'rovcp': R_D / CP, 'r': R_D, 'xlv': 2.5e6, 'svp1': 0.6112, 'svp2': 17.67,
         'svp3': 29.65, 'svpt0': 273.15, 'ep1': R_V / R_D - 1.0, 'ep2': R_D / R_V, 'karman': 0.4, 'p1000mb': 1.0e5}


def preamble(n: dict, sfc64: dict) -> dict:
    """The Noah driver's real64 SFLX arguments from the raw inputs (see module docstring)."""
    psfc = float(n['p8w1'])
    sfcprs = float(n['sfcprs'])            # p8w(2) is supplied as 2 sfcprs - p8w(1): exact
    sfctmp = float(n['t3d1'])
    qgh = float(sfc64['qgh_out'][0])
    q2sat = qgh / (1.0 + qgh)
    capa = R_D / CP
    return {
        'sfcprs': sfcprs, 'sfctmp': sfctmp,
        'q2k': float(n['qv3d1']) / (1.0 + float(n['qv3d1'])),
        'q2sat': q2sat,
        'dqsdt2': q2sat * A2 * (A3 - A4) / (sfctmp - A4) ** 2,
        'th2': sfctmp * (1.0e5 / sfcprs) ** capa / (1.0e5 / psfc) ** capa,
        'lwdn': float(n['glw']) * float(n['emissi']),
        'soldn': float(n['swdown']),
        'solnet': float(n['swdown']) * (1.0 - float(n['albedo_out'])),
        'prcp': float(n['rainbl']) / float(n['dt']),
        'ch': float(sfc64['chs_out'][0]),
    }


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    worst = {}
    for S in STEPS:
        sd = esm_dump.load(f'{DUMPS}/esm_dump_sfclayrev_{2 * S - 1}.json')
        assert int(sd['vars']['j']) == 1, 'call 2S-1 must be row j = 1'
        esm_dump.to_flat(sd, f'{OUT}/sfc_{S}.flat')
        subprocess.run([f'{K}/sfclayrev_driver', f'{OUT}/sfc_{S}.flat', f'{OUT}/sfc64_{S}.json'], check=True)
        sd['vars'].update({k: np.float64(v) for k, v in EXACT.items()})
        esm_dump.to_flat(sd, f'{OUT}/sfcx_{S}.flat')
        subprocess.run([f'{KL}/sfclayrev_driver', f'{OUT}/sfcx_{S}.flat', f'{OUT}/sfc64x_{S}.json'], check=True)
        sfc64 = esm_dump.load(f'{OUT}/sfc64x_{S}.json')['vars']

        nd = esm_dump.load(f'{DUMPS}/esm_dump_noah_{S}.json')
        n = nd['vars']
        # the coupling claim this rests on: Noah's CH and QGH are the surface layer's of the same ktau
        assert float(n['chs']) == float(sd['vars']['chs_out'][0]) and float(n['qgh']) == float(sd['vars']['qgh_out'][0])
        pre = preamble(n, sfc64)
        wrf32 = {'sfcprs': n['sfcprs'], 'sfctmp': n['sfctmp'], 'q2k': n['q2k'], 'q2sat': n['q2sat'],
                 'dqsdt2': n['dqsdt2'], 'th2': n['th2'], 'lwdn': n['lwdn'], 'soldn': n['soldn'],
                 'solnet': n['solnet'], 'prcp': n['prcp'], 'ch': n['chs']}
        for k, v in pre.items():
            ref = float(wrf32[k])
            r = abs(v - ref) / max(abs(ref), 1e-300) if ref != 0 else abs(v)
            worst[k] = max(worst.get(k, 0.0), r)
            n[k] = np.float64(v)
        esm_dump.to_flat(nd, f'{OUT}/noah_{S}.flat')
        for dt in DTS:
            subprocess.run([f'{K}/noah_driver', f'{OUT}/noah_{S}.flat', f'{OUT}/noah64_{S}_{dt}.json'],
                           check=True, env={**os.environ, 'ESM_DT': dt})
        print(f'step {S}: sfclayrev call {2 * S - 1} and noah replayed at ESM_DT = {", ".join(DTS)}')
    print('rebuilt real64 SFLX arguments vs WRF real32, worst relative difference over the steps:')
    for k, v in worst.items():
        print(f'  {k:8s} {v:.3e}')
    with open(f'{OUT}/preamble_vs_wrf32.json', 'w') as f:
        json.dump(worst, f, indent=1)


if __name__ == '__main__':
    main()
