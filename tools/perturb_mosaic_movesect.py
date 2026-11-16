#!/usr/bin/env python3
"""Write a PERTURBED move_sections entry state that reaches the branches the
chem170 column never does, for movesect_driver to replay.

The column's sections all hold their mean volume in the lower or middle third
of the section, every section is active, and none sits within 1 % of an edge
(data/eqwefic/notes/chem170_mosaic_movesect_capture.md).  Starting from one
capture step, this edits the pre-growth number (rsub_in at numptr_aer) and the
pre/after-growth dry mass and density (ms_grow) of chosen sections of chosen
levels; the components are left alone.  WRF's own move_sections, replayed on
the result, is still the reference.

    perturb_mosaic_movesect.py <dump.json> <out_prefix>
writes <out_prefix>.json (an esm dump with the edited fields) and
<out_prefix>.flat (the driver input).

Cases (level: sections 2-4 unless stated; levels 1-based):
  1  mean at 0.85 of the section width: upper-third triangle, growing
  2  mean at 0.50: full-section linear reconstruction, growing
  3  mean 0.5 % below the upper edge: inside MOSAIC's 1 % edge band
  4  mean 0.5 % above the lower edge: inside the 1 % edge band
  5  pre-growth dry density -1 (MOSAIC's no_aerosol): does not move
  6  the whole section scaled to a dry mass of 1e-21 g/mol (pre and after
     growth), below smallmassaa: does not move.  Its water is set to 0, since
     MOSAIC's initial conform zeroes the water of such a section, a projection
     this rate operator does not contain.
  7  as 1, shrinking by 1e-3 of the dry volume over the step
  8  as 2, shrinking by 1e-3
  9  section 8 (the top) at 0.5, growing: the closed upper end
 10  section 1 (the bottom) at 0.5, shrinking: the closed lower end
"""
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import esm_dump  # noqa: E402


def main():
    src, out = sys.argv[1], sys.argv[2]
    d = esm_dump.load(src)
    v = d['vars']
    npt = np.asarray(v['numptr_aer']) - 1
    wpt = np.asarray(v['waterptr_aer']) - 1
    hpt = np.asarray(v['hyswptr_aer']) - 1
    nc_ = int(v['ncomp_plustracer_aer']); nb_ = int(v['nsize_aer'])
    mpt = np.asarray(v['massptr_aer']).reshape(nc_, nb_, order='F') - 1
    vlo = np.asarray(v['volumlo_sect'], dtype=float)
    vhi = np.asarray(v['volumhi_sect'], dtype=float)
    rin = np.array(v['rsub_in'], dtype=float)
    grow = np.array(v['ms_grow'], dtype=float)

    def place(k, j, frac):
        """Set the number of section j, level k, so its mean volume sits at
        vlo + frac (vhi - vlo), from the section's pre-growth dry volume."""
        V = grow[0, j, k]/grow[1, j, k]
        rin[npt[j], k] = V/(vlo[j] + frac*(vhi[j] - vlo[j]))

    def growth(k, j, factor):
        """After-growth dry volume = factor x pre-growth, density unchanged."""
        grow[2, j, k] = grow[0, j, k]*factor
        grow[3, j, k] = grow[1, j, k]

    mid = (1, 2, 3)
    for j in mid:
        place(0, j, 0.85); growth(0, j, 1.001)
        place(1, j, 0.50); growth(1, j, 1.001)
        place(2, j, 0.995); growth(2, j, 1.001)
        place(3, j, 0.005); growth(3, j, 0.999)
        place(4, j, 0.50); growth(4, j, 1.001); grow[1, j, 4] = -1.0
        # case 6: shrink the whole section to a dry mass of 1e-21 g/mol --
        # components and hysteresis water scaled together, water zeroed -- so
        # its own volume still puts the mean mid-section (MOSAIC's initial
        # conform then has nothing to clip)
        f6 = 1.0e-21/grow[0, j, 5]
        for r_ in list(mpt[:, j]) + [hpt[j]]:
            rin[r_, 5] *= f6
        rin[wpt[j], 5] = 0.0
        grow[0, j, 5] = 1.0e-21; grow[2, j, 5] = 1.001e-21; grow[3, j, 5] = grow[1, j, 5]
        place(5, j, 0.50)
        place(6, j, 0.85); growth(6, j, 0.999)
        place(7, j, 0.50); growth(7, j, 0.999)
    place(8, 7, 0.50); growth(8, 7, 1.001)
    place(9, 0, 0.50); growth(9, 0, 0.999)
    v['rsub_in'] = rin
    v['ms_grow'] = grow
    v.pop('coag_fink', None)
    esm_dump.to_flat(d, out + '.flat')
    # an esm-dump-format copy of the edited inputs for the fill script
    import json
    doc = json.load(open(src))
    doc['vars']['rsub_in']['data'] = [float(x) for x in rin.reshape(-1, order='F')]
    doc['vars']['rsub_in']['kind'] = 'real64'
    doc['vars']['ms_grow']['data'] = [float(x) for x in grow.reshape(-1, order='F')]
    doc['vars']['ms_grow']['kind'] = 'real64'
    with open(out + '.json', 'w') as f:
        json.dump(doc, f)


if __name__ == '__main__':
    main()
