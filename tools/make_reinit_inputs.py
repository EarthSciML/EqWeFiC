#!/usr/bin/env python3
"""Extract the 80x80 input fields of LevelSetReinitWRFFire (the stage level set psi = lfn_curr and
the frozen reference psi_ref the smoothed sign is built from) from `ls_reinit` dumps of a
fire_behavior standalone run into a mount-edge template library.  Numbers only.

Usage: make_reinit_inputs.py <dump_psi.json> <dump_psiref.json> <out.esm> <name> <description>
       the two dumps are the same file for Runge-Kutta stage 1, where psi_ref == psi.
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import esm_dump

dpsi, dref, out, name, desc = sys.argv[1:6]
HALO, N = 4, 80


def field(path, key):
    a = np.asarray(esm_dump.load(path)['vars'][key])[HALO + 1:HALO + 1 + N, HALO + 1:HALO + 1 + N]
    return [[float(a[i, j]) for j in range(N)] for i in range(N)]


def tmpl(op, value):
    return {"params": [], "match": {"op": op, "args": []},
            "body": {"op": "faq", "output_idx": ["i", "j"], "args": [],
                     "ranges": {"i": {"from": "x"}, "j": {"from": "y"}},
                     "expr": {"op": "index", "args": [{"op": "const", "args": [], "value": value}, "i", "j"]}}}


doc = {"esm": "1.1.0",
       "metadata": {"name": name, "description": desc, "license": "MIT",
                    "tags": ["test-inputs", "generated"]},
       "expression_templates": {
           "input_psi": tmpl("input_psi", field(dpsi, 'lfn_curr')),
           "input_psi_ref": tmpl("input_psi_ref", field(dref, 'lfn_curr')),
       }}
with open(out, 'w') as f:
    json.dump(doc, f); f.write('\n')
print('wrote', out)
