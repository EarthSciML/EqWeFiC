#!/usr/bin/env python3
"""Extract the input fields of FireWindHinterpWRFFire -- the atmosphere-grid winds and the fractional
atmosphere-grid coordinates of each fire cell -- from the `hinterp` dumps of a fire_behavior
standalone run into a mount-edge template library.  Numbers only.

Usage: make_hinterp_inputs.py <dump_ua.json> <dump_va.json> <out.esm> <name> <description>
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import esm_dump

dua, dva, out, name, desc = sys.argv[1:6]
HALO, N = 4, 80
A = esm_dump.load(dua)['vars']
B = esm_dump.load(dva)['vars']
tile = lambda a: np.asarray(a)[HALO + 1:HALO + 1 + N, HALO + 1:HALO + 1 + N]


def grid2(a):
    return [[float(a[i, j]) for j in range(a.shape[1])] for i in range(a.shape[0])]


def tmpl(op, value):
    return {"params": [], "match": {"op": op, "args": []},
            "body": {"op": "faq", "output_idx": ["i", "j"], "args": [],
                     "ranges": {"i": {"from": IDX[op][0]}, "j": {"from": IDX[op][1]}},
                     "expr": {"op": "index", "args": [{"op": "const", "args": [], "value": value}, "i", "j"]}}}


IDX = {"input_u_src": ("src_x", "src_y"), "input_v_src": ("src_x", "src_y"),
       "input_i_real": ("x", "y"), "input_j_real": ("x", "y")}
doc = {"esm": "1.1.0",
       "metadata": {"name": name, "description": desc, "license": "MIT",
                    "tags": ["test-inputs", "generated"]},
       "expression_templates": {
           "input_u_src": tmpl("input_u_src", grid2(np.asarray(A['data_in']))),
           "input_v_src": tmpl("input_v_src", grid2(np.asarray(B['data_in']))),
           "input_i_real": tmpl("input_i_real", grid2(tile(A['i_real']))),
           "input_j_real": tmpl("input_j_real", grid2(tile(A['j_real']))),
       }}
with open(out, 'w') as f:
    json.dump(doc, f); f.write('\n')
print('wrote', out)
