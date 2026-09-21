#!/usr/bin/env python3
"""Extract the 80x80 input fields of LevelSetTendencyWRFFire (psi = lfn after Extrapol_var_at_bdys,
ros, and the mid-flame wind uf/vf) from one `ls_tend` dump of a fire_behavior standalone run into a
mount-edge template library.  Numbers only: the four `input_*` template bodies are the same shape as
the hand-authored libraries this replaces, and nothing about the component's structure is generated.

Usage: make_ls_inputs.py <esm_dump_ls_tend_N.json> <out.esm> <library-name> <description>
"""
import json, sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import esm_dump

dump, out, name, desc = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
v = esm_dump.load(dump)['vars']
HALO = 4          # ifms = -4, so fire cell (i, j) is array element [i + 4, j + 4]
N = 80


def field(key):
    a = np.asarray(v[key])[HALO + 1:HALO + 1 + N, HALO + 1:HALO + 1 + N]
    return [[float(a[i, j]) for j in range(N)] for i in range(N)]


def tmpl(op, value):
    return {
        "params": [],
        "match": {"op": op, "args": []},
        "body": {
            "op": "faq", "output_idx": ["i", "j"], "args": [],
            "ranges": {"i": {"from": "x"}, "j": {"from": "y"}},
            "expr": {"op": "index", "args": [{"op": "const", "args": [], "value": value}, "i", "j"]},
        },
    }


doc = {
    "esm": "1.1.0",
    "metadata": {"name": name, "description": desc, "license": "MIT",
                 "tags": ["test-inputs", "generated"]},
    "expression_templates": {
        "input_psi": tmpl("input_psi", field("lfn")),
        "input_ros": tmpl("input_ros", field("ros")),
        "input_wind_u": tmpl("input_wind_u", field("uf")),
        "input_wind_v": tmpl("input_wind_v", field("vf")),
    },
}
with open(out, 'w') as f:
    json.dump(doc, f)
    f.write('\n')
print('wrote', out)
