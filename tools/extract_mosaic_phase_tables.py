#!/usr/bin/env python3
"""Transcribe MESA's mutual-deliquescence tables into mosaic_phase_parameters.esm.

The tables (jsalt_index, jsulf_poor, jsulf_rich, d_mdrh) are the model's own arrays
as load_mosaic_parameters fills them, dumped by kernels/mosaic_therm_driver (WRF fork
earthsciml-instrumented-chem170) at the end of its run.  d_mdrh is real(kind=8)
filled from default-real literals, so its values are binary32 roundings.  Numbers
only; the component's equations are hand-written.

    python3 tools/extract_mosaic_phase_tables.py <driver output .json>
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from esm_dump import load  # noqa: E402

P = os.path.join(os.path.dirname(__file__), "..", "components", "aerosol", "mosaic", "mosaic_phase_parameters.esm")


def main(path):
    v = load(path)["vars"]
    doc = json.load(open(P))
    var = doc["models"]["MosaicPhaseParameters"]["variables"]
    var["jsalt_index"]["default"] = [int(x) for x in v["jsalt_index"]]
    var["jsulf_poor"]["default"] = [int(x) for x in v["jsulf_poor"]]
    var["jsulf_rich"]["default"] = [int(x) for x in v["jsulf_rich"]]
    var["d_mdrh"]["default"] = [[float(x) for x in row] for row in v["d_mdrh"]]   # (63, 4)
    assert len(var["d_mdrh"]["default"]) == 63 and len(var["d_mdrh"]["default"][0]) == 4
    with open(P, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main(sys.argv[1])
