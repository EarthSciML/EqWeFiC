#!/usr/bin/env python3
"""Prove that MadronichCloudGrid FAILS LOUDLY outside the regime it implements.

components/gaschem/madronich/cloud_grid.esm reproduces WRF's subgrid only for
idt = 1 (one inserted level per cloudy half-layer) and only on a grid of NRT
levels. Outside either bound it must stop the run with a diagnostic, never
return a silently wrong grid. Two sentinel gathers enforce that (see the
component description): `idt_supported` and `rt_level_capacity`.

An esm inline test cannot express "this run must fail" -- the test schema has no
expected-error field and no xfail marker -- and a deliberately failing inline test
would sit permanently red in the suite. So the check lives here instead. It
builds TRANSIENT documents from the current cloud_grid.esm, outside the repo so no
directory test run can ever pick them up, runs the Rust CLI on each, and asserts:

  in_range     the real reference column passes           (no false alarm)
  idt_over     a very thick cloud layer needs idt > 1     -> must fail naming 'idt_supported'
  nrt_over     cloud on all 70 input levels needs 139     -> must fail naming 'rt_level_capacity'

Usage (from the repo root, ESD_ROOT exported as for any column test):
    tools/check_cloud_grid_bounds.py [--esm ./esm]
Exit status 0 only if all three behave as required.

This script writes no model logic: it copies the component, swaps in a cloud field
and replaces the tests block, which is the test-tuple use CLAUDE.md permits.
"""
import argparse, copy, json, os, subprocess, sys, tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMP = os.path.join(REPO, "components", "gaschem", "madronich", "cloud_grid.esm")
DUMP = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/dumps/photmad/esm_dump_photmad_1410.json"


def reference_column():
    """The clear 1999-10-22 reference column (call 1410): the geometry, T, air, O3
    and aerosol every case shares. Only the cloud field is varied."""
    v = json.load(open(DUMP))["vars"]
    col = {k: [float(x) for x in v[k]["data"]] for k in ("zz", "t_col", "air_col", "o3_col", "aer_col")}
    return col


def build(cloud, case_id, desc, assertions):
    doc = json.load(open(COMP))
    m = doc["models"]["MadronichCloudGrid"]
    # resolve the component's relative subsystem ref absolutely, so the transient
    # copy works from outside the repo
    for sub in m.get("subsystems", {}).values():
        if not os.path.isabs(sub["ref"]):
            sub["ref"] = os.path.normpath(os.path.join(os.path.dirname(COMP), sub["ref"]))
    c = reference_column()
    m["tests"] = [{
        "id": case_id, "description": desc,
        "parameter_overrides": {"ze": c["zz"], "t_node": c["t_col"], "air_node": c["air_col"],
                                "o3_node": c["o3_col"], "aer_node": c["aer_col"], "cloud_node": cloud},
        "time_span": {"start": 0, "end": 1}, "tolerance": {"rel": 1e-9},
        "assertions": assertions}]
    return doc


def run(esm, doc, workdir, name):
    path = os.path.join(workdir, name + ".esm")
    json.dump(doc, open(path, "w"))
    r = subprocess.run([esm, "test", path], capture_output=True, text=True, cwd=REPO)
    return r.stdout + r.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--esm", default=os.path.join(REPO, "esm"))
    a = ap.parse_args()
    if not os.environ.get("ESD_ROOT"):
        sys.exit("export ESD_ROOT first (see CLAUDE.md)")
    n = len(reference_column()["zz"])                  # 70 input interfaces
    model_top = 60                                      # WRF w-levels; the rest are above-model
    cases = [
        # (name, cloud field, what must happen, required substrings)
        ("in_range", [0.0] * n, "pass", []),
        ("idt_over",
         [0.0] * 54 + [800.0] + [0.0] * (n - 55),      # one layer, half-layer optical depth >> 50
         "fail", ["E_TREEWALK_CONSTARRAY_OOB", "const array 'idt_supported'", "out of range 1..1"]),
        ("nrt_over",
         [2.0] * n,                                     # cloud on every input level: 139 > 131
         "fail", ["E_TREEWALK_CONSTARRAY_OOB",
                  "const array 'rt_level_capacity' index 139 out of range 1..131"]),
    ]
    ok = True
    with tempfile.TemporaryDirectory(prefix="cloud_grid_bounds_") as wd:   # a few hundred kB
        for name, cloud, want, needles in cases:
            asserts = [{"variable": "nlevel", "time": 0.0, "expected": 70.0,
                        "tolerance": {"abs": 0.0}}] if want == "pass" else \
                      [{"variable": "nlevel", "time": 0.0, "expected": 0.0, "tolerance": {"abs": 0.0}}]
            out = run(a.esm, build(cloud, name, name, asserts), wd, name)
            total = [l for l in out.splitlines() if l.strip().startswith("TOTAL")]
            counts = total[-1].split()[-3:] if total else ["?", "?", "?"]
            if want == "pass":
                good = counts == ["1", "0", "0"]
                detail = "passes, no false alarm" if good else f"expected a pass, got {counts}"
            else:
                good = counts[2] != "0" and all(s in out for s in needles)
                hit = next((l.strip() for l in out.splitlines() if "E_TREEWALK_CONSTARRAY_OOB" in l), "")
                detail = (hit[hit.find("E_TREEWALK"):hit.find(" (CONFORMANCE")] if good
                          else f"expected a stopped run naming {needles}, got {counts}")
            ok &= good
            print(f"{'OK  ' if good else 'FAIL'} {name:9s} {detail}")
    print("ALL BOUNDS BEHAVE AS REQUIRED" if ok else "BOUND CHECK BROKEN")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
