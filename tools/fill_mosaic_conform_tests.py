"""Fill the inline tests of components/aerosol/mosaic/conform_electrolytes.esm from the
conform_electrolytes replays (data/eqwefic/chem170_scm/conform64):

  col_<S>.txt / col_<S>.out   every bin mosaic_dynamic_solver conforms at model step S of the
                              chem170 MOSAIC capture (make_conform_cases.py), run through WRF's own
                              compiled conform_electrolytes by kernels/mosaic_conform_driver
  hand.txt / hand.out         hand-specified bins for the branches the column never reaches
                              (make_conform_hand.py; names in hand_names.txt)

From each column step the cells are grouped by the branches they take and up to PER_CLASS cells
of each group (spread over the levels) become tests; every hand bin is a test.

    python3 tools/fill_mosaic_conform_tests.py
"""
import collections
import json

D = "/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/chem170_scm/conform64"
ESM = "components/aerosol/mosaic/conform_electrolytes.esm"
STEPS = (364, 1026)
PER_CLASS = 3
NLEV, NBIN = 59, 8
REL = 1e-9
XT_SPLIT = 1.999899983406067          # binary32 1.9999
SO4_TINY = 1.0000000036274937e-15     # binary32 1.e-15
BUILD = ("WRF's own compiled conform_electrolytes (main/libwrflib.a of the B30-fixed chem170 build, "
         "real(kind=8)) run by kernels/mosaic_conform_driver; rel 1e-9 with an absolute floor of 1e-12 "
         "of the bin's largest amount, the domain and regime flags pinned.")


def load(p):
    rows = open(p).read().split("\n")[1:]
    return [list(map(float, r.split())) for r in rows if r.strip()]


def classify(a, xt):
    """The branches conform_electrolytes takes on clipped amounts a (so4 no3 cl nh4 msa na ca), in
    the Fortran's own order of operations; returns (poor, regime, tags)."""
    so4, no3, cl, nh4, msa, na, ca = [max(0.0, x) for x in a]
    tags = []
    if any(x < 0 for x in a):
        tags.append("negative input clipped")
    c = min(ca, so4); ca = max(0.0, ca - c); so4 = max(0.0, so4 - c)
    if c > 0: tags.append("CaSO4")
    c = min(ca, .5 * msa); ca = max(0.0, ca - c); msa = max(0.0, msa - 2 * c)
    if c > 0: tags.append("Ca(MSA)2")
    if xt >= XT_SPLIT or xt < 0:
        tags.append("sulfate-poor" + (" (X_T = -1)" if xt < 0 else " (1.9999 <= X_T < 2)" if xt < 2 else ""))
        c = min(.5 * na, so4); na = max(0.0, na - 2 * c); so4 = max(0.0, so4 - c)
        c = min(na, msa); na = max(0.0, na - c); msa = max(0.0, msa - c)
        if c > 0: tags.append("NaMSA")
        c = min(ca, .5 * no3); ca = max(0.0, ca - c); no3 = max(0.0, no3 - 2 * c)
        if c > 0: tags.append("Ca(NO3)2")
        c = min(na, no3); na = max(0.0, na - c); no3 = max(0.0, no3 - c)
        cl = cl - na
        if cl < 0: tags.append("NaCl chloride deficit")
        cl = max(0.0, cl)
        c = min(ca, .5 * cl); ca = max(0.0, ca - c); cl = max(0.0, cl - 2 * c)
        if c > 0: tags.append("CaCl2")
        if ca > 0: tags.append("CaCO3")
        c = min(.5 * nh4, so4); nh4 = max(0.0, nh4 - 2 * c)
        if so4 - c > 0: tags.append("sulfate discarded")
        c = min(nh4, msa); nh4 = max(0.0, nh4 - c); msa = max(0.0, msa - c)
        if c > 0: tags.append("NH4MSA")
        if msa > 0: tags.append("free MSA")
        c = min(nh4, no3); nh4 = max(0.0, nh4 - c); no3 = max(0.0, no3 - c)
        c = min(nh4, cl); nh4 = max(0.0, nh4 - c); cl = max(0.0, cl - c)
        for x, n in ((no3, "HNO3"), (cl, "HCl"), (nh4, "NH3")):
            if x > 0: tags.append(n + " degassed")
        return 1, None, tags
    tags.append("sulfate-rich")
    c = min(na, msa); na = max(0.0, na - c); msa = max(0.0, msa - c)
    if c > 0: tags.append("NaMSA")
    c = min(nh4, msa); nh4 = max(0.0, nh4 - c); msa = max(0.0, msa - c)
    if c > 0: tags.append("NH4MSA")
    if msa > 0: tags.append("free MSA")
    if so4 == 0.0:
        regime = 0
        tags.append("regime 0 (no sulfate after CaSO4)")
    else:
        xp, xn = (na + nh4) / so4, .5 * na / so4 + 1.
        if xp >= xn:
            regime = 1
            c = min(.5 * na, so4); s1 = max(0.0, so4 - c)
            xh = nh4 / s1 if s1 > SO4_TINY else 0.0
            tags.append("regime 1, " + ("(NH4)2SO4 + letovicite" if xh >= 1.5 else "letovicite + NH4HSO4"))
            nh4 = 0.0
        elif xp >= 1.:
            regime = 2
            tags.append("regime 2")
            nh4 = max(0.0, nh4 - min(nh4, so4))
        else:
            regime = 3
            tags.append("regime 3")
            s1 = max(0.0, so4 - min(na, so4))
            nh4 = max(0.0, nh4 - min(nh4, s1))
    for x, n in ((no3, "HNO3"), (cl, "HCl"), (nh4, "NH3")):
        if x > 0: tags.append(n + " degassed")
    return 0, regime, tags


def make_test(tid, desc, inp, g_in, out):
    """inp: so4 no3 cl nh4 msa na ca (driver order); g_in: hno3 hcl nh3; out: a driver output row."""
    xt = out[2]
    aer = out[3:11]
    elec = out[11:33]
    gas = out[33:36]
    hcl_created = out[36]
    poor, regime, tags = classify(inp[:7], xt)
    a_in = list(inp[:7])
    g5 = [0.0, g_in[0], g_in[1], g_in[2], 0.0]
    floor = 1e-12 * max([abs(x) for x in a_in + g5] + [1e-300])
    tol = {"rel": REL, "abs": floor}
    pin = {"rel": 0.0, "abs": 0.0}
    c = {"lev": 1, "bin": 1}
    A = []

    def put(var, val, t, **kw):
        A.append({"variable": var, "coords": dict(c, **kw), "time": 0.0, "expected": val, "tolerance": t})

    put("x_t", xt, tol)
    put("poor", float(poor), pin)
    if regime is not None:
        put("regime", float(regime), pin)
    for e in range(22):
        put("elec", elec[e], tol, el=e + 1)
    for s in range(7):
        put("a", aer[s], tol, sp=s + 1)
    put("co3", aer[7], tol)
    for gi, v in zip(range(5), [0.0, gas[0], gas[1], gas[2], 0.0]):
        put("g", v, tol, gas=gi + 1)
    put("hcl_created", hcl_created, tol)
    so4_0 = max(0.0, a_in[0])
    put("so4_lost", so4_0 - aer[0] if poor else 0.0, tol)
    return {
        "id": tid,
        "description": f"Covers: {'; '.join(tags)}. {desc} {BUILD}",
        "parameter_overrides": {"a_in": [[a_in]], "g_in": [[g5]]},
        "time_span": {"start": 0, "end": 1},
        "assertions": A,
    }, tuple(tags)


def column_tests(step):
    cin, cout = load(f"{D}/col_{step}.txt"), [list(map(float, l.split())) for l in open(f"{D}/col_{step}.out")]
    assert len(cin) == len(cout) == NLEV * NBIN
    cands = collections.defaultdict(list)
    for n, (r, o) in enumerate(zip(cin, cout)):
        if o[1] == 0:
            continue
        k, b = n // NBIN + 1, n % NBIN + 1
        g_in = r[11:14] if b == 1 else cout[n - 1][33:36]
        desc = (f"Source: chem170 MOSAIC capture (dumps/chem170_full_mosaic, B30-fixed build), model step "
                f"{step}, level {k}, section {b}: aer(:, jtotal) as mosaic_dynamic_solver conforms it "
                f"(map_mosaic_species of rsub0), the gases as the bins before it left them.")
        t, tags = make_test(f"col{step}_k{k:02d}_b{b}", desc, r, g_in, o)
        cands[tags].append(t)
    out = []
    for tags, ts in sorted(cands.items(), key=lambda kv: -len(kv[1])):
        pick = sorted({round(i * (len(ts) - 1) / max(1, PER_CLASS - 1)) for i in range(PER_CLASS)})
        out += [ts[i] for i in pick]
    return out, cands


def hand_tests():
    hin, hout = load(f"{D}/hand.txt"), [list(map(float, l.split())) for l in open(f"{D}/hand.out")]
    names = open(f"{D}/hand_names.txt").read().split()
    out = []
    for n, r, o in zip(names, hin, hout):
        desc = (f"Source: hand-specified bin '{n}' (chem170_scm/conform64/make_conform_hand.py), "
                f"for a branch the chem170 column never reaches.")
        out.append(make_test(f"hand_{n}", desc, r, r[11:14], o)[0])
    return out


def main():
    doc = json.load(open(ESM))
    m = doc["models"]["MosaicConformElectrolytes"]
    tests = []
    for s in STEPS:
        ts, cands = column_tests(s)
        print(s, len(ts), "tests from", sum(len(v) for v in cands.values()), "cells in", len(cands), "branch classes")
        tests += ts
    tests += hand_tests()
    m["tests"] = tests
    with open(ESM, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    print(len(tests), "tests,", sum(len(t["assertions"]) for t in tests), "assertions")


if __name__ == "__main__":
    main()
