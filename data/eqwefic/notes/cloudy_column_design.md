# Cloudy columns in the Madronich photolysis chain: the design choice, measured

Written 2026-09-21. The Madronich components are exact for a **clear** column
(`couplings/madronich_column.esm`, 1323/1323). This note is the evidence needed to
decide how — or whether — to extend them to a cloudy one. **It does not make the
choice.**

---

## 1. What `subgrid` actually does, measured

`SUBROUTINE subgrid` (`chem/module_phot_mad.F:2213-2320`) inserts extra levels inside
cloud. Per model level `i` it emits:

* **clear** (`cloud(i) <= 0`): exactly one level, at `zz(i)`;
* **cloudy**: a cloud-base mid-point level if the level below was clear, then
  `idt_below` levels up to `zz(i)`, then `idt_above` levels into the half-layer above,
  with `idt = max(ifix(cloud(i)·dzt·0.02), 1)`.

**Validated exactly.** Re-implementing that emission rule and counting reproduces
WRF's `nlevel` in **all 13 dumped calls** (70, 72, 77, 78 ×10). With `idt = 1` it
reduces to

> `nlevel = 70 + (number of cloudy levels) + (number of contiguous cloud blocks)`

and reconstructing the *values* as well gives `max |z − z_WRF| = 4.8e-7` km and
`vt` to 5.4e-8 relative.

**`idt > 1` is never reached in this SCM but is reachable in general.** It needs
`cloud·dz >= 50`; the reference column's maximum is `15.56 /km × 0.7 km ≈ 11`. A deep
convective cloud (`qll ≈ 5 g/m³` → `bext ≈ 107 /km`, `cloud ≈ 534 /km`) over 0.2 km
gives `idt = 2`.

**Bounds.** The Fortran aborts at `lev > nj = 200`. For this 60-model-level
configuration, `idt = 1` and every model level cloudy gives `nlevel = 131`; alternating
cloud/clear gives 130. So **70 ≤ nlevel ≤ 131** here, against WRF's own ceiling of 200.

## 2. The insertion is pure numerical refinement — it conserves optical depth exactly

Measured on call 1 (6 cloudy levels, peak extinction 15.56 /km):

| | total cloud optical depth |
|---|---|
| WRF's 77-level inserted grid | **54.942384** |
| plain 70-level grid, layer = mean of bounding levels | **54.942391** |

Identical. What changes is only how it is **distributed**:

| | layers spanning the cloud | per-layer optical depth |
|---|---|---|
| inserted grid | **15** | up to **5.9** |
| plain 70-level grid | **7** | up to **11.6** |

So this is not a physics difference. It is sub-gridding an optically thick layer,
and delta-Eddington's two-stream solution is strongly nonlinear in per-layer optical
depth.

## 3. What skipping it costs, measured

Full RT run both ways on call 1, compared against WRF's `phot1_s`, with the coarse
grid's cloud optical depth **conserved** (the most favourable version):

| level | worst channel | plain 70-level grid vs WRF |
|---|---|---|
| 1 (surface) | `j_o31d` | **7.2e-3** |
| 20 | `j_o31d` | **7.5e-3** |
| 40 | `j_o31d` | **7.9e-3** |
| 55 (in cloud) | `j_hcocho` | **4.6e-2** |
| 59 (model top) | `j_hcocho` | **7.1e-2** |

**0.7 % below the cloud, 4.6–7.1 % in and above it.** For scale, the clear-column
chain agrees to **8.9e-6** and the coupled test asserts at `rel 1e-3`. So skipping the
insertion is 7× to 70× outside the tolerance the clear-sky work is held to, and the
error grows with cloud optical depth.

*(An earlier pass measured 19 %; that used a non-conserving cloud placement and
overstated the case. 0.7 %/7 % is the fair number.)*

## 4. The options, and what each actually costs

### A — Fixed maximum grid, padded with zero-thickness layers

Size the RT grid at `NLEV_RT = 131` (or 199, WRF's own ceiling) and pad unused levels
with **zero-thickness** layers.

**This is exact, not an approximation.** A zero-thickness layer in `ps2str` is provably
a no-op: `taun = 0` → `expon = 1` → `cup = cuptn` and `cdn = cdntn`, so the continuity
rows reduce to an identity and the tridiagonal solution is unchanged — the system
merely gains a redundant interface. The one apparent singularity, `mu2 = (tauc(i) −
tauc(i−1))/(tausla(i) − tausla(i−1)) = 0/0`, is **already guarded in the Fortran**
(`IF(tausla(i) .EQ. tausla(i-1)) mu2(i) = SQRT(largest)`), and `divisr`'s `AMAX1(eps, …)`
handles the downstream term.

*Needs:* per-level emission counts (elementwise, `floor` for `ifix`), a prefix sum to
starting indices (a contraction — already used for `tauc`/`tausla`), and an inverse map
RT-level → (model level, sub-index) via `interp.searchsorted`, which is in the closed
function registry. Every primitive exists.

*Costs:* one substantial new component, roughly the size of `two_stream.esm`; the RT
grows from 70 to 131 levels (~1.9×) and the two-stream recurrence from 138 to 260 rows;
and `trapez` back onto the model levels becomes real work (plain piecewise-linear,
`interp.linear` per model level — cheap, and it comes free with A).

*Caveat:* my float64 reconstruction of the inserted grid currently agrees with WRF to
**2.5e-4** at most levels, not the 8.9e-6 the clear column reaches, and to 7.1e-2 at
the model top. Something in the cloudy path is still slightly off — most likely the
`vcld` placement across the `lev`/`lev-1` assignments in `subgrid`'s three branches, or
the `vo3`/`vaer` re-interpolation onto the inserted grid. **That residual has to be
chased before A is worth building**, and it is cheap to chase in the float64 reference.

### B — Adaptive-grid discretization in EarthSciDiscretizations

ESD grids size from metaparameters, which are load-time constants. A grid whose *size*
varies with state is outside the ESD model entirely, so this means designing a new
capability in another repo around one scheme's idiosyncrasy, with no second consumer to
generalise against. **High cost, speculative, and the requirement is not general.**

### C — Run on the plain 70-level grid, conserving cloud optical depth

Expressible **today** with the components that already exist; no new machinery.
Cost is §3: 0.7 % below cloud, up to 7.1 % in and above it.

*The catch:* it cannot be made to agree with WRF, so every cloudy test would carry a
7 % tolerance, which certifies nothing. It would be a *correct-ish* column tendency,
not a verified reproduction.

### D — Fixed uniform sub-division (always split each model layer into k)

Fixed size, no state dependence, trivially expressible. But it does **not** reproduce
WRF — WRF refines only inside cloud and at different boundaries. As `k → ∞` both
converge to the same continuous solution, so D can be *more accurate* than WRF while
*matching* it less. That is a different goal from everything else in this repo.

## 5. The judgement, and why it goes to a human

The options do not differ mainly in cost — they differ in **what the project is for**,
and that is not mine to decide:

* **If the goal is to reproduce WRF**, which is what every test in this repo asserts,
  then **A is the only option that can**, C and D are ruled out by construction, and B
  is the wrong place to put the complexity. A is feasible and the padding argument
  makes it exact.
* **If the goal is a correct column tendency** and a few percent in `j` under thick
  cloud is acceptable, **C is free today** and its error is now measured rather than
  guessed.

My own lean is **A, conditional on first closing the 2.5e-4 residual in the float64
reconstruction** — because a cloudy test held to 7 % would not be a test. But that is a
resourcing call about a substantial component, so it belongs to whoever is paying for
it.

**One thing is settled either way:** the clear-column chain is unaffected, and
`nlevel = 70` is not a simplification there — it is what WRF does when there is no cloud.
