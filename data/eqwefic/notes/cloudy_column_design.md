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

---

## 6. Decision, and the residual closed (2026-09-21, later)

**The user chose Option A** — the fixed maximum grid with zero-thickness padding. Recorded
here so it is not re-litigated. The deciding points were that the padding is *exact* rather
than approximate, and that A is the only option that can reproduce WRF under cloud, which is
what every test in this project asserts.

### 6.1 The 2.5e-4 residual is explained, and it was not where I guessed

§4's caveat said my float64 reconstruction of the inserted grid sat at 2.5e-4 (and 7.1e-2 at
the model top) against the clear column's 8.9e-6, and that it had to be chased before A was
built. It has been. My `vcld` hypothesis was **wrong**, and so was reading the 7.1e-2 as a
boundary case. In order:

1. **The regridding is exact.** Every quantity `subgrid` produces — `z`, `zmid`, `vt`, `vair`,
   `vo3`, `vaer`, `vcld`, `cvo2` — matches the dump to **8e-7 or better**, worst case a
   relative 7.8e-7 in `vcld`. Not one entry exceeds 1e-6. The grid and its layer means are
   reconstructed correctly.
2. **The level mapping is exact.** All 60 model levels are present in the 77-level RT grid to
   better than 1e-5 km, so picking the nearest RT level is not an approximation.
3. **The problem is well conditioned.** Perturbing every layer optical input by one real32 ulp
   moves the fluxes by **2e-8**, while the difference from WRF is 1.4e-4 to 2.8e-4 — a ratio
   of 8 000 to 18 000. So it is not input sensitivity.
4. **It is WRF's own real32 arithmetic.** Re-running the identical algorithm in emulated
   float32 moves my own answer by **1.8e-4 (clear) and 5.8e-4 (cloudy)** — *larger* than my
   disagreement with WRF. Carried through to j-values: clear column, mine-vs-WRF **6–8.5e-6**
   against a real32 envelope of **1.0–1.6e-5**; cloudy column levels 1–40, mine-vs-WRF
   **2.54e-4** against an envelope of **5.7e-4**. In both cases the agreement is *better than
   the reference is determined*.
5. **The 4.6e-2 / 7.1e-2 at cloudy levels was my comparison harness, not the physics** — and
   chasing it turned up a real Fortran bug. See below. With it fixed, **every level of the
   cloudy column collapses to a uniform 2.55e-4**, flat from the surface to the model top, and
   inside the real32 envelope.

**So: the residual is an irreducible real32 artefact of the reference, established with
numbers.** The cloudy column is intrinsically noisier than the clear one (5.7e-4 vs 1.6e-5)
because the optically thick layers make the 138-unknown tridiagonal solve less well determined
in single precision — the *inputs* are well conditioned, the *arithmetic* is not. A cloudy test
should therefore assert at about **1e-3**, and a tighter tolerance would be asserting the
reference's rounding rather than the physics.

### 6.2 What the chase turned up: FORTRAN_BUGS B11

At `chem/module_phot_mad.F:2245`, `:2274` and `:2304`, the log-linear air interpolation onto
inserted levels computes

```fortran
hlocal = 1./alog(air(i-1)/air(i))
x0     = (z(lev)-zz(i-1))/(zz(i)+zz(i-1))     ! <-- SUM, should be a difference
zair(lev) = air(i-1)*exp(-x0/hlocal)
```

`hlocal` is a *normalized* scale height, so `x0` must be the normalized offset
`(z−zz(i−1))/(zz(i)−zz(i−1))`. With the sum it is not a fraction at all and is not scale-free.
Consequently `zair` **does not return `air(i)` at `z = zz(i)`**: on the reference column, at
inserted-grid levels that coincide *exactly* with model levels, the air density is wrong by
**4.9 % at model level 55 rising to 8.9 % at level 60**, growing with altitude exactly as a
sum-denominator would. It feeds every pressure-quenched quantum yield, so it biases
`j_ch2om`, `j_ch3cho`, `j_ch3coch3`, `j_ch3coc2h5`, `j_hcocho`, `j_ch3cocho` and `j_hcochob`
inside cloud by up to ~8 %. `j_hcocho` was the worst channel at every cloudy level, which is
the fingerprint.

Clear columns are untouched: the clear branch assigns `zair(lev) = air(i)` directly.

**It must be reproduced, not corrected**, for the cloudy grid to match WRF — and it is exactly
the kind of thing that would have been silently baked into the padded grid had the residual not
been chased first.

---

## 7. Option A built — and its exactness claim, corrected

`components/gaschem/madronich/cloud_grid.esm` (100 assertions) and
`couplings/madronich_column_cloudy.esm` (315 assertions over all four dumped columns,
nlevel 70/72/77/78) are green. **§4 and §6 above describe the padding as "exact, not
approximate" because τ = 0 collapses the continuity rows to an identity. That wording was
too strong and is superseded here.**

### 7.1 What went wrong, and what it showed

The first end-to-end cloudy run returned NaN for every j-value. The cause was
`SUBROUTINE sphers`, which divides each layer's slant path by its **thickness** to form
`dsdh`; on a zero-thickness padded layer that is 0/0. The §4 argument had covered the
two-stream solve and not the geometry upstream of it.

Chasing it turned up a second, quieter gap. The argument assumed a padded layer has τ = 0.
It does not: `optics` floors `dtscat` and `dtabs` at `1/largest = 1e-36` each, so a padded
layer carries **τ = 2e-36, ω = 0.5, g = 0**.

### 7.2 What is actually true, measured

* **`dsdh` is the only divisor that can reach zero.** Every division in the chain whose
  denominator involves a height or thickness was audited. The others are `zair_rt` (a *sum*
  of heights, B11) and the ozone/aerosol interpolation (guarded, with a numerator exactly 0
  on padded layers), and both act on the 70-level *input* grid, which padding never touches.
* **The guard is new logic with no Fortran counterpart**, because WRF never produces a
  zero-thickness layer. It sets `dsdh := 0` there.
* **Its value is immaterial.** `dsdh := 0` and `dsdh :=` its true finite secant limit agree
  to 3e-15 in the fluxes, because a padded layer's `dsdh` is multiplied by `taun = 2e-36`.
  So the guard is not choosing an arbitrary value that matters.
* **The padded solve equals the unpadded one to ≤ 1.8e-13** on fluxes of order one, at every
  real level of all four columns, in float64 (1410: 1.1e-13, 1200: 1.8e-15, 1: 1.8e-13,
  60: 1.2e-14). That residual is float64 rounding from a 262-row elimination instead of a
  150-row one — **nine orders of magnitude** below the precision to which WRF's own real32
  solve determines these fluxes (1.8e-4 clear, 5.8e-4 cloudy).

**Corrected claim:** the padded computation equals WRF's unpadded one *to float64 roundoff*.
It is not an approximation in any sense that reaches a test. It is **not** bit-identical,
though, and it depends on one guard that has no counterpart in the Fortran. The component
descriptions now say this.

**Moving the padding would not have avoided it.** `sphers` loops over every layer, so a
zero-thickness layer *anywhere* makes `dsdh` 0/0. The only way to keep padding out of the
geometry is to give padded layers real thickness, and that would add physics WRF does not
have. So the design is unchanged; only the claim about it is.

### 7.3 `interp.linear` and `trapez` — tested, not read

Probed directly (scratch document, not committed):

| question | result |
|---|---|
| state-valued axis | **works** — my earlier belief that it needed a constant axis was wrong |
| shaped query, bare `fn` | **vectorizes** |
| shaped query through `faq` | **vectorizes** |
| outside the axis | **clamps** to the end value; `trapez` instead writes a `1e-12` sentinel |
| repeated knot (`[1,2,2,4]`) | **NaN, silently** — the spec says to raise `interp_non_monotonic_axis` |

Consequences:

* **subgrid's `trapez` onto layer mid-points is now `interp.linear`** over the 70-level input
  interfaces, which are strictly increasing. This replaced a hand-rolled bracketing
  contraction I had written only because of the mistaken constant-axis belief. The clamp vs
  `1e-12` difference cannot trigger, since every mid-point lies inside the input grid.
* **`interp.linear` cannot be used on the padded RT axis**, which has repeated knots by
  construction; it would return NaN.
* **`trapez`'s other use, carrying j-values back onto model levels, needs no interpolation.**
  Every model interface is itself a node of the RT grid, so the coinciding node is read
  directly.
* **The silent NaN on a repeated knot looks like an EarthSciAST conformance gap.** A
  constant axis is checked at load; a state-valued one can only be checked at evaluation,
  and it isn't. Not yet filed.

### 7.4 Bound

`70 ≤ nlevel ≤ NRT = 131` for this 60-model-level configuration with one inserted level per
emission, against WRF's own abort ceiling of `nj = 200`. If a cloud field would need more
(`idt > 1` needs `cloud·dzt ≥ 50`; the reference column peaks near 11), the excess levels
are **not emitted**: `L_eff` clamps at NRT and the top of the column is silently truncated.
The component description says a consumer in that regime must raise NRT.
