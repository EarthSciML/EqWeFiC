# EqWeather-SCM tendency error budget

**First measured 2026-09-21** on `couplings/eqweather_scm.esm`; **re-measured the same
day** after the changes this file asked for were made. EarthSciAST `main` at
`be1fafe6d` (`./esm` built 2026-09-21; carries #439 and #440), `ESD_ROOT` = the pinned
`EarthSciDiscretizations-column` worktree at `f80069a`.

## Why this document exists

A green test only means *each residual is under the tolerance that assertion was given*.
The first pass of this file measured every residual with the tolerances set to zero, and
found that the document's bounds recorded how well the column happened to agree rather
than how well it was required to agree: **15 of 134 assertions carried a tolerance larger
than the entire reference amplitude**, so they would have passed against a tendency of
zero or one of the wrong sign. This file states the actual agreement, per assertion, so
that "correct tendencies" is a measured claim rather than a pass/fail.

## What changed, and what it cost

| | before | after |
|---|---|---|
| tests on `couplings/eqweather_scm.esm` | 2 (both sunlit, both convective) | 4 (+ step 900 deep night, + step 300 dusk) |
| assertions | 134 | 278 |
| assertions with tolerance >= reference amplitude | 15 | 13 (4 x `dw_dt`, 9 x `mp_*`; all named below) |
| PBL reference | WRF's driver-level `rthblten`/... at dt = 60 s | the dt -> 0 Richardson limit of the YSU kernel |
| worst PBL relative residual | 3.1e-01 | 6.3e-03 |
| tolerances | set by eye, 4-10x a residual measured once | 4x the residual measured for that assertion |

Five more coupled documents were moved onto the same PBL reference:
`scm_physics_column{,_night}.esm`, `geometry_surface_pbl_column{,_night}.esm`,
`physics_column{,_night}.esm` and `surface_pbl_column.esm`.

## Method

The document is copied with **every tolerance set to `{abs: 0, rel: 0}`** and re-run,
which forces each assertion to report `actual` against `expected`. For a `Linf_error`
assertion the reported `actual` IS the absolute error, since `expected` is 0. The
reference amplitude is `max |reference|` taken from the assertion's own `faq` reference
array (or `|expected|` for a scalar assertion), and the relative error is the ratio of
the two. The copies are scratch files (`couplings/_z_*.esm`, gitignored); they have been
deleted. Tolerances are then set to the smallest value of {1, 2, 3, 5} x 10^k that is at
least **4x** the measured residual, in the kind (`rel` for scalars, `abs` for column
`Linf_error`) the assertion already used, and an assertion that is bit-exact keeps
`abs 0`.

## The PBL reference was wrong in kind, and that was the whole residual

WRF's driver-level `rthblten`/`rublten`/`rvblten`/`rqvblten`/`rqcblten`/`rqiblten` are the
IMPLICIT tridiagonal step's `(c_new - c_old)/dt` at dt = 60 s. That is an increment, not a
derivative, and the documents in `couplings/` compute a derivative. Eight of them are now
referenced to the **dt -> 0 limit** — the second-order Richardson extrapolation
`(8 T(dt/4) - 6 T(dt/2) + T(dt))/3` of real64 `kernels/ysu_driver` replays of the SAME
dumped column at dt = 0.04/0.02/0.01 s, which is the reference `ysu.esm`'s own tests
already used and pass at `abs 1e-9`.

**A trap, and it invalidated the plan this work started from.**
`data/eqwefic/dumps/ysu/driver_120_dt*.json` — called "the replays that already exist" for
the coupled step-60 test — are replays of `ysu_120.flat`, a DIFFERENT run with z0 = 0.15 m,
hfx 89.8 W/m², u* 0.587. The reference run `data/eqwefic/scm_ref` has z0 = 0.05 m, hfx
70.7, u* 0.506, and its column is `scm_ref_120.flat`. Using the wrong replays would have
moved every step-60 PBL reference by about 25 %. New replays were generated from
`scm_ref_{120,600,1800,2820}.flat` at dt = 0.08/0.04/0.02/0.01/0.005 s.

**Extrapolation self-consistency** — the `(0.04, 0.02, 0.01)` limit against the
`(0.08, 0.04, 0.02)` one, as a fraction of the column maximum:

| replay set | step / call | `utnp` | `ttnp` |
|---|---|---|---|
| `driver_120` (the stated check, reproduced) | — / 120 of the old run | 7.5e-10 | 4.9e-08 |
| `scmref_120` | 60 / 120 | 8.6e-10 | 4.3e-08 |
| `scmref_600` | 300 / 600 | 7.3e-10 | 3.6e-08 |
| `scmref_1800` | 900 / 1800 | 3.6e-10 | 6.6e-08 |
| `scmref_2820` | 1410 / 2820 | 1.0e-09 | 5.6e-08 |

**Measured residual, before and after, per document** (`Linf` over the column):

| document | step | `pbl_dthdt` | `pbl_dudt` | `pbl_dvdt` | `pbl_dqvdt` |
|---|---|---|---|---|---|
| `physics_column` | 60 | 2.7e-05 -> **2.86e-09** | 5.7e-05 -> **1.09e-09** | 1.6e-05 -> **2.84e-09** | 2.4e-09 -> **3.18e-12** |
| `physics_column_night` | 300 | 1.9e-06 -> **7.06e-09** | 1.5e-05 -> **6.21e-09** | 3.5e-06 -> **8.70e-09** | 5.2e-10 -> **4.93e-12** |
| `surface_pbl_column` | 60 | 2.74e-05 -> **1.33e-06** | 5.74e-05 -> **2.09e-07** | 1.65e-05 -> **3.08e-07** | 2.36e-09 -> **8.65e-11** |
| `scm_physics_column` | 1410 | 1.13e-04 -> **1.40e-06** | 6.15e-05 -> **3.79e-07** | 2.96e-05 -> **4.68e-07** | 1.41e-08 -> **2.14e-10** |
| `scm_physics_column_night` | 900 | 1.70e-06 -> **8.88e-08** | 2.51e-05 -> **6.56e-09** | 4.46e-06 -> **1.17e-08** | 2.40e-10 -> **1.05e-11** |
| `geometry_surface_pbl_column` | 1410 | 1.14e-04 -> **1.36e-06** | 6.15e-05 -> **5.27e-07** | 2.96e-05 -> **6.49e-07** | 1.41e-08 -> **2.38e-10** |
| `geometry_surface_pbl_column_night` | 900 | 1.70e-06 -> **7.31e-08** | 2.51e-05 -> **2.32e-08** | 4.46e-06 -> **3.14e-08** | 2.40e-10 -> **1.29e-11** |

Read the first row, not the last: **`physics_column`, which feeds YSU the DUMPED
geometry, lands at 1e-9 K/s and m/s², about 1e-5 of the column maximum.** That is the
component's own standalone accuracy, reproduced inside a coupled document. The documents
that build the geometry themselves sit two to three orders higher, at ~1e-6 — the same
N50 lowest-layer pressure offset (0.30 Pa) that moves the surface layer by 1e-6 to 7e-5
relative. Neither is a defect in YSU.

The size of the correction is worth stating separately, because it is what the old
tolerances were hiding: at step 1410 the dt = 60 s driver value differs from the dt -> 0
limit by 1.13e-4 K/s against a column maximum of 3.71e-4 (30 %) in theta and 6.12e-5
against 3.37e-4 (18 %) in u. At step 300 the same difference is 1.2 % and 4.0 %, because
the implicit step's error scales with the diffusivity and K_h there is 54x smaller.

## WSM6: the operator-split difference, measured

The step-60 cirrus test's three microphysical assertions were the other block of
"disagreeing" tendencies. They have now been measured rather than assumed, on the same
column (`data/eqwefic/dumps/wsm6/scm_ref_120.*`, WRF step 60 / call 120, dtcld = 60 s
because ndt = max(nint(dt/dtcldcr), 1) = 1).

**1. The whole-scheme increment does NOT converge as dt -> 0, and the reason is `pigen`,
not the saturation adjustment.** Replaying the kernel at dtcld = 0.01, 0.02 and 0.04 s
and splitting `(x_out - x_in)/dtcld` into the scheme's three stages:

| stage (q_i) | dtcld = 0.01 | 0.02 | 0.04 | 60 |
|---|---|---|---|---|
| sedimentation + melting/freezing | 3.4935e-08 | 3.4935e-08 | 3.4935e-08 | 3.6372e-08 |
| the process-rate block | 7.7819e-06 | 3.8910e-06 | 1.9455e-06 | 3.2997e-08 |
| saturation adjustment | 0 | 0 | 0 | 0 |

The fallout stage **converges** (it is dt-dependent by only 4 % between dtcld = 0.04 and
60 s). The rate block scales exactly as 1/dtcld: its *increment* is 7.78e-08 kg/kg at
every small dtcld. The term responsible is ice nucleation,
`pigen = max(0, (roqi0/den - qi)/dtcld)` (`mp_wsm6.F90:1200`) — a projection of q_i onto
the diagnosed equilibrium ice content in exactly one sub-step, so its increment is
dtcld-independent and its "rate" is whatever 1/dtcld says. **The saturation adjustment is
identically zero in this column**, so the piece that makes a dt -> 0 extrapolation
meaningless here is `pigen`, not `pcond`. This is the behaviour FORTRAN_BUGS N3 records;
what is new is which term dominates and by how much. At the model's own dtcld = 60 s
`pigen` is only 1.23e-09 of the 3.30e-08 total, i.e. 3.7 %, so at finite dtcld the scheme
is dominated by a genuine rate (`pidep`) and the ESM's formulation — which prescribes the
SAME dtcld = 60 s — is comparable with it.

**2. What is left is the sequencing, and it is measurable.** WRF applies fallout first
and evaluates the rate block at the state fallout left (x1 = x0 + dx_sed); a continuous
formulation evaluates every stage at one instantaneous state x0 and sums. To measure
`R(x1) - R(x0)` at the model's own sub-step without changing dtcld, the same kernel was
replayed at dtcld = 60 s from a synthetic input equal to x1 itself
(`data/eqwefic/dumps/wsm6/syn_postsed_120.flat`, the dt = 60 s replay's `*_rates` arrays
written back over the state), so its rate block sees x2 = x1 + dx_sed(x1):

- `pidep` at x1 (WRF's own run): **3.29972e-08** kg/kg/s
- `pidep` at x2 (the synthetic replay): **4.45250e-08**
- difference: **1.15278e-08**, i.e. 35 % of the rate at x1

The sensitivity is linear over this interval, so the backward estimate
`2 R(x1) - R(x2) = 2.1469e-08` should be `R(x0)`; the directly measured `R(x0)` — the
dtcld = 0.01 s replay, where the fallout increment is 6000x smaller — is **2.15477e-08**.
They agree to **0.4 %**, so the operator-split difference at dtcld = 60 s is
**R(x1) - R(x0) = 1.145e-08 kg/kg/s**.

**3. That is exactly the document's residual.** `mp_dqv_dt`'s measured residual is
**1.145e-08 kg/kg/s** — the same number to four significant figures. Carried into
temperature through `xl/cpm/pi` at the layer where it is largest (level 51,
xl = 2.588e6, cpm = 1004.6, pi = 0.7552) it is 3.93e-05 K/s against `mp_dtheta_dt`'s
measured 4.318e-05 K/s, 91 % of it. **So the vapour and heating residuals of the
microphysics are the operator split, in full, and are not an error in the scheme.**

**4. `mp_dqi_dt` is the one that is not explained this way**, because the ice tendency
also carries the fallout stage. See the next section: that part is a real discretisation
difference, measured like-for-like.

## WSM6 ice sedimentation: the like-for-like comparison (2026-09-22)

**Why it had to be redone.** The earlier numbers (5.304e-08 against 3.4935e-08) came from
two different computations, and the component's own 65/65 was being taken as evidence
that sedimentation was right. It is not evidence for ICE: `sedimentation.esm`'s own
cirrus test says *"Cloud-ice fluxes have no kernel counterpart above the surface (fallc
is dumped only at k = 1), so F_i and dqi_dt are checked against the same formulae
evaluated by the filler on the input state -- a wiring check, not an independent
reference."* The rain/snow/graupel fluxes ARE checked against the kernel's `fall`
arrays, and that comparison is sound, but in this column those species are identically
zero. So before today no one had compared the model's ice fallout with WRF's.

**The comparison.** One scratch copy of `couplings/eqweather_scm.esm` carrying only the
step-60 cirrus test, with `mp_sed_dqi_dt` -- the document's own ice sedimentation
tendency, evaluated by the coupled model at that state -- asserted at zero tolerance
against WRF's ice fallout tendency on the SAME column, the SAME 59 levels, the SAME
quantity (kg kg⁻¹ s⁻¹): `(qi_rates - qi_in)/dtcld` from the real64 `kernels/wsm6_driver`
replay of `scm_ref_120` at dtcld = 0.01 s. `qi_rates` is the state after the fallout and
melting/freezing block and before the process rates; nothing melts or freezes in this
column, so it is the fallout alone. That reference is converged: levels 50-57 are
identical to five significant figures at dtcld = 0.01, 0.02 and 0.04 s.

| level | model (donor-cell) | WRF, dtcld -> 0 | ratio |
|---|---|---|---|
| 50 | 3.567e-09 | 1.721e-09 | 2.07 |
| 51 | 5.304e-08 | 3.493e-08 | 1.52 |
| 52 | 7.335e-09 | 2.775e-08 | 0.26 |
| 53 | -2.918e-08 | -1.751e-08 | 1.67 |
| 54 | -1.835e-08 | -2.352e-08 | 0.78 |
| 55 | -8.916e-09 | -1.292e-08 | 0.69 |
| 56 | -2.606e-09 | -4.606e-09 | 0.57 |
| 57 | -6.679e-10 | -1.551e-09 | 0.43 |
| 58 | 0 | -1.884e-10 | (see B16) |

**Linf difference 2.041e-08 against a WRF peak of 3.493e-08: 58 %.** The ratio runs from
0.26 to 2.07, so it is not a constant factor and not a wiring error, and the disagreement
survives. The model's column is what it says it is: recomputing the donor-cell divergence
in float64 from WRF's own dumped `den`, `qi`, `delz` and WRF's own fall-speed formula
(`xni` agreeing with the dump to 1e-7) gives 5.3045e-08 at level 51, against the model's
5.3037e-08.

**What WRF's fallout is, measured.** `nislfv_rain_plm` is a forward semi-Lagrangian remap,
and its dt -> 0 limit is a second-order flux divergence with two ingredients donor-cell
lacks. Transcribing both in float64 on the same state reproduces WRF at every level from
50 to 57 **to five significant figures** (Linf 1.88e-10, 0.5 % of the peak, and all of it
at level 58 -- FORTRAN_BUGS B16):

- a THIRD-ORDER interface fall speed, `wi(k) = 9/16(ww(k)+ww(k-1)) - 1/16(ww(k+1)+ww(k-2))`,
  on its own closing 37 % of the gap (Linf difference 58.6 % -> 37.2 %);
- a monotone piecewise-linear donor-face value whose slope is the AVERAGE of the two
  one-sided differences, zeroed at extrema and reset if either face value goes negative,
  which closes the rest. A textbook minmod PLM does NOT do it: 31.6 % remains.

The exact formulae are FORTRAN_BUGS **N77**. B16 is what this found on the way: the
routine deletes the whole content of the topmost ice layers every sub-step,
4.7542e-10 kg m⁻² at dtcld = 0.01, 0.02 and 0.04 alike, and none of it reaches the
surface (`fallc(1) = 0`).

**AUTHORED 2026-09-22: EarthSciDiscretizations PR #42** (`sedimentation_jh2010_flux_D_lev`,
stacked on #36). The rule was checked against WRF's own `nislfv_rain_plm`, extracted
verbatim and Richardson-extrapolated to dt -> 0, on five synthetic columns. The relative
Linf differences were 2.4e-9 to 1.1e-7, against 32-278 % for donor-cell. On the WRF
cirrus column the like-for-like numbers are:

| comparison | donor-cell | this rule |
|---|---|---|
| component, dumped state, 40-layer window (Linf / peak) | 2.05e-8 (58.6 %) | **2.07e-13 (5.9e-6)** |
| coupled `mp_sed_dqi_dt`, below the top two layers | 2.04e-8 (58 %) | **7.7e-11 (0.22 %)** |
| coupled `mp_dqi_dt`, whole microphysics against WRF's dt = 60 s | 2.81e-8 | **1.01e-8** |

- **Only the ice tendency distinguishes the two rules.** Swapping the rule changes exactly
  one assertion; the flux assertions do not depend on the rule.
- **The top two layers are excluded, and why.** `nislfv_rain_plm` empties the top two
  layers of the column on every call (B16, re-characterised: the top two layers of the
  column, not of the hydrometeor layer). The rule keeps them conservative.
- **The 0.22 % left in the coupled column** is the coupled state's own geometry offset.
- **What is left of `mp_dqi_dt` is the operator split.** 1.01e-8 is the same pidep
  sequencing difference that makes up `mp_dqv_dt`'s 1.145e-8.
- **Not yet tightened.** The committed document still mounts the donor-cell rule, because
  no pinned ESD_ROOT carries #42 yet. Once it does, 1.01e-8 against a 2.17e-8 amplitude
  allows at best `abs 2e-8` (0.92 of the amplitude). The 4x convention (4e-8) would still
  be vacuous. Do not tighten before the rule is mounted.

**So the ESD rule IS warranted, and this is what it would have to be** (as authored in #42):
a new rule on the same face-flux structure as `sedimentation_upwind1_flux_D_lev`, i.e.
lowering `-D(W m, lev)` for a layer-centred fall velocity W and mass density m, with the
face flux `F(k) = w_face(k) m_face(k)` instead of `W(k) m(k)`:

1. `w_face` by N77's third-order interpolation, with its end-face cases
   (`wi(1) = ww(1)`, `wi(2) = (ww(2)+ww(1))/2`, `wi(N) = (ww(N)+ww(N-1))/2`,
   `wi(N+1) = ww(N)`) and the "top of group" override `if ww(k) == 0 then wi(k) = ww(k-1)`;
2. `m_face` = the donor cell's bottom-face value of N77's reconstruction, with slopes
   formed over the NON-UNIFORM layer thicknesses (`dz(k-1)+dz(k)` and `dz(k)+dz(k+1)`
   denominators), which the column grid already supplies;
3. `F(N+1) = 0` at the lid, divergence over `dz(k)` as now.

It is a five-point stencil in velocity and a three-point one in mass, and it is
non-smooth (the extremum switch and the positivity reset). The finite-dt `decfl > 0.05`
velocity clamp is NOT part of the dt -> 0 operator and should be left out. A generic PLM
or PPM limiter from the literature would not reproduce WRF, as the minmod trial shows.

**Consequence for the three `mp_*` assertions.** They keep WRF's dt = 60 s increment as
their reference, because no dt -> 0 reference exists for a scheme containing `pigen`, and
because at dtcld = 60 s the two sides ARE comparable — the residual is now a measured,
attributed quantity rather than an unexplained one. Their tolerances are set from the
measured residual like every other assertion. `mp_dqi_dt` is the single in-scope
assertion whose residual (2.812e-08) exceeds its own reference amplitude (2.173e-08), so
no tolerance can be both sufficient and non-vacuous; it is left at `abs 1e-7` and named
here rather than quietly rounded.

## The optional schemes really are off in the reference run

Checked in `data/eqwefic/scm_ref/namelist.input` and `namelist.output`, because a term
missing from the document is only correct if it is missing from the reference:
`cu_physics = 0` and `SHCU_PHYSICS = 0` (all 11 domains) — no cumulus and no shallow
cumulus; `GWD_OPT = 0` (all 11 domains) — no gravity-wave drag. Also confirmed for the
record: `ICLOUD = 1` (which is the `cal_cldfra1` path the document mounts),
`W_DAMPING = 0` (so `drw_w_damp` is correctly not mounted), `DAMP_OPT = 2` (the Rayleigh
sponge, which is), `SCM_FORCE = 0`, `DIFF_OPT = 2`, `KM_OPT = 2` (inert behind
`bl_pbl_physics /= 0`), `ISFFLX = 1`, `IFSNOW = 1`. So the mounted set is complete for
this configuration, and INVENTORY's "optional" marking on GWDO and cumulus is correct.

## What is NOT wrong

**Vertical advection is not missing.** The five mounted WRF-ARW components carry no
advection term, and `du_dt_dyn` is only Coriolis + curvature + Rayleigh. That is correct
for this configuration, confirmed three ways from `data/eqwefic/dumps/dyn_column`:
`dyn_ww` (the eta-coordinate vertical mass flux) is **identically 0.0** at all 60
interfaces at steps 60, 900 and 1410; `dyn_mu` is **-97.750000 at all three steps**, so
d(mu)/dt = 0 over 1350 steps; and the instrumented dump contains `dru_coriolis`,
`dru_curvature`, `dru_rayleigh`, `drw_pg_buoy`, `drw_rayleigh`, `drw_w_damp` and **no
advection term at all**. In eta coordinates the surfaces are mass-following and omega is
diagnosed from the column-integrated horizontal mass divergence, which is identically zero
in a horizontally homogeneous single column. `drw_w_damp` is also identically zero
(`W_DAMPING = 0`) and correctly not mounted.

**The optional schemes are off in the reference run**, so their absence is correct and
not an omission: `cu_physics = 0`, `SHCU_PHYSICS = 0`, `GWD_OPT = 0` (see the namelist
section above).

**WSM6 sedimentation is a discretisation difference, not an error in the transcription** --
but it IS a difference from WRF, 58 % Linf like-for-like, and the rule that closes it is
specified above.

## So the actual state of the column

Every tendency in `couplings/eqweather_scm.esm` that is asserted against a genuine
instantaneous derivative now agrees with WRF to **1e-9 to 1e-6 absolute** — 1e-7 to 1e-3
of its own column maximum — across four regimes spanning day and night, convective and
stable. The exceptions are named, measured and attributed:

1. **`mp_dqi_dt`** (step 60 only): residual 2.812e-08 against a reference amplitude of
   2.173e-08, dominated by donor-cell-versus-PLM sedimentation. **This is the one in-scope
   assertion whose tolerance cannot be brought below its reference amplitude.**
2. **`dw_dt`**: WRF's real32 p' rounding drives a w tendency up to half the physical one
   (N67), so w cannot be matched beyond ~3e-4 m/s².
3. The surface layer's 1e-6 to 7e-5 relative floor, which is the 0.30 Pa
   driver-versus-kernel pressure offset (N50), not a scheme error.

## Next actions implied

1. **The sedimentation flux rule specified above** (N77's interface velocity and
   face reconstruction) in EarthSciDiscretizations. It is the only thing that makes
   `mp_dqi_dt` non-vacuous, and the only known place where the SCM's tendencies differ
   from WRF's dt -> 0 operator for a reason on the model side.
2. `dw_dt`'s tolerance is set by WRF's real32 STATE, not by the model or by arithmetic,
   and a real64 kernel replay does not escape it: `pg_buoy_w` takes p' as an input, and
   replaying it in real64 on the dumped p' reproduces WRF (that is what
   `vertical_momentum.esm` already asserts, at abs 5e-8 at steps 60, 900 and 1410);
   replaying `calc_p_rho_phi` in real64 on the dumped state reproduces the MODEL's own p',
   so comparing against it is circular. N70 shows the reference w tendency at step 1410
   is almost entirely the real32 column's departure from hydrostatic balance (4.23e-4 ->
   3.3e-11 on the balanced state). The measured chain closes: the coupled `p_p` residual
   is 1.571e-2 Pa, inside N67's 1.2-1.6e-2 Pa band, and `dw_dt`'s is 3.04e-4, N67's
   "up to 3.0e-4". Only a double-precision WRF MODEL run would help, and that changes the
   reference for every assertion in the repo. The physics is already pinned three ways:
   `dw_dt_nonpg` (abs 2-3e-10), `vertical_momentum.esm` (abs 5e-8), and
   `eqweather_scm_dynamics.esm`'s hydrostatic-balance test (abs 1e-10).
3. The four regimes are all from the same 59-hour em_scm_xy run and the same sounding. A
   genuinely different column (the em_quarter_ss supercell, already used by
   `microphysics_column*.esm` and `cal_cldfra`) would test what a diurnal cycle cannot.

## Per-assertion detail, `couplings/eqweather_scm.esm`

## Result by subsystem, all four regimes

| subsystem | n | median rel err | worst rel err |
|---|---|---|---|
| WRF-ARW dynamics | 33 | 2.4e-07 | 7.2e-01 |
| surface layer (sfclayrev) | 92 | 2.7e-06 | 2.5e-04 |
| YSU PBL | 36 | 9.1e-05 | 6.3e-03 |
| RRTM longwave | 24 | 6.4e-07 | 1.5e-04 |
| Dudhia shortwave | 4 | 1.8e-06 | 5.6e-06 |
| cal_cldfra | 4 | 2.8e-05 | 7.7e-04 |
| slab LSM | 32 | 2.7e-04 | 1.1e-01 |
| WSM6 microphysics | 9 | 1.0e+00 | 9.4e+00 |
| coupled sums | 21 | 2.3e-03 | 1.2e-01 |

**Assertions whose tolerance is >= their reference amplitude: 13**

| test | variable | tolerance | ref amplitude | measured resid | tol/amp |
|---|---|---|---|---|---|
| eqweather_scm_step300_dusk | `mp_dqi_dt` | {"abs": 5e-09} | 8.204e-11 | 7.747e-10 | 60.94 |
| eqweather_scm_step900_night | `mp_dqi_dt` | {"abs": 2e-11} | 1.235e-12 | 3.335e-12 | 16.19 |
| eqweather_scm_step900_night | `mp_dtheta_dt` | {"abs": 3e-06} | 5.086e-07 | 5.086e-07 | 5.90 |
| eqweather_scm_step60_cirrus | `mp_dqi_dt` | {"abs": 1e-07} | 2.173e-08 | 2.812e-08 | 4.60 |
| eqweather_scm_step300_dusk | `mp_dqv_dt` | {"abs": 5e-09} | 1.144e-09 | 1.145e-09 | 4.37 |
| eqweather_scm_step300_dusk | `mp_dtheta_dt` | {"abs": 2e-05} | 4.578e-06 | 4.579e-06 | 4.37 |
| eqweather_scm_step60_cirrus | `dw_dt` | {"abs": 0.001} | 3.713e-04 | 2.168e-04 | 2.69 |
| eqweather_scm_step1410 | `dw_dt` | {"abs": 0.001} | 4.233e-04 | 3.042e-04 | 2.36 |
| eqweather_scm_step300_dusk | `dw_dt` | {"abs": 0.001} | 4.318e-04 | 2.456e-04 | 2.32 |
| eqweather_scm_step900_night | `dw_dt` | {"abs": 0.001} | 4.365e-04 | 2.244e-04 | 2.29 |
| eqweather_scm_step900_night | `mp_dqv_dt` | {"abs": 1e-11} | 5.821e-12 | 1.338e-12 | 1.72 |
| eqweather_scm_step60_cirrus | `mp_dtheta_dt` | {"abs": 0.0002} | 1.241e-04 | 4.318e-05 | 1.61 |
| eqweather_scm_step60_cirrus | `mp_dqv_dt` | {"abs": 5e-08} | 3.300e-08 | 1.145e-08 | 1.52 |

### `eqweather_scm_step1410`

| variable | reduce | abs resid | ref amplitude | rel err | tolerance | tol/amp |
|---|---|---|---|---|---|---|
| `dw_dt` | Linf_error | 3.042e-04 | 4.233e-04 | 7.19e-01 | {"abs": 0.001} | 2.36 |
| `dTsk_dt_wrf_step` | value | 5.788e-06 | 1.363e-04 | 4.25e-02 | {"rel": 0.15} | 0.15 |
| `lsm_dTs_dt` | Linf_error | 5.788e-06 | 2.660e-04 | 2.18e-02 | {"abs": 1e-05} | 0.04 |
| `dtheta_m_dt_phys` | Linf_error | 2.329e-06 | 3.301e-04 | 7.06e-03 | {"abs": 1e-05} | 0.03 |
| `dtheta_m_dt` | Linf_error | 2.329e-06 | 3.301e-04 | 7.06e-03 | {"abs": 1e-05} | 0.03 |
| `pbl_dthdt` | Linf_error | 2.343e-06 | 3.712e-04 | 6.31e-03 | {"abs": 1e-05} | 0.03 |
| `rth_phys_sum` | Linf_error | 2.339e-06 | 3.797e-04 | 6.16e-03 | {"abs": 1e-05} | 0.03 |
| `dv_dt` | Linf_error | 6.516e-07 | 2.390e-04 | 2.73e-03 | {"abs": 3e-06} | 0.01 |
| `pbl_dvdt` | Linf_error | 6.516e-07 | 3.248e-04 | 2.01e-03 | {"abs": 3e-06} | 0.01 |
| `du_dt` | Linf_error | 5.293e-07 | 2.655e-04 | 1.99e-03 | {"abs": 3e-06} | 0.01 |
| `pbl_dudt` | Linf_error | 5.293e-07 | 3.370e-04 | 1.57e-03 | {"abs": 3e-06} | 0.01 |
| `lsm_hfx` | value | 1.435e-01 | 1.736e+02 | 8.26e-04 | {"rel": 0.001} | 0.00 |
| `dtheta_m_conv_delta` | Linf_error | 1.169e-07 | 1.538e-04 | 7.60e-04 | {"abs": 5e-07} | 0.00 |
| `pbl_dqvdt` | Linf_error | 2.358e-10 | 3.306e-07 | 7.13e-04 | {"abs": 1e-09} | 0.00 |
| `p_p` | Linf_error | 1.571e-02 | 3.798e+01 | 4.14e-04 | {"abs": 0.05} | 0.00 |
| `lsm_qfx` | value | 2.565e-08 | 8.103e-05 | 3.17e-04 | {"rel": 0.0005} | 0.00 |
| `lsm_lh` | value | 6.412e-02 | 2.026e+02 | 3.17e-04 | {"rel": 0.0005} | 0.00 |
| `pbl_K_h` | Linf_error | 4.833e-02 | 1.768e+02 | 2.73e-04 | {"abs": 0.06} | 0.00 |
| `pbl_K_m` | Linf_error | 2.736e-02 | 1.026e+02 | 2.67e-04 | {"abs": 0.04} | 0.00 |
| `lsm_qsfc` | value | 2.059e-06 | 1.139e-02 | 1.81e-04 | {"rel": 0.0005} | 0.00 |
| `sfc_hfx` | value | 1.272e-02 | 1.735e+02 | 7.33e-05 | {"rel": 0.0003} | 0.00 |
| `sfc_mol` | value | 2.271e-05 | 3.192e-01 | 7.12e-05 | {"rel": 0.0003} | 0.00 |
| `dthdt_lw` | Linf_error | 4.380e-09 | 7.158e-05 | 6.12e-05 | {"abs": 2e-08} | 0.00 |
| `rthraten` | Linf_error | 4.325e-09 | 8.214e-05 | 5.27e-05 | {"abs": 2e-08} | 0.00 |
| `sfc_br` | value | 5.652e-06 | 1.178e-01 | 4.80e-05 | {"rel": 0.0002} | 0.00 |
| `sfc_rmol` | value | 1.296e-06 | 2.782e-02 | 4.66e-05 | {"rel": 0.0002} | 0.00 |
| `sfc_psim` | value | 2.199e-05 | 9.457e-01 | 2.33e-05 | {"rel": 7e-05} | 0.00 |
| `sfc_psih` | value | 3.233e-05 | 1.631e+00 | 1.98e-05 | {"rel": 6e-05} | 0.00 |
| `sfc_flhc` | value | 4.325e-04 | 4.708e+01 | 9.19e-06 | {"rel": 3e-05} | 0.00 |
| `sfc_fh` | value | 3.226e-05 | 4.616e+00 | 6.99e-06 | {"rel": 3e-05} | 0.00 |
| `dthdt_sw` | Linf_error | 1.536e-10 | 2.725e-05 | 5.64e-06 | {"abs": 1e-09} | 0.00 |
| `sfc_qfx` | value | 4.224e-10 | 8.095e-05 | 5.22e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_flqc` | value | 6.379e-08 | 1.225e-02 | 5.21e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_lh` | value | 1.049e-03 | 2.024e+02 | 5.18e-06 | {"rel": 2e-05} | 0.00 |
| `pbl_hpbl` | value | 4.756e-03 | 9.256e+02 | 5.14e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_fm` | value | 2.215e-05 | 5.302e+00 | 4.18e-06 | {"rel": 2e-05} | 0.00 |
| `dtheta_m_dt_dyn` | Linf_error | 3.106e-11 | 1.022e-05 | 3.04e-06 | {"abs": 2e-10} | 0.00 |
| `sfc_ust` | value | 9.593e-07 | 4.545e-01 | 2.11e-06 | {"rel": 7e-06} | 0.00 |
| `sfc_q2` | value | 9.842e-09 | 5.563e-03 | 1.77e-06 | {"rel": 6e-06} | 0.00 |
| `gsw` | value | 8.489e-04 | 5.497e+02 | 1.54e-06 | {"rel": 1e-05} | 0.00 |
| `pbl_dqidt` | Linf_error | 1.180e-22 | 1.262e-16 | 9.35e-07 | {"abs": 5e-22} | 0.00 |
| `sfc_v10` | value | 2.992e-06 | 3.773e+00 | 7.93e-07 | {"rel": 1e-06} | 0.00 |
| `sfc_u10` | value | 2.729e-06 | 3.540e+00 | 7.71e-07 | {"rel": 1e-06} | 0.00 |
| `sfc_th2` | value | 1.658e-04 | 2.882e+02 | 5.75e-07 | {"rel": 1e-06} | 0.00 |
| `sfc_t2` | value | 1.639e-04 | 2.859e+02 | 5.73e-07 | {"rel": 1e-06} | 0.00 |
| `F_up` | Linf_error | 2.148e-04 | 3.831e+02 | 5.61e-07 | {"abs": 0.001} | 0.00 |
| `F_dn` | Linf_error | 1.222e-04 | 2.689e+02 | 4.54e-07 | {"abs": 0.0005} | 0.00 |
| `glw` | value | 1.222e-04 | 2.689e+02 | 4.54e-07 | {"rel": 2e-06} | 0.00 |
| `olr` | value | 8.071e-05 | 2.560e+02 | 3.15e-07 | {"rel": 2e-06} | 0.00 |
| `du_dt_dyn` | Linf_error | 1.228e-10 | 4.254e-04 | 2.89e-07 | {"abs": 5e-10} | 0.00 |
| `dw_dt_nonpg` | Linf_error | 3.303e-11 | 1.748e-04 | 1.89e-07 | {"abs": 2e-10} | 0.00 |
| `dv_dt_dyn` | Linf_error | 2.147e-11 | 1.295e-04 | 1.66e-07 | {"abs": 1e-10} | 0.00 |
| `dtheta_m_dt_dyn` | value | 9.041e-13 | 5.957e-06 | 1.52e-07 | {"rel": 1e-06} | 0.00 |
| `du_dt_dyn` | value | 1.142e-11 | 4.254e-04 | 2.68e-08 | {"rel": 2e-07} | 0.00 |
| `sfc_gz1oz0` | value | 1.633e-07 | 6.248e+00 | 2.61e-08 | {"rel": 2e-07} | 0.00 |
| `lsm_capg` | value | 5.285e-02 | 2.365e+06 | 2.24e-08 | {"rel": 1e-07} | 0.00 |
| `sfc_wspd` | value | 9.823e-08 | 6.030e+00 | 1.63e-08 | {"rel": 1e-07} | 0.00 |
| `dphi_dt` | Linf_error | 0.000e+00 | 7.411e-04 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_qsfc` | value | 0.000e+00 | 1.138e-02 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_znt` | value | 0.000e+00 | 5.000e-02 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_regime` | value | 0.000e+00 | 4.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `pbl_kpbl` | value | 0.000e+00 | 1.600e+01 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `lsm_land` | value | 0.000e+00 | 1.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `cldfra` | Linf_error | 0.000e+00 | 1.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `pbl_dqcdt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |

### `eqweather_scm_step60_cirrus`

| variable | reduce | abs resid | ref amplitude | rel err | tolerance | tol/amp |
|---|---|---|---|---|---|---|
| `mp_dqi_dt` | Linf_error | 2.812e-08 | 2.173e-08 | 1.29e+00 | {"abs": 1e-07} | 4.60 |
| `dw_dt` | Linf_error | 2.168e-04 | 3.713e-04 | 5.84e-01 | {"abs": 0.001} | 2.69 |
| `mp_dtheta_dt` | Linf_error | 4.318e-05 | 1.241e-04 | 3.48e-01 | {"abs": 0.0002} | 1.61 |
| `mp_dqv_dt` | Linf_error | 1.145e-08 | 3.300e-08 | 3.47e-01 | {"abs": 5e-08} | 1.52 |
| `dTsk_dt_wrf_step` | value | 1.026e-05 | 9.003e-05 | 1.14e-01 | {"rel": 0.15} | 0.15 |
| `lsm_dTs_dt` | Linf_error | 9.755e-06 | 1.780e-04 | 5.48e-02 | {"abs": 5e-05} | 0.28 |
| `pbl_dthdt` | Linf_error | 9.208e-07 | 1.808e-04 | 5.09e-03 | {"abs": 5e-06} | 0.03 |
| `du_dt` | Linf_error | 1.541e-07 | 4.862e-05 | 3.17e-03 | {"abs": 1e-06} | 0.02 |
| `rth_phys_sum` | Linf_error | 9.197e-07 | 4.068e-04 | 2.26e-03 | {"abs": 5e-06} | 0.01 |
| `pbl_dvdt` | Linf_error | 4.858e-07 | 3.032e-04 | 1.60e-03 | {"abs": 2e-06} | 0.01 |
| `dv_dt` | Linf_error | 4.858e-07 | 3.192e-04 | 1.52e-03 | {"abs": 2e-06} | 0.01 |
| `pbl_dqvdt` | Linf_error | 1.145e-10 | 8.461e-08 | 1.35e-03 | {"abs": 5e-10} | 0.01 |
| `lsm_hfx` | value | 6.563e-02 | 7.073e+01 | 9.28e-04 | {"rel": 0.003} | 0.00 |
| `pbl_dudt` | Linf_error | 1.541e-07 | 2.596e-04 | 5.93e-04 | {"abs": 1e-06} | 0.00 |
| `p_p` | Linf_error | 1.830e-02 | 3.809e+01 | 4.80e-04 | {"abs": 0.05} | 0.00 |
| `pbl_K_h` | Linf_error | 2.929e-02 | 1.117e+02 | 2.62e-04 | {"abs": 0.06} | 0.00 |
| `pbl_K_m` | Linf_error | 1.866e-02 | 7.282e+01 | 2.56e-04 | {"abs": 0.04} | 0.00 |
| `lsm_qfx` | value | 1.152e-08 | 7.957e-05 | 1.45e-04 | {"rel": 0.0005} | 0.00 |
| `lsm_lh` | value | 2.880e-02 | 1.989e+02 | 1.45e-04 | {"rel": 0.0005} | 0.00 |
| `lsm_qsfc` | value | 8.953e-07 | 9.437e-03 | 9.49e-05 | {"rel": 0.0005} | 0.00 |
| `dthdt_lw` | Linf_error | 3.904e-09 | 3.646e-04 | 1.07e-05 | {"abs": 2e-08} | 0.00 |
| `rthraten` | Linf_error | 3.965e-09 | 4.042e-04 | 9.81e-06 | {"abs": 2e-08} | 0.00 |
| `cldfra` | Linf_error | 9.509e-06 | 1.000e+00 | 9.51e-06 | {"abs": 5e-05} | 0.00 |
| `pbl_dqidt` | Linf_error | 4.071e-15 | 4.594e-10 | 8.86e-06 | {"abs": 2e-14} | 0.00 |
| `sfc_br` | value | 3.403e-07 | 4.449e-02 | 7.65e-06 | {"rel": 5e-05} | 0.00 |
| `sfc_rmol` | value | 8.034e-08 | 1.071e-02 | 7.50e-06 | {"rel": 5e-05} | 0.00 |
| `sfc_psim` | value | 2.582e-06 | 5.604e-01 | 4.61e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_mol` | value | 4.891e-07 | 1.165e-01 | 4.20e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_psih` | value | 4.234e-06 | 1.017e+00 | 4.16e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_hfx` | value | 2.830e-04 | 7.067e+01 | 4.00e-06 | {"rel": 2e-05} | 0.00 |
| `pbl_hpbl` | value | 3.317e-03 | 8.588e+02 | 3.86e-06 | {"rel": 2e-05} | 0.00 |
| `dthdt_sw` | Linf_error | 1.182e-10 | 5.915e-05 | 2.00e-06 | {"abs": 5e-10} | 0.00 |
| `gsw` | value | 4.557e-04 | 3.151e+02 | 1.45e-06 | {"rel": 1e-05} | 0.00 |
| `dtheta_m_dt_dyn` | Linf_error | 2.711e-10 | 1.940e-04 | 1.40e-06 | {"abs": 2e-09} | 0.00 |
| `sfc_flhc` | value | 5.071e-05 | 4.644e+01 | 1.09e-06 | {"rel": 5e-06} | 0.00 |
| `dw_dt_nonpg` | Linf_error | 3.381e-11 | 4.030e-05 | 8.39e-07 | {"abs": 2e-10} | 0.00 |
| `sfc_fh` | value | 4.291e-06 | 5.226e+00 | 8.21e-07 | {"rel": 5e-06} | 0.00 |
| `sfc_lh` | value | 1.404e-04 | 1.988e+02 | 7.06e-07 | {"rel": 3e-06} | 0.00 |
| `sfc_qfx` | value | 5.347e-11 | 7.952e-05 | 6.72e-07 | {"rel": 3e-06} | 0.00 |
| `F_up` | Linf_error | 2.299e-04 | 3.709e+02 | 6.20e-07 | {"abs": 0.001} | 0.00 |
| `sfc_flqc` | value | 7.591e-09 | 1.282e-02 | 5.92e-07 | {"rel": 3e-06} | 0.00 |
| `du_dt_dyn` | Linf_error | 1.176e-10 | 2.110e-04 | 5.57e-07 | {"abs": 5e-10} | 0.00 |
| `sfc_q2` | value | 2.130e-09 | 4.156e-03 | 5.12e-07 | {"rel": 3e-06} | 0.00 |
| `dv_dt_dyn` | Linf_error | 2.491e-11 | 4.874e-05 | 5.11e-07 | {"abs": 1e-10} | 0.00 |
| `F_dn` | Linf_error | 1.325e-04 | 2.928e+02 | 4.52e-07 | {"abs": 0.001} | 0.00 |
| `sfc_fm` | value | 2.341e-06 | 5.682e+00 | 4.12e-07 | {"rel": 2e-06} | 0.00 |
| `sfc_ust` | value | 1.012e-07 | 5.065e-01 | 2.00e-07 | {"rel": 1e-06} | 0.00 |
| `olr` | value | 2.535e-05 | 1.297e+02 | 1.95e-07 | {"rel": 1e-06} | 0.00 |
| `glw` | value | 5.332e-05 | 2.928e+02 | 1.82e-07 | {"rel": 1e-06} | 0.00 |
| `sfc_t2` | value | 4.000e-05 | 2.847e+02 | 1.41e-07 | {"rel": 1e-06} | 0.00 |
| `sfc_u10` | value | 2.360e-07 | 2.160e+00 | 1.09e-07 | {"rel": 5e-07} | 0.00 |
| `sfc_th2` | value | 2.238e-05 | 2.870e+02 | 7.80e-08 | {"rel": 5e-07} | 0.00 |
| `sfc_wspd` | value | 2.170e-07 | 7.182e+00 | 3.02e-08 | {"rel": 2e-07} | 0.00 |
| `lsm_capg` | value | 5.285e-02 | 2.365e+06 | 2.24e-08 | {"rel": 1e-07} | 0.00 |
| `sfc_v10` | value | 1.296e-07 | 5.840e+00 | 2.22e-08 | {"rel": 1e-07} | 0.00 |
| `sfc_gz1oz0` | value | 6.176e-08 | 6.243e+00 | 9.89e-09 | {"rel": 5e-08} | 0.00 |
| `sfc_qsfc` | value | 0.000e+00 | 9.434e-03 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_znt` | value | 0.000e+00 | 5.000e-02 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_regime` | value | 0.000e+00 | 4.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `pbl_kpbl` | value | 0.000e+00 | 1.500e+01 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `lsm_land` | value | 0.000e+00 | 1.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `pbl_dqcdt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqc_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqr_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqs_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqg_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_melt_dqi_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_rate_dqi_dt` | Linf_error | 2.154e-08 | 0.000e+00 | n/a | {"abs": 2.5e-08} | — |
| `mp_sed_dqi_dt` | Linf_error | 5.304e-08 | 0.000e+00 | n/a | {"abs": 6e-08} | — |

### `eqweather_scm_step900_night`

| variable | reduce | abs resid | ref amplitude | rel err | tolerance | tol/amp |
|---|---|---|---|---|---|---|
| `mp_dqi_dt` | Linf_error | 3.335e-12 | 1.235e-12 | 2.70e+00 | {"abs": 2e-11} | 16.19 |
| `mp_dtheta_dt` | Linf_error | 5.086e-07 | 5.086e-07 | 1.00e+00 | {"abs": 3e-06} | 5.90 |
| `dw_dt` | Linf_error | 2.244e-04 | 4.365e-04 | 5.14e-01 | {"abs": 0.001} | 2.29 |
| `mp_dqv_dt` | Linf_error | 1.338e-12 | 5.821e-12 | 2.30e-01 | {"abs": 1e-11} | 1.72 |
| `lsm_dTs_dt` | Linf_error | 4.578e-07 | 5.646e-05 | 8.11e-03 | {"abs": 2e-06} | 0.04 |
| `dtheta_m_conv_delta` | Linf_error | 1.895e-08 | 6.263e-06 | 3.03e-03 | {"abs": 1e-07} | 0.02 |
| `dTsk_dt_wrf_step` | value | 1.696e-07 | 5.646e-05 | 3.00e-03 | {"rel": 0.02} | 0.02 |
| `pbl_dthdt` | Linf_error | 7.106e-08 | 7.890e-05 | 9.01e-04 | {"abs": 3e-07} | 0.00 |
| `pbl_dqvdt` | Linf_error | 1.287e-11 | 1.516e-08 | 8.49e-04 | {"abs": 1e-10} | 0.01 |
| `dtheta_m_dt_phys` | Linf_error | 7.634e-08 | 9.085e-05 | 8.40e-04 | {"abs": 5e-07} | 0.01 |
| `dtheta_m_dt` | Linf_error | 7.634e-08 | 9.085e-05 | 8.40e-04 | {"abs": 5e-07} | 0.01 |
| `cldfra` | Linf_error | 7.698e-04 | 1.000e+00 | 7.70e-04 | {"abs": 0.005} | 0.01 |
| `rth_phys_sum` | Linf_error | 6.998e-08 | 9.521e-05 | 7.35e-04 | {"abs": 3e-07} | 0.00 |
| `p_p` | Linf_error | 1.156e-02 | 3.802e+01 | 3.04e-04 | {"abs": 0.05} | 0.00 |
| `lsm_hfx` | value | 5.029e-03 | 1.908e+01 | 2.64e-04 | {"rel": 0.002} | 0.00 |
| `lsm_qfx` | value | 6.965e-10 | 2.963e-06 | 2.35e-04 | {"rel": 0.001} | 0.00 |
| `lsm_lh` | value | 1.741e-03 | 7.407e+00 | 2.35e-04 | {"rel": 0.001} | 0.00 |
| `dv_dt` | Linf_error | 3.009e-08 | 1.799e-04 | 1.67e-04 | {"abs": 2e-07} | 0.00 |
| `dthdt_lw` | Linf_error | 6.268e-09 | 4.057e-05 | 1.54e-04 | {"abs": 3e-08} | 0.00 |
| `rthraten` | Linf_error | 6.268e-09 | 4.057e-05 | 1.54e-04 | {"abs": 3e-08} | 0.00 |
| `pbl_dvdt` | Linf_error | 3.008e-08 | 2.584e-04 | 1.16e-04 | {"abs": 2e-07} | 0.00 |
| `du_dt` | Linf_error | 2.229e-08 | 1.928e-04 | 1.16e-04 | {"abs": 1e-07} | 0.00 |
| `pbl_dudt` | Linf_error | 2.234e-08 | 5.332e-04 | 4.19e-05 | {"abs": 1e-07} | 0.00 |
| `lsm_qsfc` | value | 1.215e-07 | 5.455e-03 | 2.23e-05 | {"rel": 0.0001} | 0.00 |
| `sfc_rmol` | value | 1.870e-07 | 1.069e-02 | 1.75e-05 | {"rel": 0.0001} | 0.00 |
| `sfc_psim` | value | 2.522e-05 | 1.529e+00 | 1.65e-05 | {"rel": 0.0001} | 0.00 |
| `sfc_psih` | value | 3.112e-05 | 2.088e+00 | 1.49e-05 | {"rel": 0.0001} | 0.00 |
| `sfc_br` | value | 5.547e-07 | 3.732e-02 | 1.49e-05 | {"rel": 0.0001} | 0.00 |
| `sfc_hfx` | value | 2.300e-04 | 1.907e+01 | 1.21e-05 | {"rel": 5e-05} | 0.00 |
| `sfc_mol` | value | 6.110e-07 | 5.858e-02 | 1.04e-05 | {"rel": 5e-05} | 0.00 |
| `pbl_K_m` | Linf_error | 1.570e-05 | 1.757e+00 | 8.94e-06 | {"abs": 0.0001} | 0.00 |
| `pbl_K_h` | Linf_error | 1.574e-05 | 2.205e+00 | 7.14e-06 | {"abs": 0.0001} | 0.00 |
| `sfc_flhc` | value | 8.464e-05 | 1.566e+01 | 5.40e-06 | {"rel": 3e-05} | 0.00 |
| `sfc_flqc` | value | 2.180e-08 | 5.640e-03 | 3.87e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_qfx` | value | 1.141e-11 | 2.970e-06 | 3.84e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_lh` | value | 2.844e-05 | 7.425e+00 | 3.83e-06 | {"rel": 2e-05} | 0.00 |
| `sfc_fh` | value | 3.151e-05 | 8.315e+00 | 3.79e-06 | {"rel": 2e-05} | 0.00 |
| `pbl_hpbl` | value | 9.262e-04 | 2.625e+02 | 3.53e-06 | {"rel": 2e-05} | 0.00 |
| `pbl_dqidt` | Linf_error | 1.049e-20 | 3.159e-15 | 3.32e-06 | {"abs": 5e-20} | 0.00 |
| `sfc_fm` | value | 2.525e-05 | 7.756e+00 | 3.26e-06 | {"rel": 2e-05} | 0.00 |
| `F_dn` | Linf_error | 6.055e-04 | 2.519e+02 | 2.40e-06 | {"abs": 0.003} | 0.00 |
| `dtheta_m_dt_dyn` | Linf_error | 4.892e-11 | 2.854e-05 | 1.71e-06 | {"abs": 2e-10} | 0.00 |
| `sfc_u10` | value | 4.747e-06 | 2.952e+00 | 1.61e-06 | {"rel": 1e-05} | 0.00 |
| `sfc_ust` | value | 4.145e-07 | 2.666e-01 | 1.55e-06 | {"rel": 1e-05} | 0.00 |
| `sfc_v10` | value | 4.075e-06 | 2.629e+00 | 1.55e-06 | {"rel": 1e-05} | 0.00 |
| `F_up` | Linf_error | 2.220e-04 | 3.366e+02 | 6.60e-07 | {"abs": 0.001} | 0.00 |
| `glw` | value | 8.249e-05 | 2.519e+02 | 3.27e-07 | {"rel": 2e-06} | 0.00 |
| `olr` | value | 6.283e-05 | 2.401e+02 | 2.62e-07 | {"rel": 2e-06} | 0.00 |
| `du_dt_dyn` | Linf_error | 1.209e-10 | 4.950e-04 | 2.44e-07 | {"abs": 5e-10} | 0.00 |
| `sfc_q2` | value | 7.254e-10 | 5.100e-03 | 1.42e-07 | {"rel": 1e-06} | 0.00 |
| `dv_dt_dyn` | Linf_error | 2.512e-11 | 2.140e-04 | 1.17e-07 | {"abs": 2e-10} | 0.00 |
| `sfc_t2` | value | 2.536e-05 | 2.783e+02 | 9.11e-08 | {"rel": 5e-07} | 0.00 |
| `du_dt_dyn` | value | 4.461e-11 | 4.950e-04 | 9.01e-08 | {"rel": 5e-07} | 0.00 |
| `sfc_th2` | value | 2.217e-05 | 2.806e+02 | 7.90e-08 | {"rel": 5e-07} | 0.00 |
| `dw_dt_nonpg` | Linf_error | 2.109e-11 | 2.728e-04 | 7.73e-08 | {"abs": 1e-10} | 0.00 |
| `dtheta_m_dt_dyn` | value | 2.413e-13 | 4.513e-06 | 5.35e-08 | {"rel": 3e-07} | 0.00 |
| `sfc_wspd` | value | 2.123e-07 | 5.169e+00 | 4.11e-08 | {"rel": 2e-07} | 0.00 |
| `lsm_capg` | value | 5.285e-02 | 2.365e+06 | 2.24e-08 | {"rel": 1e-07} | 0.00 |
| `sfc_gz1oz0` | value | 8.677e-08 | 6.227e+00 | 1.39e-08 | {"rel": 1e-07} | 0.00 |
| `dphi_dt` | Linf_error | 0.000e+00 | 6.451e-04 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_qsfc` | value | 0.000e+00 | 5.457e-03 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_znt` | value | 0.000e+00 | 5.000e-02 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_regime` | value | 0.000e+00 | 1.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `pbl_kpbl` | value | 0.000e+00 | 6.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `lsm_land` | value | 0.000e+00 | 1.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `pbl_dqcdt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `dthdt_sw` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `gsw` | value | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqc_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqr_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqs_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqg_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |

### `eqweather_scm_step300_dusk`

| variable | reduce | abs resid | ref amplitude | rel err | tolerance | tol/amp |
|---|---|---|---|---|---|---|
| `mp_dqi_dt` | Linf_error | 7.747e-10 | 8.204e-11 | 9.44e+00 | {"abs": 5e-09} | 60.94 |
| `mp_dqv_dt` | Linf_error | 1.145e-09 | 1.144e-09 | 1.00e+00 | {"abs": 5e-09} | 4.37 |
| `mp_dtheta_dt` | Linf_error | 4.579e-06 | 4.578e-06 | 1.00e+00 | {"abs": 2e-05} | 4.37 |
| `dw_dt` | Linf_error | 2.456e-04 | 4.318e-04 | 5.69e-01 | {"abs": 0.001} | 2.32 |
| `dtheta_m_conv_delta` | Linf_error | 3.294e-06 | 2.829e-05 | 1.16e-01 | {"abs": 2e-05} | 0.71 |
| `dtheta_m_dt` | Linf_error | 3.293e-06 | 1.948e-04 | 1.69e-02 | {"abs": 2e-05} | 0.10 |
| `dtheta_m_dt_phys` | Linf_error | 3.293e-06 | 1.948e-04 | 1.69e-02 | {"abs": 2e-05} | 0.10 |
| `pbl_dthdt` | Linf_error | 6.601e-07 | 1.591e-04 | 4.15e-03 | {"abs": 3e-06} | 0.02 |
| `dTsk_dt_wrf_step` | value | 1.089e-06 | 2.680e-04 | 4.06e-03 | {"rel": 0.02} | 0.02 |
| `lsm_dTs_dt` | Linf_error | 1.089e-06 | 3.006e-04 | 3.62e-03 | {"abs": 5e-06} | 0.02 |
| `rth_phys_sum` | Linf_error | 6.641e-07 | 2.230e-04 | 2.98e-03 | {"abs": 3e-06} | 0.01 |
| `lsm_hfx` | value | 4.221e-02 | 1.913e+01 | 2.21e-03 | {"rel": 0.01} | 0.01 |
| `pbl_dqvdt` | Linf_error | 1.399e-10 | 6.510e-08 | 2.15e-03 | {"abs": 1e-09} | 0.02 |
| `lsm_lh` | value | 1.772e-02 | 2.624e+01 | 6.75e-04 | {"rel": 0.003} | 0.00 |
| `lsm_qfx` | value | 7.088e-09 | 1.049e-05 | 6.75e-04 | {"rel": 0.003} | 0.00 |
| `p_p` | Linf_error | 1.512e-02 | 3.812e+01 | 3.97e-04 | {"abs": 0.1} | 0.00 |
| `lsm_qsfc` | value | 1.957e-06 | 7.117e-03 | 2.75e-04 | {"rel": 0.002} | 0.00 |
| `sfc_rmol` | value | 4.634e-06 | 1.839e-02 | 2.52e-04 | {"rel": 0.002} | 0.00 |
| `sfc_psim` | value | 6.137e-04 | 2.593e+00 | 2.37e-04 | {"rel": 0.001} | 0.00 |
| `pbl_dvdt` | Linf_error | 3.792e-08 | 1.629e-04 | 2.33e-04 | {"abs": 2e-07} | 0.00 |
| `dv_dt` | Linf_error | 3.792e-08 | 1.719e-04 | 2.21e-04 | {"abs": 2e-07} | 0.00 |
| `sfc_psih` | value | 6.387e-04 | 3.292e+00 | 1.94e-04 | {"rel": 0.001} | 0.00 |
| `sfc_br` | value | 1.037e-05 | 5.759e-02 | 1.80e-04 | {"rel": 0.001} | 0.00 |
| `sfc_flhc` | value | 1.067e-03 | 1.048e+01 | 1.02e-04 | {"rel": 0.0005} | 0.00 |
| `du_dt` | Linf_error | 2.267e-08 | 2.643e-04 | 8.58e-05 | {"abs": 1e-07} | 0.00 |
| `sfc_qfx` | value | 8.219e-10 | 1.053e-05 | 7.81e-05 | {"rel": 0.0005} | 0.00 |
| `sfc_flqc` | value | 3.154e-07 | 4.041e-03 | 7.81e-05 | {"rel": 0.0005} | 0.00 |
| `sfc_lh` | value | 2.054e-03 | 2.632e+01 | 7.80e-05 | {"rel": 0.0005} | 0.00 |
| `sfc_fm` | value | 6.140e-04 | 8.834e+00 | 6.95e-05 | {"rel": 0.0003} | 0.00 |
| `sfc_fh` | value | 6.392e-04 | 9.534e+00 | 6.70e-05 | {"rel": 0.0003} | 0.00 |
| `pbl_K_m` | Linf_error | 1.657e-04 | 2.551e+00 | 6.50e-05 | {"abs": 0.001} | 0.00 |
| `dthdt_lw` | Linf_error | 4.007e-09 | 6.395e-05 | 6.27e-05 | {"abs": 2e-08} | 0.00 |
| `rthraten` | Linf_error | 4.007e-09 | 6.395e-05 | 6.27e-05 | {"abs": 2e-08} | 0.00 |
| `pbl_dudt` | Linf_error | 2.269e-08 | 3.739e-04 | 6.07e-05 | {"abs": 1e-07} | 0.00 |
| `sfc_mol` | value | 3.969e-06 | 7.642e-02 | 5.19e-05 | {"rel": 0.0003} | 0.00 |
| `pbl_K_h` | Linf_error | 1.657e-04 | 3.256e+00 | 5.09e-05 | {"abs": 0.001} | 0.00 |
| `cldfra` | Linf_error | 4.694e-05 | 1.000e+00 | 4.69e-05 | {"abs": 0.0002} | 0.00 |
| `sfc_ust` | value | 7.213e-06 | 2.078e-01 | 3.47e-05 | {"rel": 0.0002} | 0.00 |
| `sfc_v10` | value | 6.418e-05 | 2.193e+00 | 2.93e-05 | {"rel": 0.0002} | 0.00 |
| `sfc_u10` | value | 7.212e-05 | 2.471e+00 | 2.92e-05 | {"rel": 0.0002} | 0.00 |
| `pbl_hpbl` | value | 1.355e-02 | 5.558e+02 | 2.44e-05 | {"rel": 0.0001} | 0.00 |
| `sfc_hfx` | value | 3.289e-04 | 1.909e+01 | 1.72e-05 | {"rel": 0.0001} | 0.00 |
| `sfc_q2` | value | 5.636e-08 | 5.500e-03 | 1.02e-05 | {"rel": 5e-05} | 0.00 |
| `pbl_dqidt` | Linf_error | 2.055e-17 | 8.911e-12 | 2.31e-06 | {"abs": 1e-16} | 0.00 |
| `F_dn` | Linf_error | 4.581e-04 | 2.644e+02 | 1.73e-06 | {"abs": 0.002} | 0.00 |
| `F_up` | Linf_error | 2.592e-04 | 3.521e+02 | 7.36e-07 | {"abs": 0.002} | 0.00 |
| `dtheta_m_dt_dyn` | Linf_error | 1.458e-11 | 2.568e-05 | 5.68e-07 | {"abs": 1e-10} | 0.00 |
| `glw` | value | 7.917e-05 | 2.644e+02 | 2.99e-07 | {"rel": 2e-06} | 0.00 |
| `du_dt_dyn` | Linf_error | 1.219e-10 | 5.298e-04 | 2.30e-07 | {"abs": 5e-10} | 0.00 |
| `dv_dt_dyn` | Linf_error | 2.703e-11 | 1.729e-04 | 1.56e-07 | {"abs": 2e-10} | 0.00 |
| `sfc_th2` | value | 3.713e-05 | 2.846e+02 | 1.30e-07 | {"rel": 1e-06} | 0.00 |
| `sfc_t2` | value | 3.592e-05 | 2.823e+02 | 1.27e-07 | {"rel": 1e-06} | 0.00 |
| `dw_dt_nonpg` | Linf_error | 2.929e-11 | 2.354e-04 | 1.24e-07 | {"abs": 2e-10} | 0.00 |
| `olr` | value | 2.191e-05 | 2.392e+02 | 9.16e-08 | {"rel": 5e-07} | 0.00 |
| `dtheta_m_dt_dyn` | value | 2.403e-13 | 4.513e-06 | 5.33e-08 | {"rel": 3e-07} | 0.00 |
| `lsm_capg` | value | 5.285e-02 | 2.365e+06 | 2.24e-08 | {"rel": 1e-07} | 0.00 |
| `du_dt_dyn` | value | 9.115e-12 | 5.298e-04 | 1.72e-08 | {"rel": 1e-07} | 0.00 |
| `sfc_gz1oz0` | value | 6.725e-08 | 6.242e+00 | 1.08e-08 | {"rel": 5e-08} | 0.00 |
| `sfc_wspd` | value | 4.319e-08 | 4.588e+00 | 9.41e-09 | {"rel": 5e-08} | 0.00 |
| `dphi_dt` | Linf_error | 0.000e+00 | 6.645e-04 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_qsfc` | value | 0.000e+00 | 7.125e-03 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_znt` | value | 0.000e+00 | 5.000e-02 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `sfc_regime` | value | 0.000e+00 | 1.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `pbl_kpbl` | value | 0.000e+00 | 1.000e+01 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `lsm_land` | value | 0.000e+00 | 1.000e+00 | 0.00e+00 | {"abs": 0.0} | 0.00 |
| `pbl_dqcdt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `dthdt_sw` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `gsw` | value | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqc_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqr_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqs_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
| `mp_dqg_dt` | Linf_error | 0.000e+00 | 0.000e+00 | n/a | {"abs": 0.0} | — |
