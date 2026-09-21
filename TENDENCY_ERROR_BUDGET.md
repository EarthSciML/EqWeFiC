# EqWeather-SCM tendency error budget

**Measured 2026-09-21** on `couplings/eqweather_scm.esm`, EarthSciAST `main` at `be1fafe6d`
(`./esm` built 2026-09-21; carries #439 and #440), `ESD_ROOT` = the pinned
`EarthSciDiscretizations-column` worktree at `f80069a`.

## Why this document exists

The document's 134 assertions are green, but green only means *each residual is under the
tolerance that assertion was given*. Those tolerances were set by measuring the residual and
backing off a factor of 4-10, which records how well the column happened to agree, not how
well it is required to agree. This file states the actual agreement, per assertion, so that
"correct tendencies" is a measured claim rather than a pass/fail.

## Method

The document was copied with **every tolerance set to `{abs: 0, rel: 0}`** and re-run, which
forces each assertion to report `actual` against `expected`. For a `Linf_error` assertion the
reported `actual` IS the absolute error, since `expected` is 0. The reference amplitude is
`max |reference|` taken from the assertion's own `faq` reference array (or `|expected|` for a
scalar assertion), and the relative error is the ratio of the two. Cost: 432.97 s, peak 495 MB,
115 of 134 assertions reporting a non-zero residual and 19 bit-exact.

The copy was a scratch file (`couplings/_errbudget.esm`); it has been deleted. Nothing in the
committed document changed, and no tolerance in it was altered by this exercise.

## Result by subsystem

| subsystem | n | median rel err | worst rel err |
|---|---|---|---|
| Dudhia shortwave | 4 | 1.8e-06 | 5.6e-06 |
| RRTM longwave | 12 | 5.1e-07 | 6.1e-05 |
| surface layer (sfclayrev) | 40 | 1.4e-06 | 7.3e-05 |
| WRF-ARW dynamics | 12 | 5.3e-07 | 4.8e-04 |
| slab LSM | 14 | 3.2e-04 | 1.1e-01 |
| **YSU PBL** | 16 | 2.0e-02 | **3.1e-01** |
| **WSM6 microphysics** | 5 | 1.3e+00 | **1.3e+00** |
| coupled sums (dyn+phys) | 12 | 2.7e-01 | 6.7e+00 |

The column splits cleanly in two. **Dynamics, radiation and the surface layer agree with WRF
to 1e-7 - 1e-4** — that is at or near single-precision round-off of the reference build, and
those tendencies can be called correct. **The PBL and microphysics tendencies agree only to
tens of percent**, and every entry in the "coupled sums" row is large only because it contains
a PBL or microphysics term.

## The finding that matters: 15 assertions constrain nothing

The last column below is tolerance divided by reference amplitude. When it is **>= 1 the
tolerance is larger than the entire signal being asserted**, so the assertion would pass against
a tendency of zero — or against one of the wrong sign.

| tol/amp | test | variable | tolerance | ref amplitude | measured rel err |
|---|---|---|---|---|---|
| 23.43 | step60_cirrus | `du_dt` | abs 2e-4 | 8.54e-06 | 6.71e+00 |
| 5.53 | step60_cirrus | `pbl_dthdt` | abs 1e-3 | 1.81e-04 | 1.55e-01 |
| 4.60 | step60_cirrus | `mp_dqi_dt` | abs 1e-7 | 2.17e-08 | 1.29e+00 |
| 3.47 | step1410 | `dtheta_m_dt_phys` | abs 1e-3 | 2.88e-04 | 3.91e-01 |
| 2.70 | step1410 | `pbl_dthdt` | abs 1e-3 | 3.71e-04 | 3.06e-01 |
| 2.69 | step60_cirrus | `dw_dt` | abs 1e-3 | 3.71e-04 | 5.84e-01 |
| 2.63 | step1410 | `rth_phys_sum` | abs 1e-3 | 3.80e-04 | 2.99e-01 |
| 2.47 | step60_cirrus | `pbl_dudt` | abs 5e-4 | 2.02e-04 | 2.83e-01 |
| 2.46 | step60_cirrus | `rth_phys_sum` | abs 1e-3 | 4.07e-04 | 6.90e-02 |
| 2.36 | step1410 | `dw_dt` | abs 1e-3 | 4.23e-04 | 7.19e-01 |
| 1.81 | step1410 | `pbl_dudt` | abs 5e-4 | 2.76e-04 | 2.24e-01 |
| 1.61 | step60_cirrus | `mp_dtheta_dt` | abs 2e-4 | 1.24e-04 | 3.48e-01 |
| 1.52 | step60_cirrus | `mp_dqv_dt` | abs 5e-8 | 3.30e-08 | 3.47e-01 |
| 1.39 | step1410 | `dtheta_m_dt` | abs 4e-4 | 2.88e-04 | 3.91e-01 |
| 1.19 | step60_cirrus | `pbl_dqvdt` | abs 1e-7 | 8.37e-08 | 2.88e-02 |

Three more are weak (tol/amp 0.6-0.8): `step1410 du_dt`, `step60 pbl_dvdt`, `step1410 pbl_dvdt`.

`du_dt` at step 60 is the extreme case and worth reading carefully: its reference amplitude is
8.5e-06 because the dynamics and PBL contributions very nearly cancel there, while the absolute
error (5.73e-05) is exactly `pbl_dudt`'s absolute error. The assertion is not testing a small
tendency accurately; it is testing a large error against a tolerance 23x the signal.

## What is NOT wrong

**Vertical advection is not missing.** The five mounted WRF-ARW components carry no advection
term, and `du_dt_dyn` is only Coriolis + curvature + Rayleigh. That is correct for this
configuration, confirmed three ways from `data/eqwefic/dumps/dyn_column`:

- `dyn_ww` (the eta-coordinate vertical mass flux, i.e. omega) is **identically 0.0** at all 60
  interfaces at steps 60, 900 and 1410;
- `dyn_mu` is **-97.750000 at all three steps**, so d(mu)/dt = 0 over 1350 steps;
- the instrumented dump contains `dru_coriolis`, `dru_curvature`, `dru_rayleigh`,
  `drw_pg_buoy`, `drw_rayleigh`, `drw_w_damp` and **no advection term at all**.

In eta coordinates the surfaces are mass-following, and omega is diagnosed from the
column-integrated horizontal mass divergence, which is identically zero in a horizontally
homogeneous single column. So `w` is non-zero (a buoyancy oscillation) but transports nothing.
`drw_w_damp` is also identically zero (w-damping off), and correctly not mounted.

**WSM6 sedimentation is not wrong.** An earlier note in this project compared the ESM
sedimentation ice tendency (5.30e-08 kg/kg/s) against WRF's *whole* ice tendency (2.17e-08) and
called the 2.4x ratio an error. That comparison was invalid in two ways: it set one process
against the sum of all microphysical processes, which partly cancel; and WRF's fallout is
`nislfv_rain_plm` (Juang & Hong 2010), a **fully discrete forward semi-Lagrangian scheme whose
Courant number lives inside the operator**, so its increment over a step is not `dt` times any
derivative. `components/atmospheric_dynamics/wsm6/sedimentation.esm` already pins the correct
target — the dt -> 0 limit of the kernel's own `fall` fluxes — and is **65/65 green**. No PLM
discretization rule is needed for tendency correctness.

## Where the remaining error actually is: in the REFERENCE, not the model

The PBL and microphysics residuals are **not** ESM errors. Both are being compared against a
WRF quantity that is not an instantaneous derivative.

### YSU: the model is right to 1e-9; the reference is a dt = 60 s implicit step

Three measurements, in order:

1. **YSU standalone is accurate to ~1e-7 relative** on `dudt`/`dvdt`/`dthdt`/`dqvdt` across all
   five regimes (convective, stable, shallow stable, cloud-top, first step), on full-column
   `Linf_error` assertions at `abs 1e-9` — including `scm_call120_convective`, which is the same
   physics call as the coupled step-60 test. So the scheme is not wrong.
2. **The coupled residual does not come from the coupling.** `scm_physics_column` (YSU fed the
   DUMPED geometry) and `geometry_surface_pbl_column` (YSU fed the `Geo`-derived geometry) give
   residuals identical to four significant figures — `pbl_dthdt` 3.07e-01, `pbl_dudt` 2.23e-01
   vs 2.24e-01, `pbl_dvdt` 9.14e-02, `pbl_dqvdt` 4.26e-02. Geometry, radiation and the surface
   layer are eliminated as causes.
3. **The residual IS WRF's own implicit-step error.** Richardson-extrapolating the YSU kernel
   replays at dt = 0.04/0.02/0.01 s, `(8T(dt/4) - 6T(dt/2) + T(dt))/3`, and comparing the dt -> 0
   limit against the dt = 60 s step that the coupled document asserts against:

   | | dt->0 vs dt=60 s | coupled document residual |
   |---|---|---|
   | `utnp` | 24.6% | `pbl_dudt` 22.4% |
   | `ttnp` | 41.1% | `pbl_dthdt` 30.7% |
   | `vtnp` | 5.9% | `pbl_dvdt` 9.1% |
   | `qvtnp` | 2.2% | `pbl_dqvdt` 4.3% |

   Same magnitudes, same ordering. The extrapolation is self-consistent to 7.5e-10 (`utnp`) and
   4.9e-08 (`ttnp`) when recomputed from dt = 0.08/0.04/0.02 instead.

The standalone YSU tests already use the extrapolated dt -> 0 reference and pass at 1e-9. The
coupled document uses the driver-level `rthblten`/`rublten`/... at dt = 60 s, and its own
description says so: *"The PBL TENDENCIES are asserted only LOOSELY and document a known gap
rather than test the physics."* The budget above quantifies that gap; it does not find a defect.

### WSM6: the same category of mismatch

`mp_dqv_dt`, `mp_dqi_dt` and `mp_dtheta_dt` are compared against WRF increments over a 60 s
sub-step from a scheme that stages its processes sequentially. The step-60 test description
already states this in terms: *"THE THREE MICROPHYSICAL TENDENCIES DO NOT AGREE WITH WRF AND
THEIR TOLERANCES RECORD THAT, they do not certify it."* The O(dt_substep) operator-split
difference has still never been measured.

### So the actual state of the column

Every tendency in the document that is asserted against a genuine instantaneous derivative
agrees with WRF to 1e-7 - 1e-4. Every tendency that looks wrong is asserted against a
finite-dt increment. **There is currently no evidence of an incorrect tendency anywhere in
EqWeather-SCM.** What there is, is 15 assertions that constrain nothing and 21 that are pinned
to the wrong kind of reference.

## Next actions implied

1. **Re-reference the coupled PBL assertions to the dt -> 0 limit.** For step 60 / call 120 the
   replays already exist (`data/eqwefic/dumps/ysu/driver_120_dt{0.08,0.04,0.02,0.01,0.005,0.001}.json`),
   so this needs no new Fortran run. Step 1410 is physics call 2820, for which **no driver replay
   exists** — it needs one from `dumps/ysu/bin`. Once re-referenced, the six PBL tolerances can
   go from `abs 1e-3`/`5e-4`/`2e-4`/`1e-7` down to roughly 1e-9, and the 15 vacuous assertions
   mostly stop being vacuous.
2. **Measure the WSM6 operator-split O(dt_substep) difference** and re-reference `mp_*` the same
   way. Note the whole-scheme increment does NOT converge as dt -> 0 — the saturation adjustment
   is a projection to equilibrium, not a rate, and `(x_out - x_in)/dt` from the wsm6 driver
   replays scales exactly as 1/dt (2.3063e-02, 1.1529e-02, 5.7628e-03 at dt = 0.01/0.02/0.04).
   So the extrapolation must be done per process stage, not on the scheme as a whole.
3. Only then tighten tolerances, deriving each from its measured residual.
4. Add the night (step 900) and step 300 columns.

## Per-assertion detail


### Test 1 — `eqweather_scm_step1410` (clear sky, 12:02 local solar)

| variable | reduce | abs err | ref amplitude | **rel err** | tolerance | tol/amp |
|---|---|---|---|---|---|---|
| `dw_dt` | Linf_error | 3.042e-04 | 4.233e-04 | **7.19e-01** | abs 0.001 | 2.36 |
| `dtheta_m_dt_phys` | Linf_error | 1.128e-04 | 2.884e-04 | **3.91e-01** | abs 0.001 | 3.47 |
| `dtheta_m_dt` | Linf_error | 1.128e-04 | 2.884e-04 | **3.91e-01** | abs 0.0004 | 1.39 |
| `pbl_dthdt` | Linf_error | 1.134e-04 | 3.710e-04 | **3.06e-01** | abs 0.001 | 2.70 |
| `rth_phys_sum` | Linf_error | 1.134e-04 | 3.796e-04 | **2.99e-01** | abs 0.001 | 2.63 |
| `du_dt` | Linf_error | 6.178e-05 | 2.607e-04 | **2.37e-01** | abs 0.0002 | 0.77 |
| `pbl_dudt` | Linf_error | 6.178e-05 | 2.757e-04 | **2.24e-01** | abs 0.0005 | 1.81 |
| `dv_dt` | Linf_error | 2.958e-05 | 2.375e-04 | **1.25e-01** | abs 0.0001 | 0.42 |
| `pbl_dvdt` | Linf_error | 2.958e-05 | 3.233e-04 | **9.15e-02** | abs 0.0002 | 0.62 |
| `dtheta_m_conv_delta` | Linf_error | 6.602e-06 | 1.537e-04 | **4.29e-02** | abs 2e-05 | 0.13 |
| `pbl_dqvdt` | Linf_error | 1.410e-08 | 3.305e-07 | **4.27e-02** | abs 1e-07 | 0.30 |
| `dTsk_dt_wrf_step` | value | 5.788e-06 | 1.363e-04 | **4.25e-02** | rel 0.15 | 0.15 |
| `lsm_dTs_dt` | Linf_error | 5.788e-06 | 2.660e-04 | **2.18e-02** | abs 1e-05 | 0.04 |
| `pbl_dqidt` | Linf_error | 1.524e-18 | 1.277e-16 | **1.19e-02** | abs 1e-17 | 0.08 |
| `lsm_hfx` | value | 1.435e-01 | 1.736e+02 | 8.26e-04 | rel 0.001 | 0.00 |
| `p_p` | Linf_error | 1.571e-02 | 3.798e+01 | 4.14e-04 | abs 0.05 | 0.00 |
| `lsm_qfx` | value | 2.565e-08 | 8.103e-05 | 3.17e-04 | rel 0.0005 | 0.00 |
| `lsm_lh` | value | 6.412e-02 | 2.026e+02 | 3.17e-04 | rel 0.0005 | 0.00 |
| `pbl_K_h` | Linf_error | 4.833e-02 | 1.768e+02 | 2.73e-04 | abs 0.06 | 0.00 |
| `pbl_K_m` | Linf_error | 2.736e-02 | 1.026e+02 | 2.67e-04 | abs 0.04 | 0.00 |
| `lsm_qsfc` | value | 2.059e-06 | 1.139e-02 | 1.81e-04 | rel 0.0005 | 0.00 |
| `sfc_hfx` | value | 1.272e-02 | 1.735e+02 | 7.33e-05 | rel 0.0003 | 0.00 |
| `sfc_mol` | value | 2.271e-05 | 3.192e-01 | 7.12e-05 | rel 0.0003 | 0.00 |
| `dthdt_lw` | Linf_error | 4.380e-09 | 7.158e-05 | 6.12e-05 | abs 5e-08 | 0.00 |
| `rthraten` | Linf_error | 4.325e-09 | 8.214e-05 | 5.27e-05 | abs 5e-08 | 0.00 |
| `sfc_br` | value | 5.652e-06 | 1.178e-01 | 4.80e-05 | rel 0.0002 | 0.00 |
| `sfc_rmol` | value | 1.296e-06 | 2.782e-02 | 4.66e-05 | rel 0.0002 | 0.00 |
| `sfc_psim` | value | 2.199e-05 | 9.457e-01 | 2.33e-05 | rel 7e-05 | 0.00 |
| `sfc_psih` | value | 3.233e-05 | 1.631e+00 | 1.98e-05 | rel 6e-05 | 0.00 |
| `sfc_flhc` | value | 4.325e-04 | 4.708e+01 | 9.19e-06 | rel 3e-05 | 0.00 |
| `sfc_fh` | value | 3.226e-05 | 4.616e+00 | 6.99e-06 | rel 3e-05 | 0.00 |
| `dthdt_sw` | Linf_error | 1.536e-10 | 2.725e-05 | 5.64e-06 | abs 1e-09 | 0.00 |
| `sfc_qfx` | value | 4.224e-10 | 8.095e-05 | 5.22e-06 | rel 2e-05 | 0.00 |
| `sfc_flqc` | value | 6.379e-08 | 1.225e-02 | 5.21e-06 | rel 2e-05 | 0.00 |
| `sfc_lh` | value | 1.049e-03 | 2.024e+02 | 5.18e-06 | rel 2e-05 | 0.00 |
| `pbl_hpbl` | value | 4.756e-03 | 9.256e+02 | 5.14e-06 | rel 2e-05 | 0.00 |
| `sfc_fm` | value | 2.215e-05 | 5.302e+00 | 4.18e-06 | rel 2e-05 | 0.00 |
| `dtheta_m_dt_dyn` | Linf_error | 3.106e-11 | 1.022e-05 | 3.04e-06 | abs 2e-09 | 0.00 |
| `sfc_ust` | value | 9.593e-07 | 4.545e-01 | 2.11e-06 | rel 7e-06 | 0.00 |
| `sfc_q2` | value | 9.842e-09 | 5.563e-03 | 1.77e-06 | rel 6e-06 | 0.00 |
| `gsw` | value | 8.489e-04 | 5.497e+02 | 1.54e-06 | rel 1e-05 | 0.00 |
| `sfc_v10` | value | 2.992e-06 | 3.773e+00 | 7.93e-07 | rel 1e-06 | 0.00 |
| `sfc_u10` | value | 2.729e-06 | 3.540e+00 | 7.71e-07 | rel 1e-06 | 0.00 |
| `sfc_th2` | value | 1.658e-04 | 2.882e+02 | 5.75e-07 | rel 1e-06 | 0.00 |
| `sfc_t2` | value | 1.639e-04 | 2.859e+02 | 5.73e-07 | rel 1e-06 | 0.00 |
| `F_up` | Linf_error | 2.148e-04 | 3.831e+02 | 5.61e-07 | abs 0.002 | 0.00 |
| `F_dn` | Linf_error | 1.222e-04 | 2.689e+02 | 4.54e-07 | abs 0.001 | 0.00 |
| `glw` | value | 1.222e-04 | 2.689e+02 | 4.54e-07 | rel 2e-06 | 0.00 |
| `olr` | value | 8.071e-05 | 2.560e+02 | 3.15e-07 | rel 2e-06 | 0.00 |
| `du_dt_dyn` | Linf_error | 1.228e-10 | 4.254e-04 | 2.89e-07 | abs 1e-09 | 0.00 |
| `dw_dt_nonpg` | Linf_error | 3.303e-11 | 1.748e-04 | 1.89e-07 | abs 3e-10 | 0.00 |
| `dv_dt_dyn` | Linf_error | 2.147e-11 | 1.295e-04 | 1.66e-07 | abs 2e-10 | 0.00 |
| `dtheta_m_dt_dyn` | value | 9.041e-13 | 5.957e-06 | 1.52e-07 | rel 1e-05 | 0.00 |
| `du_dt_dyn` | value | 1.142e-11 | 4.254e-04 | 2.68e-08 | rel 1e-05 | 0.00 |
| `sfc_gz1oz0` | value | 1.633e-07 | 6.248e+00 | 2.61e-08 | rel 1e-06 | 0.00 |
| `lsm_capg` | value | 5.285e-02 | 2.365e+06 | 2.24e-08 | rel 1e-06 | 0.00 |
| `sfc_wspd` | value | 9.823e-08 | 6.030e+00 | 1.63e-08 | rel 1e-06 | 0.00 |
| `dphi_dt` | Linf_error | **0 (bit-exact)** | 7.411e-04 | 0 | abs 1e-12 | 0.00 |
| `sfc_qsfc` | value | **0 (bit-exact)** | 1.138e-02 | 0 | abs 0 | — |
| `sfc_znt` | value | **0 (bit-exact)** | 5.000e-02 | 0 | abs 0 | — |
| `sfc_regime` | value | **0 (bit-exact)** | 4.000e+00 | 0 | abs 0 | — |
| `pbl_kpbl` | value | **0 (bit-exact)** | 1.600e+01 | 0 | abs 0 | — |
| `pbl_dqcdt` | Linf_error | **0 (bit-exact)** | 0.000e+00 | 0 | abs 0 | — |
| `lsm_land` | value | **0 (bit-exact)** | 1.000e+00 | 0 | abs 0 | — |
| `cldfra` | Linf_error | **0 (bit-exact)** | 1.000e+00 | 0 | abs 0 | — |

### Test 2 — `eqweather_scm_step60_cirrus` (cirrus deck aloft)

| variable | reduce | abs err | ref amplitude | **rel err** | tolerance | tol/amp |
|---|---|---|---|---|---|---|
| `du_dt` | Linf_error | 5.732e-05 | 8.537e-06 | **6.71e+00** | abs 0.0002 | 23.43 |
| `mp_dqi_dt` | Linf_error | 2.812e-08 | 2.173e-08 | **1.29e+00** | abs 1e-07 | 4.60 |
| `dw_dt` | Linf_error | 2.168e-04 | 3.713e-04 | **5.84e-01** | abs 0.001 | 2.69 |
| `mp_dtheta_dt` | Linf_error | 4.318e-05 | 1.241e-04 | **3.48e-01** | abs 0.0002 | 1.61 |
| `mp_dqv_dt` | Linf_error | 1.145e-08 | 3.300e-08 | **3.47e-01** | abs 5e-08 | 1.52 |
| `pbl_dudt` | Linf_error | 5.732e-05 | 2.025e-04 | **2.83e-01** | abs 0.0005 | 2.47 |
| `pbl_dthdt` | Linf_error | 2.809e-05 | 1.808e-04 | **1.55e-01** | abs 0.001 | 5.53 |
| `dTsk_dt_wrf_step` | value | 1.026e-05 | 9.003e-05 | **1.14e-01** | rel 0.15 | 0.15 |
| `rth_phys_sum` | Linf_error | 2.809e-05 | 4.068e-04 | **6.90e-02** | abs 0.001 | 2.46 |
| `lsm_dTs_dt` | Linf_error | 9.755e-06 | 1.780e-04 | **5.48e-02** | abs 5e-05 | 0.28 |
| `pbl_dvdt` | Linf_error | 1.601e-05 | 2.963e-04 | **5.40e-02** | abs 0.0002 | 0.67 |
| `dv_dt` | Linf_error | 1.601e-05 | 3.111e-04 | **5.15e-02** | abs 0.0001 | 0.32 |
| `pbl_dqvdt` | Linf_error | 2.415e-09 | 8.375e-08 | **2.88e-02** | abs 1e-07 | 1.19 |
| `pbl_dqidt` | Linf_error | 2.441e-12 | 4.570e-10 | 5.34e-03 | abs 2e-11 | 0.04 |
| `lsm_hfx` | value | 6.563e-02 | 7.073e+01 | 9.28e-04 | rel 0.003 | 0.00 |
| `p_p` | Linf_error | 1.830e-02 | 3.809e+01 | 4.80e-04 | abs 0.05 | 0.00 |
| `pbl_K_h` | Linf_error | 2.929e-02 | 1.117e+02 | 2.62e-04 | abs 0.06 | 0.00 |
| `pbl_K_m` | Linf_error | 1.866e-02 | 7.282e+01 | 2.56e-04 | abs 0.04 | 0.00 |
| `lsm_qfx` | value | 1.152e-08 | 7.957e-05 | 1.45e-04 | rel 0.0005 | 0.00 |
| `lsm_lh` | value | 2.880e-02 | 1.989e+02 | 1.45e-04 | rel 0.0005 | 0.00 |
| `lsm_qsfc` | value | 8.953e-07 | 9.437e-03 | 9.49e-05 | rel 0.0005 | 0.00 |
| `dthdt_lw` | Linf_error | 3.904e-09 | 3.646e-04 | 1.07e-05 | abs 5e-08 | 0.00 |
| `rthraten` | Linf_error | 3.965e-09 | 4.042e-04 | 9.81e-06 | abs 5e-08 | 0.00 |
| `cldfra` | Linf_error | 9.509e-06 | 1.000e+00 | 9.51e-06 | abs 0.0001 | 0.00 |
| `sfc_br` | value | 3.403e-07 | 4.449e-02 | 7.65e-06 | rel 0.0002 | 0.00 |
| `sfc_rmol` | value | 8.034e-08 | 1.071e-02 | 7.50e-06 | rel 0.0002 | 0.00 |
| `sfc_psim` | value | 2.582e-06 | 5.604e-01 | 4.61e-06 | rel 7e-05 | 0.00 |
| `sfc_mol` | value | 4.891e-07 | 1.165e-01 | 4.20e-06 | rel 0.0003 | 0.00 |
| `sfc_psih` | value | 4.234e-06 | 1.017e+00 | 4.16e-06 | rel 6e-05 | 0.00 |
| `sfc_hfx` | value | 2.830e-04 | 7.067e+01 | 4.00e-06 | rel 0.0003 | 0.00 |
| `pbl_hpbl` | value | 3.317e-03 | 8.588e+02 | 3.86e-06 | rel 2e-05 | 0.00 |
| `dthdt_sw` | Linf_error | 1.182e-10 | 5.915e-05 | 2.00e-06 | abs 1e-09 | 0.00 |
| `gsw` | value | 4.557e-04 | 3.151e+02 | 1.45e-06 | rel 1e-05 | 0.00 |
| `dtheta_m_dt_dyn` | Linf_error | 2.711e-10 | 1.940e-04 | 1.40e-06 | abs 2e-09 | 0.00 |
| `sfc_flhc` | value | 5.071e-05 | 4.644e+01 | 1.09e-06 | rel 3e-05 | 0.00 |
| `dw_dt_nonpg` | Linf_error | 3.381e-11 | 4.030e-05 | 8.39e-07 | abs 3e-10 | 0.00 |
| `sfc_fh` | value | 4.291e-06 | 5.226e+00 | 8.21e-07 | rel 3e-05 | 0.00 |
| `sfc_lh` | value | 1.404e-04 | 1.988e+02 | 7.06e-07 | rel 2e-05 | 0.00 |
| `sfc_qfx` | value | 5.347e-11 | 7.952e-05 | 6.72e-07 | rel 2e-05 | 0.00 |
| `F_up` | Linf_error | 2.299e-04 | 3.709e+02 | 6.20e-07 | abs 0.002 | 0.00 |
| `sfc_flqc` | value | 7.591e-09 | 1.282e-02 | 5.92e-07 | rel 2e-05 | 0.00 |
| `du_dt_dyn` | Linf_error | 1.176e-10 | 2.110e-04 | 5.57e-07 | abs 1e-09 | 0.00 |
| `sfc_q2` | value | 2.130e-09 | 4.156e-03 | 5.12e-07 | rel 6e-06 | 0.00 |
| `dv_dt_dyn` | Linf_error | 2.491e-11 | 4.874e-05 | 5.11e-07 | abs 2e-10 | 0.00 |
| `F_dn` | Linf_error | 1.325e-04 | 2.928e+02 | 4.52e-07 | abs 0.001 | 0.00 |
| `sfc_fm` | value | 2.341e-06 | 5.682e+00 | 4.12e-07 | rel 2e-05 | 0.00 |
| `sfc_ust` | value | 1.012e-07 | 5.065e-01 | 2.00e-07 | rel 7e-06 | 0.00 |
| `olr` | value | 2.535e-05 | 1.297e+02 | 1.95e-07 | rel 2e-06 | 0.00 |
| `glw` | value | 5.332e-05 | 2.928e+02 | 1.82e-07 | rel 2e-06 | 0.00 |
| `sfc_t2` | value | 4.000e-05 | 2.847e+02 | 1.41e-07 | rel 1e-06 | 0.00 |
| `sfc_u10` | value | 2.360e-07 | 2.160e+00 | 1.09e-07 | rel 1e-06 | 0.00 |
| `sfc_th2` | value | 2.238e-05 | 2.870e+02 | 7.80e-08 | rel 1e-06 | 0.00 |
| `sfc_wspd` | value | 2.170e-07 | 7.182e+00 | 3.02e-08 | rel 1e-06 | 0.00 |
| `lsm_capg` | value | 5.285e-02 | 2.365e+06 | 2.24e-08 | rel 1e-06 | 0.00 |
| `sfc_v10` | value | 1.296e-07 | 5.840e+00 | 2.22e-08 | rel 1e-06 | 0.00 |
| `sfc_gz1oz0` | value | 6.176e-08 | 6.243e+00 | 9.89e-09 | rel 1e-06 | 0.00 |
| `sfc_qsfc` | value | **0 (bit-exact)** | 9.434e-03 | 0 | abs 0 | — |
| `sfc_znt` | value | **0 (bit-exact)** | 5.000e-02 | 0 | abs 0 | — |
| `sfc_regime` | value | **0 (bit-exact)** | 4.000e+00 | 0 | abs 0 | — |
| `pbl_kpbl` | value | **0 (bit-exact)** | 1.500e+01 | 0 | abs 0 | — |
| `pbl_dqcdt` | Linf_error | **0 (bit-exact)** | 0.000e+00 | 0 | abs 0 | — |
| `mp_dqc_dt` | Linf_error | **0 (bit-exact)** | 0.000e+00 | 0 | abs 0 | — |
| `mp_dqr_dt` | Linf_error | **0 (bit-exact)** | 0.000e+00 | 0 | abs 0 | — |
| `mp_dqs_dt` | Linf_error | **0 (bit-exact)** | 0.000e+00 | 0 | abs 0 | — |
| `mp_dqg_dt` | Linf_error | **0 (bit-exact)** | 0.000e+00 | 0 | abs 0 | — |
| `mp_melt_dqi_dt` | Linf_error | **0 (bit-exact)** | 0.000e+00 | 0 | abs 0 | — |
| `mp_rate_dqi_dt` | Linf_error | 2.154e-08 | 0.000e+00 | nan | abs 2.5e-08 | — |
| `mp_sed_dqi_dt` | Linf_error | 5.304e-08 | 0.000e+00 | nan | abs 6e-08 | — |
| `lsm_land` | value | **0 (bit-exact)** | 1.000e+00 | 0 | abs 0 | — |
