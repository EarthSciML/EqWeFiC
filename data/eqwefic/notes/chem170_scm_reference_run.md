# chem_opt = 170 SCM reference run (stage 1), 2026-09-22

Stage-1 reference data for **EqAtmChem at `chem_opt = 170`** (CBM-Z gas phase through
KPP + 8-bin MOSAIC aerosol), run on the idealized `em_scm_xy` column together with
**Noah** (`sf_surface_physics = 2`) and **New Tiedtke** (`cu_physics = 16`), both of
which became required schemes on 2026-09-22. GWDO does nothing on this column; it is
covered by the real-terrain run in `data/eqwefic/notes/gwdo_real_reference_run.md`,
which also carries real-terrain Tiedtke and Noah dumps.

`data/eqwefic/...` below means the disk-backed data tree
`/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/`, not this repo.

## Where everything is

| what | where |
|---|---|
| WRF source | fork `ctessum-claude/WRF`, branch `earthsciml-instrumented-chem170` @ `eae782d` (off `earthsciml-instrumented-chem` @ cc541d4); worktree `data/eqwefic/wrf-chem170` |
| build | `WRF_CHEM=1 WRF_KPP=1 ./compile -j 8 em_scm_xy` inside the existing `apptainer/wrf-chem-build.sif` (recipe unchanged — every KPP mechanism, `cbmz_mosaic` included, is already compiled by that build), 11–16 min per pass |
| run directory | `data/eqwefic/chem170_scm` (sounding, soil, tables, namelist) and its `main/` (`wrfout_d01_1997-06-20_12:00:00`, logs) |
| launcher | `data/eqwefic/chem170_scm/run.sh` (`MODE=probe` for the step-selection pass, `MODE=main` for the dumped run) |
| dumps | `data/eqwefic/dumps/chem170_wrf/esm_dump_<scheme>_<step>.json`, 204 files, 135 MB |
| kernel replays | `data/eqwefic/wrf-chem170/kernels/build/{ntiedtke,noah}_driver` (real64); scratch and `replay_summary.txt` in `chem170_scm/replay` |

## The case

Weisman & Klemp (1982) analytic sounding (θ₀ = 300 K, θ_tr = 343 K at z_tr = 12 km,
RH = 1 − 0.75 (z/z_tr)^1.25 capped at q_v = 14 g/kg, westerly shear 0 → 5 m/s over the
lowest 3 km, p_sfc = 1000 hPa), Kansas (37.6 N, 96.7 W), **1997-06-20 12 UTC + 48 h**,
`dt` = 60 s (2880 steps), 59 layers to `ztop` = 20 km, `dx = dy` = 20 km. Noah soil:
silt loam (`scm_isltyp = 4`), dryland cropland (`scm_lu_index = 2`), `vegfra` 50 %,
θ_soil = 0.30 through the column, T_skin = 300 K. `scm_force = 0`; the retimed
`force_ideal.nc` supplies U_g = 5 m/s, V_g = 0 and no advective forcing.

Why this sounding: the CASES-99 profile of the existing `chem_scm` run is stable and
dry (q_v = 2.5 g/kg), so New Tiedtke never fires there. Why `dx` = 20 km: the scheme's
scale-dependency factor switches at `dxref` = 15 km (`cu_ntiedtke.F90:232`), and 20 km
is on the coarse-grid branch the scheme was tuned on, which is also where a cumulus
parameterization belongs.

Chemistry options: `chem_opt = 170`, `phot_opt = 1` (Madronich), `chemdt = photdt` =
1 min, `gaschem_onoff = aerchem_onoff = vertmix_onoff = 1`, `gas_drydep_opt = 1`
(Wesely) and `aer_drydep_opt = 1` (MOSAIC), `gas_ic_opt = aer_ic_opt = 1`,
`chem_in_opt = 0`, `emiss_opt = 0`, `wetscav_onoff = cldchem_onoff = 0`,
`aer_ra_feedback = 0`, `chem_conv_tr = 0`.

- **Emissions are off**, as in the RADM2 run. The default profiles plus 48 h of
  photochemistry exercise every rate; emissions would need an inventory this
  idealized column has no geography for.
- **Wet scavenging and cloud chemistry are off**: `chem_opt = 170` is the non-aqueous
  MOSAIC package (`cbmz_mosaic_kpp`), and those are separate `_aq` options.
- **Convective tracer transport is off.** WRF-Chem's `chem_conv_tr` is wired to the
  Grell cumulus schemes, not to New Tiedtke, so with `cu_physics = 16` there is no
  convective chemistry transport to dump.

## Coverage actually achieved

16 steps are dumped: `1, 60, 230, 364, 601, 639, 689, 918, 1080, 1621, 1730, 1801,
2101, 2134, 2400, 2521`. They were chosen from a probe pass of the same run that dumped
New Tiedtke at **every** step (1772 of 2880 steps have non-zero convective tendencies),
then classified by replaying the kernel in real64, which reports `ktype`. Verified on
the dumps themselves (`chem170_scm/replay/replay_summary.txt`):

| step | local | Tiedtke `ktype` | cloud base → top (kernel levels) | max \|rthcuten\| (K/s) | conv. precip over dt (mm) | MOSAIC nucleation |
|---|---|---|---|---|---|---|
| 1 | 06 LST | 1 deep | 44 → 5 | 2.6e-2 | 0.489 | active |
| 60 | 07 LST | 0 none | — | 0 | 0 | — |
| 230 | 09 LST | 0 none | — | 0 | 0 | active |
| 364 | 12 LST | 1 deep | 46 → 5 | 8.6e-5 | 1.5e-4 | — |
| 601 | 16 LST | 1 deep | 46 → 6 | 3.7e-3 | 0.112 | — |
| 639 | 16 LST | 2 shallow | 45 → 32 | 2.4e-4 | 2.5e-3 | active |
| 689 | 17 LST | 2 shallow | 45 → 32 | 2.5e-4 | 1.6e-3 | — |
| 918 | 21 LST | 3 mid-level | 10 → 10 | 9.1e-6 | 0 | — |
| 1080 | 00 LST | 3 mid-level | 13 → 13 | 0 | 0 | — |
| 1621 | 09 LST | 3 mid-level | 14 → 14 | 1.0e-6 | 0 | — |
| 1730 | 10 LST | 1 deep | 48 → 10 | 8.4e-5 | 5.5e-4 | active |
| 1801 | 12 LST | 1 deep | 48 → 9 | 8.2e-5 | 3.9e-4 | — |
| 2101 | 17 LST | 1 deep | 46 → 10 | 1.1e-3 | 2.5e-2 | — |
| 2134 | 17 LST | 2 shallow | 46 → 30 | 1.1e-3 | 1.9e-3 | — |
| 2400 | 22 LST | 3 mid-level | 10 → 10 | 0 | 0 | — |
| 2521 | 00 LST | 3 mid-level | 13 → 13 | 0 | 0 | — |

So all three New Tiedtke convection types fire, plus the inactive branch; Noah runs
day and night, dry and precipitating; and the MOSAIC gas-particle, coagulation and
nucleation sub-processes are all exercised (coagulation changes `rsub` at every dumped
step, nucleation only in bursts — steps 1, 230, 639 and 1730 of those dumped).

Madronich (`photmad`) has 12 dumps, not 16: `photolysis_driver` is not called at the
four night steps (918, 1080, 2400 and 2521). The corresponding
`cbmz_kpp` dumps have `jv = 0` throughout, so the dark chemistry is still covered.

## Instrumentation (branch `earthsciml-instrumented-chem170` @ `eae782d`)

Five new dump schemes, all driven by `ESM_DUMP` / `ESM_DUMP_CALLS` / `ESM_DUMP_DIR`:

| scheme | where | what |
|---|---|---|
| `ntiedtke` | `phys/module_cu_ntiedtke.F`, row `j = jts` | driver inputs (`dz8w`, `pi3d`, `pcps`, `p8w`, `w`, `t3d`, `qv3d`, `qc3d`, `qi3d`, `rho3d`, `u3d`, `v3d`, `qvften`, `thften`, `dx`, `hfx`, `qfx`, `xland`, `dt`, `stepcu`, `itimestep`, the seven constants); the kernel-level state `cu_ntiedtke_pre_run` builds in the scheme's own **top-down** ordering (`pu`, `pv`, `pt`, `pqv`, `pqc`, `pqi`, `pqvf`, `ptf`, `poz`, `pzz`, `pomg`, `pap`, `paph`, `slimsk`, `delt`); the updated kernel state, `zprecc`, and the `rthcuten`/`rqvcuten`/`rqccuten`/`rqicuten`/`rucuten`/`rvcuten` tendencies |
| `noah` | `phys/module_sf_noahdrv.F` around `SFLX`, point `(its, jts)` | every `SFLX` input, the state it advances (`cmc`, `t1`, `stc`, `smc`, `sh2o`, `snowh`, `sneqv`) and all its fluxes and diagnostics, **plus the VEGPARM/SOILPARM/GENPARM tables `REDPRM` reads** so the call replays without a table file |
| `cbmz_kpp` | `chem/KPP/inc/cbmz_mosaic/kpp_mechd_{u,l,b,ib,ia,a}_cbmz_mosaic.inc` | the same six KPP hook points used for RADM2: per level `var`, `fix`, `RCONST` after `Update_Rconst`, `jv`, `TEMP`, `C_M`, `C_H2O`, `p`, `rho`, `qv` before `INTEGRATE`, and `var` after (66 species, 142 reactions, 52 photolysis rates) |
| `mosaic` | `chem/module_mosaic_driver.F`, column `(its, jts)` | the WRF `chem` and `moist` columns in and out, the MOSAIC column state `rsub` with its 235 species names, and `rsub` **after each stage separately**: `aerchemistry` (gas–particle transfer, water uptake, section moving), `mosaic_newnuc_1clm`, `mosaic_coag_1clm` |
| `mosaic_drydep` | `chem/module_mosaic_drydep.F` | surface-layer inputs (`u*`, aerodynamic resistance, T, p, ρ, the surface `rsub` row), the section geometry and composition pointers, the per-bin deposition velocities and the `ddvel` row handed to vertical mixing |

The existing `photmad`, `wesely` and `chem_vertmx` hooks work unchanged at
`chem_opt = 170`, and the `wesely` dump records `chem_opt`. The commit also carries the
`chem_vertmx` and `ESM_SWEEP` hooks that had been left uncommitted in the
`earthsciml-instrumented-chem` worktree.

Two real64 kernel drivers, both built by `kernels/Makefile`:

- `ntiedtke_driver` replays `cu_ntiedtke_run` and re-forms the WRF tendencies exactly as
  `cu_ntiedtke_post_run` does. The Makefile compiles the kernel from a **build-time
  copy** of `physics_mmm/cu_ntiedtke.F90` with `ntiedtke_esm_inner.inc` spliced in after
  the `cumastrn` call, so `ktype`, cloud base and top, and `cumastrn`'s own tendencies
  and mass fluxes come out too. The submodule source is untouched.
- `noah_driver` replays `SFLX` with the dumped parameter tables;
  `kernels/stub_wrf_error.F90` stands in for `frame/module_wrf_error`.

Replay agreement against the in-model real32 dumps (all 16 steps, `replay_summary.txt`):

- **Noah**: `t1` within 3.3e-5 K, sensible heat flux within 4.1e-4 W m⁻² of 130 W m⁻²
  (3e-6 relative) — i.e. at real32 round-off. The driver is validated.
- **New Tiedtke**: `|r64 − r32|` ≤ 8.9e-7 K/s on `rthcuten`. That is **not** a
  transcription gap but cancellation: WRF forms the tendency as
  `(T_out − T_in)/π/dt` in real32, and one ulp of a ~300 K temperature over dt = 60 s
  is already 2.4e-7 K/s (FORTRAN_BUGS N84). **Use the real64 replay, not the in-model
  `rthcuten`, as the reference**, and keep the in-model dump only as a cross-check.

## Tendency form

- **New Tiedtke** updates the state in place as `x + (dx/dt)·dt` and the driver divides
  the increment by `dt` again, so `rthcuten` is an explicit-Euler rate evaluated at the
  input state. It is not exactly dt-independent (the deep closure carries a CAPE
  adjustment time scale and `nonequil`), so the driver honours `ESM_DT`.
- **Noah `SFLX`** is an in-place step: `(x_out − x_in)/dt` converges to the
  instantaneous rate as dt shrinks, and `ESM_DT` overrides dt in the replay.
- **MOSAIC** is operator-split by construction: each dumped stage boundary is an
  increment over `dtchem` = 60 s, not a rate. The stage dumps are what makes each
  sub-process separable.
- **CBM-Z** is the KPP Rosenbrock integration over `dtstepc` = 60 s, exactly like
  RADM2: `var_in`/`var_out` are a trajectory, and `RCONST` with `var_in` gives the
  instantaneous rates.

## Things found on the way

- **B18: MOSAIC initial aerosol was identically zero.** `module_mosaic_initmixrats.F`
  hard-wires `CASENAME 4`, and the initializer body is `#if`'d to CASENAME 0–3, so with
  `chem_in_opt = 0` every bin starts at exactly 0 and nothing warns. The branch admits
  CASENAME 4, so `aer_ic_opt = 1` gives `mosaic_init_wrf_mixrats_opt1`'s four lognormal
  modes (accumulation so4/nh4/oc/bc, Aitken, coarse dust, coarse sea salt). All eight
  bins are then populated: number 9.2e5 → 12 cm⁻³ and so4 6.9e-13 → 1.9e-13 at level 1.
- **N83: a hard-coded component index in that initializer.** The coarse dust mode sets
  its hysteresis water from `aaprof_nsm(9,nsm)`, meaning `oin`, but the component order
  for `cbmz_mosaic` is so4, no3, cl, co3, nh4, na, ca, oin, oc, bc — index 9 is `oc`,
  which is zero in that mode, so the dust mode starts with no hysteresis water.
- **N84**: the real32 cancellation in `rthcuten` described above.
- **The chemistry changes the meteorology — resolved 2026-09-22, and it is legitimate
  physics.** Running the identical case with `chem_opt = 0` is bitwise identical at step 1
  but diverges by step 60 (0.04 K, 0.17 W m⁻² in `hfx`) and by ~2 K within hours, enough
  to change which convection type fires. The cause is **Dudhia shortwave aerosol
  scattering**, which `aer_ra_feedback = 0` does not switch off (FORTRAN_BUGS N86): when
  `pm2_5_dry`/`pm2_5_water` are PRESENT — i.e. in any WRF-Chem run with an aerosol package
  — `module_ra_sw.F:465` adds them to the layer scattering `XSCA`, and `sum_pm_driver`
  refills them every step. See "Bisecting the divergence" below.

- The run prints `Warning: refi is larger than lookup table range ... SW band 1`
  repeatedly. It comes from the optical driver, does not touch radiation at
  `aer_ra_feedback = 0`, and is not investigated here.
- In the container, `wrf.exe` segfaults in OpenSSL cleanup *after* printing
  `SUCCESS COMPLETE WRF`; `run.sh` therefore checks the log rather than the exit code.



## Kernel replays of the chemistry (2026-09-22)

Four real64 drivers now exist beside the two weather ones; all are built by
`wrf-chem170/kernels/Makefile` (`make` for real64, `make PREC= B=build32` for real32).

| driver | what it replays | agreement with the in-model dump |
|---|---|---|
| `cbmz_driver` | the whole KPP CBM-Z interface: `Update_RCONST` from the dumped T, C_M, C_H2O and j values; `Fun` (the instantaneous production-minus-loss rates); and the Rosenbrock step | **exact (Linf = 0)** for both `RCONST` and `var_out`, at all 16 dumped steps. WRF's KPP is already double precision, so this is an independent replay of the same arithmetic rather than a precision upgrade |
| `mosaic_drydep_driver` | `mosaic_drydep_1clm` + `aerosol_depvel_2`, extracted verbatim at build time | per-bin deposition velocities within **2.3e-7 to 1.6e-6 relative** over the 16 steps, i.e. real32 round-off |
| `mosaic_subproc_driver` | `mosaic_newnuc_1clm` and `mosaic_coag_1clm`, WRF's own modules compiled standalone | **real32 build: bit for bit** (nucleation exactly 0; coagulation 1e-22..1e-20, the decimal round-off of the dump itself) at all tested steps. See the caveat below |
| `ntiedtke_driver`, `noah_driver` | as before | see above |

**What this buys stage 2.** CBM-Z can be pinned as tightly as the .esm evaluation allows,
like RADM2. The MOSAIC deposition velocities support `rel 1e-5`. And the replay settles
what the *instantaneous* rates are: at the dumped daytime steps `Fun(var_in)` reaches
1.1e6 to 3.6e7 molec cm⁻³ s⁻¹ while the finite-step increment `(var_out - var_in)/60 s`
reaches only 2.7e5 to 8.0e6 — a factor of 4 to 35. A reaction-system component asserting
instantaneous rates must be referenced against `Fun`, not against the KPP step. At night
the two agree, because the chemistry is then slow.

**The MOSAIC stage increments are cancellation-limited in the dumps** (FORTRAN_BUGS N87).
Measured as (|rsub_before| + |rsub_after|)·2⁻²⁴ / |increment| over six steps:
coagulation has a **median of 0.04-0.09** (one to two significant digits survive) with
worst entries at 2 (the increment is below the representable resolution); gas-particle
transfer 5e-5 to 3e-4; nucleation 6e-8 to 1e-4. So the stage dumps cannot be used
directly as references for coagulation.

**Open item: the real64 build of `mosaic_subproc_driver` is not yet trustworthy.** Its
real32 build reproduces WRF exactly, so the dumped state is complete and the driver's
wiring is right, but the real64 build returns *no nucleation at all* and a different
coagulation state. Every module variable was verified correct at the call site
(`rsub(ktemp)`, `rsub0`, `cairclm`, `relhumclm`, `afracsubarea`, `nsubareas` and the row
indices all print correctly immediately before the call), the symbol table shows one
shared copy of each data module, and the behaviour is identical at -O0 and -O2 and with
`-fdefault-double-8` added. Inside the routine the same reads come back zero or shifted.
That points at the `-fdefault-real-8` promotion of the MOSAIC data modules rather than at
the physics, and it must be resolved before real64 references for nucleation and
coagulation can be produced. Until then the honest stage-2 position is: deposition and
CBM-Z are ready, the two aerosol dynamics stages are not.


## MOSAIC writes its own rates now (2026-09-22, WRF fork `a40edfe`)

Differencing two real32 stage states was never going to give a usable reference for
coagulation (N87), so each sub-process now writes its **instantaneous rate where it forms
it**, at the state the stage starts from (`chem/module_esm_mosaic_rates.F`, inert unless
the driver switches it on for the dumped column). What each one records:

| process | what is recorded | is it a true rate? |
|---|---|---|
| nucleation | the projection increment itself (`rate_newnuc`, mol/mol-air per call), the same divided by `dtnuc` (`rate_newnuc_finitedt`), and the composition and thresholds B20 corrupted (`newnuc_composition`) | **no rate exists**: with B20 fixed the scheme relaxes the vapour to critical in one call, the same increment whatever `dtnuc` is |
| coagulation | the same production and loss sums `coagsolv` uses, evaluated with the **input** distribution and divided by `deltat` (`rate_coag_num`, `rate_coag_vol`) | yes |
| gas-particle transfer, non-volatile | `gas(iv)*kg(iv,ibin)` per bin — the **exact** dt → 0 limit of ASTEM's analytical exponential decay (`rate_cond_nonvol`) | yes, exactly |
| gas-particle transfer, semi-volatile | `flux_s + flux_l` on ASTEM's **first sub-step**, i.e. kg·(gas − surface value) at the input state (`rate_cond_semivol`) | yes |
| aerosol dry deposition | unchanged: the per-bin velocities, cross-checked by `mosaic_drydep_driver` | yes (a velocity, not an increment) |

### Rate against finite increment, over the 16 dumped steps

- **Coagulation.** Where the increment is resolved (more than 10× the real32 resolution)
  the rate reproduces it to **0.1 % median, 4.4 % at the 90th percentile, 16 % worst**;
  the residual is the semi-implicit solver's own nonlinearity over 60 s. At every step
  **90–113 of the 472 number entries have a genuine rate but an increment at or below the
  real32 resolution** — those bins simply cannot be referenced from the state dumps.
- **Gas-particle transfer.** H2SO4 1.00–1.01 and HNO3 1.00 — linear over the step, so
  either quantity would do. HCl 0.67–1.00. **NH3 ranges from −856 to +438**: the ammonia
  flux reverses sign inside the step, so its increment is not a tendency at all and only
  the rate can be referenced.
- **Nucleation.** The increment WRF applies is **not** the rate: sulfate differs by 1×–56×,
  number by 1.4×–79×, and ammonium by **363×–6.4e4×**. The cause is recorded as
  FORTRAN_BUGS B20 — `qh2so4_avail` is an excess rate times `dtnuc`, and the composition
  partition then divides `qnh3_cur` by that dt-scaled amount, so both the size and the
  neutralisation of the new particles move with the step. A component must say which of
  the two it transcribes; both are dumped.

### Tolerances this supports

| quantity | reference | tolerance |
|---|---|---|
| CBM-Z rate coefficients and integrated step | `cbmz_driver`, exact | as tight as the .esm arithmetic allows (RADM2 uses `rel 1e-9`) |
| CBM-Z instantaneous rates | `Fun` from the same driver | same |
| MOSAIC aerosol deposition velocities | `mosaic_drydep_driver`, real64 | `rel 1e-5` (real32 round-off, measured 2.3e-7–1.6e-6) |
| MOSAIC coagulation and condensation rates | the in-model rate hooks, real32 | `rel 1e-5`; the rates are single-precision values from the model, so the bound is WRF's own precision, not the differencing noise |
| MOSAIC nucleation | the in-model projection increment (B20-corrected) | `rel 1e-5`, and asserted as an increment per call, not a rate |
| MOSAIC stage increments | the stage states | **not usable for coagulation** (N87) and not a tendency for NH3 or nucleation |

`tools/mosaic_rate_check.py` prints the comparison for any dumped column; the run over
all 16 steps is in `chem170_scm/replay/rate_summary.txt`.

**This takes the `-fdefault-real-8` promotion off the critical path.** The rates are
computed inside the model in its own precision, so nucleation and coagulation no longer
need the real64 replay whose build was unresolved. That driver's real32 build still
reproduces WRF bit for bit and stays as a cross-check.


## The nucleation bug (B20) is fixed in the fork, and the reference recaptured

`wexler_nuc_mosaic_1box` multiplied an excess mixing ratio by the time step and then used
the product as a mixing ratio. The fix is one line, `qh2so4_avail = qh2so4_cur -
qh2so4_crit`, which is what **WRF's own newer copy of the routine already has**
(`chem/module_mosaic_newnucb.F:1668`). The other nucleation route in the same file (the
ternary/Napari one at `:612-626`, unreachable because `newnuc_method = 2`) is correct:
`ratenuclt*dtnuc*mass_part` is a genuine rate times a step, and its composition partition
divides by `qh2so4_cur`, a mixing ratio.

**The clamp analysis was verified against the stock dumps before anything was changed.**
At step 1 all 46 firing levels have `|qh2so4_del| = 0.9999*qh2so4_cur` exactly, and the
measured rate/increment ratio matches `60*(1 - qh2so4_crit/qh2so4_cur)` to 1e-4 across
5.1x-59.6x. The three later single-level events sit 0.5-0.9 % above critical — inside the
~1.7 % band where the clamp does not bind — and there the increment is `(qcur -
qcrit)*dtnuc`, the other half of the same bug.

**What the fix changes, on identical input columns** (same column through both schemes,
`kernels/mosaic_subproc_driver`; `chem170_scm/replay/b20_kernel_compare.txt`):

| quantity | stock WRF | fixed |
|---|---|---|
| H2SO4 vapour left / vapour in, step 1 | 1.0e-4 | 0.063 (= the critical concentration) |
| same, at the three near-critical events | 0.44-0.68 | 0.991-0.995 |
| (mole NH4)/(mole SO4) of the new particles | 0.0009-0.0055 | 0.054-0.33 |
| particle number made | — | -4.5 % at step 1; 20x-60x fewer at the near-critical events |
| levels the guard lets through | 46 / 1 / 0 | 46 / 1 / 0 — unchanged on identical inputs |

**The corrected scheme is exactly dt-independent.** Replayed at dtnuc = 60, 6 and 0.6 s it
returns bit-identical increments (`b20_dt_independence.txt`), which is what a projection
should do. That is why the instrumentation now records the projection increment itself and
does not manufacture a rate for it.

**It moves the meteorology**, through the N86 path (aerosol mass feeds Dudhia's
scattering): identical at step 1, 0.0013 K by step 60, and up to 1.9 K, 2.6e-3 kg/kg in
q_v and 24 W/m2 in `hfx` by step 601 (`b20_met_impact.txt`). The recaptured reference is
therefore a slightly different trajectory from the stock one; both are internally
consistent, and each dumped call is still self-contained.

**Both dump sets are kept**: `dumps/chem170_wrf` is the corrected reference,
`dumps/chem170_wrf_stockB20` the stock cross-check. Nothing else in the configuration
changed between them.

## Bisecting the chem-on / chem-off divergence

One binary (`wrf-chem170/main/wrf.exe`), one namelist switch at a time
(`chem170_scm/bisect.sh`, dumps in `dumps/chem170_bisect`), comparing the New Tiedtke
driver inputs (which carry `t3d`, `qv3d`, `qc3d`, `hfx`, `qfx`) and every other scheme's
dump at steps 1-5:

| variant | vs `chem_opt = 0` |
|---|---|
| `chem_opt = 1` (RADM2 through the hand-coded path, **no aerosol package**) | **bitwise identical** on every dumped field at steps 1-5 |
| `chem_opt = 170` with gas chemistry, aerosol chemistry, photolysis, dry deposition and vertical mixing **all off** (only the aerosol initial condition present) | differs from the first radiation call onwards |
| `chem_opt = 170` full | same, slightly larger |

So it is not the chem arrays in advection, not a heap-layout or uninitialised-memory
effect (which would not spare `chem_opt = 1`), and not compile flags (one binary).

Tracing it through the dumps: the first differing quantity is Dudhia's own output at the
model's second radiation call, with **identical inputs** (`t3d`, `qv3d`, `qc3d`, `qi3d`,
`albedo`, `coszen` all bitwise equal) — layer scattering `xsca` 1.21x the clear-sky value
at the surface rising to 3.11x at the model top, `sdown` differing by 17 W m⁻² and net
surface shortwave `gsw` 13.13 -> 8.49 W m⁻². The extra scattering is the PM2.5 term in
`XSCA`. Everything after that is the SCM responding to a different surface energy budget,
amplified by convection.

Two consequences worth stating plainly:

- **It does not threaten the stage-1/stage-2 use of this run.** Every dump is a
  self-contained kernel call: the inputs and outputs of one scheme at one step. A
  component test replays that call and never depends on the trajectory that produced it.
- **It does constrain stage 3.** This run's meteorology is aerosol-coupled, so an
  EqWeather-only model cannot reproduce its trajectory. A stage-3 weather comparison must
  either use a `chem_opt = 0`/`1` run (which are the same trajectory, as shown above) or
  mount the Dudhia PM2.5 scattering term in the coupled model.
- The initial aerosol is a **constant mass mixing ratio through the whole 20 km column**
  (`mosaic_init_wrf_mixrats_opt1`, `iiprof_nsm = 1`), so there is aerosol in the
  stratosphere and the scattering enhancement grows with height. That is an artefact of
  the initial-condition path, not of MOSAIC.

## Cost

Each 48 h dumped run is ~4.5 min and the build 11–16 min, so the chemistry was cheap;
the expensive part was rebuilding WRF-Chem three times (the KPP interface's `.inc`
dependencies are not in the makefile, so the hooks needed a forced re-`cpp`).
