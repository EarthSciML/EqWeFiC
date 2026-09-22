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
| WRF source | fork `ctessum-claude/WRF`, branch `earthsciml-instrumented-chem170` @ `ac6d4cc` (off `earthsciml-instrumented-chem` @ cc541d4); worktree `data/eqwefic/wrf-chem170` |
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

## Instrumentation (branch `earthsciml-instrumented-chem170` @ `ac6d4cc`)

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
- **The chemistry changes the meteorology.** Running the identical case with
  `chem_opt = 0` gives a bitwise-identical step 1 but diverges by step 60 (0.04 K,
  0.17 W m⁻² in `hfx`) and by hours 2.2 K, enough to move which convection type fires.
  `aer_ra_feedback = 0`, and the `mosaic` dumps show `moist` unchanged across the
  aerosol driver, so it is not aerosol–radiation or water-uptake feedback; it is most
  likely a scalar-array-dependent path in the dynamics amplified by a convecting
  column. **Consequence for this reference set: every dumped step, including the step
  selection, comes from the chemistry run itself**, and the physics-only probe was used
  only to shortlist candidates (all of which were re-classified against the chemistry
  run's own dumps).
- The run prints `Warning: refi is larger than lookup table range ... SW band 1`
  repeatedly. It comes from the optical driver, does not touch radiation at
  `aer_ra_feedback = 0`, and is not investigated here.
- In the container, `wrf.exe` segfaults in OpenSSL cleanup *after* printing
  `SUCCESS COMPLETE WRF`; `run.sh` therefore checks the log rather than the exit code.

## Cost

Each 48 h dumped run is ~4.5 min and the build 11–16 min, so the chemistry was cheap;
the expensive part was rebuilding WRF-Chem three times (the KPP interface's `.inc`
dependencies are not in the makefile, so the hooks needed a forced re-`cpp`).
