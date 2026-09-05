# EqWeFiC implementation plan

Written 2026-09-04 after a de-risking review of the CLAUDE.md plan against the
current EarthSciAST spec (esm 1.0.0, EarthSciAST @ a1dc9bb), EarthSciModels,
EarthSciDiscretizations, WRF v4.8.0, and WRF-SFIRE (W4.4-S0.1 base).
Sibling checkouts live in `../`. Status tracking is in `INVENTORY.md`.

## 0. What the de-risking established

Facts that shape the plan (each was checked, not assumed):

1. **Derivative tests are expressible today.** esm §6.6 tests are simulation
   assertions `(variable, time, expected)`, but an *observed* (algebraic)
   variable asserted at `time: 0.0` evaluates the right-hand side at the given
   state. Every existing EarthSciModels component (e.g.
   `atmospheric_dynamics/holtslag_boville/*.esm`) uses exactly this form. So
   "instantaneous derivative" tests need no new format feature: tendencies are
   observeds.
2. **Column physics is expressible.** Variables carry `shape: ["lev"]` over an
   `index_sets` axis; PDE tests select a scalar with `coords` (index-space) or
   `reduce`. `aggregate` covers sums, prefix sums, and (§4.3.1.1) causal
   recurrences along one axis (needed for two-stream up/down sweeps). The Rust
   CLI runs recurrence fixtures and a latlon3d column-advection problem
   (23 s wall for the latter).
3. **`integral` is an open-tier op with no lowering anywhere yet** (esm §4.2;
   EarthSciDiscretizations tracks it as esd-yfj). A column integral written as
   `integral` will load but not simulate until a rule exists. That rule is the
   first EarthSciDiscretizations PR this project needs (Phase 0).
4. **Vertical diffusion already has rules.** EarthSciDiscretizations
   `grids/latlon3d/rules/varcoeff_laplacian_lev_{noflux,robin_surface}_bc.esm`
   is the conservative `d/dz(K du/dz)` operator with a surface-exchange BC;
   `grids/cartesian_nonuniform_1d/` is a 1-D non-uniform cell grid with
   consumer-supplied edge coordinates. A single-column WRF grid can be built
   from the nonuniform-1d pattern (Phase 0).
5. **WRF 4.8 ships several schemes in CCPP-style, self-contained F90**
   (`phys/physics_mmm/`, shared with MPAS): `bl_ysu`, `bl_mynn`, `bl_shinhong`,
   `bl_gwdo`, `sf_sfclayrev`, `sf_mynn`, `mp_wsm6`, `cu_ntiedtke`. Each has a
   `*_run(...)` entry taking plain arrays and scalar constants and
   `use`s only `ccpp_kind_types`. These are ideal kernel-harness targets:
   no WRF build, no netCDF, single `gfortran` compile, and they can be built
   in double precision by setting `kind_phys`.
6. **WRF-SFIRE already has a standalone driver** (`standalone/fire_ros.exe`,
   `fire.exe`) that runs the fuel/ROS subsystem and the full fire model
   without WRF. That is the fire harness.
7. **WRF physics call order** (dyn_em/module_first_rk_step_part1.F, solve_em.F):
   `radiation_driver → surface_driver → pbl_driver → cumulus_driver →
   shallowcu_driver` inside RK step 1, then dynamics, then
   `microphysics_driver` after the RK loop. Tendencies (`RTHBLTEN`,
   `RQVBLTEN`, `RTHRATEN`, `RTHCUTEN`, …) are Registry state with io `r`
   (restart only); adding `h` puts them in history output — that is the whole
   "instrumentation" needed at the SCM/subassembly level.
8. **The WRF single-column case** (`test/em_scm_xy`) runs on a 3×3 periodic
   stencil with prescribed forcing and the default suite
   `mp=2 (Lin), ra_lw=1 (RRTM), ra_sw=1 (Dudhia), sfclay=1 (sfclayrev),
   sf_surface=2 (Noah), pbl=1 (YSU), cu=0`, 60 levels, dt = 60 s.
9. **EarthSciModels CI has no skip/xfail.** Every inline test runs on every
   push (`tools/run_esm_inline_tests.py`, Python runner). Stubs with red tests
   cannot be merged; they must stay on a branch until stage 2 makes them pass.
10. **Existing overlap in EarthSciModels**: Monin-Obukhov surface layer
    (`local_scale/surface_layer_profile.esm`), Holtslag-Boville PBL and surface
    flux, saturation vapor pressure and thermodynamics (`sp_ch1/*`), Wesely
    dry deposition, Fast-JX photolysis, Rothermel ROS, level-set fire spread,
    fuel-model lookup, midflame wind, fire heat flux. None of the WRF-Chem gas
    mechanisms (RADM2, CBMZ, MOZART, SAPRC, CB05, CRI) exist there.
11. **Toolchain on this cluster**: system `gfortran` 11.5; modules `gcc/12.4.0`,
    `gcc/13.3.0`, `openmpi/5.0.1-gcc-13.3.0`; `apptainer` 1.5.2. No netCDF
    module. The spack netCDF/OpenMPI install under `libs/spack` is stale
    (rhel7 paths, gcc 12.2.0 that no longer exists). Kernel harnesses do not
    need netCDF; SCM and 3-D WRF runs do (Phase 0 task).
12. **GitHub**: `gh` is authenticated as `ctessum-claude` with push on
    EarthSciModels, EarthSciDiscretizations, and EarthSciAST.

## 1. Scope and reference physics suite

Target one suite first; add schemes later as separate inventory rows.

| Category | WRF option | Scheme | Source | Size | Harness route |
|---|---|---|---|---|---|
| Surface layer | `sf_sfclay_physics=1` | Revised MM5 (sfclayrev) | `phys/physics_mmm/sf_sfclayrev.F90` | 1.1k | CCPP kernel |
| PBL | `bl_pbl_physics=1` | YSU | `phys/physics_mmm/bl_ysu.F90` | 1.7k | CCPP kernel |
| Microphysics | `mp_physics=6` | WSM6 | `phys/physics_mmm/mp_wsm6.F90` (+`mp_radar`) | 2.4k | CCPP kernel |
| SW radiation | `ra_sw_physics=1` | Dudhia | `phys/module_ra_sw.F` | 0.5k | legacy kernel |
| LW radiation | `ra_lw_physics=1` | RRTM | `phys/module_ra_rrtm.F` | 7.6k | legacy kernel |
| Land surface | `sf_surface_physics=1` then `2` | 5-layer slab, then Noah | `module_sf_slab.F`, `module_sf_noahlsm.F` | 0.6k / 4.8k | legacy kernel |
| Cumulus | `cu_physics=16` (optional) | New Tiedtke | `phys/physics_mmm/cu_ntiedtke.F90` | 3.6k | CCPP kernel |
| Gravity-wave drag | `gwd_opt=1` (optional) | GWDO | `phys/physics_mmm/bl_gwdo.F90` | 0.7k | CCPP kernel |

Rationale: this is the SCM default suite except WSM6 replaces Lin (WSM6 is
in the portable set and is the more common choice). RRTMG (27k lines,
k-distribution tables) and Noah-MP (25k lines) are deferred; both are
"more of the same" once RRTM and Noah are done, and they are separate
inventory rows. MYNN is deferred for the same reason.

WRF-Chem (EqAtmChem): `chem_opt=1` (RADM2 gas phase, no aerosol) first, then
`chem_opt=300` (GOCART simple aerosols). Components: anthropogenic emissions
(`emissions_driver.F`), biogenic emissions, dry deposition
(`dry_dep_driver.F`, Wesely: reuse EarthSciModels), photolysis
(`module_phot_fastj.F` / TUV: reuse Fast-JX where it matches), gas mechanism
(`chem/KPP/mechanisms/radm2` → `reaction_systems`), wet scavenging.

Fire (EqAtmFire): the fire code is NCAR's **Community Fire Behavior Model**
(`NCAR/fire_behavior`, the `phys/fire_behavior` submodule of WRF 4.8, pushed
2026-09-04), not openwfm/WRF-SFIRE (WRF 4.4 base, last commit 2026-05-01).
Decided 2026-09-04: it is the more recent and actively maintained line, it is
the fire module already wired into the WRF version we instrument, and it has
a standalone driver (`driver/fire_behavior.F90`) with CI tests. Its physics
(`physics/`, 4k lines): Anderson fuel categories (`fuel_anderson_mod`),
fuel-moisture model (`fmc_wrffire_mod`), Rothermel ROS with wind/slope limits
(`ros_wrffire_mod`), level-set front propagation (`level_set_mod`), and the
fire→atmosphere flux in `fire_model_mod`/`fire_driver_mod`. Most of these
exist in EarthSciModels `wildland_fire/`; the fire work is mainly *adding
fire_behavior-derived tests* to them plus the fuel-moisture model and the
coupling file. WRF-SFIRE-only physics (Balbi ROS, firebrand spotting) are
optional later rows.

## 2. Phase 0 — foundations (before any component)

0.1 **Repo scaffolding (this repo).** `INVENTORY.md`, `FORTRAN_BUGS.md`,
    `lib/wrf_constants.esm` (transcribed from `share/module_model_constants.F`
    with the exact WRF values: `g=9.81`, `r_d=287`, `cp=7 r_d/2`, `r_v=461.6`,
    `xlv=2.5e6`, `svp1..3`, `ep_1/ep_2`, `karman=0.4`, `p1000mb`, `t0`, …),
    `lib/wrf_thermo.esm` (Exner function, potential/virtual temperature,
    saturation mixing ratio, moist static quantities, factored as templates
    and reused by every scheme), `.gitignore` for the `esm` binary.
0.2 **Column grid.** Add to EarthSciDiscretizations (PR) a
    `grids/column_nonuniform_1d/` (or reuse `cartesian_nonuniform_1d` under a
    `lev` rename) with: cell axis `lev` (N, metaparameter `NLEV`), edge axis
    `lev_nodes` (N+1), consumer-supplied edge heights `z_edge` and mid-level
    `z`, `dz`, plus rules: `varcoeff_laplacian_lev_{noflux,robin_surface}_bc`,
    `upwind1_flux_D_lev` (sedimentation, downward only), and the first
    `integral` lowering (whole-column and cumulative, mass-weighted). Verify
    with MMS problems, per that repo's AGENTS.md.
0.3 **Fortran harness repo.** Forks under the `ctessum-claude` account
    (done 2026-09-04): `ctessum-claude/WRF` with branch
    `earthsciml-instrumented` (pushed; `../WRF` has it checked out with remote
    `fork`) and `ctessum-claude/fire_behavior`. Directory `kernels/` on that
    branch holds: one small driver per scheme that reads a JSON/namelist column
    state, calls `<scheme>_run` once (or the legacy entry
    `ysu/sfclayrev/mp_gt_driver/…` with 1×1 horizontal extent), and writes
    inputs and outputs as JSON. Build with `gfortran` only, `kind_phys=real64`
    where the module allows. Sub-process rates (autoconversion, accretion,
    evaporation, K-profiles, φ functions, per-band fluxes) are exposed by
    adding `intent(out)` diagnostic arrays or a `debug_dump` module rather
    than `print` statements, so the dump is structured.
0.4 **Extraction tool (this repo, `tools/`).** A small script that reads a
    kernel JSON dump plus a hand-written test skeleton and fills in
    `parameter_overrides` and `assertions` with the selected regimes. It never
    writes equations.
0.5 **netCDF toolchain** (needed from Phase 1.4 on): apptainer (decided
    2026-09-04). Unprivileged image *builds* fail on this cluster (no
    subuid/subgid for fakeroot), but *pulls* work, so the approach is to pull
    a prebuilt WRF toolchain image (DTC `dtcenter/wps_wrf`) and compile the
    fork inside it with the source bind-mounted; image and cache live under
    `data/eqwefic/apptainer/`. Pulled 2026-09-04: `wps_wrf_latest.sif` (1.0 GB;
    gfortran 8 via devtoolset-8, netCDF-Fortran 4.4.6, MPI, a prebuilt WRF 4.3
    for reference). Next: build the 4.8 fork inside it with `../WRF` bind-mounted,
    then `em_scm_xy` (serial) and confirm the SCM case runs.
0.6 **Three de-risking spikes**, each ending in a passing `./esm test` file
    on this repo's `main` (not yet a PR):
    - **Spike A — YSU as a PDE.** Column diffusion `∂θ/∂t = ∂/∂z(K ∂θ/∂z) −
      ∂/∂z(K γ)` with the YSU K-profile and nonlocal term, lowered by the 0.2
      rule, tested against `bl_ysu_run` tendencies at `dt = 0.01 s`. Proves the
      "implicit Fortran vs explicit derivative" tolerance story.
    - **Spike B — Dudhia SW as an integral.** Column optical depth via
      `integral` (cumulative) and the surface flux; tested against
      `module_ra_sw.F`. Proves the `integral` lowering end to end.
    - **Spike C — WSM6 warm-rain process rates.** `praut`, `pracw`, `prevp`
      as pointwise observeds with WRF constants; tested against instrumented
      WSM6. Proves the sub-process factoring pattern and the constants
      library.
    Any spike that fails on a format gap becomes an EarthSciAST/EarthSciDiscretizations
    issue before Phase 1 starts.

## 3. Phase 1 — instrument and stub

For each scheme in the suite (order: sfclayrev → YSU → slab → WSM6 → Dudhia
SW → RRTM LW → Noah → Tiedtke → GWDO):

1.1 Write the kernel driver (0.3) and dump inputs/outputs for a set of
    regimes chosen to cover every branch: stable / unstable / neutral surface
    layer; night / day; clear / cloudy; land / water; PBL top inside / above
    cloud; mixed-phase temperatures; etc. Record regime, build precision, and
    WRF commit in the dump metadata.
1.2 Decompose the scheme into sub-components at the granularity of a
    named physical process or a textbook equation (φ-functions, bulk Richardson
    PBL-height search, K-profile, countergradient term, entrainment flux; per
    hydrometeor process rates; per-band optical depth; soil heat diffusion,
    Penman potential evaporation, canopy resistance, …). Write the stub `.esm`
    for each (variables with units, description, references, empty
    `equations`, full `tests`) and the assembly stub that mounts them as
    subsystems.
1.3 For components that already exist in EarthSciModels, add the WRF-derived
    tests to the existing file and open that PR immediately (these are green
    by construction if the physics matches; if not, that is a finding for
    `FORTRAN_BUGS.md` or the EarthSciModels file).
1.4 Subassembly data: with the SCM build (0.5) and history output of the
    `R*TEN` tendencies, run `em_scm_xy` with the suite and dump per-time-step
    physics inputs and tendencies for the surface+PBL, radiation, and
    full-physics subassemblies. Select a handful of time steps as the
    subassembly tests.
1.5 Stubs live on this repo's `main` under `components/<domain>/…` mirroring
    the EarthSciModels layout. They are *not* opened as EarthSciModels PRs
    until their tests pass (finding 9).

## 4. Phase 2 — physics and subassemblies

2.1 Fill in each stub so `./esm test` passes, then run the EarthSciModels
    Python gate, then open the PR (one PR per scheme, sub-components included,
    with the WRF commit and kernel-dump provenance in `metadata`).
2.2 Wire subassemblies by reference (`subsystems: {ref: ...}` and
    `coupling` blocks): surface-layer + land + PBL column; SW + LW radiation
    column; microphysics column; then the full SCM physics suite with the
    WRF call order. Test against the 1.4 SCM dumps.
2.3 Each PR is human-reviewed and merged before the next subassembly
    depends on it.

## 5. Phase 3 — top-level models

3.1 **EqWeather.esm**: SCM mode first (physics suite + prescribed forcing,
    compared with `em_scm_xy` history over 24 h; acceptance: tendencies within
    1e-3 relative RMS, state within 0.5 K / 5 % qv at 24 h). Then 3-D
    idealized cases (`em_hill2d_x`, `em_quarter_ss`).

    **Dynamics split (decided 2026-09-04).** "Dynamical core" means the part of
    WRF that is not a physics parameterization: the governing equations for
    wind, pressure, and temperature plus the numerics that step them. Per the
    project rule, these are separated:
    - The *equations* — WRF-ARW's flux-form compressible non-hydrostatic Euler
      equations in the terrain-following hybrid sigma-pressure (dry-mass)
      vertical coordinate, with the perturbation base-state split and the
      moist thermodynamic closure (Skamarock et al. 2021, Ch. 2) — become an
      EarthSciModels component `atmospheric_dynamics/wrf_arw/` written as
      PDEs with `D` wrt `x`, `y`, `eta`, `t`, factored into mass, momentum,
      thermodynamic, geopotential, and diagnostic-pressure sub-components (the
      existing `clark1977/` anelastic set is the layout precedent).
    - The *discretization* — Arakawa-C staggering on the mass coordinate,
      5th-order horizontal / 3rd-order vertical upwind advection with
      positive-definite limiting, divergence damping, and the vertical
      coordinate metric — becomes EarthSciDiscretizations grid + rule files
      (`grids/wrf_arw_c/`), imported by EqWeather at the coupling edge.
    - The *time integration* (RK3 with acoustic sub-stepping) is the solver's
      job and is not written in `.esm`; the runner picks the integrator. This
      is why 3-D agreement is statistical (domain-mean profiles, precipitation
      totals, front position), not bitwise. The SCM comparison has no
      dynamics and is where tight agreement is required.
    `../simpleclimate.esm` (Held-Suarez, hydrostatic primitive equations on
    latlon3d with PPM rules) is the working precedent for "PDEs in
    EarthSciModels, numerics in EarthSciDiscretizations", but its equations
    are not WRF's, so it is reused only as a pattern, not as EqWeather's
    dynamics.
3.2 **EqAtmChem.esm**: EqWeather + emissions + RADM2 + deposition +
    photolysis, compared with a WRF-Chem `chem_opt=1` SCM-like run.
3.3 **EqAtmFire.esm**: EqWeather + SFIRE fire components via the
    `wildlandfire.esm`-style coupling, compared with `test/em_fire/hill`.

## 6. Test and authoring conventions (normative; summarized in CLAUDE.md)

- Tendencies are observeds; assert at `time: 0.0`. Trajectory assertions only
  for genuinely in-place steps.
- Implicit-solve schemes: reference tendency at small `dt`; document it.
- Column components shaped over `lev`; geometry supplied by the consumer.
- Constants only from `lib/wrf_constants.esm`.
- Tolerances: `rel 1e-9` for real64 kernel references, `rel 1e-5` for real32.
- Every test `description` states regime, WRF commit, kernel build precision,
  and which Fortran routine produced the values.
- Metadata `references` cite the scheme paper and the Fortran file:line.

## 7. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Implicit Fortran tendencies do not equal instantaneous derivatives | Small-`dt` references (Spike A); state tolerance explicitly; a one-step trajectory test at WRF `dt` as a secondary check |
| `integral` has no lowering | Phase 0.2 rule in EarthSciDiscretizations; fall back to `aggregate` prefix sums only if the rule PR stalls (and record the deviation from CLAUDE.md) |
| Radiation (RRTM) is large and table-driven | Absorption data as `function_tables`/`data_sources`; two-stream sweeps as recurrences; do Dudhia first |
| Rust CLI cannot integrate stiff column models | Most tests are algebraic at t=0 so stiffness is irrelevant; for trajectory tests use the Python gate (LSODA) as the reference runner |
| Physics in WRF-SFIRE (4.4 base) differs from WRF 4.8 | Fire components are tested against SFIRE's standalone driver; atmosphere physics against WRF 4.8; document the version split in EqAtmFire |
| Stubs with failing tests break EarthSciModels CI | Keep stubs here until green (finding 9) |
| Precision mismatch (WRF real32) | Build kernels real64 where possible; otherwise `rel 1e-5` |
| Volume of dumps | Disk-backed `data/eqwefic/`; never in git or /tmp |

## 8. Decisions (answered 2026-09-04)

1. Physics suite: the §1 table. RRTMG, Noah-MP, MYNN, Thompson are later rows.
2. Instrumented Fortran: forks under `ctessum-claude` (`WRF`, `fire_behavior`).
3. Dynamics: PDEs in EarthSciModels, discretization rules in
   EarthSciDiscretizations, time integration in the runner (see §5, 3.1).
4. WRF-Chem: RADM2 gas phase (`chem_opt=1`) first, then GOCART; scripts may
   translate KPP `.eqn` tables into `reaction_systems`.
5. Fire code: NCAR `fire_behavior` (Community Fire Behavior Model), the more
   recent and actively developed line; WRF-SFIRE extras are optional rows.
6. netCDF/WRF builds: apptainer with a pulled toolchain image (see 0.5).
