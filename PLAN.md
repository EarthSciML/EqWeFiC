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
   Decided 2026-09-05: that gate is a merge condition, not a PR-creation
   condition, and is not run locally; local testing is the Rust CLI only
   (Julia/Python only to debug a suspected Rust CLI fault). The Python
   binding's slowness on array operations is a separate performance item.
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
0.2 **Column grid.** Done: EarthSciDiscretizations PR #34
    (branch `column-nonuniform-1d`, 2026-09-05; full julia/python/rust gate
    green on all 13 new cases): a native
    `grids/column_nonuniform_1d/` (cell axis `lev`, edge axis `lev_nodes`,
    metaparameter `NLEV`, consumer-supplied interface array `ze`, grid-derived
    `zc`/`dz`; layer 1 = surface, WRF/CCPP order) rather than a rename of
    `cartesian_nonuniform_1d`, because esm-spec §9.7.7 renaming rewrites only
    `wrt`/`dim` and not an `integral` node's `var`/bounds (verified in the Rust
    CLI; filed upstream). Rules: `face_flux_D_lev_supplied_faces`
    (`D(F, lev)` on an interface flux), `varcoeff_face_laplacian_lev_flux_bc`
    (`D(K·D(u,lev),lev)` with K on interfaces exactly as YSU's `xkzh`, free names
    `kdudz_bot`/`kdudz_top` for the surface/top values of K ∂u/∂z, so a WRF
    surface flux enters as `kdudz_bot = -hfx/(rho cp)`), and five `integral`
    lowerings (`integral_lev_whole`, `integral_lev_cumulative_{from_bottom,to_top}`
    to layer centres, `integral_lev_nodes_cumulative_{from_bottom,to_top}` to
    interfaces; the bound literal `lev`/`lev_nodes` selects the form). Verified
    numerically in the Rust CLI (27 probe assertions, three MMS problems all
    green; the constant-K column problem reproduces the cartesian
    `heat_1d_nonuniform_neumann` error at N=64). Sedimentation done 2026-09-06: ESD PR #36
    (stacked on #35) `sedimentation_upwind1_flux_D_lev` — donor-cell
    `D(W·q, lev)` for a layer-centred W = −v_t ≤ 0, no inflow at the top,
    surface outflow = layer-1 flux; MMS order 0.87–0.95; rust+julia gates
    green, python left to CI.
    **Format gaps found (EarthSciAST PR #177 fixes (a)–(c) in all three
    bindings — (b)/(c) turned out to be one Python bug, a multi-model cell-name
    collision in the inline runner, not an aggregate bug; issues #174 (rename
    gap (d)), #175 (Julia gather gap (e)), #176 (Julia dead observeds)):**
    (a) shaped state-dependent observeds (e.g. a column tendency `dudt[lev]`)
    cannot be asserted in §6.6 inline tests in the Rust or Python runners
    ("array state has no cells in var_map"); only ODE states and state-free
    array observeds are supported, so the derivative-test shape of §6 needs
    this fix before Phase 1 tests can be green; (b) the Python binding drops
    the elementwise term of `aggregate + elementwise` array sums (breaks the
    centre-cumulative integral rules in Python only; manifests carry
    `blocked_upstream_bindings`); (c) Python disagrees with Rust on the
    node-cumulative forms (cause under investigation upstream); (d) the
    §9.7.7 rename gap above; (e) the Julia reference runner faults
    (`E_TREEWALK_UNBOUND_VARIABLE`) on an elementwise-defined array observed
    (`f = 1 + cos(π zc)`) that is consumed only through an aggregate gather —
    spell such integrands as explicit gathers `aggregate(i; 1 + cos(π zc[i]))`
    until fixed. With that spelling Python also passes all five integral
    forms, so (b)/(c) may be the same root cause.
0.3 **Fortran harness repo.** Forks under the `ctessum-claude` account:
    `ctessum-claude/WRF` branch `earthsciml-instrumented` (pushed; `../WRF` has
    it checked out with remote `fork`) and `ctessum-claude/fire_behavior`.
    Done 2026-09-04 (commit cf3aed7): `phys/module_esm_dump.F` writes one
    JSON file per dumped kernel call (`ESM_DUMP=<schemes|all>`,
    `ESM_DUMP_CALLS=<1-based call indices>`, `ESM_DUMP_DIR`; in the SCM every
    j-row of the tile is one call, so call 2n is step n); `module_bl_ysu.F`
    dumps every `bl_ysu_run` input/output; `kernels/` holds the flat-input
    reader and `ysu_driver` (built real64 with `-DDOUBLE_PRECISION`; `ESM_DT`
    overrides the step), which replays a dump and writes the same JSON.
    Verified: at dt = 60 s the real64 driver reproduces WRF's real32 tendencies
    to ≈1e-5 relative (≈1e-3 for θ, single-precision cancellation in
    (θ_new − θ_old)/dt), and the diagnostic outputs (`exch_hx`, `hpbl`,
    `wstar`) to ≈1e-6. At dt = 0.01 s the implicit tendency differs from the
    dt = 60 s one by up to 41 % (θ) and 33 % (u), so the small-dt replay is
    REQUIRED for derivative tests (as §6 prescribes) and stage-3 comparisons
    must expect the implicit-vs-explicit gap at the WRF step. Dumps and driver
    outputs live in `data/eqwefic/dumps/<scheme>/`. Remaining: the same hook
    for sfclayrev, WSM6, Dudhia SW, RRTM LW, slab (one driver each), and
    exposing sub-process diagnostics (e.g. YSU's `hgamt`), which needs a fork
    of NCAR/MMM-physics (the `phys/physics_mmm` submodule) — deferred until a
    spike needs it.
    Done 2026-09-05 (WRF fork e42a269, MMM-physics fork
    ctessum-claude/MMM-physics@52cfbd7, branch `earthsciml-instrumented`,
    submodule pointer updated): hooks for sfclayrev, WSM6, Dudhia SW, RRTM LW
    and slab; sub-process diagnostics via `esm_dump_inner` (WSM6 rates
    raw/final + evaluation state, SWPARA per-layer optics, RRTM level fluxes,
    slab energy budget); drivers `sfclayrev_driver`, `wsm6_driver`,
    `sw_driver`, `slab_driver` (real64; `make PREC= B=build32` gives real32
    drivers that reproduce the in-model dumps bit-for-bit — the definitive
    check of the dump→flat→driver chain; RRTM has no driver). Reference runs:
    `data/eqwefic/scm_ref/` (SCM with mp=6, sfclay=1 (= sfclayrev), slab,
    num_soil_layers=5; calls 1,120,600,1440,1800,2400,2820) and
    `data/eqwefic/qss_ref/` (em_quarter_ss supercell, WSM6 warm-rain rows;
    call 9983 = step 250 row 23 is the Spike C reference). Real64 replays
    agree with the real32 dumps to 1e-5 on outputs (larger only in
    threshold-sensitive internals). WSM6 is not rate×dt: limiters scale with
    1/dt and pigen/pcond/freezing are adjustments, so Spike C must test the
    dumped `<rate>_raw` rates pointwise against the dumped `*_rates` state
    rather than replay at small dt. Note `sf_sfclay_physics=1` is sfclayrev
    (91 = the old MM5 scheme).
0.4 **Extraction tool (this repo, `tools/`).** A small script that reads a
    kernel JSON dump plus a hand-written test skeleton and fills in
    `parameter_overrides` and `assertions` with the selected regimes. It never
    writes equations.
    Done 2026-09-05: `tools/fill_tests.py` (sidecar `<esm>.fill.json`; fills
    scalar overrides, per-test input-field libraries for the `input_<x>`
    rewrite targets, and Richardson-extrapolated reference columns as inline
    `Linf_error` references). The input-library route is needed because
    esm-spec §6.6 admits scalar overrides only (format gap, to be filed
    upstream).
0.5 **netCDF toolchain.** Done 2026-09-04. Apptainer image builds work here in
    setuid mode (never pass `--fakeroot`; same route as
    `../moves.rs/characterization/apptainer/build-sif.sh`). Recipe
    `tools/apptainer/wrf-build.def` (Ubuntu 24.04, gfortran 13.3, netCDF-C/F
    4.5.4, HDF5, OpenMPI, python3-netCDF4) → `data/eqwefic/apptainer/wrf-build.sif`
    (334 MB); set `APPTAINER_CACHEDIR`/`APPTAINER_TMPDIR` under
    `data/eqwefic/apptainer/`, never `/tmp`. The 4.8 fork compiles inside it
    (`printf "32\n\n" | ./configure` = GNU serial, `./compile -j 8 em_scm_xy`,
    ~25 min): `main/ideal.exe` and `main/wrf.exe` built 2026-09-04; the
    `em_scm_xy` case (59 h, dt = 60 s, e_vert = 60, mp=2 lw=1 sw=1 sfclay=1
    sf_surface=2 pbl=1) is the SCM reference run (ran 2026-09-04 inside the container: `ideal.exe` then `wrf.exe`, "SUCCESS COMPLETE WRF", 59 h in seconds, output `wrfout_d01_1999-10-22_19:00:00`). The DTC image pulled earlier
    is superseded and can be deleted.
0.6 **Three de-risking spikes**, each ending in a passing `./esm test` file
    on this repo's `main` (not yet a PR):
    - **Spike A — YSU as a PDE.** Column diffusion `∂θ/∂t = ∂/∂z(K ∂θ/∂z) −
      ∂/∂z(K γ)` with the YSU K-profile and nonlocal term, lowered by the 0.2
      rule, tested against `bl_ysu_run` tendencies at `dt = 0.01 s`. Proves the
      "implicit Fortran vs explicit derivative" tolerance story.
      Done 2026-09-05: `components/atmospheric_dynamics/ysu/` (YSU as a
      density-weighted flux-form PDE on `column_nonuniform_1d`, PBL-height
      scans as aggregates, `lib/wrf_thermo.esm` templates), 4 SCM regimes,
      117 assertions green in Rust and Python. Measured implicit-vs-explicit
      gap up to 9e-4 relative at dt = 0.01 s, so references are the
      (0.04, 0.02, 0.01 s) second-order Richardson limit; residual mismatch
      1e-7–1e-6 relative comes from real32 constants/literals inside the
      real64 kernel. Gaps found: (f) shaped `parameter_overrides`/`initial_conditions`
      are number-only in the schema (worked around with test-injected
      `input_<x>` template libraries); (g) the EarthSciModels gate
      `run_esm_inline_tests.py` has its own sampler without shaped-observed
      support (101 errors); (h) Rust rejects the dimension name `lev` as a free
      variable in an inline `reference` and qualified subsystem overrides
      (`P.wrf.g`); (i) the ESD `varcoeff_face_laplacian_lev_flux_bc` rule's
      fixed free names `kdudz_bot/top` give one BC pair per model (proposed
      follow-up rule `D(K·D(u) − F, lev)`). Repros in
      `data/eqwefic/esm-repro/spikeA/`. Cloud-top radiative-entrainment
      branch transcribed 2026-09-05 (`radsum` as the from-bottom cumulative
      integral of cp ρ max(0, −rthraten·exner); `wstar3_2`, `we_rad`,
      `hgamt2`, `wscalek2`, `wscale_c`) with a Fortran-generated cloud-topped
      regime (qc 3e-4 and rthraten −1e-4 K/s in layers 12–14 of call 120,
      replayed at dt 0.04/0.02/0.01): 157 assertions green; residuals
      unchanged (K 5e-6 m²/s, scalars 4e-8 rel). Lesson: WRF's in-branch
      `wscale` is post-diagnosis; folding it into the pre-diagnosis `wscale`
      makes a cycle that the Rust CLI reports as an unrelated
      `E_TREEWALK_UNBOUND_NAME` (worth an EarthSciAST issue).
      Upstream status 2026-09-05: (f) EarthSciAST issue #178 (array-valued
      `default`/`parameter_overrides`/`initial_conditions`, observed-valued
      `ic`); (g) EarthSciModels PR #2 (`run_esm_inline_tests.py` samples
      shaped observeds; the YSU file reaches 117 pass, the other 113 files
      are unchanged at 224 OK / 3 pre-existing load errors); (h) EarthSciAST
      PR #179 (an inline `reference` binds the field's dimension names;
      override keys resolve to the longest dotted suffix; `P.wrf.g` in a
      Rust single-model equation), independent of #177; (i) EarthSciDiscretizations
      PR #35 (branch `column-diffusion-face-flux-rule`, stacked on PR #34;
      rule `varcoeff_face_flux_laplacian_lev`: `D(K·D(u,lev) − F, lev)` with
      no free names, MMS problem `heat_column_varcoeff_faceflux_forced`,
      observed order 2.00; Rust and Julia gates green, Python convergence
      left to CI).
    - **Spike B — Dudhia SW as an integral.** Column optical depth via
      `integral` (cumulative) and the surface flux; tested against
      `module_ra_sw.F`. Proves the `integral` lowering end to end.
      Done 2026-09-05: `components/atmospheric_radiation/dudhia_sw/` — column
      vapour/aerosol/cloud paths as `integral_lev_nodes_cumulative_to_top`
      integrals (WRF's top-down SWPARA arrays reversed onto `lev`/`lev_nodes`),
      Lacis–Hansen vapour absorptance and the Stephens ALBTAB/ABSTAB cloud
      tables as `function_tables` (bilinear), downward beam with the
      cumulative depletions, `gsw`, and the heating rate as
      `face_flux_D_lev_supplied_faces` divided by ρ cp and Exner; 66 equations,
      174 assertions green in Rust. Reference = the real64 `sw_driver` replay
      (diagnostic scheme, no dt); measured mismatch is binary64 roundoff
      (≤ 1e-15 relative), tolerances at 1e-9 of column maxima. Regimes: SCM
      calls 1, 120 (cirrus), 2400 (low sun), 2820, 600 (night, all zero) and a
      synthetic stratus deck (qc = 3e-4 inserted into the call-120 column and
      replayed through the Fortran driver) because the SCM has no liquid cloud.
      Renormalisation (`ff`, `tau_min`) and beam floor (`S_min`) transcribed
      2026-09-05 with a Fortran-generated surface-fog regime (coszen 0.05,
      swrad_scat 3, qc 8.2e-3 in layer 1: exactly the lowest layer
      renormalised; 222 assertions green). Gap (k): below a renormalised layer
      SWPARA's bookkeeping is a coupled 4-component nonlinear sweep (raw
      fractions accumulated, rescaled ones applied), which esm §4.3.1.1 cannot
      express (vector recurrence rejected; repro
      `esm-repro/spikeB/probe_vector_recurrence.esm`), so the component is
      exact only when no depleting layer lies below a renormalised one; the
      floor `S_min` is therefore never reached in a valid regime.
    - **Spike C — WSM6 warm-rain process rates.** `praut`, `pracw`, `prevp`
      as pointwise observeds with WRF constants; tested against instrumented
      WSM6. Proves the sub-process factoring pattern and the constants
      library.
      Done 2026-09-05: WSM6 warm-rain rates as pointwise observeds (unlimited
      physics + separate `dtcld` limiters) over `lev`, saturation state as its
      own model, shared pieces as templates (`lib/wrf_thermo.esm`
      fpvs/qsat/L(T)/cpm, `lib/wrf_air_properties.esm` ν/Dv/ka,
      `wsm6/microphysics_templates.esm`), constants from `wsm6_parameters.esm`
      + `lib/wrf_constants.esm`; 4 warm-rain regimes (supercell mature, second
      column, early storm with the `prevp` limiter active, SCM zero) + 2
      saturation regimes, 194 assertions green in Rust and Python. Reference =
      real64 driver replay at the WRF step, pointwise (no small dt); residual
      ≤ 4e-7 relative from real32 constants/literals; tolerances 1e-6 rel or
      1e-6 of column max abs. Gap (j): §6.6 tests cannot rebind a metaparameter
      (`NLEV`), so the SCM tests are truncated to 40 layers.
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
    Done 2026-09-06 (fork 592f693): driver-level dumps `subassembly_rad`,
    `subassembly_pbl`, `subassembly_mp`, `subassembly_all` for steps 1, 60,
    300, 720, 900, 1200, 1410 (`data/eqwefic/dumps/subassembly`). First
    assemblies live in `couplings/` (EarthSciModels' directory for coupled
    documents; kept out of `components/` because the Rust CLI re-runs the
    mounted components' inline tests under the coupling and a whole-file run
    goes red — run them with `./esm test --model <Assembly> <file>`):
    `surface_pbl_column.esm` (SfclayRev + YSU) and `surface_soil_column.esm`
    (SfclayRev + slab) as top-level ref mounts coupled by `variable_map
    param_to_var`, column state injected per mount (§9.7.10), tests in the
    assembly model (31 + 20 assertions green, step 60). Findings: YSU's
    `psim/psih` = sfclayrev's `fm/fh`; WRF's PBL sees the slab's sub-stepped
    hfx; implicit-vs-instantaneous PBL tendency gap 15 % at dt = 60 (asserted
    loosely, documented). Gaps (EarthSciAST issue #198): (n) document-scoped
    index sets forbid two column components with different NLEV in one
    document; coupling endpoints in models with subsystems → NaN or silently
    dropped; mounted components' inline tests re-run under coupling; top-level
    mounts don't merge leaf index sets; subsystem parameters not settable from
    a parent test (PR #179). Repros: `data/eqwefic/esm-repro/assemblies/`.
    Remaining: the full-physics-sum assembly with microphysics (dumps exist);
    the other six dumped steps as additional regimes; tighten the fork's RRTM
    hook gating (N39). Microphysics-step assembly done 2026-09-07 (below).
1.5 Stubs live on this repo's `main` under `components/<domain>/…` mirroring
    the EarthSciModels layout. They are *not* opened as EarthSciModels PRs
    until their tests pass (finding 9).


### Phase 1 progress (2026-09-05)

- **EarthSciModels PRs opened 2026-09-06 (user decision):** #15 lib
  (`wrf_constants`/`wrf_thermo`/`wrf_air_properties`), #16 sfclayrev, #17 slab
  (new domain `land_surface/`), #18 ysu, #19 wsm6, #20 dudhia_sw and #21
  rrtm_lw stages 1–3 (new domain `atmospheric_radiation/`; refreshed 2026-09-06
  after stage 3, 888 assertions); #18 and #21 refreshed again with the additive
  `<x>_in` coupling-target parameters (EqWeFiC 0496f20); #19 refreshed
  2026-09-07 to 790 assertions (mass conservation), head 7faecd0; each
  component PR carries the lib files and is self-contained. EarthSciModels'
  own CLAUDE.md (upstream b7912bc) forbids `Co-Authored-By: Claude ...`
  trailers, so branches there must drop them (keep `Claude-Session:`);
  branches opened before 2026-09-07 still carry one. ESD rules are referenced by
  relative sibling path because the Rust CLI does not expand `${ESD_ROOT}` in
  template imports. Merge conditions (not ours to address): ESD #34/#35/#36
  merged and reachable from CI, EarthSciAST #177 + EarthSciModels #2 for the
  Python gate.
- **Physics-column assembly done 2026-09-06.** `couplings/physics_column{,_night}.esm`
  mount SfclayRev, YSU, the five RRTM stages and DudhiaSW as top-level systems
  (131 `variable_map` edges); YSU gained the additive `rthraten_in`
  coupling-target parameter (157/157 unchanged), Dudhia needed none; radiation →
  rthraten → YSU and surface layer → YSU couplings reproduce radiation_driver,
  surface_driver (sfclayrev part) and pbl_driver at steps 60 and 300 (40
  assertions each): rthraten Linf ≤ 4.4e-9 K/s, fluxes ≤ 2e-4 W/m², surface-
  layer outputs at the real32 kernel's own floor (1.6e-5 day, 4.9e-5 stable
  night: the stable regime amplifies real32 roundoff in θ_v1 − θ_vs;
  real32-vs-real64 replay 5.3e-5), h 9e-7 / 3.3e-5, PBL tendencies at the
  implicit-explicit gap (15 % of the column max by day, 1.2 % at night;
  asserted loosely). Inputs to the surface layer and YSU come from the KERNEL
  dumps because the driver-level dump's p/psfc differ from the kernel
  arguments by up to 0.28 Pa (N50). Remaining: slab (gap (n)) and WSM6
  (in-place step) in the same document; the other five dumped steps as
  regimes; tighten the `subassembly_pbl` hook's pressure fields in the fork.
- **EarthSciModels fire PRs opened 2026-09-06:** #22 (4 fire_behavior-derived
  tests appended to the existing `wildland_fire/level_set/fire_heat_flux.esm`,
  PLAN 1.3) and #23 (`components/wildland_fire/fire_behavior/`, 464
  assertions).
- **EarthSciModels RADM2 PR opened 2026-09-06:** #24
  (`components/gaschem/radm2/`, 540 assertions), base `main`. RADM2 is
  self-contained — its only refs are `./radm2_ratelaws.esm` and the four
  `./tests/*_inputs.esm`, it mounts no `lib/` subsystem and no
  EarthSciDiscretizations rule — so it does not stack on #15. Finding: none
  of #15-#23 is stacked in git either; each has the then-current `main`
  (7259fd8) as its single parent and the physics branches duplicate the
  `lib/` content inside their own commit (it drops out on rebase). A
  cross-repo PR base must be a branch in `EarthSciML/EarthSciModels`, and the
  `eqwefic/*` branches exist only on the fork, so `main` is the only possible
  base for any of them.

- **Slab done.** `SlabLandSurface` written as instantaneous tendencies: surface
  budget from consumer-supplied `FLHC/FLQC`, soil heat equation
  `D(K·D(T) − F, lev)/capg` with the surface flux G and the fixed deepest layer
  as prescribed interface fluxes (`lev` = soil layer, `ze` = depth,
  NLEV = num_soil_layers − 1). SLAB1D's small-dt replay is the exact derivative
  (nsoil = 1 ⇒ `dthgdt` bit-identical for dt = 0.04…0.005 s; the WRF-step
  change differs by 0.5–11 %). Residuals: real32 `svpt0` → 4e-7 in `es_g`,
  amplified by the saturation deficit and by G = rnet − qs to 5e-5 relative in
  `dTsk_dt`; tolerances 5–10× above, per assertion. 84 assertions green. The
  Bolton saturation template lives in `slab_templates.esm` pending a second
  consumer (sfclayrev) before moving to `lib/wrf_thermo.esm`.
- **sfclayrev done.** Full `sf_sfclayrev_run` as scalar observeds (regimes
  1/3/4, first-guess branch, Fairall/Garratt/Zilitinkevich roughness,
  Charnock/AHW/shallow-water z0, `isfflx`/`scm_force_flux` options); the
  Ri_b → z/L secant iteration is transcribed step by step as a 33-cell
  §4.3.1.1 recurrence because the kernel stops it at |x1−x2| ≤ 0.01 (up to
  8e-5 relative from the root; the residual equation is exposed for a consumer
  who wants the converged physics); psi lookup tables reproduced; 610
  assertions in 18 regimes (12 synthetic, Fortran-generated water/option
  regimes since the SCM is land-only with default options); residual ≤ 4e-6
  relative (real32 `svpt0` through θ_v1 − θ_vs), tolerances 2e-5/1e-6.
  EarthSciModels `surface_layer_profile.esm` (Businger-Dyer) shares no
  function with sfclayrev (Cheng-Brutsaert / Kansas-free-convection blend), so
  nothing was reused. Lessons: a template library may layer on another via
  top-level `expression_template_imports`; the units checker rejects `^0.33`
  of a dimensioned base (unit-valued parameters); coupled iterations must be
  interleaved into one recurrence array.
- **WSM6 done to the process-rate level.** `hydrometeor_slopes.esm`
  (distributions, N_i, D_i, fall speeds; mounted as subsystem `sd` by the rate
  models — the first use of a shaped model as a subsystem, which the Rust CLI
  resolves through the test-injected input libraries), `cold_accretion.esm`
  (13 collection terms), `ice_deposition.esm` (deposition chain with WSM6's
  sequential vapour budget, nucleation/aggregation as target + adjustment-rate
  pairs, evaporation of melting precipitation), `melting_freezing.esm`
  (in-place block as increments + post-block state; post-sedimentation inputs
  reconstructed from the dumps' water and energy budgets); 6 regimes
  (supercell 9983 i=24/25, 9982 i=25, 7980 i=20; SCM 120 lowest-40 and layers
  20–59), 504 assertions green. Residuals 1e-8–1e-7 except where the kernel's
  real32 `t0c`/`pfrz2`/`dimax` are amplified by cancellation (documented per
  test, tolerances 5e-5–2e-4 there). Remaining for WSM6: sedimentation (needs
  the downward upwind/PLM flux rule in EarthSciDiscretizations), the
  saturation adjustment `pcond` (an adjustment formulation), the
  mass-conservation rescaling (time-step artefact). pcond done 2026-09-05:
  `WSM6SaturationAdjustment` (`saturation.esm` mounted as subsystem `sat` for
  qsat at the pre-adjustment state; Asai/RH83 one-step Newton increment, WSM6
  bounds, target state T*, qv*, qc*, kernel rate pcond = dq/dtcld and
  relaxation form dq_adj/τ); 4 regimes (condensing updraft, evaporating
  downdraft with the −qc bound, mixed early storm, SCM zero), 52 assertions
  green against real64 replays; residual ≤ 1.7e-10 kg/kg/s from real32
  constants in qsat. Sedimentation done 2026-09-06 (`sedimentation.esm`, 65
  assertions): the kernel's `fall` arrays are post-fallout density × speed /
  delz, so references are the Richardson dt→0 limit of real64 replays at
  0.04/0.02/0.01 s, which equals the donor-cell flux ρ v_t q; WSM6's PLM
  increment at dt = 12 s differs from the donor-cell divergence by factors
  0.4–1.6 (rain) with sign flips for snow/graupel — the instantaneous-
  derivative form is scheme-independent only at t = 0. Gap (m): subsystem
  index sets fold metaparameters before merging, so grid imports must be bound
  to the literal `NLEV` when a shaped model is mounted. Format notes: a test cannot
  assert a mounted subsystem's variable; a subsystem's `index_sets`/
  `metaparameters` must not be redeclared by the importer; `-` is strictly
  unary/binary while `+`/`*`/`min`/`max` are n-ary.
  Mass-conservation rescaling done 2026-09-07 (`mass_conservation.esm`, 169
  assertions): the block "check mass conservation of generation terms and
  feedback to the large scale" as one cumulative factor per species budget
  (`mass_conservation_factor(q, q_floor, S) = max(q_floor,q)/max(S,max(q_floor,q))`,
  the branch-free form of `if (source > value) factor = value/source`) times the
  raw rate — praut f_qc f_qr, paacw f_qc f_qs f_qg, piacr f_qr f_qs f_qg — plus
  the state update it feeds, both branches blended by a `cold` indicator; new
  scheme constant `par.q_delta = 1e-4` (WSM6's unnamed delta2/delta3 literal).
  Reproduces the real64 replay EXACTLY (0.0) in all 22 `<rate>_final` columns,
  the switches and the `*_upd` state; residual only dqv_dt 2e-19 kg/kg/s,
  dT_dt 2e-15 K/s, one ulp in a factor. Activation over the whole dumped set:
  no SCM column binds at all; in the supercell the rain budget binds hardest
  (f_qr = 0.5877 at 9983 col 30 layer 16), then ice (0.9289) and snow (0.9736);
  the cloud-water and graupel budgets and the ENTIRE warm branch never bind
  materially (N59), so the warm-branch factors are transcribed but untested at
  the binding level — their update path is exercised in every regime.
- **Microphysics-step assembly done 2026-09-07.** `couplings/microphysics_column{,_mixed}.esm`
  mount seven WSM6 stages as top-level systems (106 `variable_map` edges) in the
  order mp_wsm6_run calls them: Saturation -> MeltingFreezing -> {WarmRain,
  ColdAccretion, IceDeposition} -> MassConservation -> SaturationAdjustment,
  with WarmRain's prevp opening the ice deposition budget and Saturation coupled
  forward into every stage from the START-of-sub-step state (N4). Each stage
  gained additive `<x>_in` coupling-target parameters (standalone 790/790
  unchanged) that mount-edge libraries `couplings/tests/wsm6_couple_*_inputs.esm`
  lower `input_<x>` to. NEW: coupling INTO a mounted subsystem's parameter
  (`Cold.sd.qr_in`, `IceDep.sd.qs_in`, `SatAdj.sat.T_in`) works in the current
  Rust CLI — the `probe_wrap_couple` NaN row of `esm-repro/assemblies/README.md`
  no longer holds — which is what makes the full chain expressible. Sedimentation
  is NOT in the chain: WSM6's PLM semi-Lagrangian fallout is not a rate (its
  sub-step increment differs from the donor-cell divergence by 0.4-1.6x with sign
  flips), so the post-sedimentation qi/qr/qs/qg entering MeltingFreezing are
  dumped profiles and the chain covers six of the seven stages. 2 x 40/40 with
  `./esm test --model MicrophysicsColumn <file>`; end-to-end residual (Linf as a
  fraction of the column maximum) <= 2.4e-7 in the melting/freezing outputs,
  6.1e-6 in the rates, 3.5e-7 in the updated state and 1.2e-5 in pcond, one to
  two orders above the standalone 1e-7 from accumulated real32 constants (N2);
  tolerances 5e-8 (T, latent heats), 5e-6 (mixing ratios), 5e-5 (rates), 1e-4
  (pcond). Gap (m) extended: an importing DOCUMENT's own `index_sets` are also
  compared against the folded subsystem declaration, so `lev` must be sized by
  the literal 40 and not by the document's `NLEV` metaparameter, or the load
  fails with `[subsystem_index_set_conflict]` (EarthSciAST #198).
- **Thermo consolidation done.** `exner_function`,
  `bolton_saturation_vapor_pressure`, `bolton_saturation_mixing_ratio` and
  `dry_air_density` moved into `lib/wrf_thermo.esm` from slab, sfclayrev, ysu,
  dudhia_sw and rrtm_column (bodies verbatim, add-only in lib); suite identical
  at 1845/1845. Remaining single-consumer helpers stay beside their components
  (Beljaars w_c in `sfclayrev_thermo`; bulk flux / black-body / net radiation
  in `slab_templates`).
- **fire_behavior done to the component level (2026-09-06).** Fork
  ctessum-claude/fire_behavior@c9dda02 (`-DESM_DUMP=ON` CMake option, dumps of
  ROS/level-set/fuel/flux/moisture internals; standalone CMake build in the
  container; `tests/test7` rerun with `num_tiles = 1`, plus `fire_upwinding = 4`
  and `fmoist_run` variants; dumps `data/eqwefic/dumps/fire/`). Nine components
  under `components/wildland_fire/fire_behavior/` (parameters, Anderson13 table
  incl. waf, Rothermel-WRFFire fuel-bed parameters, ROS with projection/caps/
  chaparral, fuel-fraction cell with 2×2 submesh, five-class fuel moisture ODE,
  ignition time, Byram flame length, level-set PDE) = 464 assertions, plus 4
  PR-ready tests on a copy of `level_set/fire_heat_flux.esm` (22); 486/486
  Rust. Only `fire_heat_flux` matches EarthSciModels algebraically;
  `fuel_model_lookup`, `rothermel/fire_spread` and `level_set_fire_spread`
  differ (see component descriptions). Missing EarthSciDiscretizations rules
  for the default `fire_upwinding = 9`: WENO5/ENO1 hybrid |∇ψ| and the
  Godunov-upwind first derivative; boundary = linear-extrapolation halo.
  fire_behavior is real32 only (rel 1e-5). Gap (o): an observed named `t`
  collides with the time variable without a validate error
  (`esm-repro/fire/probe_observed_named_t.esm`). Not started: atmosphere→fire
  wind (`Interp_profile`), the coupling file, reinitialisation, ignition lines,
  smoke; the rain/wetting moisture branch is transcribed but untested.
  Atmosphere–fire coupling 2026-09-06: fire_behavior fork 8d519c6 adds
  ESM_DUMP hooks to `Calc_fire_wind` / `Interp_wrfwinds_to_cfbm` /
  `Provide_atm_feedback`; WRF fork branch `earthsciml-instrumented-fire`
  (c01746f) builds WRF 4.8 with CMake + ENABLE_CFBM + ESM_DUMP in the container
  (25 min; the classic `./compile` has no CFBM), but no shipped case can run
  CFBM (Lambert-only; ideal em_fire is map_proj 0) and test7/wrf.nc is an
  ifire = 2 (SFIRE) output, so the coupled references are fire_behavior's own
  WRF-side routines replayed by `data/eqwefic/fire_build/feedback_driver` on
  wrf.nc columns (`dumps/fire_wrf/`). Components: `FireWindWRFFire` (21
  assertions, all Interp_profile branches via fire-wind-height overrides) and
  `FireFluxToAtmosphereWRFFire` (54; grnhfx aggregation identical to SFIRE's,
  Fire_tendency differs from SFIRE's prop_heat profile by 0.7–1.5 %);
  `couplings/fire_atmosphere_column.esm` chains FireHeatFlux → flux/tendency (9
  assertions). Remaining: a real-data Lambert WRF+CFBM case for a true in-model
  coupled reference; the horizontal atm→fire mapping (nearest/bilinear,
  projection) and the fire → smoke tracer path. Authoring note: an `aggregate`
  with `output_idx: ["k"]` is a shaped map; a scalar reduction needs
  `output_idx: []`, otherwise the observed is silently unmaterialised.
- **RADM2 (EqAtmChem) 2026-09-06.** KPP 2.1 built on the host (flex/bison;
  two Makefile patches), radm2 generated from `chem/KPP/mechanisms/radm2`,
  standalone real64 box driver in the container (`data/eqwefic/chem_build/
  driver`: `Update_RCONST` with WRF's include lists expanded, `radm2_Fun`/
  `IRRFun`, Rodas3 trajectories) dumping RCONST/A/Vdot/trajectories for three
  documented typical states (no WRF-Chem run). `components/gaschem/radm2/`:
  11 hand-authored WRF rate-law templates; `radm2.esm` translated from the
  .eqn by `radm2_from_kpp.py` (`data/eqwefic/chem_build`) with the reaction
  rates as scoped references to the `RADM2RateConstants` model so every
  coefficient lives once; negative KPP product coefficient (CSL+OH, −0.9 OH)
  as a shadow reaction. 405 assertions green in Rust at rel 1e-9 (rc_n2o5
  1e-6, WRF evaluates it in real32). Gaps (EarthSciAST issue #206): (p) the
  Rust CLI runs no `reaction_systems` tests; (q) an observed `D(x,t)`
  evaluates to 0, so Vdot cannot be asserted at t=0; bare
  `parameter_overrides` resolve per component in Rust but document-wide in
  Python (model parameter renamed `M_air` after a one-off gate check).
  WRF-Chem reference run done (2026-09-06): chem-enabled WRF built in
  `apptainer/wrf-chem-build.sif` (recipe `tools/apptainer/wrf-chem-build.def`;
  flex/bison, two KPP-2.1 Makefile patches, WRF_CHEM=1 WRF_KPP=1, compile twice
  — N57), idealized `em_scm_xy` column `data/eqwefic/chem_scm` with
  `chem_opt = 101` (the KPP RADM2; `chem_opt = 1` is the hand-coded
  `module_radm.F`), `phot_opt = 1`, `emiss_opt = 0`, `chem_in_opt = 0`, 59 h,
  59 levels. Dumps of the KPP interface (`var`, `fix`, `RCONST(1:156)`,
  `jv(1:52)`, TEMP, C_M, C_H2O, rc_n2o5, p, rho, qv, pre/post-INTEGRATE `var`)
  through the six `kpp_mechd_*_radm2.inc` coupler hooks, fork
  `ctessum-claude/WRF` branch `earthsciml-instrumented-chem` @ 2ce0388 —
  the only insertion points the KPP coupler offers inside the k loop.
  `radm2.esm` gained a fourth regime `wrfchem_scm_midtrop_rate_constants`
  (135 assertions, level 45, T = 254.7 K): 540/540 green at rel 1e-9. The KPP
  interface promotes t_phy/rho_phy/ph_* to real64 before Update_Rconst, so the
  single-precision WRF build does NOT force rel 1e-5 on the rate coefficients;
  only rc_n2o5 is real32 (rel 1e-6). Cross-check: WRF-Chem and box-driver
  RCONST agree bit-for-bit at 5 of 6 dumped levels. The run's SPECIES are not
  usable as a reference: bug B8 (CO2 missing from the interface copy loops)
  poisons the KPP vector at call 3 and the Rosenbrock fails for the rest of
  the run. Remaining: Vdot/stoichiometry verification once (p)/(q) are fixed;
  a WRF-Chem species trajectory once B8 is patched; RADM2SORG/aerosol;
  emissions, Wesely and FastJX couplings (FastJX covers 10 of 21 j inputs;
  Wesely ≈20 species).
- **RRTM LW stage 1 done.** Heating-rate convention proved:
  `HTR(L−1) = HEATFAC (FNET(L−1) − FNET(L))/(PZ(L−1) − PZ(L))` is the heating
  of layer L; RRTM indexes bottom-up so `TOTUFLUX/TOTDFLUX(0..kte)` map onto
  `lev_nodes` with no reversal; WRF's OLR is the model-top flux. Written as
  `dT/dt = −(g/cp) D(F_net, lev)` with the column grid's `ze = p_sfc − p_e`
  (pressure coordinate), 56 assertions; the MM5ATM/SETCOEF column mapping
  (NBUF = nint(p_top/(100 deltap)) buffer layers, standard-atmosphere
  temperatures, coldry and column amounts, cloud optical depth) as
  `RRTMColumn`, 200 assertions; both real32 references at rel 1e-5. Stage 2
  (RTRN sweep, per-g-point recurrences on `rlay`) is blocked on
  instrumentation: dump `TAUG, PFRAC, ITR` (NGPT × NLAYERS),
  `TOTUCLFL/TOTDCLFL`, band Planck integrals; stage 3 (TAUGB1–16 k-tables,
  ≈1e5 coefficients) needs `data_sources`. `colo3` deferred (O3DATA as a
  pressure integral). Gap (l): Rust `esm test` does not lower `table_lookup`
  (esm-spec §9.5.3; repro `esm-repro/rrtm/probe_table_lookup.esm`; EarthSciAST
  issue #188), so tables
  are spelled as `fn interp.linear` on `const` arrays; also `ifelse` branches
  evaluate eagerly, so index gathers in inactive branches must be clamped.
  Stage 2 done 2026-09-06 (fork 0e52608 dumps TAUG/PFRAC/ITR, SETCOEF
  coefficients, RTRN Planck integrals, clear-sky and per-band fluxes; dumps
  `scm_ref2_<call>.json`): `rtrn_sweep.esm` transcribes RTRN as two esm
  §4.3.1.1 rank-2 recurrences over the level index (one cell frame
  `[gpts, rlev]`; the downward sweep on a top-counted index since a self-read
  must be strictly earlier), TF/TAU as closed forms of the quantised ITR,
  TOTPLNK/DELWAVE/NGB as a script-transcribed table file; 112 assertions in 4
  regimes, residuals ≤ 7e-7 relative on fluxes and 4e-5 on HTR. Stage 3 (gas
  optics) needs: the 16 TAUGBn k-tables ABSA/ABSB/SELFREF/FORREF/FRACREFA/B
  (≈1e5 coefficients, `data_sources`) and the SETCOEF interpolation (dumped
  fac00..fac11, forfac, selffac, selffrac, jp, jt, jt1, indself, laytrop,
  layswtch, laylow, co2mult now available as inputs/references), plus the
  GASABS quantisation `ITR = INT(5000 od/(BPADE+od) + 0.5)`. Tooling: the
  filler has no rank-2 fields (flattened `gl` inputs) or rank-2 reductions.
  Stage 3 done 2026-09-06: `rrtm_setcoef.esm` (SETCOEF, 140 assertions) and
  `rrtm_gas_optics.esm` (TAUGB1–16 + GASABS ITR, 340 assertions, all 16 bands,
  50 s in the Rust CLI); k-tables transcribed by script from RRTM_DATA
  (big-endian real32 sequential file) + module DATA statements after the CMBGB
  256→140 g-point reduction into 16 `const` template libraries (1.5 MB;
  `data_sources`/`from_file` are not usable in inline tests, gap (f)); `colo3`
  added to `rrtm_column.esm` as the O3DATA layer integral (240 assertions) and
  WRF's off-by-one ozone layer confirmed (B5). Residuals: taug ≤ 2e-7 of band
  maxima, pfrac ≤ 9e-8, itr exact up to quantisation flips. Remaining:
  nothing at the stage level. (The 2026-09-06 note that "all SCM clouds are
  overcast" was wrong: calls 600 and 1440 already carried fractional layers
  of 0.015-0.037.)
  End-to-end assembly done 2026-09-06: `couplings/rrtm_lw_column*.esm` chain
  the five stages as top-level ref mounts; each stage keeps its `input_<x>`
  rewrite targets and gained shaped `<x>_in` coupling-target parameters
  (additive, standalone 888/888 unchanged) that the mount-edge libraries
  `couplings/tests/rrtm_couple_*_inputs.esm` make the targets resolve to,
  filled by `variable_map param_to_var` from the upstream stage — the general
  coupling path for every `input_<x>` column component; rank-2 fields
  flattened (gas optics gained `pfrac_flat`). 4 regimes × 25 assertions green
  (18 s each); end-to-end residuals: fluxes ≤ 2e-4 W/m² except 3.3e-3 at call
  1440 (±1 ITR quantisation flips propagated through the sweep), htr ≤ 1.3e-3
  K/day, dTdt ≤ 1.5e-8 K/s, glw/olr ≤ 1.3e-5 rel; tolerances 5e-3 W/m² /
  5e-3 K/day / 2e-8 K/s / rel 2e-5. `couplings/radiation_column{,_night}.esm`
  add DudhiaSW and reproduce radiation_driver's rthraten, gsw, glw, olr at
  steps 60 and 300 (30 assertions each). Gap (r) (EarthSciAST issue #207):
  the Rust CLI folds a gather from a bare `const` array at interpreter-build
  time before document couplings are applied (`E_TREEWALK_CONSTARRAY_OOB` with
  the coupled temperature = 0); wrapping the table in an aggregate defers it
  (`rtrn_sweep.esm` `totplnk`; repro `esm-repro/assemblies/rrtm/`).
  Cloud-overlap regimes done 2026-09-07. RTRN's cloudy branch and MM5ATM's
  cloud block were already transcribed; what was missing was coverage. Four
  regimes added (324 component assertions, rrtm_lw 888 -> 1212, plus one coupled
  document): optically thick overcast (call 120, 7 layers, TAUCLOUD <= 0.760,
  ABSCLD <= 0.72, HTR - HTRC = 25 K/day), fractional cloud (call 240, layer 55
  CLDFRAC = 0.1470 at TAUCLOUD = 0.009966), the zero-fraction override (call
  396, CLDFRA = 0 with TAUCLOUD = 0.01072 handed to RTRN as overcast) and a
  three-layer fractional overlap (call 648). Key finding (N58): MM5ATM's
  `IF (TAUCLOUD > 0.01) CLDFRC = 1` (`module_ra_rrtm.F:4254`) means RTRN can
  only ever see 0 < CLDFRAC < 1 where the cloud optical depth is at or below
  0.01, so the fractional branch is reachable only with ABSCLD <= 0.0165; call
  240 sits exactly at that ceiling and is the strongest fractional case the
  scheme admits. References came from re-running the existing em_scm_xy case
  for 8 h with 71 calls dumped (`data/eqwefic/scm_cloud`, dumps
  `dumps/scm_cloud_raw`, 97 MB); the re-run reproduces `scm_ref2_120.json`
  bit-for-bit, so it is the same trajectory and no fork change was needed.
  Residuals (real32, rel 1e-5 contract): fluxes Linf <= 2.5e-4 W/m^2, band
  fluxes <= 1.5e-5, Planck <= 2.2e-10, HTR/HTRC <= 4.3e-4 K/day, TAUCLOUD
  <= 7.8e-8 on a column maximum of 0.760 (so the thick regime's taucloud
  tolerance is abs 1e-6, not the 1e-8 of the thin ones), CLDFRAC exact; the HTR
  spot checks of the fractional regimes use rel 2e-4 rather than 1e-4 because
  |HTR| there is only 0.6-1.8 K/day. `couplings/rrtm_lw_column_call240.esm`
  adds the fractional regime end-to-end (25/25, residuals <= 1.5e-4 W/m^2,
  2.8e-4 K/day, 3.4e-9 K/s); the overcast regime was already covered end-to-end
  by `couplings/radiation_column.esm`/`physics_column.esm` (call 120) and the
  weak fractional one by their `_night` variants (call 600), so no other coupled
  document needed a new regime. Mutation-checked: replacing EFCLFRAC =
  ABSCLD*CLDFRAC by ABSCLD fails 16/28 assertions in each new fractional regime
  and leaves clear and pure-overcast green; removing the TAUCLOUD override fails
  all four new RRTMColumn regimes. No esm expressiveness gap: there is no
  CLDPROP and no RTRNMR in `phys/module_ra_rrtm.F` (RTRNMR is RRTMG's, a
  deferred row), and every term of loops 220/2000/4000 is expressible. Still
  untested: the liquid, rain and snow optical-depth terms (ABCW = 0.144,
  ABRN = 0.330e-3, ABSN = 2.34e-3) -- the SCM reference case is ice-only
  (QCLOUD = QRAIN = QSNOW = QGRAUP = 0 for the whole 59 h), so a warm- or
  mixed-phase case with radiation on (e.g. em_quarter_ss with
  ra_lw_physics = 1) would be needed to exercise them.

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
