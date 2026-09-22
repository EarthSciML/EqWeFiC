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
**`chem_opt=170` (CBMZ_MOSAIC_KPP: CBM-Z gas phase + 8-bin MOSAIC aerosol)**
— revised 2026-09-22, replacing `chem_opt=300` (GOCART), because WRF-Chem has
no RADM2 + MOSAIC pairing; see section 1.4b. Components: anthropogenic emissions
(`emissions_driver.F`), biogenic emissions, dry deposition
(`dry_dep_driver.F`, Wesely: reuse EarthSciModels), photolysis
(`module_phot_mad.F`, Madronich -- the SCM runs `phot_opt = 1`; Fast-JX is NOT
a WRF-Chem option at any `phot_opt`, so the EarthSciModels `fastjx/*` components
match nothing here and are NOT reused -- see `data/eqwefic/notes/photolysis_diagnosis.md`),
gas mechanism
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
    Remaining: the other four dumped steps as additional regimes; tighten the
    fork's RRTM hook gating (N39). Microphysics-step assembly done 2026-09-07
    and the full SCM physics-suite assembly done 2026-09-07 (both below); the
    full-physics-sum assembly WITH microphysics in one document is blocked, not
    outstanding work — gap (m) below.
1.4a GWDO reference run on real terrain (done 2026-09-22; details in
    `data/eqwefic/notes/gwdo_real_reference_run.md`). GWDO does nothing on the
    idealized SCM column, which has no sub-grid orography.
    - **Run.** A real-data run supplies the reference: d01 of the Last Chance
      Gulch WPS domain (8.1 km, 150 × 180 × 51, Colorado Rockies), 18 h from
      2012-06-25 12Z, with `gwd_opt = 1`, `cu_physics = 16` and
      `sf_surface_physics = 2` over the EqWeather suite. It uses a dmpar build
      of WRF fork branch `earthsciml-instrumented-gwdo` (790bb4f).
    - **Dumps.** New positional dump selection by step and global row/column
      produced 140 GWDO rows, 140 New Tiedtke rows and 210 Noah `SFLX` points.
    - **Kernel replay.** A real32 replay of `bl_gwdo_run` reproduces WRF bit for
      bit on all 140 rows; the real64 driver is `kernels/gwdo_driver`.
    - **`ELVMAX`.** WRF 4.8's GWDO needs it, and WPS 4.3 does not write it.
      Without it the scheme is silently off everywhere (FORTRAN_BUGS B17). It
      was supplied by WPS `util/compute_gwdo.py`.
    - **Coarser grid.** A 12–30 km run needs fresh ERA5 intermediate files.
      It would add regimes, not code paths.
1.4b chem_opt = 170 reference run on the idealized SCM column (done 2026-09-22;
    details in `data/eqwefic/notes/chem170_scm_reference_run.md`). This is the
    stage-1 data for the chemistry model's new target and for the two newly
    required weather schemes on the column where they belong.
    - **Decision (2026-09-22).** EqAtmChem targets **`chem_opt = 170`,
      `CBMZ_MOSAIC_KPP`** — CBM-Z gas phase through KPP plus 8-bin MOSAIC —
      instead of the GOCART plan in section 1 (`chem_opt = 300`). WRF-Chem has
      **no RADM2 + MOSAIC option**: `Registry/registry.chem` pairs MOSAIC only
      with CBM-Z, SAPRC-99, MOZART and CRI, and pairs RADM2 only with SORGAM
      (2, 11, 41, 106) or GOCART (303). 170 keeps the KPP route the RADM2 work
      already uses, so Madronich photolysis, Wesely deposition and the chem
      vertical mixing carry over unchanged and only the reaction list is new.
      The finished RADM2 components stay as they are; they are simply not the
      mechanism the aerosol model runs on.
    - **Cumulus, GWDO and Noah are now required, not optional.** They stay out
      of the EqWeather SCM reference suite, which has them off, but the model
      cannot run at coarse resolution or over real terrain without them.
    - **Run.** `em_scm_xy`, Weisman-Klemp sounding, Kansas, 1997-06-20 12 UTC
      + 48 h, dt = 60 s, 59 levels to 20 km, dx = 20 km, with New Tiedtke
      (`cu_physics = 16`) and Noah (`sf_surface_physics = 2`). The existing
      CASES-99 profile is stable and dry and never triggers convection.
      WRF fork branch `earthsciml-instrumented-chem170` (eae782d); the
      unchanged `wrf-chem-build.sif` already compiles every KPP mechanism.
    - **Dumps.** 16 steps × 13 schemes (`data/eqwefic/dumps/chem170_wrf`),
      including five new hooks: `ntiedtke`, `noah`, `cbmz_kpp`, `mosaic` (the
      column state after EACH of gas-particle transfer, nucleation and
      coagulation, separately) and `mosaic_drydep`. All three Tiedtke
      convection types plus the inactive branch, Noah day and night, and all
      three MOSAIC sub-processes are covered.

1.4c Noah stage 2, first tranche (done 2026-09-22). Five hand-authored
    components under `components/land_surface/noah/`, 448/448 green at
    rel 1e-9 against the real64 `noah_driver` replay of the 16 idealized-column
    dumps of 1.4b:

    - `noah_parameters.esm` -- Noah's own thermodynamic constants and the
      GENPARM scalars. Deliberately separate from `lib/wrf_constants.esm`:
      Noah shadows rd and sigma, and PENMAN shadows cp again, so one call
      evaluates three different values of two constants (FORTRAN_BUGS N88).
      The GENPARM defaults are the REAL32-WIDENED table values, not the
      table's decimal text: those numbers reach the scheme through a real32
      array, and using the decimal for `cmcmax` puts a 4.7e-8 error in a
      canopy-wetness ratio near one, which the transpiration amplifies to
      1.5e-5.
    - `canopy_resistance.esm` (CANRES), `potential_evaporation.esm` (PENMAN),
      `evapotranspiration.esm` (EVAPO/DEVAP/TRANSP, including the dew branch)
      and `soil_heat.esm` (TDFCND + HRT + the yy/zz1 surface closure).
    - The soil heat equation is the flux-form PDE D(K D(T) - F) lowered by the
      existing ESD rule `varcoeff_face_flux_laplacian_lev` -- the same rule the
      slab scheme uses; no new ESD rule was needed. Both boundary fluxes are
      prescribed (the linearised surface energy balance at the top, conduction
      to `tbot` at the bottom), which is exactly the form the rule expects.
      The interface conductivity is UPWIND, not a face average, because that is
      what the Fortran does (FORTRAN_BUGS N90).
    - The root-zone stress factors are written as column integrals: `rcsoil` as
      a thickness-weighted mean over the root zone, and TRANSP's `sgx` as the
      UNWEIGHTED layer mean, which is the integral of gx/dz over the masked
      root zone divided by nroot.
    - **Two dt scales are needed for the references and the distinction is
      load-bearing.** SFLX is an in-place step: the surface closure (yy, zz1,
      t1, ssoil, sheat) is referenced at ESM_DT = 1e-6 s, where the O(dt)
      contamination from SMFLX having already moved the soil moisture is below
      1e-10 relative; the soil-temperature tendency must use ESM_DT = 0.01 s,
      because at 1e-6 s the increment of two ~300 K temperatures is lost to
      cancellation. WRF's own dt = 60 s increment is off the dt -> 0 limit by
      up to 2.7e-7 K/s -- the implicit solve, not the physics.
    - **Instrumentation added.** `df1`, `yy`, `zz1` and the pre-step ground
      heat flux are not Registry variables. `kernels/noah_esm_{state,penman,
      nopac}.inc` are spliced by `kernels/Makefile` into a BUILD-TIME COPY of
      `phys/module_sf_noahlsm.F` (the WRF source is untouched), following the
      pattern the New Tiedtke driver already uses.
    - **Still to do:** soil moisture (SMFLX/SRT/SSTEP, Richards plus the three
      runoff terms), the frozen-soil sink (SNKSRC/FRH2O/TMPAVG/TBND), the snow
      pack (SNOPAC and friends), REDPRM's table lookups as a component, and the
      urban and `opt_thcnd = 2` branches. The reference set is above freezing
      and snow-free, so those branches are UNTESTED, not merely unused -- a
      consumer must not mount these components where they matter.
    - **Kernel replays.** New real64 drivers `kernels/{ntiedtke,noah}_driver`.
      Noah agrees with WRF to real32 round-off. Tiedtke's in-model `rthcuten`
      is cancellation-limited at ~3e-7 K/s (N84), so the real64 replay is the
      reference and the in-model dump only a cross-check.
    - **Kernel replays for the chemistry (added 2026-09-22).** `cbmz_driver`
      reproduces the KPP CBM-Z rate coefficients and integrated step EXACTLY,
      and supplies `Fun` — the instantaneous rates, 4x to 35x the finite-step
      increment by day, which is what a reaction-system component must be
      referenced against. `mosaic_drydep_driver` agrees to real32 round-off.
      `mosaic_subproc_driver` (nucleation, coagulation) reproduces WRF bit for
      bit in its real32 build; its real64 build is an open item, so the real64
      reference for those two stages is still owed (FORTRAN_BUGS N87).
    - **The chem-on / chem-off divergence is explained (N86):** Dudhia
      shortwave adds the chem PM2.5 mass to its layer scattering whatever
      `aer_ra_feedback` says, so this run's meteorology is aerosol-coupled. It
      does not affect stage-1/stage-2 component references; it does mean a
      stage-3 EqWeather comparison must use the `chem_opt = 0`/`1` trajectory
      (bitwise identical to each other) or mount that scattering term.
      **Decided 2026-09-22 (user): stage-3 EqWeather references the
      `chem_opt = 0` trajectory.** The chemistry run is then a reference for
      chemistry components and for a future aerosol-coupled EqAtmChem
      assembly, never for a weather-only comparison. Mounting the Dudhia
      aerosol scattering term is deferred until MOSAIC exists in .esm; when
      it does, it is the coupling edge that lets a coupled model reproduce
      this run, and it belongs in `couplings/lib/` beside the other pairwise
      libraries.
    - **Emissions, wet scavenging, cloud chemistry and convective tracer
      transport stay off** (the last has no New Tiedtke path in WRF-Chem), so
      they remain stage-1 gaps, as they were for RADM2.
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
  branches opened before 2026-09-07 still carry one. ESD rules were referenced by
  relative sibling path because the Rust CLI did not expand `${ESD_ROOT}` in
  template imports; **EarthSciAST #400 fixed that (merged 2026-09-17,
  `b136a9b92`) and the refs have since moved to `${ESD_ROOT}`**, which is how
  they resolve in CI (EarthSciModels `.github/workflows/test-esm.yml:62-67`
  clones ESD and exports the variable). Merge conditions (not ours to address): ESD #34/#35/#36
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

  **Slab co-mounted 2026-09-12 (gap (n) closed).** EarthSciAST #213 + #298
  removed the index-set obstruction, and `couplings/scm_physics_column{,_night}.esm`
  now mount nine components — SfclayRev, the slab, YSU, the five RRTM stages and
  DudhiaSW — at 2 x 50/50 (167 coupling edges each). The soil column enters
  through a mount-edge `index_set_rename` (`lev` -> `soil`, `lev_nodes` ->
  `soil_nodes`, `NSOIL = 4` bound at the edge), radiation feeds the surface
  energy balance directly (the dumped `pbl_gsw`/`pbl_glw` are bit-identical to
  `sw_gsw`/`rrtm_glw`, so radiation_driver's own output is what surface_driver
  saw), and the slab's `hfx`/`qfx` now DRIVE YSU instead of being prescribed
  from the dump. Mutation-checked: cutting `Sfc.flhc -> Lsm.flhc` fails 6
  assertions including `pbl_hpbl`, `K_h` and `K_m`; cutting `Sw.gsw -> Lsm.gsw`
  fails both soil-temperature tendencies. The residual against WRF is the N12
  sub-stepping gap (the esm slab is an instantaneous flux calculator, WRF's SLAB
  reports the flux at the end of its internal sub-step): hfx 7.5e-4 day /
  2.7e-4 night, qfx and lh 3.1e-4 / 2.4e-4, qsfc 1.8e-4, `lsm_dTs_dt` Linf
  5.2e-6 K/s; `lsm_capg` and `lsm_land` exact. Because YSU now sees the
  instantaneous flux, the day column's `K_h`/`K_m` moved by Linf 0.042 / 0.025
  m^2/s (2.4e-4 of their column maxima) and those two assertions loosened from
  `abs` 0.002/0.001 to 0.06/0.04 — the one place this work weakened an existing
  assertion. Gap (m)'s other half is also gone: a document must no longer
  redeclare a mounted leaf's index sets, and doing so now breaks the load.
  **WSM6 is absent for physics reasons, not format reasons**: the SCM reference
  column is cloud-free (`qc = qr = qs = qg = 0` in every layer at both steps,
  `qi` a 2.1e-9 / 5.2e-8 kg/kg trace), the merged suite dump carries no `mp_`
  variables at all, and microphysics is a separate driver call at the end of
  WRF's step — mounting it would add ~7 components and ~100 edges asserting
  nothing. It is covered on the supercell column by
  `couplings/microphysics_column{,_mixed}.esm` instead.
  Three hazards found while doing this: (i) `tools/fill_tests.py` is
  DESTRUCTIVE on this document family — regenerating the day document reproduced
  all 50 assertions and 47 overrides bit-identically (a strong check on the hand
  transcription) but wiped `tests/rrtm_couple_heating_call2820_inputs.esm`,
  which is only half generated: its two coupling-target lowerings
  `input_F_up -> F_up_in` and `input_F_dn -> F_dn_in` cannot be expressed in a
  sidecar `libraries` entry, and without them every assertion fails with
  `unlowered_operator`. The file was restored and both sidecars carry a warning
  `_note`; the `libraries` schema needs a way to carry hand-authored lowerings
  before that script is safe to re-run here. (ii) A parent test still cannot set
  a mounted component's parameter (`"Lsm.thc": 0.02` -> `Invalid parameter
  'Lsm.thc'`, probed on current main), which is what forces the
  prescribe-and-couple idiom. (iii) `.esm` refs reach
  `../EarthSciDiscretizations` by relative sibling path, which does not exist
  from a nested git worktree; a worktree agent must symlink it.
- **EarthSciModels fire PRs opened 2026-09-06:** #22 (4 fire_behavior-derived
  tests appended to the existing `wildland_fire/level_set/fire_heat_flux.esm`,
  PLAN 1.3) and #23 (`components/wildland_fire/fire_behavior/`, 464
  assertions).
- **PR refresh 2026-09-15.** #15-#22 rebased onto EarthSciModels main b7912bc
  with their `Co-Authored-By` trailers dropped; no file changes were needed and
  every branch verifies at its full count in the EarthSciModels tree. #23 went
  red on the current CLI (425/0/39, all in `godunov_step20_stage1`) from an
  EarthSciAST loader defect, filed as **EarthSciAST #358** (Rust only; Python
  passes, Julia's injection step does too): `aggregate` normalisation raises the
  version only on the file that spelled it, so when a test injects a library the
  re-serialised leaf, still at 1.0.0 but now holding the rule's `faq` body, fails
  `faq_version_too_old`. Worked around by declaring `"esm": "1.1.0"` on
  `level_set_tendency_wrffire.esm` in the PR (head 623fd30, 464/0/0); EqWeFiC's
  copy is unaffected because it imports the rules at test scope. The newer fire
  work (`fire_wind_wrffire`, `fire_flux_to_atmosphere_wrffire`, the WENO5/ENO1
  and driver regimes) will follow as a separate PR once ESD #34 and #37 merge,
  since it adds those dependencies.
- **Upstream unblocking, 2026-09-16.** All four PRs opened from here merged:
  EarthSciAST **#362** (Rust `maxiters` uncapped by default — the 24 h
  single-column test is no longer capped, retiring the ≤6 h window workaround),
  **#369** (a `coupling_import` ref resolves against the importing document),
  **#370** (a coupling library's role-bound edges are validated, merged as "in
  all five bindings" after review extended the Rust+TS change), and **#371**
  (the `Cargo.lock` for the crate that ships `esm`, which unbreaks
  EarthSciModels' `rust-cli-inline-tests` job). #369 + #370 were the two
  prerequisites for the pairwise-coupling-library composition method, so that
  refactor is now unblocked. Verified after the merges: the CLI rebuilt from
  main (699d03d16) runs `components/` + `lib/` at **5152/0/0** — the dozen
  other PRs that merged alongside, several tightening validation, break
  nothing here. **Both of the EarthSciAST PRs opened from here have since merged
  (2026-09-17): #400** (`${VAR}` expansion in every §4.7 ref, in all five
  bindings, `b136a9b92`), which let this repo's EarthSciDiscretizations refs
  move from relative sibling paths to `${ESD_ROOT}`, **and #401** (a coupling
  library's refs name roles, not systems, `8ac536bf0`), which cleared the six
  structural errors the `couplings/lib/` libraries reported. #401 picked up
  three review commits (`cb2e955f2` pinning the full §10.10.2 occurrence
  surface against false positives, plus a docs and a Go style commit) that were
  not in the version opened from here. Nothing in EarthSciAST now blocks this
  repo; the only open PR there is #361, which is not ours.
- **EarthSciModels CI repairs, 2026-09-16.** Three failures that hit every
  EqWeFiC PR were that repo's own, failing on its `main` too. **PR #26**
  renames `era5.esm`'s temperature parameter, which was named `t` while the
  same model's equations use `t` as time, so the file failed to load in every
  Python gate run. **PR #27** makes its Julia gate call EarthSciAST's runner
  instead of its own 565-line copy, which compiled each container alone and so
  never received upstream #315 — `radm2.esm` goes **0 pass / 115 errors ->
  825/0/0**, matching the Rust CLI; the same PR loads `OrdinaryDiffEqRosenbrock`
  in `runtests.jl`, without which a document declaring `solver.stiffness =
  "high"` silently stayed on the non-stiff solver and overflowed to -1.7e35
  instead of failing loudly. The third, the `--locked` Rust job, is
  EarthSciAST #371 above.
- **Geometry refresh 2026-09-15.** #15–#21 each gained one commit syncing the
  geometry changes from EqWeFiC 9fca477: `lib/wrf_thermo.esm`'s three new 0D
  templates (all seven PRs carry it), YSU's `p_in`/`p_int_in`/`exner_in`/`ze_in`
  (#18), Dudhia's `p_in`/`exner_in`/`ze_in` (#20), and the RRTM column and
  heating coupling targets (#21). Every branch was drift-free beforehand and
  verifies at its unchanged count (12, 622, 96, 169, 802, 234, 1463); new heads
  24c8bb8, 75ecef3, 2cf75ee, db3c15c, 8ea2241, c294fcd, 9a819d3. The PR bodies
  note that the column-geometry component that fills those targets lives in
  EqWeFiC and is not part of the PRs.
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

- **Couplings repaired for EarthSciAST main @ 3564d09b5 (2026-09-12).** The
  update brought three behaviour changes that the assemblies in `couplings/`
  were written around. (1) A mounted component's inline tests no longer run
  under the coupling (#214), so the whole directory runs as
  `./esm test couplings` and the `--model` selector CLAUDE.md prescribed is
  obsolete. (2) A top-level `models.<k>` `{ref}` mount now runs the full §4.7
  edge pipeline (#261/#298): the leaf resolves in its own scope, folds its
  metaparameters and merges its index sets into the importing document, so the
  document-level redeclarations every RRTM-bearing assembly carried now collide
  (folded literal 103 against the unfolded `NLEV + NBUF`) and were deleted —
  the metaparameters stay, they still close each mount edge's `bindings`.
  (3) A mount-edge `expression_template_imports` injection does not reach a
  component the mounted leaf itself mounts through a nested `subsystems.<k>`
  ref: the injection is consumed by `resolve_template_machinery` -> `lower` ->
  `expand` before the leaf's nested refs are loaded, and `expand` then seals the
  leaf. That blocks `couplings/microphysics_column{,_mixed}.esm` (2 x 40),
  whose Melt/Cold/IceDep/SatAdj stages each carry the shared `sd`/`sat`
  size-distribution subsystem — filed as **EarthSciAST #311** with a four-case
  repro in `data/eqwefic/esm-repro/mount-nested-subsystem/` (a TEST-scope
  injection into the same nested subsystem passes, and both mount forms fail
  alike, which is what makes it a defect rather than a design limit). There is
  no fix inside the coupling document: a document-level and a test-level import
  both fail to reach the sealed leaf, and `TemplateImport` has no way to name a
  scope inside it. **Fixed upstream the same day by EarthSciAST PR #313**
  (opened from here; 21/21 CI checks green): the nested `{ref}` walk now runs
  BEFORE the injection and the §9.6.3 fixpoint at both mount edges, so the
  grandchild no longer arrives after `expand` has sealed the leaf, and nested
  edge `bindings` fold against the leaf's own closed metaparameter environment
  rather than an empty map. The semantics were settled by probe, not argument:
  §9.7.10 defines all three injection forms as extending the target's effective
  scope "as if the target had added those entries to the END of its own
  `expression_template_imports`", and a leaf that DECLARES the library itself
  does lower a rewrite target in a nested `subsystems.<k>` mount — so the
  injected form must too. The contrary sentence ("a parent still cannot
  discretize a grandchild it does not directly mount") is amended by that PR;
  what stays deferred is TARGETING, not reach. Rust only, per the #298 -> #307
  precedent; Julia and Python have the same inversion by inspection.
  **Verified independently in this session**: a CLI built from the PR branch
  runs `./esm test couplings` at **505/0/0 over 52 files**, the microphysics
  pair included. #313 merged 2026-09-14 (with #312 and #315); **the tracked
  `./esm`, rebuilt from main @ 58d21a1a1 on 2026-09-15, reproduces 505/0/0.**
  EarthSciAST #198, #236 and #239 are also closed upstream; #274
  (`table_lookup` on the `esm_problem` carrier) is the one filed gap still open. The alternative to the upstream fix is for
  `hydrometeor_slopes.esm` to import the `input_<x>()` -> `<x>_in` library
  itself (it already declares those parameters), which is a WSM6 component
  decision that collides with the const-array library its own standalone tests
  inject at test scope. All other assemblies are green at their historical
  counts: 325 assertions over 14 files.

- **RADM2 reaction-system tests executed for the first time, 2026-09-12.** The
  three box trajectory tests written 2026-09-06 had never run: the Rust CLI
  reported "no inline tests found" (gap (p), EarthSciAST #206, closed by #251).
  With the runner fixed they failed to LOAD — `H2O` and `M` carry
  `"constant": true`, and esm-spec §7.4 lowers a reservoir species to a
  PARAMETER, so their values belong in `parameter_overrides`, not
  `initial_conditions`. With that fixed only the night regime ran: both daylight
  regimes collapsed the step size to zero at t ~ 6.5e-15 s under the runner's
  default `abstol` of 1e-6, which is meaningless for a state carried in
  molec/cm3 (1e-3 ... 1e19). Isolated by flag (`--abstol 1e-3` alone fixes it,
  `--reltol 1e-9` alone does not) and fixed durably by the document-scoped
  `solver` block (esm-spec §2.2) carrying the box driver's own RTOL 1e-9 /
  ATOL 1e-3. That block forces `"esm": "1.1.0"` (below it,
  `solver_version_too_old`); the bump is required ONLY by the two daylight
  trajectory tests — the 540 rate-constant and 177 tendency assertions are green
  at 1.0.0 with no block. **It does NOT block EarthSciModels PR #24** (checked
  2026-09-12): that repo's CI installs `earthsci_ast` from
  `git+https://github.com/EarthSciML/EarthSciAST.git@main`
  (`.github/workflows/test-esm.yml`, both jobs), never from a release, and main
  has carried esm 1.1.0 since PR #299. No version needs minting; the only real
  condition is that the PYTHON binding on main honours the `solver` block the
  same way the Rust CLI does, which the gate itself will show.
  **The mechanism reproduces the KPP Rosenbrock reference at better than
  rel 1e-7 on every species at every output time** (all 108 pass at 1e-7; 6 fail
  at 1e-8 — the fast NO3/N2O5 night pair; 34 at 1e-9), so the declared rel 1e-5
  has ~100x margin. Gap (q) closed too: a right-hand-side `D(state, t)` now
  resolves, so `models.RADM2Tendencies` exposes one observed
  `d<SPC>_dt = D(RADM2.<SPC>, t)` per variable species — the derivative-test
  shape of §6 — asserted against `Vdot` from `radm2_Fun` at all three box states
  (177 assertions). `constraint_equations` on the reaction system is NOT a
  usable home for these (the observed never materialises); the sibling-model
  form documented in `flatten.rs` phase 5b' is. 102 of the 177 are bit-identical
  to KPP, max residual 9.5e-7 (`dALD_dt`, rural morning), every other species
  within 3.9e-8; the outlier is cancellation, not stoichiometry — ALD's net is
  1.1 % of its largest of 35 terms, and reconstructing the budget from this
  document's stoichiometry times the dump's `A` reproduces the esm value, not
  the Fortran one. Tolerance held at rel 1e-5 because evaluation ORDER differs
  between bindings. Mutation-checked: a 10 % error in `k_R040` fails 34
  assertions across all three surfaces. Note for any index-mapped comparison
  against KPP `A`/`RCONST`: the document has 157 reactions to KPP's 156,
  `R060x` being a shadow reaction of rate `0.9 k_R060` standing in for
  radm2.eqn's negative product `-0.9 OH`.

- **B8 fixed in the fork, 2026-09-15 (55d1c08).** The missing CO2 copy was a
  generator bug in the KPP coupler, not in the Registry: a 2018 fixed-species
  override clobbered CO2's match to the transported `chem(P_co2)` in all 16
  mechanisms where CO2 is `#DEFVAR`. With the fix, the idealized `em_scm_xy`
  `chem_opt = 101` column runs 59 h with vertical mixing on and no Rosenbrock
  aborts (`data/eqwefic/chem_scm_b8`; the buggy run logged 529,584 and every
  species went NaN, worse than the "decays to ~0" first recorded), and call-1
  rate constants are bit-identical to the buggy build. New dumps
  `data/eqwefic/dumps/radm2_wrf_b8` (13 calls × 59 levels) give 767 day and
  night WRF-Chem box states. Each is one 60 s chemistry-only window integrated by
  WRF at RTOL 1e-3, so direct `var_out` assertions need about rel 1e-3; tight
  trajectory and tendency references come from replaying the dumped `var_in` +
  RCONST through the real64 box driver (`chem_build/wrf_dump_to_box.py`).
  Operator splitting with vertical mixing means there is still no multi-step WRF
  species trajectory.

- **Julia could not run any of it, and one bug was why (2026-09-12).** The
  refreshed EarthSciModels PR #24 came back with `radm2.esm` at 0 pass / 115
  errors in `julia-inline-tests` while the Rust CLI ran the same file at
  825/825. Two symptoms — `Unsupported operator: apply_expression_template` (112
  errors, present since the branch was opened) and `Variable 'RADM2.SULF' not
  found in variable dictionary` (3, on the new tendency model) — turned out to
  be **one root cause**: `run_file_tests!` compiled each container through
  `MTK.System(model::Model)`, which wraps that one model in a SYNTHETIC
  single-model `EsmFile` (`flatten(::Model)`, `src/flatten.jl:1313`) and
  flattens that, dropping everything the document supplies around the
  container. `Model` has no `expression_templates` field — the registry is
  document-scoped — so `expand_flattened_refs` returned early on an empty
  registry and a §9.6.4 Option-B `apply_expression_template` reference (the form
  the spec says SURVIVES the load fixpoint for the build to resolve) reached
  `_esm_to_symbolic` with no arm; and sibling components' variables were absent,
  so any equation reading another component died in `_resolve_lowering_var`.
  §6.6 selects which tests RUN, not what the system CONTAINS: Rust and the
  Python gate both build over the whole flattened document, and Julia was the
  1-of-3 outlier (Go and TS execute no inline tests). **Two of this session's
  working assumptions were wrong**: the `SULF` failure has nothing to do with
  product-only species — the minimal repro fails on an ordinary substrate, and a
  plain read without `D` fails identically — and the template gap was not the
  larger of the two, it was the same fix. Filed as **EarthSciAST #314**, fixed
  by **PR #315** (`fix/julia-gaps`, CI 20/20 green incl. julia 1.10/1.12/1.13;
  local `Pkg.test()` 25571 pass / 7 pre-existing broken; the new
  `container_in_document_test.jl` is 16/16 with the fix and 7/5 without).
  Effect on EarthSciModels: `gaschem/radm2/radm2.esm` 0/115 -> **825/0** and
  `urban_canopy/urban_radiation.esm` 0/13 -> **31/0**, both matching Rust. The
  other five files failing in that Julia shard are unrelated and unchanged. Two
  further Julia gaps surfaced and are NOT addressed: `Unsupported operator: and`
  (`fuel_model_lookup.esm`) and an unqualified subsystem-scoped
  `Variable 'T_ww' not found` (`urban_canopy_model.esm`). Repro fixtures kept in
  `data/eqwefic/esm-repro/julia-gaps/`.

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
- **Gap (m) is only HALF resolved (measured 2026-09-07).** The literal-size
  workaround above works only in a document where EVERY contributor of an axis
  arrives FOLDED. An index set that a mounted component obtained through an
  `expression_template_imports` edge (YSU, Dudhia and slab all take `lev` from
  the EarthSciDiscretizations `column_nonuniform_1d` grid that way) is
  re-exported to the mounting document UNFOLDED, as `size: "NLEV"` under the
  CHILD's metaparameter name — so the document must declare a metaparameter
  literally called `NLEV` and must NOT declare `lev` itself. A folded literal
  and that symbolic form are not deep-equal, so the two cannot meet: mounting
  any WSM6 stage that carries the `sd`/`sat` subsystem (melting_freezing,
  cold_accretion, ice_deposition, saturation_adjustment) beside YSU fails with
  `[template_import_index_set_conflict] models.Pbl` EVEN AT THE SAME SIZE
  (probed at NLEV = 40 and 59; also with `NLEV` threaded into the nested
  subsystem edge, with the document declaring `lev` symbolically or literally,
  and with a renamed grid injected at the mount edge). Per esm-spec §4.7 a
  subsystem ref's registry should be "fully concrete when it splices in", so
  this is a Rust-CLI deviation from the spec, not an authoring mistake — filed
  2026-09-07 as EarthSciAST issue #236, with the probe documents kept in
  `data/eqwefic/phase3_probes/index_set_merge/`. The discriminator is HOW the
  contribution arrives, not its size: YSU beside `hydrometeor_slopes.esm`
  mounted at TOP level loads, while YSU beside `cold_accretion.esm`, which
  mounts that same file as a NESTED subsystem, does not. Gap (n) — two
  components needing the SAME axis name at DIFFERENT sizes, the 4-cell soil
  column and the 59-cell atmospheric column — is NOT filed separately: it is
  already #198 item 4 and EarthSciAST **PR #213** (`index_set_rename` on the
  mount edge) is open against it, quoting this exact soil/atmosphere case.
  A THIRD symptom of the same root was found while checking that and filed as
  EarthSciAST issue #239: mounting a component in a document that merely
  declares a metaparameter of the same name silently overrides the component's
  own default, with NO `bindings` on the mount edge and no diagnostic. Slab
  (4 soil layers, `NLEV` default 4) mounted alone under a document `NLEV` of 59
  loads and then fails 8 of its own 84 assertions with WRONG NUMBERS — the
  deep-reservoir assertion that should be exactly 0 reads 0.0514, and the
  tendency is off by ~600x — with nothing pointing at the grid. Declaring
  `bindings` explicitly and omitting it are indistinguishable. Every column
  component here names its size metaparameter `NLEV`, so any assembly mounting
  two at different resolutions silently resolves both to one number.
  PR #213 grants `prefix`/`rename` at a §4.7 subsystem edge (esm-spec §9.7.7 grants
  those three fields to `expression_template_imports` only;
  esm-schema.json `$defs/SubsystemRef` has `ref`/`model`/`reaction_system`/
  `bindings`/`expression_template_imports` and nothing else).
- **SCM physics-suite assembly done 2026-09-07.**
  `couplings/scm_physics_column{,_night}.esm` mount SfclayRev, YSU, the five
  RRTM stages and DudhiaSW (8 mounts, 133 `variable_map` edges) at two regimes
  the earlier `physics_column` pair does not cover — step 1410 (midday, the
  run's strongest forcing: gsw 550 W/m², hfx 174, h 926 m) and step 900 (deep
  stable night: gsw exactly 0, hfx −19, h 263 m) — and assert WRF's
  DRIVER-LEVEL tendency vector, including `rqcblten`/`rqiblten` for the first
  time. 2 × 42/42 with `./esm test --model ScmPhysicsColumn <file>`. The seven
  dumps of a step are merged into one namespace by `tools/merge_dumps.py`
  because `tools/fill_tests.py` admits two per test. Residuals (zero-tolerance
  run): surface-layer outputs ≤ 1.9e-6 relative by day and ≤ 2.5e-5 at night
  (the stable regime's real32 floor), h 6.4e-7 / 3.5e-6, K_h/K_m Linf 2.3e-4 /
  2.3e-5 m²/s, longwave fluxes ≤ 6.2e-4 W/m², dthdt_lw and rthraten ≤ 5.9e-9
  K/s, gsw/glw/olr ≤ 1.4e-6 relative, PBL tendencies at the
  implicit-vs-instantaneous gap (30 % of the column max by day, 2 % at night;
  asserted loosely as in `physics_column`). Every tolerance is the smallest of
  {1, 2, 5}×10^k ≥ 5× the measured residual. Slab and WSM6 are not in the
  document — see gap (m) above, not a physics limitation.
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
  differ (see component descriptions). The ESD rules for the DEFAULT
  `fire_upwinding = 9` landed as EarthSciDiscretizations PR #37 (2026-09-07,
  branch `cartesian2d-levelset-weno` off ESD main e7defb2, not stacked on
  #34/#35/#36): 7 rules on `cartesian_uniform_2d` — the compound hybrid
  WENO5/ENO1 ∣∇ψ∣ (`weno5_eno1_norm_D2_extrapolate_bc`, priority 20), the two
  per-axis hybrid derivatives, the two Godunov-upwind per-axis quantities
  max(D⁻,0) − min(D⁺,0) the front normal is built from, and the two centred
  second derivatives of the viscosity term — plus 9 stencils and the new
  boundary tag `bc:extrapolate_linear` for fire_behavior's one-cell guarded
  linear-extrapolation halo max(2a − b, a, b). Two MMS problems gate the
  order: observed L2/Linf 1.00 for the Godunov case and observed L2 1.47
  (expected 1.5) for the hybrid, the latter necessarily in ACCUMULATOR form
  (u frozen at the manufactured g, one state w integrating
  w_t = S(∣∇g∣ − ∣∇u∣_h)) because the hybrid's ENO1 branch is a minmod, not an
  upwind discretisation for the transported error, so the coupled forced-steady
  form the Godunov sibling uses does not converge and stalls Tsit5.
  `level_set_tendency_wrffire.esm` now covers both regimes, 203/203 Rust.
  Still missing: `fire_upwinding` 0–3, 5–8 and 10 (the blending zone), the
  `tbound` CFL reduction (a solver step limit, not a discretisation) and
  reinitialisation.
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
  emissions and the Madronich photolysis / WRF-Wesely couplings. (The earlier
  note here -- "FastJX covers 10 of 21 j inputs; Wesely ~20 species" -- described
  reusing the EarthSciModels `fastjx/*` and `wesley_dry_gas.esm` components, and
  both reuses are now withdrawn: Fast-JX is not WRF-Chem's scheme at all, and the
  existing Wesely is an AtmosphericDeposition.jl migration that differs from
  `module_dep_simple.F` on every one of the 19 overlapping species. Both were
  replaced by fresh transcriptions: `gaschem/madronich/*` 409/409 and
  `atmospheric_deposition/wrf_wesely/*` 2048/2048.)
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
  deferred row), and every term of loops 220/2000/4000 is expressible.
  Warm and mixed-phase cloud optics done 2026-09-07, closing the last gap:
  the liquid, rain and snow coefficients (ABCW = 0.144, ABRN = 0.330e-3,
  ABSN = 2.34e-3) were transcribed but unverified because the em_scm_xy
  reference case is ice-only. New case `data/eqwefic/qss_rad`: the idealized
  supercell em_quarter_ss (WSM6, Weisman-Klemp sounding, 41 x 41, dx = 2 km,
  dt = 6 s, 1 h) with `ra_lw_physics = 1` and `ra_sw_physics = 1` (WRF aborts
  on LW without SW) -- but NOT on its stock grid: esm inline tests have no
  metaparameter overrides (`esm-schema.json` `$defs/Test`), so a new regime must
  land on the component's NLEV = 59 / NBUF = 44 column, which pins e_vert = 60
  and p_top in [174, 178) hPa; `ztop = 12800 m` gives p_top = 177.07 hPa
  (NLAYERS = 103, the SCM's exactly), found by running `ideal.exe` at four
  ztop values. Second obstacle: the hooks dumped only tile column `its`, the
  western domain edge, never the storm; fork commit 5cb62da adds
  `esm_dump_col`/`ESM_DUMP_I` to `module_esm_dump` and an `i_col` record to the
  RRTM hook. Three RRTMColumn regimes (rrtm_column 480 -> 691 assertions) and
  one RTRNSweep regime (196 -> 224), rrtm_lw 1212 -> 1451: a pure WARM LIQUID
  cloud (call 61 = radiation call 2, j = 20, i = 20; qi = qr = qs = 0
  everywhere, TAUCLOUD 1.26-6.43 in layers 10-15 = ABCW CLWP alone), a
  MIXED-PHASE storm column (call 268, j = 22, i = 20; 57 cloudy layers,
  TAUCLOUD <= 18.74, with layers 1-10 pure rain, 11-19 liquid + rain, 22-33 all
  four terms, 39-41 pure ice and 51-54 snow-dominated), and an ICE-AND-SNOW
  ANVIL (call 396, j = 27, i = 40; no liquid and no rain, layers 17-30 pure
  snow with TAUCLOUD 0.047-0.465 = ABSN PIWP alone, layers 51-53 pure ice).
  The RTRNSweep regime is the opaque limit the SCM never reaches:
  ABSCLD = 1 - exp(-1.66 TAUCLOUD) is 1 to 3e-14 in 25 layers (SCM maximum
  0.72), glw 451.1 against a clear-sky 405.7 W/m^2, olr 160.5 against 259.5,
  HTR - HTRC = 41.0 K/day. Residuals (real32, rel 1e-5 contract): TAUCLOUD Linf
  8.5e-7 / 3.1e-6 / 6.4e-8 on column maxima 6.43 / 18.74 / 1.04 (<= 1.4e-7
  relative at every spot-check layer) with tolerances abs 1e-5 / 3e-5 / 1e-6,
  10-16x those floors and 1e-6 of each maximum; colh2o needed abs 1.5e-3
  instead of the SCM's 2e-4 because the tropical column holds 131 cm^-2; the
  RTRN fluxes agree to Linf 2.2e-4 W/m^2 and HTR/HTRC to 2.6e-4 K/day.
  Mutation-checked at +10 % on each coefficient: ABCW fails 10 assertions (7
  liquid, 3 mixed), ABRN 6 (5 mixed, 1 anvil), ABSN 12 (8 anvil, 4 mixed) and
  ABICE 22 (8 anvil, 5 mixed, 9 spread over the eight SCM regimes); the SCM
  regimes stay green under ABCW/ABRN/ABSN and the RTRNSweep tests never move,
  because they take TAUCLOUD as a dumped input.
  Cross-checks: re-running the archived SCM case with the rebuilt binary
  reproduces 119 of the 122 variables of `scm_ref2_120.json` BIT-FOR-BIT (the
  three that differ are uninitialised above LAYTROP, B9); the four dump passes
  with different `ESM_DUMP_I` give bit-identical row arrays and a wrfout whose
  md5 matches the no-dump scan run; and each dumped 1-D column equals the tile
  array at the selected i exactly (p3d only to 2.9e-5 hPa, the real32 Pa -> hPa
  division). New Fortran findings: B9 (SETCOEF leaves SELFFAC/SELFFRAC/INDSELF
  undefined above LAYTROP; unread by TAUGB, but it is why three dump variables
  are not build-reproducible) and N62 (MM5ATM takes QG and never uses it -- no
  graupel term and no ABGR in TAUCLOUD, while the supercell carries up to
  1.3e-2 kg/kg of graupel). Nothing in the cloud-optics block is now untested
  except the graupel that WRF itself discards. The nine coupled documents that
  mount RRTM were re-run unchanged after the additions: `rrtm_lw_column*` 5 x
  25/25, `radiation_column{,_night}` 2 x 30/30, `physics_column{,_night}` 2 x
  40/40. No coupled supercell regime was added -- the RRTM coupling documents
  have no `.fill.json` sidecar and their mount-edge libraries would have to be
  rebuilt by hand; the end-to-end chain is still covered by the SCM regimes.
  EarthSciModels PR #21 refreshed the same day: branch `eqwefic/rrtm-lw` kept
  its three commits off 7259fd8 (base `main`, duplicated `lib/`), the tip
  amended with the new content, all three messages rewritten to drop the
  `Co-Authored-By: Claude` trailer (upstream b7912bc forbids it) and keep
  `Claude-Session:`, force-pushed with `--force-with-lease`; verified 1463/1463
  from a clean worktree at `code/EarthSciModels-rrtmlw`, a sibling of
  `EarthSciDiscretizations` so the `../../../../` refs resolve.

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

    **What the reference case prescribes vs. computes (established 2026-09-07).**
    `data/eqwefic/scm_ref` runs with `scm_force = 0`, and `force_scm`'s third
    executable statement is `IF (scm_force .EQ. 0) return`
    (`dyn_em/module_force_scm.F:183`), so **no** large-scale advective tendency,
    subsidence, geostrophic-wind forcing or nudging is applied, and the
    `scm_th_adv` / `scm_qv_adv` / `scm_wind_adv` / `scm_vert_adv` flags are read
    only after that early return. `force_ideal.nc` is opened and "processed"
    (auxinput3 has no `scm_force` guard in `share/mediation_integrate.F:350`)
    but nothing is read from it: the forcing fields sit behind
    `package scmopt scm_force==1` (`Registry/Registry.EM_COMMON:3322`), so they
    are never appended to `grid%tail_statevars`, which is the list
    `share/input_wrf.F:1296` walks. The column is therefore free-running under
    its physics — but it is *not* physics-only. A faithful single-column
    EqWeather must take as **prescribed input**:
    - the initial sounding and the stretched eta grid
      `eta(k) = 1 - (e^((k-1)/40) - 1)/(e^(59/40) - 1)`
      (`dyn_em/module_initialize_scm_xy.F:322-331`), `p_top = 17761.31 Pa`,
      `mu = 79438.7 Pa` constant for the whole run. Note theta below 200 m is
      the constant-extrapolated 286 K (`module_init_utilities.F:66-71`), not
      the sounding's 288 K surface line;
    - lat 37.6 / lon -96.7 / 1999-10-22 19:00 UTC, giving `f = 8.899e-5`,
      `e = 1.1554e-4` and a real diurnal zenith angle (`radt = 0`, so radiation
      runs every step; `module_radiation_driver.F:1119-1123`, `:3514-3541`);
    - **the geostrophic wind (3, -9) at every level**, supplied implicitly by
      `pert_coriolis` through `u_base`/`v_base`
      (`dyn_em/module_em.F:747`, `module_big_step_utilities_em.F:3854-4171`).
      Plain Coriolis on the total wind gives the wrong trajectory;
    - **the Rayleigh sponge** above `ztop - zdamp` ~ 7 km, relaxing u, v, w and
      theta to the initial sounding with `dampcoef = 0.003`
      (`module_em.F:933-943`, `module_big_step_utilities_em.F:5836-6110`).
      This is the largest non-physics term and cannot be dropped;
    - constant surface properties from the USGS cat-2 **WINTER** row of
      `LANDUSE.TBL` (Oct 22 = day 295, `module_physics_init.F:1959-1960`):
      albedo 0.20, emiss 0.92, z0 0.05 m, mavail 0.60, thc 0.04, so
      `capg = 5.9114e7 * thc` (`phys/module_sf_slab.F:387`);
    - soil grid `dzs = 0.01/0.02/0.04/0.08/0.16`, `zs = 0.005/0.02/0.05/0.11/0.23`
      (`share/module_soil_pre.F:1078-1126`, which *discards* `input_soil`'s
      profile under `sf_surface_physics = 1`), with `TSLB(5)` frozen as the
      lower BC (`module_sf_slab.F:470`) and `TMN` a pure IC (its `HM` term is
      used only in the `num_soil_layers == 1` branch, `:455-457`).
    Everything the dynamical core contributes to *transport* is identically
    zero: `ww == 0` exactly for a horizontally uniform periodic column
    (`calc_ww_cp`, `module_big_step_utilities_em.F:640-782`), so there is no
    coordinate vertical advection; horizontal advection and diffusion vanish on
    the 2x2 periodic mass grid; and the `km_opt = 2` TKE closure runs but is
    **inert**, because `vertical_diffusion_2` sits behind
    `IF (config_flags%bl_pbl_physics .eq. 0)`
    (`dyn_em/module_first_rk_step_part2.F:1023-1025`) and YSU is 1. So
    EqWeather-SCM = the six physics schemes + a Rayleigh-sponge component + a
    perturbation-Coriolis component + a 1-D fixed-mu column-geometry component
    (w != 0 even though ww == 0). No advection operator, no horizontal
    diffusion, no TKE closure.

    **Reference trajectory.** `data/eqwefic/scm_ref_traj/` (`.npz` + CSVs), read
    by `tools/scm_ref_traj.py`. **60 hourly frames, t = 0 to 59 h**; the column
    is bitwise horizontally uniform at every level and time, so a single column
    is exact, not approximate. Dry and cloud-free in the PBL (no precipitation
    for 59 h). Signal against which the acceptance criteria must be read:
    lowest-level Tk swings 7.1 K over the first day with a +1.65 K net drift,
    TSK 11.5 K; qv at the lowest level rises **+93 %** (0.00250 -> 0.00483) in
    24 h, almost all of it in the first 6 h. So 0.5 K is ~7 % of the diurnal
    amplitude and 5 % qv is ~10 % of the change actually produced — meaningful
    but forgiving. Caveats: every 2 m / 10 m diagnostic is zero in the t = 0
    frame (compare at the lowest model level, z = 25.59 m); `SMOIS`/`SH2O` are
    identically zero (slab has no prognostic soil moisture); `GRDFLX` is
    identically zero for all 59 h, so the surface energy budget cannot be closed
    from the history.

    **The test shape works (settled 2026-09-07 by probe, not by reading the
    schema).** A 24 h trajectory assertion **is** expressible as an esm 6.6
    inline test and the Rust CLI integrates it accurately. `time_span` may be
    `{start: 0, end: 86400}` and assertions may name any `time` in that span;
    the runner integrates once per *test* (not per assertion, measured) and
    interpolates. Cross-validated against SciPy LSODA at rtol 1e-12 on a
    59-level column: BDF agrees to 1.5e-7 K at 24 h, ERK to 3e-8 K, SDIRK to
    1e-6 K — six orders below the 0.5 K criterion, so the integrator is not the
    limiting factor. Confirmed working at a late time: `coords` point sampling;
    `reduce: mean/integral/max/min`; `reduce: L2_error` / `Linf_error` against
    both an inline `const` array and `reference: {type: "from_file"}`; `t` as
    simulation time in an equation (so the diurnal cycle is expressible); an
    `ic` op equation for a non-uniform 59-level initial profile (array ICs
    cannot go through `initial_conditions`, which the schema restricts to
    scalars). Negative controls fail as they must.
    Therefore the plan's criteria map onto assertions directly:
    **state within 0.5 K at 24 h** = `reduce: "Linf_error"` with a `from_file`
    profile and `tolerance: {abs: 0.5}` (verified: perturbing one level of the
    reference by 0.6 K flips it to FAIL); **5 % qv** = `reduce: "L2_error"`
    (which is *relative* L2) with `tolerance: {abs: 0.05}`. The one criterion
    that does **not** map is "tendencies within 1e-3 relative RMS" *over the
    trajectory*: `reduce` is spatial-only at a single `time`, and there is no
    temporal reduction. Write it as N per-hour `L2_error` assertions on the
    tendency (relative RMS over the column at each hour) instead of one.

    **Prototype (2026-09-07).** `components/land_surface/slab/slab.esm` turned
    out **not** to be steppable: it has no `D(T_s, t)` equation — `T_s` is
    algebraic (`T_s = input_T_s()`) and `dTs_dt` is a pure observed, so the
    whole component is an instantaneous-tendency calculator. (Only `ysu.esm`,
    `rrtm_lw_heating.esm` and `dudhia_sw.esm` carry `D(..., t)` today.) The
    prototype was therefore built as a standalone steppable soil column
    (`data/eqwefic/phase3_probes/proto_soil.esm`, not committed: it duplicates
    slab's physics and so must not enter `couplings/`): WRF's soil heat equation
    over layers 2-4, top node prescribed from the reference hourly TSK,
    bottom pinned at TMN, stepped 24 h and asserted against the run's own TSLB.
    It is **green at 0.5 K with ~12x margin**; the achievable tolerance is
    **abs 0.05 K** (passes at 0.05, fails at 0.02), and the residual is
    dominated by hourly sampling of the prescribed forcing, not by the physics
    or the integrator — max error is 0.04 K at 24 h but 0.45 K mid-window during
    rapid transitions. **So assert at fixed checkpoints, not at every hour**, or
    dump the forcing at model resolution.

    **Three upstream blockers found (repros in `data/eqwefic/phase3_probes/`).**
    (s) **DEFERRED 2026-09-12 (user decision): 3-D simulations will be solved
    by a different method, so the inline-test runner's iteration cap is not on
    the critical path and no upstream issue was filed.**
    **Update 2026-09-15:** the 24 h single-column test runs in the inline-test
    runner, so the cap does matter there. EarthSciAST **PR #362** makes Rust's
    `SolveOptions::maxiters` an `Option` defaulting to no cap (matching Python);
    an explicit cap keeps its old meaning. Verified from source for the PR:
    SciML's `maxiters` caps integrator time-loop iterations for explicit and
    implicit methods alike (SciMLBase `integrator_interface.jl:587`; default
    1,000,000 for adaptive methods, OrdinaryDiffEqCore `solve.jl:52`), and the
    Newton iteration inside implicit steps has its own `max_iter`; Rust already
    counted steps the SciML way, only its 10,000 default was wrong. Python counts
    right-hand-side evaluations instead, a spec gap noted in the PR. Original finding:
    `SolveOptions::maxiters` is pinned at its 10 000 default by the
    inline-test runner (`pkg/earthsci-ast-rs/src/bin/esm.rs:3719`, `..Default::default()`)
    and there is no CLI flag and no document field for it. A *scalar* ODE with
    600 s structure already exceeds it between a 54 000 s span (passes) and a
    57 600 s span (`solver retcode MaxIters`) — `repro_maxiters.esm`. Stiffness
    alone is fine (BDF absorbs a 59-level column at surface-layer stiffness over
    24 h), so the risk is *non-smooth* RHS structure — WSM6 saturation
    adjustment, radiation switching at sunrise — not stiffness.
    (t) `table_lookup` over a top-level `function_tables` entry does not
    evaluate: `unevaluable_operator` in the array interpreter, and **silently
    NaN** in the scalar one, although esm-spec 4.2 says it lowers at load time
    to `interp.linear`. The primitive itself works, so prescribed time-series
    forcing must be spelled `{"op": "fn", "name": "interp.linear",
    "args": [values, axis, "t"]}` for now — `repro_table_lookup.esm`.
    (u) A semantically inert alias between two `unknown`s breaks the solve
    once the system passes a size threshold — EarthSciAST issue #234, filed
    2026-09-07. Two documents identical except for one extra bare-alias
    equation (`Tsfc = Tsfc_raw`): `delta_direct.esm` is green, `delta_split.esm`
    errors with `diffsol: Exceeded maximum number of nonlinear solver failures`
    on all three solvers. Characterised by sweeping `lev`: both pass at sizes
    3-58, the split form fails from 59 upward, and the direct form is still
    green at 70 (141 unknowns, MORE than the failing split system's 120), so
    neither the alias nor the size is sufficient alone. NOT a long-horizon
    effect: the failure time is identical (t = 39.763 s) at every size and
    every span, including a 60 s span — the earlier "fragile over long spans"
    reading was wrong. This is still the biggest risk to 3.1, because both of
    this repo's own conventions produce the shape: components are factored and
    imported by reference, and the `<x>_in` + mount-edge + `variable_map`
    coupling pattern is a chain of pass-throughs over a 59-level column. Every
    inline test here asserts instantaneous derivatives over a ~1 s span, so
    none of them exercise it; the first long single-column integration would.

    **(u) RESOLVED 2026-09-12.** Fixed upstream by PR #255 (`8e115ce44`, merged
    at `60da85338`), which was aimed at #207 and never linked to #234. Bisected:
    the probe errors at `ebc432873` with the reported message and failure time
    (t = 39.763 s) and passes at the next commit. The cause was tape export
    ordering — an observed whose body is a bare `Expr::Variable` emits nothing
    into its home chunk, so its `Export` ran after the `Fallback` of a later
    rule reading it and the reader took the preallocated `0.0`; in the probe
    that ghost zero lands in the surface relaxation term and drags cell 1 toward
    0 K until Newton gives up. On current main (`3564d09b5`)
    `delta_direct`/`delta_split` are green at all three spans, a `lev` sweep of
    55-70 is green in both forms on all three solvers, and `ESS_TAPE_CHECK=1`
    reports no tape-vs-oracle divergence on any probe. **The pre-fix corruption
    was size-independent**: `lev = 58` "passed" only because the corrupted ODE
    happened not to defeat Newton, and the N = 4 minimal case diverges on the
    tape check too — so any long integration run on a CLI older than
    `60da85338` is suspect, not merely the 59-level ones. An earlier attempt on
    the abandoned branch `claude/issue-234-alias-elimination` diagnosed this
    correctly (naming `compute_exports` in `tape/lower.rs`, the function #255
    rewrote) but landed no fix. Nothing upstream guards the shape: #255's tests
    are `coupled_const_array_fold.rs` and `scalar_param_array_default.rs`,
    neither a bare alias feeding a per-cell rule — the silent form this repo's
    factoring and `<x>_in`/`variable_map` conventions generate. **EarthSciAST
    PR #312** (opened 2026-09-12 from here) closes that hole:
    `tests/alias_observed_export_order.rs`, four cells and four tests, pinning
    the analytic zero tendency, the taped-vs-oracle path bitwise (the
    `ESS_TAPE_CHECK=1` invariant without the env var), and that removing the
    inert alias moves no bit. Verified to FAIL at `ebc432873` (4/4 failed, cell
    0 tendency -0.288 K/s against 0) and pass on main. No cross-binding fixture:
    the defect is in the Rust array runtime's tape lowering, which no other
    binding has. Phase 3.1 is
    unblocked on this axis; the remaining Phase 3 blockers are (s) `maxiters`
    (still unfiled upstream) and (t) `table_lookup` on the `esm_problem`
    carrier (EarthSciAST #274).

    **Column dynamics done 2026-09-15 (derivative level).** EqWeather-SCM's
    non-physics terms are FIVE components, not the three listed above: the
    Rayleigh sponge, perturbation Coriolis, vertical momentum + geopotential
    (`pg_buoy_w`, dφ/dt = g w), the column geometry / equation of state
    (`calc_p_rho_phi` + `phy_prep`), and Earth curvature, which `rk_tendency`
    applies even to this Cartesian column (N66, 1.4e-5 m/s² in w). All live in
    `components/atmospheric_dynamics/wrf_arw/` with a shared
    `column_staggering.esm` template library: 185 assertions against the
    `dyn_column` dumps of fork commit 83347b2, each mutation-checked except the
    sponge's above-top base-state extrapolation (u, v, θ are constant above 4 km
    in this sounding). Coupled as `couplings/scm_dynamics_column{,_night}.esm`
    (2 × 8). `w_damp` is inert (`w_damping = 0`). WRF's real32 p′ rounding drives
    a w tendency of up to half the physical one (N67), so w cannot be matched to
    WRF beyond ~3e-4 m/s²; θm, u and v tendencies match to ≤1e-10. The sponge
    relaxes θm toward a DRY reference (N64). The physics suite cannot be mounted
    as a unit (a top-level `{ref}` to a multi-model document needs a model
    selector), so EqWeather-SCM must re-mount all nine physics components and
    their ~167 edges. **24 h acceptance still needs:** (1) a state document with
    D(u, v, w, φ′, θm, qv, soil T, TSK) = physics + dynamics tendencies, including
    WRF's `use_theta_m` conversion of physics θ/qv tendencies into θm; (2)
    geometry → physics coupling-target parameters (p, π, z, dz8w, ρ, p8w, t8w,
    p_hyd) in sfclayrev, YSU, slab, Dudhia and the RRTM stages, replacing the
    dumped profile libraries — EarthSciModels PR #16–#21 territory; (3)
    time-dependent solar geometry for Dudhia and RRTM; (4) a float64-balanced
    initial state, since WRF's real32 state starts an acoustic adjustment in the
    esm column; (5) checkpointed windows of ≤6 h restarted from wrfout frames,
    because the pinned `maxiters` (gap (s)) is likely to be exceeded by
    undamped vertical acoustic modes (WRF damps them with `epssm`) plus the
    physics' regime switching over 24 h.

    **Solar geometry and the θm conversion done 2026-09-15 (derivative level).**
    Two decision-independent pieces of the 24 h list now exist as components.
    `WRFMoistThetaTendency` (`components/atmospheric_dynamics/wrf_arw/moist_theta_tendency.esm`)
    transcribes `conv_t_tendf_to_moist` as the product rule on θm = θ(1 + ε qv),
    with 28 assertions over 7 steps; WRF converts all 59 layers, and the
    conversion changes the physics θ tendency by up to 5.8e-4 K/s, almost
    entirely through the vapour term. `WRFSolarGeometry`
    (`components/atmospheric_radiation/wrf_solar/`) transcribes radconst +
    calc_coszen as algebraic functions of xtime and julian, with 63 assertions (13
    SCM regimes including the sunrise and sunset crossings and the 1440-minute
    wrap, plus 6 synthetic regimes from a kernel driver that extracts both
    subroutines verbatim), each physical term mutation-checked. Only DudhiaSW
    takes solar geometry (`csza`, `solcon`); RRTM LW takes none. `lib/solar.esm`
    cannot stand in for it: its formulation differs from WRF's by 0.022 rad in
    declination. Neither component is wired into the physics or the SCM coupling
    yet. Two REAL*4 precision notes set the crossing tolerances (N68, N69). Fork
    1b4b556 adds the `theta_m_conv` hook and `kernels/solar_driver.F90`.

    **Geometry inputs wired, 2026-09-15.** The physics components now have
    additive coupling-target parameters for geometry: YSU
    `p_in`/`p_int_in`/`exner_in`/`ze_in`, Dudhia `p_in`/`exner_in`/`ze_in`, RRTM
    column `p_in`/`p_e_in`/`T_in`/`T_e_in`/`dz_m_in`, RRTM heating
    `p_e_in`/`exner_in`; every standalone count is unchanged. sfclayrev and slab
    needed none, because their lowest-layer inputs are plain parameters that a
    coupling fills by indexing the geometry column (`Geo.p[1]`), and Dudhia's
    `csza`/`solcon` are mapped directly from `WRFSolarGeometry`. **WRF passes YSU
    (p2di) and RRTM (p8w) the hydrostatic interface pressure** — `psfc` at the
    ground and `p_top` at the top — not the extrapolated `p8w` of `phy_prep`, so
    those inputs come from `p_hyd_w`; `p8w` puts RRTM's top interface off by
    64.7 Pa. `column_geometry.esm` gained `z_w_agl` and was re-factored onto 0D
    templates in `lib/wrf_thermo.esm` (reused `exner_function`,
    `temperature_from_theta`; new `wrf_equation_of_state`,
    `moist_theta_to_dry_theta`, `moist_density_from_inverse_dry_density`), per
    the user's 2026-09-15 rule that 0D processes are authored separately and
    composed into columns; `moist_theta_tendency.esm` still inlines
    θ = θm/(1 + ε qv) and should apply the new template. The wiring is
    demonstrated by `couplings/geometry_radiation_column{,_night}.esm` (2 × 14)
    and `geometry_surface_pbl_column{,_night}.esm` (2 × 41). The latter's
    sfclayrev → slab → YSU edges and assertions were copied from
    `scm_physics_column{,_night}.esm`, the duplication the assembly-composition
    analysis (`data/eqwefic/design/assembly_composition.md`) proposes to remove
    with pairwise coupling libraries. By day, geometry's p[1] is 0.30 Pa below
    the surface-layer kernel's `p1` (N50 class), which moves surface fluxes by up
    to 7e-5 relative; 14 day-time surface-layer tolerances are 3× measured. State
    inputs (θ, T, u, v, q) and `psfc` remain prescribed pending that design.
    `wrf_solar_geometry.esm` must be mounted under its own model name: its
    model-qualified constant references (`WRFSolarGeometry.wrf.degrad`) are not
    rewritten at a mount edge with a different key, and it is the only component
    written that way even though CLAUDE.md's Constants bullet prescribes the
    qualified form.

    **Composition method decided 2026-09-15: pairwise coupling libraries.**
    Assemblies are NOT mounted as units. What crosses a slot boundary depends on
    the specific pair of schemes (MYNN exchanges TKE, boundary-layer clouds and
    species mixing that YSU does not; Noah takes `chs`/`chs2`/`cqs2` where slab
    takes `flhc`/`flqc`; RRTMG reads effective radii, aerosol optics and PBL
    clouds that RRTM LW does not; Thompson with fire requires an active PBL and a
    fire-emissions aerosol source), so a unit drawn around "physics" or
    "dynamics" has an interface that changes with its contents. Instead, the
    edges between each pair of components live in a hand-authored coupling
    library (`coupling_roles` + `coupling_import`, esm-spec §10.9–§10.10),
    imported by every document that uses that pair: today's subassembly
    documents become the per-pair tests, `scm_physics_column.esm`'s ~167 edges
    are split into pair libraries, and EqWeather-SCM becomes ~14 mounts plus
    ~15 imports and a few inline edges; swapping a scheme swaps one mount and
    its pair imports. INVENTORY records which pair libraries exist. Full
    analysis: `data/eqwefic/design/assembly_composition.md`. Prerequisite
    upstream fixes, both approved 2026-09-15: (i) `coupling_import` refs
    resolved against the working directory rather than the importing document
    (`flatten.rs:633` used `CouplingImportOptions::default()`, base path ".") —
    **fixed by EarthSciAST PR #369** (all five bindings record the document's
    base at load time and prefer it, the authored `ref` still round-trips
    verbatim per §10.10.3; 18/18 Rust coupling tests including a new
    working-directory regression test, and a CLI built from the branch runs the
    probe document from any directory where it previously failed); and (ii) the
    validation gaps — **fixed by EarthSciAST PR #370** (stacked on #369):
    `validate` now expands each import and checks the role-bound edges,
    reporting each finding at the import's own pointer, and the spec states the
    idiom for a REQUIRED import — declare the coupling-target parameter without
    a `default`, so an omitted or mis-bound import fails loudly when the problem
    is built instead of running with a placeholder (measured: dropping the
    import left the run green with a silently wrong 0.0; dropping the default
    too gave `Invalid parameter`). Rust and TypeScript; the other three
    bindings iterate the un-expanded coupling the same way and are a follow-up.
    20 coupling tests and 733 library tests pass.

    **EqWeather-SCM state document, 2026-09-17: the dynamics half runs, and the
    binding constraint is cost, not `maxiters`.** Of the five prerequisites
    listed above, (3) is done, (2) is confirmed already done, (4) is built and
    measured and turns out not to be the lever, (1) is half built, and (5) is
    deliberately not built because the reason given for it was wrong.

    - **(1) State document — dynamics half only.**
      `couplings/eqweather_scm_dynamics.esm` (45 edges, 31 equations, 15
      assertions) closes D(u), D(v), D(w), D(φ′), D(θm), D(qv) over the five
      WRF-ARW column-dynamics components and feeds the ODE state BACK into
      geometry, vertical momentum, Coriolis, curvature and the sponge, so
      heights, eta weights, the equation of state and the sponge reference
      profiles are recomputed at every RHS evaluation. Mechanism: 14 additive
      coupling-target parameters on the five components, lowered from the
      `input_<x>` rewrite targets by five mount-edge libraries
      `couplings/tests/scmstate_couple_*_inputs.esm`. Base-state fields (eta,
      phb, alb, pb, u/v/t/z_base) stay on the dumped library — they are
      genuinely time-invariant in this run. Closing the loop changes no bit of
      the RHS: the derivative test reproduces `scm_dynamics_column.esm`'s
      references at its tolerances, plus dφ/dt = g w to 1e-12.
      `WRFMoistThetaTendency` is now mounted in `couplings/scm_physics_column.esm`
      and driven by that assembly's own `rth_phys_sum` and YSU `rqvblten`; the
      assembly exposes `dtheta_m_dt_phys`, which is what D(θm) must take,
      because the physics suite's own output is a DRY θ tendency.
      **Not done:** the two halves are not mounted in one document.
    - **(2) Geometry → physics: verified complete.** Audited both geometry
      couplings; nothing geometric still comes from a dumped profile library.
      What is still dumped is *state* (θ, u, v, q, cloud, psfc), which is item
      (1)'s remaining work.
    - **(3) Solar geometry: done.** `couplings/solar_diurnal.esm` supplies the
      clock (`xtime = xtime0 + t/60`, `julian = julian0 + t/86400`), mounts
      `WRFSolarGeometry` behind it and reads coszen/solcon off ONE 58 h
      integration at the 13 WRF steps the component's regime tests pin,
      including the three sunrise/sunset crossings: 26/26 at the component's
      own tolerances, not loosened. With every tolerance set to abs 0, 20 of 26
      agree with WRF to better than 1e-6 relative; the six that do not are the
      crossings, largest absolute residual 8.3e-7 (xtime 2510), i.e. WRF's
      REAL*4 zenith cosine (N69, 7.7e-7), not the clock — the REAL*4 `julian`
      differs from `julian0 + t/86400` by at most 2.0e-5 d over 58 h, worth
      2.7e-7 in coszen.
    - **(4) Float64-balanced initial state: built, and it is not the lever.**
      Recipe and numbers in **N70**. It annihilates the vertical
      pressure-gradient term (4.23e-4 → 3.3e-11 m/s²) but cuts the acoustic
      ringing only 27 %, and by t = 600 s the balanced and unbalanced columns
      have the same max|dw/dt| to 2 %. The mode is permanent and equation-level.
    - **(5) Checkpointed ≤6 h windows: not built, and the stated reason was
      wrong.** The pinned `maxiters` never triggered in any run: gap (s) is not
      what stops a 24 h span. The dynamics half alone is stable and expensive —
      18.0 s wall for a 60 s span, 164.9 s for 600 s, ≈0.27 s of wall per
      simulated second, so 24 h of the *dynamics half alone* is ≈6.6 h — with
      max|w| going 7.41e-4 → 9.60e-4 → 9.35e-4 m/s over 0/60/600 s and
      max|dw/dt| FALLING 4.23e-4 → 3.13e-4. Window length is therefore a cost
      decision, not a solver-failure one.

    **The unmeasured number that decides 3.1's feasibility.** One single-step
    evaluation of the nine-component physics column costs 86.35 s end to end
    against 0.15 s for the five-component dynamics column (`./esm validate` of
    the physics column is 3.18 s, so ≈83 s is interpreter build + solve +
    assertions). **One-time interpreter build has not been separated from
    per-RHS-evaluation cost**, and until it is, nothing about a 24 h full-column
    integration is decidable: if per-evaluation cost is near 86 s the run is out
    of reach by orders of magnitude and the acceptance test needs a different
    execution path, not smaller windows. This is the next measurement to make.

    **Acceptance criteria not yet measured** — no 24 h run exists. For scale,
    the reference signal from `scm_ref_traj`: θ changes by 2.811 K (Linf) over
    24 h, so the 0.5 K criterion is 18 % of the signal; q_v's relative L2 change
    over the column is 0.377, so 5 % is 13 % of it; the lowest level goes
    u 3.00 → 4.23 and v −9.00 → −4.29 m/s. At step 1410 the physics D(θm) is
    3.80e-4 K/s against the dynamics' 1.02e-5, i.e. **the dynamics half alone
    reproduces ~2 % of the θm tendency**; D(u) runs the other way, dynamics
    4.25e-4 against PBL 2.76e-4.

    **Both of these were BUILT on 2026-09-20; the physics half is now mountable.**

    - **20 state coupling-target parameters (the estimate said ~19; it missed
      Dudhia's `T_in`) and five mount-edge libraries.** YSU `theta/u/v/qv/qc/qi_in`
      (`scmstate_couple_ysu_inputs.esm`); Dudhia `T/qv/qc/qr/qi/qs/qg_in`
      (`..._sw_inputs.esm`); RRTM column `qv/qc/qr/qi/qs_in` + `cldfra_in`
      (`..._rrtmcol_inputs.esm`); slab `T_s_in` (`..._slab_inputs.esm`); and
      CalCldfra's own six (`..._cldfra_inputs.esm`). Slab's soil geometry `ze`
      is deliberately NOT a target — the layer depths come from SOILPARM and are
      time-invariant, like the dynamics half's base state.

      **Two components had to stop owning their state for these edges to mean
      anything.** YSU carried `ic(x) = input_x` with `D(x) = dxdt` for all six of
      theta, u, v, qv, qc, qi, and Dudhia the same for `T`. Lowering `input_x` to
      a coupling target on those would have bound only the INITIAL CONDITION: the
      component would then integrate its own copy under the PBL (or radiative)
      tendency alone, while the document integrated the shared state under the
      sum, and the two would silently diverge over a 24 h run. Both are tendency
      diagnostics — YSU's outputs are dthdt/dudt/dvdt/dqvdt/dqcdt/dqidt and
      Dudhia's is dTdt — so the seven `ic`/`D` pairs became plain algebraic reads
      `x = input_x`, which is the instantaneous-derivative form CLAUDE.md
      prescribes and the form every other component here already had (slab and
      the RRTM stages were already algebraic). Both files now have no structural
      equation at all, which post-EarthSciAST #412 (`static_evaluation_assertions`)
      is a supported document kind in Rust, Julia and Python. Every assertion in
      both is at t = 0 on an observed, so the counts are unchanged: YSU 157/157,
      Dudhia 222/222.

    - **`cal_cldfra` exists**, as
      `components/atmospheric_radiation/cloud_fraction/{cal_cldfra,cal_cldfra_parameters}.esm`,
      **63/63**. It is WRF's `cal_cldfra1` (icloud = 1): the condensate threshold,
      the ice-weighted blend of the Murray over-water and over-ice saturation
      fits, and the Randall (1994) fractional fit with its −6.9 exponential floor
      and 0.01 snap, pointwise in `lev`. It observes `branch`, the Fortran's own
      `cldfra1_flag`, so a test pins which arm fired and not only the number.

      **The reference had to be a new real64 kernel replay, and that is the
      interesting part.** The in-model real32 CLDFRA is NOT an adequate reference
      here: the Randall branch divides the condensate by
      (RHGRID·QVS_WEIGHT − QV)^GAMMA, a cancellation of two near-equal numbers,
      and near saturation that costs WRF's single precision up to 5.1e-5
      absolute — 1.2e-4 relative, twelve times looser than the rel 1e-5 contract
      for real32 references. Measured: re-evaluating this component's own
      expressions in binary32 reproduces the in-model CLDFRA to 6.0e-8, so the
      gap is WRF's arithmetic and not the transcription. `kernels/cldfra_driver`
      (+ `cldfra_state_stub.F90` and a Makefile rule) therefore compiles
      cal_cldfra1 EXTRACTED VERBATIM from `module_radiation_driver.F` at build
      time, the same awk-extraction trick the solar kernels use, so the reference
      cannot drift from the model. Residual against it is binary64 roundoff
      (≤ 4.3e-16), so tolerances are abs 1e-12 on the columns and rel 1e-9 on the
      spot checks.

      Six regimes, from two reference runs, chosen to cover branches the SCM
      alone cannot: em_scm_xy calls 1 (fully clear — and not trivially so, RHUM
      reaches 27 in the stratosphere, so the condensate guard is the only thing
      holding those five layers at 0), 120 and 360 (pure-ice high cloud); and the
      em_quarter_ss supercell columns 20/40 at calls 61 (pure LIQUID, the only
      case that exercises the over-water arm alone), 268 (the full ice-weight
      sweep 0 → 1, and the only case that reaches the −6.9 floor, 8 layers, and
      the 0.01 snap, 10 layers) and 396 (deep ice anvil, 30 fractional layers).
      Every scm_cloud dump in the repo is pure ice, which is why the supercell
      case is here. Mutation-checked: replacing the ice fit with the water fit
      fails 44/63, dropping the phase blend 44/63, dropping the −6.9 floor 4, the
      0.01 snap 2, the condensate guard 1 — the last of those is why the clear
      column is a test at all, since the first five regimes leave that guard
      unexercised.

      The dumps are built by `tools/prep_cldfra_dumps.py` (dump extraction only,
      no equations) into `data/eqwefic/dumps/cldfra/`. A Fortran defect found
      while transcribing is **FORTRAN_BUGS.md B13**: cal_cldfra1's `ELSE` arm
      writes a CLDFRA that the block below unconditionally overwrites, and in that
      arm `QCLD` and `weight` are never assigned, so the code reads an
      uninitialised local on the first iteration and the previous cell's value
      afterwards. Harmless for every microphysics option WRF ships through this
      path; not transcribed.

    **The full state document was BUILT on 2026-09-20: `couplings/eqweather_scm.esm`,
    65/65.** Sixteen components co-mounted over ten ODE states (u, v, w, φ′, θm,
    q_v, q_c, q_i, T_s, T_sk) through 217 coupling edges, with every tendency a
    sum: D(u)/D(v) = Coriolis + curvature + sponge + RUBLTEN/RVBLTEN, D(θm) =
    sponge + WRFMoistThetaTendency's conversion of the radiative and PBL
    tendencies, D(q_v)/D(q_c)/D(q_i) = YSU's, D(T_s)/D(T_sk) = the slab's.
    Nothing about the column is read from a dump the state could supply; only
    q_r/q_s/q_g stay frozen, because no microphysics is mounted, and they are
    read by reference from the existing hand-authored libraries rather than
    transcribed. One more component had to give up private state first:
    `rrtm_lw_heating` carried `ic(T) = input_T` with `D(T,t) = dTdt` for a `T`
    that nothing in it ever reads, i.e. 59 spurious ODE states in every document
    mounting it; it is now `T = input_T` with a `T_in` target fed from `Geo.T`
    (56/56 unchanged, and all thirteen documents that mount it still green).

    **What the test showed.** The fourteen surface-layer residuals come out at
    exactly `couplings/geometry_surface_pbl_column.esm`'s N50 numbers to two
    significant figures (sfc_hfx 7.3e-5, sfc_mol 7.1e-5, … sfc_ust 2.1e-6
    relative) — that agreement is the result, because it says co-mounting the two
    halves introduced nothing of its own; the gap is the 0.30 Pa
    driver-versus-kernel offset in the lowest layer, already isolated. Everything
    else meets the two halves' own tolerances unchanged. The three SUMS land at
    Linf 6.2e-5 (du_dt), 3.0e-5 (dv_dt) and 1.1e-4 K/s (dtheta_m_dt), each inside
    the tolerance its own PBL part already carries, so no term is dropped or
    doubled. `cldfra` matches WRF's CLDFRA exactly at abs 0 — but step 1410 is a
    clear column (trace cirrus only, all below the 0.01 snap), so that assertion
    pins the wiring and the zero, not the fractional fit; cal_cldfra's own
    six-regime set is where the fit is exercised.

    **(5) REVISITED, and the answer is not windowing. MEASURED: one right-hand-side
    evaluation of the full document costs 1 h 47 m** (106 m 56 s wall, 99 % of one
    core, 481 MB) on an EMPTY time span — esm-spec §6.6.2's instantaneous-derivative
    shape, so nothing is integrated and no Jacobian is formed. Against 86 s for the
    nine-component physics column. It is not the component count and not the
    assertion count: a copy of the document carrying a SINGLE assertion costs the
    same order. It is that no other document in `couplings/` drives the FULL RRTM
    longwave chain (Col → Setcoef → GasOpt → Rtrn) from LIVE geometry —
    `radiation_column.esm` runs the whole chain in seconds from dumped inputs, and
    `geometry_radiation_column.esm` drives Col and Heat from live geometry in
    seconds but stops before SETCOEF. Here the geometry expressions stop being
    constant leaves and propagate through 103 layers × 140 g-points. BDF over a 1 s
    span, tried first, had not finished a single step in 14 minutes, which is the
    ~300-state Jacobian on top of that. **So the unmeasured number that decided
    3.1's feasibility is now measured, and it was an interpreter defect, not the
    model: EarthSciAST #438, fixed by PR #439 (see (aa) below). The driver built
    an implicit solver even for a run that cannot advance, and materialising its
    finite-difference Jacobian cost 2·n_states + 1 full right-hand-side
    evaluations to return the initial state unchanged. With the fix the same
    65/65 takes 2 m 31 s. THE FEASIBILITY QUESTION IS THEREFORE REOPENED, NOT
    CLOSED: one RHS evaluation is ~2.5 s, not ~107 min, so a 24 h run is worth
    re-costing from scratch once #439 merges — the earlier conclusion that it was
    four orders of magnitude out of reach was measuring the defect, not the
    model.**
    The practical consequence for this repo is recorded in CLAUDE.md: `./esm test .`
    is now ~2 hours, and that one file should be excluded from a fast whole-repo
    check and run on its own.

    **WSM6 MOUNTED 2026-09-20, and the document now has a CLOUDY test.** Eight more
    mounts (Sat, Sed, Melt, Warm, Cold, IceDep, Rescale, SatAdj) take it to 24
    components, 339 coupling edges and 656 ODE states; q_r, q_s and q_g became
    states, and with all six moisture species prognostic `q_t` stopped being a
    frozen input and is now computed in the document and fed to the geometry and
    the vertical-momentum buoyancy term (bit-identical to the old frozen profile,
    which was exactly q_v + q_c + q_i). The step-1410 test is **65/65 unchanged**.
    The new step-60 / call-120 test is the cloudy one: cirrus ice to 4.64e-5 kg/kg,
    cal_cldfra = 1.0 in seven layers, OLR 129.7 W/m² against 256.0 in the clear
    column — so cal_cldfra's Randall fit, the RRTM cloud optics and WSM6's ice
    branch are all exercised for the first time in a coupled document. **69/69,
    but read the tolerances before believing it**: the three microphysical
    tendencies DISAGREE with WRF and their bounds record the disagreement rather
    than certify agreement (mp_dqv_dt 1.1e-8 against a reference maximum of
    3.3e-8, mp_dqi_dt 2.8e-8 against 2.2e-8, mp_dtheta_dt 4.3e-5 K/s against
    1.2e-4). The cause is measured, not assumed: three diagnostics split D(q_i)|mp
    into its stages at exactly 0 / 2.15e-8 / 5.30e-8 kg/kg/s, so **sedimentation
    alone is 2.4× WRF's whole ice tendency**. That is the PLM-vs-donor-cell gap
    `sedimentation.esm` already documents about itself, appearing in a coupled
    document for the first time because step 1410 has no falling ice. Closing it
    needs a PLM-equivalent flux rule in EarthSciDiscretizations, or a dt→0
    extrapolated reference of the kind the component's own sedimentation tests
    use. It is the next thing to fix on this document. The vapour and heating
    residuals are a different, expected effect: WRF stages fallout → melting →
    rates → adjustment sequentially, each on the state the last one left, while a
    continuous formulation evaluates every stage at one instantaneous state and
    sums — O(dtcld) by construction, and not removable without abandoning the ODE
    form.
    Cost after WSM6: **7 m 58 s for both tests on the PR #439 binary**, 3 m 47 s
    for the clear-sky test alone (up from 2 m 31 s). On `main` it got worse, not
    better: 479 → 656 states means #438 now charges 1313 full RHS evaluations per
    run instead of 959, and a run was killed unfinished at 15 minutes.

    **The PBL reference was WRONG IN KIND, and fixing it is what 2026-09-21 did
    (TENDENCY_ERROR_BUDGET.md).** Every coupled document in `couplings/` that
    mounts YSU was asserting its PBL tendencies against WRF's driver-level
    `rthblten`/`rublten`/`rvblten`/`rqvblten`/`rqcblten`/`rqiblten`. Those are
    the IMPLICIT tridiagonal step's `(c_new - c_old)/dt` at dt = 60 s — an
    increment, not a derivative — while the documents evaluate an instantaneous
    right-hand side. The gap that produced (2–31 % of the column maximum) was
    being carried in the tolerances, which is why six of them sat at `abs 1e-3`
    to `1e-7` and why 15 of `eqweather_scm.esm`'s 134 assertions had a tolerance
    larger than the whole signal. Eight documents are now referenced to the
    **dt → 0 Richardson limit** `(8T(dt/4) − 6T(dt/2) + T(dt))/3` of real64
    `kernels/ysu_driver` replays of the SAME dumped column at dt = 0.04/0.02/0.01 s
    (`data/eqwefic/dumps/ysu/scmref_<call>_dt*.json`, transcribed by
    `tools/reref_pbl_dt0.py`) — the reference `ysu.esm`'s own tests already used.
    **A trap found on the way:** `driver_120_dt*.json`, which an earlier note
    called "the replays that already exist" for the coupled step-60 test, are
    replays of a DIFFERENT run (`ysu_120.flat`, z0 = 0.15 m); the reference run
    `data/eqwefic/scm_ref` has z0 = 0.05 m and its column is `scm_ref_120.flat`.
    Using the wrong one would have moved every PBL reference by ~25 %.
    The result is the measurement that matters: with the dumped geometry
    (`physics_column.esm`) the residual is **1e-9 K/s and m/s²**, about 1e-5 of
    the column maximum, so the ESM's YSU tendency IS the kernel's dt → 0 limit;
    with a computed geometry the residual rises to ~1e-6, which is the same N50
    lowest-layer pressure offset the surface layer shows. Tolerances everywhere
    are now 4× the measured residual.

    **Two more regimes on `couplings/eqweather_scm.esm` (2026-09-21).** The two
    existing tests were both sunlit and both convective, so they could not
    separate a state-dependent tendency error from a constant one and never
    entered YSU's stable branch or the longwave-only radiation path. Added:
    **step 900** (physics call 1800, 03:33 local, cos zenith −0.552, Ri_b +0.037,
    surface-layer regime 1, h = 262 m over 6 layers, K_h max 2.20 m²/s against
    176.8 at step 1410) and **step 300** (call 600, 17:33 local, just after
    sunset, cos zenith −0.088, h = 556 m over 10 layers, and the most exercised
    cloud fraction of the four — `cal_cldfra` returns 0.0272/0.0373/0.0670/0.1356
    and 1.0, so the Randall fit rather than only its snap is evaluated inside a
    coupled document). `gsw` and `dthdt_sw` are asserted at **abs 0** in both, so
    the shortwave path is pinned OFF rather than small. The step-300 dumps had
    to be merged first (`esm_dump_dyn_suite_300.json`, `esm_dump_scm_suite_300.json`
    via `tools/merge_dumps.py`; the recipe was validated by rebuilding the
    step-1410 merges byte for byte). The numbers are transcribed by
    `tools/add_scm_regime.py`, whose expression table re-derives BOTH pre-existing
    tests as a check.

    **The WSM6 operator-split difference is now MEASURED, not assumed.** See
    TENDENCY_ERROR_BUDGET.md for the numbers. Summary: the whole-scheme increment
    has no dt → 0 limit, and the term responsible is **`pigen`**, ice nucleation
    written as `max(0, (roqi0/den − q_i)/dtcld)` — a projection onto the diagnosed
    equilibrium ice content in one sub-step, whose increment is dtcld-independent
    (7.78e-8 kg/kg at dt = 0.01, 0.02 and 0.04 s); the saturation adjustment is
    identically zero in this column. At the model's own dtcld = 60 s `pigen` is
    only 3.7 % of the total, so the ESM (which prescribes the same dtcld) IS
    comparable, and the remaining difference is the sequencing. Measured by
    replaying the kernel at dtcld = 60 s from a synthetic post-sedimentation state
    (`data/eqwefic/dumps/wsm6/syn_postsed_120.flat`): `pidep` moves
    3.29972e-8 → 4.45250e-8, so `R(x1) − R(x0) = 1.145e-8` kg/kg/s — **exactly the
    document's `mp_dqv_dt` residual, to four significant figures**, and 91 % of its
    `mp_dtheta_dt` residual once carried through `xl/cpm/pi`. `mp_dqi_dt` is the
    exception: its residual is dominated by the fallout stage, where the document's
    donor-cell flux divergence gives 5.304e-8 against WRF's own PLM increment
    3.637e-8 (1.46×, and 1.52× against the dtcld → 0 fallout 3.4935e-8). The
    earlier "sedimentation alone is 2.4× WRF's whole ice tendency" was comparing
    one process against the sum of all of them; the like-for-like ratio is 1.46.

    **The optional schemes are confirmed off in the reference run**
    (`data/eqwefic/scm_ref/namelist.{input,output}`): `cu_physics = 0`,
    `SHCU_PHYSICS = 0`, `GWD_OPT = 0`, `W_DAMPING = 0`, `ICLOUD = 1`,
    `DAMP_OPT = 2`, `SCM_FORCE = 0`, `DIFF_OPT = KM_OPT = 2`. So the absence of
    cumulus, shallow cumulus, gravity-wave drag and w-damping from the mounted set
    is correct rather than an omission.

    Easy win still not taken: the same θm wiring for
    `couplings/scm_physics_column_night.esm` (conversion dumps exist at steps
    1, 60, 300, 900, 1410, 2160, 3000).

    **New EarthSciAST gaps found, repros under `data/eqwefic/phase3_probes/`:**
    (v) array `initial_conditions` **do** work in the Rust CLI for a shaped
    unknown (row-major nested array, esm-spec §11.4 run-time overrides) — the
    earlier claim that the schema restricts them to scalars is **stale**, and
    the state document's three tests rely on it (`repro_wrt_default_shaped_ok.esm`);
    (w) **FIXED upstream, EarthSciAST #406 → #412 (merged 2026-09-19).** The
    inline runner integrated a document that has nothing to integrate instead of
    evaluating it once the way `simulate` does, so an **algebraic-only** model
    failed with `Exceeded maximum number of nonlinear solver failures` at t = 0
    even for `y ~ a·t` (`repro_algebraic_only_no_integration.esm`). The fix adds
    the conformance category `static_evaluation_assertions` and makes Julia
    evaluate when `isempty(prob.u0)`; a §6.6 assertion is now evaluated at its
    own `time`, not at t = 0. `solar_diurnal.esm` still carries a TOA-insolation
    integral as its ODE state — that is now a modelling choice, not a
    workaround, and could be dropped;
    (x) **FIXED upstream, EarthSciAST #407 → #411 (merged 2026-09-19).** `D`
    with `wrt` omitted (spec §4.2: absent means `t`) worked for a **shaped**
    state but failed for a **scalar** one with `State variable 'X' has no
    D(X, t) = ... equation in flat.equations`
    (`repro_wrt_default_scalar.esm` / `_fixed.esm`). The root cause was the
    §4.2 default being applied at one `D` consumer and not the rest; it now
    lives in `op_registry` (`STRUCTURAL_DERIVATIVE_WRT`/`derivative_wrt`) and is
    applied by every consumer in Rust, Julia, Python and Go. Writing `"wrt": "t"`
    explicitly is still the clearer habit but is no longer required;
    (y) **FIXED upstream, EarthSciAST #408 → #410 (merged 2026-09-19).** A bare
    subsystem mount name (`wrf.g`) did not resolve inside an assertion
    `reference` expression, only the model-qualified form did. The fix rebinds
    the assertion parameter scope over flattened name → owner-relative remainder
    → globally unambiguous dotted suffixes. **The deviation is retired:**
    `couplings/eqweather_scm_dynamics.esm` is back to the bare `wrf.g` and
    passes 15/15, so CLAUDE.md's bare-mount-name rule now holds with no
    exception anywhere in the repo;
    (z) an `Assertion` admits no `description` field (schema
    `additionalProperties: false`), so per-assertion provenance has to go in the
    test description.
    (aa) **FILED 2026-09-20 as EarthSciAST #438, fixed by PR #439** (branch
    `perf/issue-438-nonadvancing-run`, 5 files, +395/-37), the cost defect behind
    `couplings/eqweather_scm.esm`'s 1 h 51 m. **The title I filed it under states
    the wrong mechanism and the issue comment corrects it.** What I measured —
    a static document flat in N and an ODE one quadratic — was real, but the
    reason is not the `faq` recurrence. The Rust driver built a diffsol problem
    and an implicit solver for EVERY run, including one that cannot advance: an
    empty `time_span`, or a non-empty one whose whole output grid sits at `t0`.
    Constructing an implicit solver materialises the Jacobian, and that crate's
    is a matrix-free finite-difference JVP, so diffsol calls it once per state
    column and each call evaluates the whole right-hand side twice —
    **2·n_states + 1 full evaluations to return the untouched initial state.**
    In the reproducer NL is BOTH the state count and the sweep depth, which is
    what made it look quadratic in the sweep.

    **The control that settles it** (`data/eqwefic/phase3_probes/`, old binary):
    a 103-deep self-referential sweep with NO state costs 0.03 s; the SAME sweep
    with ONE scalar state that it DEPENDS ON costs 0.09 s; the same sweep over
    103 states costs 1.84 s. A state-dependent recurrence is cheap — it is the
    state COUNT that multiplies. Two further hypotheses were ruled out earlier:
    not the `param_to_var` edge or a const-node-versus-parameter distinction
    (19.27 → 19.85 s, nothing), and not ODE states as such in a cheap document
    (an unused scalar `D(s) = 0` leaves the static radiation column at 19.06 s).
    Corpus-scale before the fix: one input of `couplings/radiation_column.esm`
    changed from an observed to a state, same numbers and same 30/30, is
    **19.27 s → 647.86 s**. **After the fix: `eqweather_scm.esm` is 2 m 31 s for
    the same 65/65, verified here with the PR's binary** (110 m 36 s on `main`).
    The fix also covers a non-empty span whose output grid never leaves `t0`,
    which is this repo's standard `{start: 0, end: 1}` algebraic-test shape, so
    **the whole suite should be re-timed when #439 merges.** Also recorded on
    #438, not fixed: the **Julia** tree-walk runner refuses the §4.3.1.1 causal
    self-reference outright (`E_TREEWALK_UNSUPPORTED_RECURRENCE`), so the
    discrete-recurrence shape the spec sanctions is evaluable in Rust and Python
    only — an honest refusal rather than a wrong answer, and worth its own issue.

    **Stale caveat dropped:** `wrf_solar_geometry.esm` mounts fine under a key
    that is not its own model name — the file carries no model-qualified
    constant reference, so the note above requiring `WRFSolarGeometry` as the
    mount key no longer applies.

    **Trig now honours a unit's SCALE (EarthSciAST #409 → #413, merged
    2026-09-19) — and it caught a real mis-declaration here.** Before #413 the
    transcendentals tested only the DIMENSION of their argument, so
    `sin(90 deg)` evaluated `sin(90 rad)`. #413 makes them apply the declared
    scale. Rebuilding the CLI turned 17 assertions red across
    `wrf_solar_geometry.esm` (63 assertions), `geometry_radiation_column.esm`,
    `geometry_radiation_column_night.esm` and `solar_diurnal.esm`, all of them
    `declin`, `coszen` or a SW flux downstream of one.

    The regression was **ours, not upstream's**. `lib/wrf_constants.esm`
    declared `degrad` as `units: "1"`, transcribing WRF's DEGRAD as the bare
    number it is in Fortran. That made `sin(sun_longitude·wrf.degrad)` type as
    **deg** while its value was already in radians, so #413 applied pi/180 a
    second time — to both sine arguments of
    `declin = arcsin(sin(obliquity·degrad)·sin(sun_longitude·degrad))`,
    collapsing the declination to ~2e-4 rad everywhere. The other trig calls
    read variables *declared* `rad` (`lat_rad`, `hour_angle`, `orbit_angle`,
    `day_angle`) and were never affected: only a `deg`-typed **expression**
    passed directly to a trig function double-converts.

    Fixed by declaring `degrad` with its honest units, `rad/deg` (one line, no
    equation changes): the stored value stays pi/180, `deg · rad/deg` types as
    `rad` at scale 1, and Eq 2 went from unchecked to `consistent [angle]`.
    63/63 green. The same declaration also makes Eqs 3, 11 and 12
    (`orbit_angle`, `hour_angle`, `lat_rad`, each declared `rad` but previously
    assigned a `deg`-typed right-hand side) type-correct, closing a latent
    inconsistency that had been invisible only because scale was unenforced.

    **Rule for the rest of the project:** a conversion factor gets the units of
    the conversion (`rad/deg`), never `"1"`, even when the Fortran it is
    transcribed from carries no units. Any `deg`-valued quantity handed to a
    trig function must either be declared `deg` (and left unscaled) or be
    converted through a variable declared `rad`. This will matter again in the
    fire components, where CFBM carries slope and wind directions in degrees.

    **The §6.6.5 `reference` scope disagrees between `validate` and `test`
    (filed 2026-09-19, repros under `data/eqwefic/phase3_probes/`).** Chasing
    the two loose ends left after #410/#412 turned up one defect with two faces,
    both in the STRUCTURAL VALIDATOR rather than the runner:
    - **EarthSciAST #421** — the validator's assertion-`reference` scope is too
      NARROW. #410 gave the runner three spellings for a mounted parameter
      (flattened, owner-relative, and any globally unambiguous dotted suffix);
      the validator only learned the first two. `./esm validate` reports
      `Variable 'g' referenced in equation is not declared` on a document that
      `./esm test` passes 1/1 — the same binding contradicting itself
      (`repro_reference_bare_tail.esm`). Rust, Julia, Go and TypeScript all
      reject; Python accepts, but indiscriminately — it accepts an AMBIGUOUS
      bare tail too (`repro_reference_bare_tail_ambiguous.esm`), where the Rust
      runner correctly errors, which is the validator-side face of the already
      open **#419**. Commented there rather than filing a third issue.
    - **EarthSciAST #422** — the same scope is too PERMISSIVE about `t`. #412
      made all three executing bindings refuse a reference mentioning the
      independent variable (`REFERENCE_MENTIONS_TIME`, commit `556675570`,
      because the build-time evaluator's time slot holds 0.0 and such a
      reference silently answered with the start-of-span value). No validator
      knows: all five pass `repro_reference_mentions_time.esm`, which no runner
      will run. It is a purely syntactic property — does `t` occur free — so it
      is decidable in the structural pass, and for Go and TypeScript, which
      never execute a test, the validator is the only place it can ever be
      caught.

    **REGRESSION on EarthSciAST `main`: commit `556675570` (PR #412), filed as
    EarthSciAST #432 (2026-09-19).** Rebuilding the CLI to pick up #410-#413
    broke three documents that were green on the parent commit. Bisected in a
    clean worktree (`/scratch.local/ctessum/bisect-eqwefic`), one commit, three
    symptoms that appear and disappear together:

    | document | `d7ad059d3` (last good) | `556675570` -> current `main` |
    |---|---|---|
    | `wildland_fire/.../fire_wind_wrffire.esm` | 21 / 0 | **0 / 21** |
    | `atmospheric_dynamics/wrf_arw/perturbation_coriolis.esm` | 213 / 0 | **198 / 15** |
    | `atmospheric_dynamics/sfclayrev/` | 610 / 0, **22 s** | **no result at 900 s** |

    - `kstar`, a `faq` with `reduce: "min"` over a 44-element index set, answers
      **`inf`** -- the identity of a `min` over ZERO iterations -- and `uf`/`vf`
      go `NaN` downstream of it. The same document's other reductions over the
      same index set are correct on the same run.
    - `perturbation_coriolis` tendencies that are identically zero come back at
      1e-8, and the non-zero ones agree only to ~6e-5 against a 1e-5 tolerance:
      the numbers of a slightly different state, not of a different formula.
    - `sfclayrev` goes from 22 s to not finishing in 40x that. **This is what
      ate the 6 h 17 m whole-repo run** -- `/proc` showed 99 % CPU on one thread
      with 54 read syscalls and 353 kB read in the whole six hours, i.e. spinning
      inside one document rather than working through files. No stack sample was
      possible: `kernel.yama.ptrace_scope=3` on the cluster blocks both `gdb -p`
      and `eu-stack` (the same setting the fire sbatch template works around for
      MPI shared memory).

    The mechanism is that commit's "build once more with the build pipeline on
    and read `observed_field`" retry. All three documents are algebraic and
    state-free, exactly the class the retry targets. The `inf` is the decisive
    one: a reduction returning its identity element was never evaluated, so the
    retry produces a SILENT WRONG NUMBER -- the failure mode #406 was filed to
    remove.

    **Not our corpus.** Checked before blaming upstream: the pre-`faq`-rename
    `aggregate` spelling of `fire_wind_wrffire.esm`, recovered with
    `git show 52b9321^:`, fails identically on the new binary. The `faq`
    migration is exonerated.

    **RESOLVED 2026-09-20 by EarthSciAST #434 and #433, both merged; the pin on
    `./esm` is lifted and the CLI is back on `main` (`3cf6f6f92`).** The fix
    separated two defects that the single `inf` had been hiding:

    - **The deeper one PREDATES #412 and was never about the range.**
      `prepare::eval_observed` boxed a scalar observed as a rank-1, ONE-CELL
      array, but `lookup_variable` answers `Value::Scalar` only for a rank-0
      entry -- so a one-cell vector came back as a FIELD. A `faq` body reading
      such a name therefore evaluated to an array, `reduce_contraction` collapsed
      each term through `as_scalar()` (`None` above rank 0) and produced `NaN`
      per iteration, and then the reducer decided what you saw: `+` accumulated
      `NaN`, while `min` seeded at `+inf` and kept it, because IEEE-754
      `f64::min(inf, NaN) == inf` DROPS the NaN operand. An unreduced identity is
      indistinguishable from a reduction over zero iterations -- which is why
      `kstar` looked like "44 elements, answered like zero". In
      `fire_wind_wrffire.esm` the body reads the scalar observed `kdmax`: the
      READ was empty, not the range. **Our reading in #432 was wrong on this
      point**, and the correction matters because it means `esm simulate` had
      been silently answering `NaN`/`inf` for this whole class on every binary,
      including the `d7ad059d3` we pinned to -- verified directly:

          simulate @ d7ad059d3   above[4] = NaN   kstar = inf   n_above = NaN
          simulate @ 3cf6f6f92   above[4] = 1     kstar = 3     n_above = 2

      Filed and fixed separately as EarthSciAST **#431 / #433** ("a build-time
      field's rank is its declaration's"), with `tests/valid/faq/` and
      `tests/scalar_operand_in_faq.rs` pinning it.
    - **The #412 regression proper** was `with_build_pipeline_if_needed` running
      eagerly inside `build_for_test`, ahead of any solve, so the runner answered
      from `observed_field` and never called `solve` -- replacing a correct
      answer rather than supplying a missing one, and paying a second whole
      document build per `BuildKey` (`sfclayrev`'s 22 s -> no result at 900 s).
      #434 makes it a LAST RESORT, reached only after `solve` refuses a problem
      with `has_nothing_to_integrate`, and a failed retry now leaves `solve`'s
      own diagnostic standing.

    Verified on `3cf6f6f92`: `fire_wind_wrffire.esm` 21/0, `wrf_arw` 213/0,
    `sfclayrev` 610/0 back at 22 s.

    **Measured per-family cost** (pinned binary, for future reference): couplings
    674 assertions / 565 s; rrtm_lw 1451 / 120 s; wildland_fire 725 / 363 s;
    gaschem 825 / 2 s; wsm6 790 / 3 s; sfclayrev 610 / 22 s; dudhia_sw 222 / 5 s;
    wrf_arw 213 / 1 s; ysu 157 / 3 s; land_surface 84 / <1 s; wrf_solar 63 / <1 s;
    lib 12 / 3 s. The whole-repo run is ~20 min, not hours; a run that exceeds
    that is a symptom, not a big suite.

    **This retires the item carried since the #408 work** ("Go/TypeScript reject
    a bare-tail name that the three executing bindings accept"). The measured
    shape is different: it is not Go/TS versus the executing bindings, it is
    every validator versus every runner, and Rust disagrees with itself. The
    second carried item ("Julia's analytic `reference` path still evaluates at
    `t = 0`") is **already fixed** by #412 — by refusing `t` outright rather
    than by re-timing the evaluation — so nothing was filed for it.

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
3.3 **EqAtmFire.esm**: EqWeather + the fire_behavior (CFBM) components via the
    `wildlandfire.esm`-style coupling, compared with a **real-data WRF 4.8 + CFBM
    simulation of the Last Chance Fire (eastern Colorado, 2012-06-25)** with the
    reference physics suite on (decided 2026-09-15; replaces `test/em_fire/hill`).

    **Why not an idealized case.** CFBM refuses any non-Lambert domain
    (`dyn_em/start_em.F:2261`) and locates each fire cell's atmospheric column by
    inverse-projecting its lat/lon (`phys/fire_behavior/io/wrf_mod.F90`), while
    the ideal `em_fire` init writes Cartesian metres into those arrays
    (`set_ideal_coord`, "fake coordinates, in m"), so no shipped ideal case runs
    CFBM. The alternatives considered were WRF-SFIRE `em_fire/hill` (runs, but
    SFIRE on a WRF 4.4 base, contradicting decision 5), WRF 4.8's in-tree
    WRF-Fire `em_fire` (not CFBM), patching the ideal init to Lambert (cheap and
    CFBM-true), and the existing offline replay (one-way only). Every idealized
    fire case also switches the whole physics suite OFF and runs a 3-D TKE
    closure, so it would test none of the six schemes; the real-data case is the
    only option that exercises them.

    **Why this fire.** It is CFBM's own test7 case (Lambert 100 m domain,
    Anderson-13 fuels, ignition line at 39.685 N -103.585 W, 18:00 UTC), and the
    level-set algorithm transcribed here was evaluated on it (Munoz-Esparza et al.
    2018, JAMES), so the existing 486 fire assertions and the offline replay are
    already on this domain. Flat grass plains on a late-June midday exercise
    sfclayrev, YSU, slab, Dudhia and RRTM without steep-terrain grey-zone issues.
    About 45,000 acres (~182 km2) burned by the next morning. Runner-up: the
    Marshall Fire (2021-12-30; HRRR forcing, but a winter downslope windstorm in
    Front Range terrain).

    **Data path (verified reachable 2026-09-15).** HRRR's AWS archive starts
    2014-07-30 and GFS's 2021-02, so forcing is ERA5 (0.25 deg hourly, NCAR AWS
    bucket, netCDF needing conversion) or NARR (32 km 3-hourly, NCEI). WPS 4.3 in
    `data/eqwefic/apptainer/wps_wrf_latest.sif` carries `Vtable.NARR` and a
    `GEOGRID.TBL.FIRE` reading NFUEL_CAT from LANDFIRE; geog static data and
    LANDFIRE are downloadable. test7's 1.9 km domain is far too small for the
    burn scar, so only its projection, location and ignition are reused; its
    `wrf.nc` came from WRF 4.3.3 with SFIRE, MYNN and RRTMG (N52) and is not a
    reference. Setup, builds and runs live under `data/eqwefic/lastchance/`.

    **The run is one continuous job, not a restart chain (decided 2026-09-16).**
    The chain was built (19 x 45 simulated minutes on `secondary`, whose wall
    limit is 4 h) and held pending a burned-area check. That check passed --- a
    5-minute run from the 16 Z `real.exe` output (job 10593079) put 7,011 fire
    subgrid cells at non-zero `FIRE_AREA`, drew `FUEL_FRAC` down to 1.9e-16,
    raised `GRNHFX` to 620 kW/m2 over 478 atmospheric cells and drove `RQVFRTEN`
    in 10,415 of them, so fire, fire-to-atmosphere coupling and CFBM's own
    `fire_output_*.nc` all work --- but reading the restart path for the chain
    turned up **FORTRAN_BUGS.md B12: CFBM ignores `config_flags%restart`**. The
    registry flags the whole fire state for restart and WRF does restore it into
    `grid%*`, but the only copy runs `fire_state -> grid`, and `start_em.F:2290`
    re-initializes `fire_state` unconditionally, so every segment boundary would
    have reset the fire to unburned and re-fired the ignition (which is anchored
    to the *segment* start). The chain would have produced 19 disjoint fires.

    So the case runs as a single 14 h job on the `ctessum` partition (3 d wall
    limit): `runs/templates/run_wrf_exe_ctessum.sbatch` with
    `namelists/namelist.input.cont16` (`run_hours = 14`, no restarts), submitted
    2026-09-16 as job 10593449 in `data/eqwefic/lastchance/runs/full`. Cost from
    the 10593079 timings: the nest integrates at 164 s wall per simulated minute
    (d01's recursive "Timing for main" covers all five domains), so 840 simulated
    minutes is **~38.5 h wall on 40 cores, ~1,600 core-hours**, plus ~32 GB of
    `wrfout` and `fire_output` files. Because CFBM writes `fire_output_*.nc`
    every 600 s, a job that dies early still leaves every completed hour usable,
    so the no-restart exposure degrades gracefully. Fixing B12 properly (export
    `lfn_hist`, `fmc_g` and the `fire_*_old` fields, then seed `fire_state` from
    the restored arrays when `config_flags%restart`) is worth doing in the
    `ctessum-claude/WRF` fork, but is not on the critical path for this run.

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

## 6.9 Pairwise coupling libraries (2026-09-17)

Strategy (e) of `data/eqwefic/design/assembly_composition.md` is implemented.
`couplings/lib/` holds six libraries — `sfclayrev_ysu` (8 edges), `sfclayrev_slab`
(2), `slab_ysu` (2), `dudhia_slab` (1), `rrtm_lw_slab` (1) and `rrtm_lw_chain`
(35, four roles) — and the 15 assemblies import them instead of repeating the
edges: **393 inline edges replaced by 25 imports, a net -2913 lines**. The
duplication was larger than the design doc estimated (the RRTM longwave chain's
35 edges were copy-pasted into nine assemblies). `INVENTORY.md` carries the
catalogue. Still inline: the 118 column-state edges, a second family
(`column_state_<scheme>`) that pins a naming convention for the state carrier
rather than wiring two schemes.

Two verification notes worth keeping. First, neither `esm diff` nor
`esm coupling-analysis` can show that an import and the edges it expands to are
equivalent — both compare the DOCUMENT, so an import always reads as different.
The usable proof is a negative control: delete one edge from inside a library and
confirm the importing assembly's Fortran-dump tests fail. Second, the refactor
and the CLI fix below landed together, so the refactor was re-measured with the
UNCHANGED binary (631/0/0, identical to baseline) before the new binary was
installed, keeping the two claims separable.

**EarthSciAST PR #401** came out of this. A coupling library's `coupling[].from`
/`.to` prefixes name `coupling_roles`, but four of the five bindings resolved them
against a symbol table of models — which a library has none of — so every
well-formed library was rejected, including EarthSciModels' own
`fastjx_superfast.esm`, `fastjx_geoschem.esm` and `wildlandfire_behavior.esm`.
That is why EarthSciModels CI sweeps `components lib registered_functions` and
never `couplings/`. Fixed in Rust, Python, Julia and TypeScript with a role-based
check (not a skip, so role typos are still caught, and unlike the system path it
also checks `to`). **Go was not affected** — `isLibraryDocument` already
short-circuits — so the Go patch was reverted rather than left as dead code.

## 6.10 The `faq` node tag (2026-09-17)

The unified query node was renamed `aggregate` → `faq` (Functional Aggregate
Query) at esm 1.1.0; `aggregate` is a deprecated alias, normalized on load with
a warning, and **removed at 2.0.0**
(`../EarthSciAST/docs/content/rfcs/faq-node-rename.md`). This repo migrated ahead
of the removal: **193 files, 2771 occurrences**, plus **183 version bumps**
(`"esm": "1.0.0"` → `"1.1.0"`), because a document that spells `faq` while
declaring a version below 1.1.0 is a hard structural error
(`faq_version_too_old`), the gate reading the document as authored. Zero
transitive bumps were needed: no document lowers to `faq` only through an
imported template. `arrayop` — removed outright at 0.8.0, and rejected by name
rather than falling through to the open operator tier — was already absent.
Assertion counts did not move.

`"op"` is the only key in this repo that ever took `"aggregate"` as a value, so
the rename was textually unambiguous. Seven prose uses of the word survive in
`description` fields ("ice aggregates to snow", "written as an aggregate minimum
over the level index", …); several now describe a node the file spells `faq`,
which is stale wording to fix when those files are next edited for other reasons
— rewriting them is authoring, not a wire-tag rename.

**The deprecation warnings do not go to zero from here, and the reason matters.**
267 remain (70 components+lib, 197 couplings) with **zero originating in
EqWeFiC**. They come from `{ref}`-loaded siblings — `../EarthSciDiscretizations`
carries **2022 files / 5753 occurrences** still on the alias, `../EarthSciModels`
2 files / 3 — and the alias check fires at each binding's single wire boundary,
which every document a load touches passes through. So ESD, not this repo, is the
real exposure at esm 2.0.0, and clearing it needs a migration PR there.

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
