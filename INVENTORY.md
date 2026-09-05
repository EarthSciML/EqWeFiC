# Component inventory

Status values: `not started` · `instrumented` · `stub` (tests written, on this repo) · `physics` (tests pass locally) · `PR` · `merged`.
`Existing` names an EarthSciModels file that already implements (part of) the component.

## EqWeather — reference suite (WRF v4.8.0)

| Component | Sub-components (planned) | Fortran source | Existing in EarthSciModels | Status |
|---|---|---|---|---|
| `lib/wrf_constants.esm` | — | `share/module_model_constants.F` | — | physics (tests pass, 2026-09-04) |
| `lib/wrf_thermo.esm` | T↔θ, θv, θli, moist density (qsat/Exner still to add) | `module_model_constants.F`, per-scheme | `atmospheric_dynamics/sp_ch1/*` | physics (templates used by YSU, 2026-09-05) |
| Column grid + rules | `column_nonuniform_1d` grid, face-flux divergence, interface-K diffusion with prescribed fluxes, 5 `integral` forms (downward sedimentation flux deferred) | EarthSciDiscretizations PR #34 | `grids/column_nonuniform_1d/` | physics PR (julia/python/rust gate green, 2026-09-05) |
| Surface layer (sfclayrev) | φm/φh, bulk Richardson regimes, z0 over water, u*, exchange coefficients, 2 m/10 m diagnostics | `phys/physics_mmm/sf_sfclayrev.F90` | `local_scale/surface_layer_profile.esm` | not started |
| PBL (YSU) | PBL height (bulk Ri), K-profile, countergradient term, entrainment flux, top-down mixing, diffusion PDE | `phys/physics_mmm/bl_ysu.F90` | `holtslag_boville/*` (related) | physics (Spike A, 2026-09-05): `components/atmospheric_dynamics/ysu/ysu.esm`, 117/117 `./esm test` and Python binding; EarthSciModels gate blocked by its shaped-observed sampler; cloud-top (`ysu_topdown_pblmix`) branch not transcribed |
| Land surface (slab) | 5-layer soil heat diffusion, surface energy balance | `phys/module_sf_slab.F` | — | not started |
| Land surface (Noah) | soil heat/moisture, canopy resistance, Penman, snow | `phys/module_sf_noahlsm.F` | `urban_canopy/hydro/*` (related) | not started |
| Microphysics (WSM6) | saturation adjustment; warm rain (praut, pracw, prevp); ice processes; melting/freezing; sedimentation | `phys/physics_mmm/mp_wsm6.F90` | — | not started |
| SW radiation (Dudhia) | clear-sky absorption/scattering, cloud albedo/absorption, column integrals | `phys/module_ra_sw.F` | — | not started |
| LW radiation (RRTM) | per-band optical depth, Planck, up/down sweeps, heating rate | `phys/module_ra_rrtm.F` | — | not started |
| Cumulus (New Tiedtke) | trigger, updraft/downdraft, closure | `phys/physics_mmm/cu_ntiedtke.F90` | — | not started (optional) |
| Gravity-wave drag (GWDO) | — | `phys/physics_mmm/bl_gwdo.F90` | — | not started (optional) |
| Subassembly: surface+land+PBL column | — | SCM dumps | — | not started |
| Subassembly: radiation column | — | SCM dumps | — | not started |
| Subassembly: full SCM physics suite | — | SCM dumps | — | not started |
| Deferred: RRTMG LW/SW, Noah-MP, MYNN, Thompson, KF | | | | — |

## EqAtmChem — WRF-Chem `chem_opt=1` then `300`

| Component | Fortran source | Existing | Status |
|---|---|---|---|
| Anthropogenic emissions | `chem/emissions_driver.F` | `earthsci_data/nei2016_monthly.esm` (pattern) | not started |
| Biogenic emissions | `chem/module_bioemi_*.F` | — | not started |
| Dry deposition (Wesely) | `chem/dry_dep_driver.F`, `module_dep_simple.F` | `atmospheric_deposition/wesley_dry_gas.esm` | not started |
| Photolysis | `chem/module_phot_fastj.F` / TUV | `gaschem/fastjx/*` | not started |
| RADM2 gas mechanism | `chem/KPP/mechanisms/radm2` | — | not started |
| Wet scavenging | `chem/module_wetscav_driver.F` | `atmospheric_deposition/wet_deposition.esm` | not started |
| GOCART aerosols | `chem/module_gocart_*.F` | — | not started |

## EqAtmFire — NCAR fire_behavior (Community Fire Behavior Model, WRF 4.8 `phys/fire_behavior`)

| Component | Fortran source | Existing | Status |
|---|---|---|---|
| Fuel categories table (Anderson 13) | `physics/fuel_anderson_mod.F90` | `wildland_fire/fuel_model_lookup.esm` | not started |
| Fuel moisture model | `physics/fmc_wrffire_mod.F90` | `wildland_fire/nfdrs/*` (related) | not started |
| Rate of spread (Rothermel + limits) | `physics/ros_wrffire_mod.F90` | `wildland_fire/rothermel/*` | not started |
| Level-set front propagation | `physics/level_set_mod.F90` | `wildland_fire/level_set_fire_spread.esm` | not started |
| Fire → atmosphere heat/moisture flux | `physics/fire_model_mod.F90`, `fire_driver_mod.F90` | `wildland_fire/level_set/fire_heat_flux.esm` | not started |
| Atmosphere → fire wind | `physics/fire_driver_mod.F90` | `wildland_fire/midflame_wind.esm` | not started |
| Fire–atmosphere coupling file | `driver/fire_behavior.F90`, `tests/` | `couplings/wildlandfire_behavior.esm`, `../wildlandfire.esm` | not started |

## Dynamics (stage 3)

| Component | Content | Fortran source | Where | Status |
|---|---|---|---|---|
| `atmospheric_dynamics/wrf_arw/*` | flux-form Euler equations, hybrid sigma-pressure coordinate, base-state split, diagnostic pressure | `dyn_em/module_big_step_utilities_em.F`, `module_small_step_em.F`, Skamarock et al. 2021 Ch. 2 | EarthSciModels | not started |
| `grids/wrf_arw_c/*` | C-grid mass-coordinate grid, 5th/3rd-order upwind advection with PD limiter, divergence damping, metric terms | `dyn_em/module_advect_em.F`, Ch. 3 | EarthSciDiscretizations | not started |
