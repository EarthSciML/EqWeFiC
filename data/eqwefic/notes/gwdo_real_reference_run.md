# GWDO real-terrain reference run (stage 1), 2026-09-22

This run supplies stage-1 reference data for **GWDO** (`gwd_opt = 1`, `phys/physics_mmm/bl_gwdo.F90`). The idealized SCM column has no sub-grid orography, so GWDO does nothing there. The same run also dumps **New Tiedtke** (`cu_physics = 16`) and **Noah** (`sf_surface_physics = 2`) over real terrain. Those dumps complement the idealized-SCM reference (`chem170_scm/`), which another run produces.

## Where everything is

| what | where |
|---|---|
| WRF source | fork `ctessum-claude/WRF`, branch `earthsciml-instrumented-gwdo` (off `earthsciml-instrumented` @ 3eccae3); worktree `data/eqwefic/wrf-gwdo` |
| build | `printf "34\n1\n" \| ./configure` (GNU dmpar), `./compile -j 8 em_real` inside `apptainer/wrf-build.sif`, 18 min |
| run directory | `data/eqwefic/gwdo_real/run` (namelist, wrfinput/wrfbdy, 19 hourly `wrfout_d01_*`) |
| launcher | `data/eqwefic/gwdo_real/run_wrf_gwdo.sbatch` (Slurm job 10689122, `secondary` partition, 20 MPI ranks) |
| dumps | `data/eqwefic/gwdo_real/dumps/esm_dump_{gwdo,ntiedtke,noah}_s<step>_j<j>_i<i>.json` (490 files, 356 MB) |
| kernel replay | `data/eqwefic/gwdo_real/replay_all.py` → `replay_summary.txt`; drivers `wrf-gwdo/kernels/build/gwdo_driver` (real64) and `build32/gwdo_driver` (real32) |

(`data/eqwefic/...` above means the disk-backed data tree `/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/`, not this repo.)

## Domain and configuration

- **Domain.** d01 of the Last Chance Gulch WPS setup: Lambert, 8.1 km, 150 × 180 cells, 51 levels (`auto_levels_opt = 2`, `p_top` 50 hPa). It covers Colorado and the Front Range, 34–46 N and 96.5–110.6 W. Terrain reaches 3711 m and the sub-grid standard deviation `VAR` reaches 606 m.
  - The plan asked for 12–30 km. This run uses 8.1 km because the ERA5 intermediate files behind the existing `met_em` files have been deleted. A coarser grid would need a fresh ERA5 download and a new `ungrib`/`metgrid` pass.
  - GWDO is fully active at 8.1 km (see below). `dxmeter` and `dx_factor` enter the scheme explicitly, so a coarser run adds regimes but no new code paths.
- **Period.** 2012-06-25 12Z to 06-26 06Z (18 h), `dt` = 45 s, so there are 1440 steps.
- **Forcing.** ERA5 through the existing `met_em.d01.*` files, with hourly lateral boundaries.
- **Physics.**

  | option | value |
  |---|---|
  | `mp_physics` (WSM6) | 6 |
  | `ra_lw_physics` (RRTM) | 1 |
  | `ra_sw_physics` (Dudhia) | 1 |
  | `radt` | 10 |
  | `sf_sfclay_physics` (sfclayrev) | 1 |
  | `sf_surface_physics` (Noah) | 2 |
  | `num_soil_layers` | 4 |
  | `bl_pbl_physics` (YSU) | 1 |
  | `cu_physics` (New Tiedtke) | 16 |
  | `cudt`, `bldt` | 0 |
  | `gwd_opt` | 1 |
  | `icloud` | 1 |

  - `gwd_opt` lives in `&dynamics`, not `&physics` (Registry `namelist,dynamics`). The first submission in `&physics` failed with a namelist read error.
  - Dynamics are copied from the Last Chance d01 run.
- **Chemistry: off.** Adding `chem_opt = 170` here was optional. It would need a chem-enabled build, CBM-Z/MOSAIC initial and boundary conditions, and anthropogenic and biogenic emissions, none of which exist for this domain. It would also have put the meteorological reference behind a much longer chain. The idealized-SCM run covers `chem_opt = 170`.

## GWDO needs ELVMAX, and the existing WPS data did not have it

WRF 4.8's GWDO (the Hong et al. 2025 revision) reads `ELVMAX`, the maximum sub-grid orographic height (`omax`). It sets `ldrag = omax <= hmt_min (50 m)`, meaning no drag (`bl_gwdo.F90:416`).
- WPS 4.3 (the version in `wps_wrf_latest.sif`) writes no `ELVMAX`, so `real.exe` leaves it at 0. **GWDO is then silently off at every point with `gwd_opt = 1`.** Nothing in `module_check_a_mundo.F` or `real.exe` warns about this. It is recorded as FORTRAN_BUGS B17.
- Fix used here: WPS's `util/compute_gwdo.py` (WPS master @ 5feccec, copied into `gwdo_real/`) recomputed `CON`, `VAR`, `OA1-4`, `OL1-4` and `ELVMAX` from `topo_gmted2010_30s` on a copy of `geo_em.d01.nc`. `ELVMAX` spans 231–4227 m.
  - `SourceData.py` imports `cartopy` without using it, so the script ran with a stub module on `PYTHONPATH`.
- `metgrid` copies static fields verbatim, and the ERA5 intermediate files were gone. So the new fields were copied into copies of the 19 `met_em.d01.*` files (`gwdo_real/patch_met_em.py`). That produces what a `metgrid` rerun would have written.
- The recomputed fields replace the legacy `orogwd_10m` values. `VAR` differs by up to 394 m, because `compute_gwdo.py` computes the statistics over a 2·dx box from 30″ terrain.

## Instrumentation (branch `earthsciml-instrumented-gwdo`, commit 790bb4f)

`module_esm_dump` gains positional selection for 3-D MPI runs:
- `esm_dump_want_at(scheme, itimestep, j)` selects by `ESM_DUMP`, `ESM_DUMP_STEPS` and `ESM_DUMP_J`.
- `esm_dump_want_i(i)` selects by `ESM_DUMP_I_LIST`.
- `esm_dump_open_at` names each file by step, global row and `its` (or point `i`), so ranks do not collide.

The existing per-scheme call counter is left alone.

| scheme | where | contents |
|---|---|---|
| `gwdo` | `phys/module_bl_gwdo.F`, per tile row | Every `bl_gwdo_run` input, including the orography statistics and `omax`, plus every output. `rublten`/`rvblten` are INOUT: the YSU tendency goes in, and YSU + GWDO comes out. The GWDO part is `rublten − rublten_in`, which equals `dtaux3d` up to the real32 rounding of the sum. Assert on `dtaux3d`/`dtauy3d`. |
| `ntiedtke` | `phys/module_cu_ntiedtke.F`, per tile row | Driver-level inputs, the vertically flipped kernel-boundary arrays from `cu_ntiedtke_pre_run` (`k_*`), the kernel outputs (`k_*_out`, `k_zprecc`), and the `post_run` tendencies `r{th,qv,qc,qi,u,v}cuten`. |
| `noah` | `phys/module_sf_noahdrv.F` (`lsm`), per point | Every `SFLX` input and inout state before the call, and every output after it. |

**Selection used.**
- Steps 2, 400, 560, 720, 880, 1040 and 1440. These are 12:01Z, 17Z, 19Z, 21Z, 23Z, 01Z and 06Z; local time is UTC − 6.
- Rows `j` = 61, 81, 91, 101 and 141, which cross the Sangre de Cristo, Sawatch, Front Range and Wyoming ranges.
- Noah points `i` = 26, 39, 51, 75, 110 and 140, spanning the mountains and the plains.

With the 4 × 5 decomposition this gives:
- 140 GWDO row files and 140 New Tiedtke row files, each 750 columns per step;
- 210 Noah point files.

## Results

**Kernel replay (GWDO).** All 140 dumped rows were replayed through a standalone `bl_gwdo_run` driver.
- The real32 build (`PREC=`, WRF's own flags) reproduces WRF **bit for bit in every output on 140/140 rows**, so the dump captures every input.
- The real64 build (`kind_phys = real64`) agrees with the real32 in-model values to a median of 1.4e-5 and a worst case of 1.5e-3 of each field's row maximum. The worst case is `dtaux3d` at s2_j81_i76.
  - The mechanism was isolated afterwards and is **not** a threshold flip: see "Where the real64 replay differs from WRF" below.
  - Tests should use the real64 replay and state that.

**GWDO is active, strongly diurnal:**

| step (UTC) | max \|dtau\| m s⁻² (dumped rows) | max surface stress Pa | columns with stress (of 750) |
|---|---|---|---|
| 2 (12:01Z) | 8.7e-3 | 1.77 | 446 |
| 400 (17Z) | 0 | 0 | 0 |
| 560 (19Z) | 0 | 0 | 0 |
| 720 (21Z) | 0 | 0 | 0 |
| 880 (23Z) | 2.4e-2 | 2.06 | 8 |
| 1040 (01Z) | 2.0e-2 | 1.35 | 383 |
| 1440 (06Z) | 2.8e-2 | 3.22 | 636 |

The domain-wide `wrfout` fields agree with this:
- At night the stress is non-zero in 18–24 k of 27 k columns, with surface stress up to 7.8 Pa and |DTAUX3D| up to 4.8e-2 m s⁻².
- In the afternoon only 31–53 columns have any stress.

The scheme switches drag off where `bnv2 < 0` or `velco < 0` below the reference level, so the daytime mixed layer turns it off. That is the scheme's design, not a fault. The daytime dumps therefore exercise the all-off branch, and the night and transition dumps exercise the wave, blocking and limiter branches. The wind-reversal limiter holds in every dumped column: |dtau|·dt / |wind| ≤ 0.97.

**New Tiedtke is active.**
- Convective rain falls in 4–23 of 750 dumped columns per step from 12Z to 01Z, and in none at 06Z.
- The maximum `pratec` is 9.3e-4 mm s⁻¹.
- |`rthcuten`| reaches 0.165 K s⁻¹ and |`rucuten`| reaches 5.8e-2 m s⁻².
- Domain `RAINC` reaches 42 mm by 02Z.
- At step 2 the scheme runs with `qvften = thften = 0` (`itimestep == 1` branch). Only step 2 and later are dumped, so that branch is not exercised.

**Noah** is dumped at 30 points per step, across the diurnal cycle:

| quantity | range, W m⁻² |
|---|---|
| sensible heat (`sheat`) | −84 to 442 |
| latent heat (`eta`) | −3.5 to 491 |
| ground flux (`ssoil`) | −175 to 81 |

## Not done here

- No `.esm` authoring, by design.
- No real64 driver for New Tiedtke or Noah yet. The dumps carry everything a driver needs; Noah additionally needs the `VEGPARM`/`SOILPARM`/`GENPARM` tables through `SOIL_VEG_GEN_PARM`.
- No coarser (12–30 km) run, for lack of ERA5 intermediate files.

## Where the real64 replay differs from WRF, and what tolerance that supports (2026-09-22)

The first pass reported that the real64 GWDO replay differs from the in-model real32 values by a median 1.4e-5 and a worst 1.5e-3 of each row maximum, and guessed "single-precision threshold tests". That guess was tested and is **wrong in its mechanism**: no threshold flips. The difference is single-precision round-off amplified by two near-cancellations, one of which sits exactly on the critical Richardson number.

**Method.** `kernels/gwdo_diag_patch.py` generates a copy of `bl_gwdo.F90` that also exposes each column's branch record (`d_ldrag`, `d_kref`, `d_komax`, `d_kbomax`, `d_kblk`, `d_icrilv`, the per-level saturation flag `d_sat`, the limiter `d_dtfac`, and `d_taup`, `d_velco`, `d_usqj`, `d_fr`, `d_taub`, `d_ulow`, `d_bnv2low`). The physics is untouched: the real32 diagnostic build is bit-identical to WRF on the worst row. `kernels/gwdo_diag_driver.F90` runs it, and `gwdo_real/diag_compare.py` replays all 140 rows at both precisions.

**Result: no branch flips at all.** Over 140 rows (52 500 columns of branch record and 262 500 cells of limiter state):

| branch decision | real32 vs real64 disagreements |
|---|---|
| `ldrag` (the ELVMAX/50 m, `var`, critical-level and stability gates) | 0 |
| `kref`, `komax`, `kbomax` (reference and blocking levels) | 0 |
| `kblk` (flow-blocking layer) | 0 |
| `icrilv` (critical level above the reference level) | 0 |
| saturation hypothesis vs no wave breaking, per level | 0 |
| limiter active vs inactive, per cell | 0 of 262 500 |

So the mechanism is **(b), accumulation**, not (a) a threshold flip; and it is not (c) a transcription difference, since the real32 driver reproduces WRF bit for bit on all 140 rows.

**The amplifier.** Traced on the worst cell (row `s2_j81_i76`, column 31, level 11, where the saturation branch is taken in both precisions):

| step | real32 | real64 | relative difference |
|---|---|---|---|
| `usqj` (the Richardson number) | 0.25660196 | 0.25664472 | 1.7e-4 |
| `2*sqrt(temc) - temc`, `temc = 2 + 1/sqrt(Ri)` | 1.2927e-2 | 1.3009e-2 | 6.3e-3 |
| `taup(k+1)` (goes as that factor squared) | 7.3528e-3 | 7.4454e-3 | 1.2e-2 |
| `taup(k+1) - taup(k)` (the tendency is this difference) | −4.6500e-3 | −4.5575e-3 | 2.0e-2 |
| `dtaux3d` at the cell | −2.3632e-5 | −2.3161e-5 | 2.0e-2 (1.5e-3 of the row maximum) |

Two cancellations, in series:
1. `2*sqrt(temc) - temc` is **zero at `temc = 4`, that is at `Ri = 0.25`, which is exactly `ric`**, the value the branch test compares against. Columns that take the saturation branch marginally therefore have the smallest factor and the largest amplification. Here it turns 1.7e-4 into 6.3e-3, a factor of 37. Recorded as FORTRAN_BUGS N85.
2. The tendency is `g*(taup(k+1) - taup(k))/del`, a difference of stresses that are within a factor of 2.6 of each other at this cell.

The row-level residual tracks the first mechanism: taking each row's smallest `|2*sqrt(temc) - temc|` over the levels where the saturation test was evaluated,

| rows | count | residual/amplitude, `dtaux3d` |
|---|---|---|
| margin < 0.05 (`Ri` within ~2 % of 0.25) | 27 | median 9.2e-5, max 1.5e-3 |
| margin ≥ 0.05 | 19 | median 1.1e-5, max 8.7e-4 |

The rank correlation between the row residual and the margin is 0.58, so the marginal-`Ri` conditioning explains much but not all of it; the rest is the second cancellation, which is present in every column.

**Tolerance that the dumps support.** Two cases, and they differ by six orders of magnitude, so a test must say which reference it uses.

- **Against the real64 replay** (recommended, per the CLAUDE.md rule for double-precision kernel references): the reference is exact to the `.esm` formulation, so `rel: 1e-9`. The GWDO scheme itself carries no real32 module constants of the kind noted for YSU and WSM6 in N2.
- **Against the in-model real32 values**, measured over the 57 rows with GWDO signal, as residual / row maximum:

  | field | median | p90 | max | suggested bound |
  |---|---|---|---|---|
  | `dtaux3d` | 5.7e-5 | 4.2e-4 | 1.5e-3 | 2e-3 of the row maximum |
  | `dtauy3d` | 6.2e-5 | 6.3e-4 | 1.1e-3 | 2e-3 |
  | `dusfcg` | 5.2e-6 | 5.3e-5 | 1.2e-4 | 3e-4 |
  | `dvsfcg` | 4.6e-6 | 5.5e-5 | 1.8e-4 | 3e-4 |

  Per-cell relative tolerances are not supportable at all: at a cell where the stress divergence nearly cancels, the relative difference reaches 2e-2. Normalise by the row (or column) maximum.

**Rows that cannot be pinned tightly** are those with a near-critical Richardson number in the saturating levels. The eight worst, by `dtaux3d` residual / row amplitude: `s2_j81_i76` (1.5e-3, margin 1.6e-3), `s2_j61_i76` (1.3e-3, 5.3e-3), `s2_j61_i114` (1.1e-3, 3.7e-2), `s2_j141_i1` (8.7e-4), `s2_j61_i1` (6.9e-4), `s1440_j61_i1` (4.3e-4), `s2_j141_i39` (4.1e-4), `s1440_j91_i39` (3.2e-4). Rows with a comfortable margin, such as `s2_j91_i39` and `s2_j101_i39`, sit at 3.5e-7 to 1.2e-6 and can carry a much tighter bound if a test wants one.

## New Tiedtke: the chem170 driver runs on these rows too (2026-09-22)

No new driver was written. `ntiedtke_driver.F90` from fork branch `earthsciml-instrumented-chem170` (@ ac6d4cc) builds against this worktree unchanged; `gwdo_real/nt_adapt.py` renames this run's kernel-boundary keys (`k_pu` → `pu` and so on, plus the three sizes `im`, `kx`, `kx1`) into the flat input it expects. All 140 real-terrain rows replay (`gwdo_real/nt_compare.py`).

- **Trigger decisions are precision-independent here too:** `ktype`, `kcbot` and `kctop` agree between the real32 and real64 replays in every column of every row, 0 disagreements.
- **The real32 replay is bit-identical to WRF on 57 of 140 rows**, and on the other 83 it differs by at most 1.2e-7 of the row amplitude, i.e. float32 round-off. `cu_ntiedtke_post_run` forms the tendencies as `(xf - x) * rdelt` with `rdelt = 1/delt` precomputed, while the driver divides by `delt`; that is a one-ulp difference and it explains the size, though it was not isolated further.
- **Precision sensitivity is about 100x GWDO's.** Over the 24 rows with convective precipitation, real64 versus in-model real32, residual / row maximum:

  | field | median | p90 | max |
  |---|---|---|---|
  | `rthcuten` | 5.5e-4 | 4.3e-3 | 4.1e-2 |
  | `rqvcuten` | 4.8e-4 | 4.4e-3 | 1.5e-2 |
  | `rqccuten` | 1.8e-4 | 3.4e-3 | 1.9e-2 |
  | `rqicuten` | 4.1e-4 | 1.7e-3 | 1.9e-2 |
  | `rucuten` | 3.0e-4 | 1.6e-3 | 6.5e-3 |
  | `rvcuten` | 4.3e-4 | 2.7e-3 | 1.1e-2 |

  So a stage-2 New Tiedtke test on these real-terrain rows should reference the real64 replay, not the in-model values; against in-model values 5e-2 of the row maximum is the honest bound. Rows whose amplitude is dominated by a single triggered column are the loose ones.

**Noah was not replayed on these rows.** The chem170 `noah_driver` reads the Noah parameter tables from the dump (`tbl_*`, as `REDPRM` uses them), and this run's `SFLX` dump does not carry them; its inout state also uses `_in` suffixes. Reusing it here needs that branch's table-dump hook added to this branch's wrapper and the run repeated (about 25 min of build plus run), not a naming adapter. Left for whoever authors the Noah component, who has the chem170 SCM dumps already.
