# Bugs and suspected errors found in the Fortran sources

One entry per finding. Include file:line, WRF/WRF-SFIRE commit, what the code does, what it should do, and how it was noticed (which test).

| # | Repo @ commit | File:line | Description | Impact on .esm | Reported upstream |
|---|---|---|---|---|---|
| N1 (note) | NCAR/MMM-physics @ 550b5b4 | `bl_ysu.F90:1` | `#define NEED_B4B_DURING_CCPP_TESTING 1` makes `bl_ysu_run`'s `ttnp` a temperature tendency (θ tendency × `pi2d`) and `dtsfc` use `del/pi2d`; the WRF wrapper compensates. | Standalone consumers and tests must divide by `pi2d` (Spike A reference is `ttnp/pi2d`). | no (documented behaviour) |
| N2 (note) | NCAR/MMM-physics @ 550b5b4 | `bl_ysu.F90` parameters and literals | `kind_phys` parameters (`h1=0.33333335`, `h2=0.6666667`, `afac/bfac=6.8`, `d1..d3`, `tmin`) and in-line literals (15.9, 0.15, 1.746, 1.286, 2.1, 0.16, 1.e-7, −0.18) are default-real, so a real64 build runs with single-precision roundings (`h1 = 0.3333333432674408 ≠ 1/3`); with the real32 module constants this limits real64 agreement to ≈1e-7 relative. | Spike A tolerances 1e-6 rather than 1e-9. | not yet |
| B1 | NCAR/MMM-physics @ 550b5b4 | `sf_sfclayrev.F90` | Not kind-generic: `alog/amax1/amin1` and default-real literals passed to `kind_phys` dummies; fails to compile with `kind_phys` = real64. | Fixed in the ctessum-claude fork (52cfbd7) so a real64 driver can be built. | not yet |
| B2 | wrf-model/WRF v4.8.0 | `phys/Makefile` target `submodules` | Tests `run/NoahmpTable.TBL` relative to `phys/` (should be `../run/`), so the test always fails and every `./compile` runs `git submodule update --init --recursive`, silently resetting instrumented submodule checkouts. | Commit the submodule pointer before compiling. | not yet |
| N3 (note) | NCAR/MMM-physics @ 550b5b4 | `mp_wsm6.F90` | Process rates are dt-dependent by construction: limiters scale with `1/dtcld` (`min(rate, q/dtcld)`, `satdt`), and pigen/pcond/melting/freezing are instantaneous adjustments whose increment is dt-independent. | WSM6 tests compare the dumped `<rate>_raw` rates pointwise against the dumped evaluation state instead of a small-dt replay. | n/a |
