# Bugs and suspected errors found in the Fortran sources

One entry per finding. Include file:line, WRF/WRF-SFIRE commit, what the code does, what it should do, and how it was noticed (which test).

| # | Repo @ commit | File:line | Description | Impact on .esm | Reported upstream |
|---|---|---|---|---|---|
| N1 (note) | NCAR/MMM-physics @ 550b5b4 | `bl_ysu.F90:1` | `#define NEED_B4B_DURING_CCPP_TESTING 1` makes `bl_ysu_run`'s `ttnp` a temperature tendency (θ tendency × `pi2d`) and `dtsfc` use `del/pi2d`; the WRF wrapper compensates. | Standalone consumers and tests must divide by `pi2d` (Spike A reference is `ttnp/pi2d`). | no (documented behaviour) |
| N2 (note) | NCAR/MMM-physics @ 550b5b4 | `bl_ysu.F90` parameters and literals | `kind_phys` parameters (`h1=0.33333335`, `h2=0.6666667`, `afac/bfac=6.8`, `d1..d3`, `tmin`) and in-line literals (15.9, 0.15, 1.746, 1.286, 2.1, 0.16, 1.e-7, −0.18) are default-real, so a real64 build runs with single-precision roundings (`h1 = 0.3333333432674408 ≠ 1/3`); with the real32 module constants this limits real64 agreement to ≈1e-7 relative. | Spike A tolerances 1e-6 rather than 1e-9. | not yet |
