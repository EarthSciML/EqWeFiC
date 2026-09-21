# Photolysis in the WRF-Chem single-column reference run: what runs, and why `gaschem/fastjx/*` does not reproduce it

Written 2026-09-21, before any rewrite, as required by the task brief.
All Fortran line numbers are `../WRF` (WRF 4.8.0) unless stated otherwise.

---

## 1. The headline: the SCM did not run Fast-J, and it certainly did not run Fast-JX

`/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/chem_scm/namelist.input`:

```
 &chem
 chem_opt    = 101,
 phot_opt    = 1,
 photdt      = 1,
 chemdt      = 1,
 gas_drydep_opt = 0,
 aerchem_onoff  = 0,
 ...
```

`Registry/registry.chem:4098-4101` fixes what `phot_opt` means:

| `phot_opt` | package | driver | source |
|---|---|---|---|
| **1** | **`photmad`** | **`madronich1_driver`** | **`chem/module_phot_mad.F`** |
| 2 | `photfastj` | `fastj_driver` | `chem/module_phot_fastj.F` |
| 3 | `ftuv` | `ftuv_driver` | `chem/module_ftuv_driver.F` |
| 4 | `tuv` | `tuv_driver` | `chem/module_phot_tuv.F` |

and `chem/photolysis_driver.F:160-165` dispatches `CASE (PHOTMAD)` to `madronich1_driver`.

**The reference run used the Madronich scheme** (`chem/module_phot_mad.F`, 3673 lines):
S. Madronich's 1987 spectral code as modified by W. R. Stockwell (IFU) and
R. Forkel (1995), with `ps2str`/`sphers` added by S. Massie in Aug 2008 and the
O3 cross sections / O(1D) quantum yields updated to NASA-JPL-2011 by S. McKeen
in Feb 2013.

Two further points worth stating plainly:

* **Fast-JX is not a WRF-Chem option at all.** WRF 4.8 offers Fast-J
  (`phot_opt=2`, Wild/Prather v3-era, `module_phot_fastj.F` + `module_fastj_mie.F`),
  not Fast-JX (the Prather v6/v7 successor). So
  `../EarthSciModels/components/gaschem/fastjx/*` does not correspond to *any*
  WRF-Chem photolysis option, let alone the one that ran.
* Those two `.esm` files are self-described migrations of
  `GasChem.jl/src/interpolations_FastJX.jl` (`FastJX_interpolation_troposphere`)
  out of `polecat/esm-lnl6`. They are a GEOS-Chem-lineage tropospheric lookup,
  not a transcription of any Fortran in this project's reference model.

---

## 2. What RADM2 actually needs

`components/gaschem/radm2/radm2.esm` declares exactly **21** photolysis
parameters on `RADM2RateConstants`, each `units: "1/s"` and documented as
`ph_* (min^-1) / 60`. The WRF side confirms this exactly.

`chem/KPP/mechanisms/radm2/radm2.eqn:2-22` has 21 photolysis reactions, whose
rates are the KPP tokens `j(Pj_*)`. The WRF-KPP coupler
(`chem/KPP/util/wkc/gen_kpp_interf_utils.c:111-123`) generates, into
`chem/module_kpp_radm2_interface.F`:

```fortran
INTEGER, PARAMETER :: njv=52
jv(Pj_o31d) = REAL(ph_o31d(i,k,j)/60., KIND=dp)
...
```

verified against the actually-generated file in the build tree
(`data/eqwefic/wrf-chem/chem/module_kpp_radm2_interface.F:28-49, 160, 245`):

```
Pj_o31d=1  Pj_o33p=2  Pj_no2=3   Pj_no3o2=4   Pj_no3o=5   Pj_hno2=6
Pj_hno3=7  Pj_hno4=8  Pj_h2o2=9  Pj_ch2or=10  Pj_ch2om=11 Pj_ch3cho=12
Pj_ch3coch3=13 Pj_ch3coc2h5=14 Pj_hcocho=15 Pj_ch3cocho=16 Pj_hcochest=17
Pj_ch3o2h=18 Pj_ch3coo2h=19 Pj_ch3ono2=20 Pj_hcochob=21 Pj_macr=22 ...
```

`jv` is 52 long because the coupler emits one entry per `PHOTR*` registry
field for every mechanism. **RADM2 reads only `jv(1:21)`.** `jv(22)` (`ph_macr`)
is computed by Madronich but unused; `jv(23:52)` are never assigned by
`madronich1_driver` at all.

### The identity that makes this testable

`madronich1_driver:1595-1599` multiplies by 60 (s⁻¹ → min⁻¹) and
`:1604-1625` assigns `ph_<x> = phot1(n)`; the KPP interface then divides by 60.
Therefore, for n = 1..22:

> **`jv(n)` (as dumped in `dumps/radm2_wrf/*.json`) == `phot1(n)` in s⁻¹,
> exactly as `photolysis_mad` returns it, bit for bit up to the real32 round trip.**

So the 52x59 `jv` block already sitting in the RADM2 dumps is a ready-made,
independent check on any Madronich reimplementation, at seven model times.

### phot1 index vs. cross-section index — NOT an off-by-one bug

`photolysis_mad:1987-2005` sets `phot1(nr,·)` from `d(nr+1,·)`. The comment
block at `module_phot_mad.F:1501-1524` numbers the **cross-section table**
(`xs`/`xqy` second index, `nr` = 1..22+), whose entry 1 is **O2 absorption**.
`phot1` drops that O2 entry, so `phot1(n) = d(n+1)`:

| `d`/`xs` index | process | `phot1` index | `ph_*` |
|---|---|---|---|
| 1 | O2 absorption (Schumann-Runge) | — | *(not exported)* |
| 2 | O3 -> O(1D) | 1 | `ph_o31d` |
| 3 | O3 -> O(3P) | 2 | `ph_o33p` |
| 4 | NO2 | 3 | `ph_no2` |
| 5 | NO3 -> NO+O2 | 4 | `ph_no3o2` |
| 6 | NO3 -> NO2+O | 5 | `ph_no3o` |
| 7 | HNO2 | 6 | `ph_hno2` |
| 8 | HNO3 | 7 | `ph_hno3` |
| 9 | HNO4 | 8 | `ph_hno4` |
| 10 | H2O2 | 9 | `ph_h2o2` |
| 11 | CH2O -> radical | 10 | `ph_ch2or` |
| 12 | CH2O -> molecular | 11 | `ph_ch2om` |
| 13 | CH3CHO | 12 | `ph_ch3cho` |
| 14 | CH3COCH3 | 13 | `ph_ch3coch3` |
| 15 | CH3COC2H5 | 14 | `ph_ch3coc2h5` |
| 16 | HCOCHO process a | 15 | `ph_hcocho` |
| 17 | CH3COCHO | 16 | `ph_ch3cocho` |
| 18 | HCOCH=CHCHO (*estimate, no reliable measurement*) | 17 | `ph_hcochest` |
| 19 | CH3O2H | 18 | `ph_ch3o2h` |
| 20 | CH3COO2H (*"actually use 0.28*(h2o2 value)"*) | 19 | `ph_ch3coo2h` |
| 21 | CH3ONO2 | 20 | `ph_ch3ono2` |
| 22 | HCOCHO process b | 21 | `ph_hcochob` |
| 23 | MACR | 22 | `ph_macr` (unused by RADM2) |

The comment block and the assignments agree once the `d(nr+1)` shift is applied.
(An earlier pass through this code mistook the shift for a stale comment; it is not.)

---

## 3. What the Madronich scheme actually computes

Per column, once per `photdt` (= 1 min here, i.e. every step):

1. **Solar geometry.** `calc_zenith(xlat, -xlong, julday, gmtp, azimuth, zenith)`
   (`:2653-2759`) — its own ephemeris (geometric mean longitude
   `ml = 279.2801988 + .9856473354*d + 2.267e-13*d*d`, equation of time,
   `pi = 3.1415926535590`). Note WRF passes **`-xlong`**, a deliberate sign flip.
2. **Nighttime cutoff.** `zenith == 90` -> 89.9; `zenith >= 90` -> skip RT entirely,
   all `phot1 = 0`.
3. **Effective cosine.** `zenita = cos(zenith)`, but for `zenith > 75` it is
   replaced by **`1/chap(zenith)`**, the reciprocal Chapman function, a 22-point
   table (`:2138-2165`) linearly interpolated per whole degree.
4. **Column assembly** (`:1526-1563`): `tt` (T), `rhoa` (rho), `o33` = model O3 in ppmv
   **floored at 1e-3**, `aerext`, `qll` = 1e3*(qc+qi)*rho in g/m3 **zeroed below 1e-5**,
   `phizz` = height above ground in km.
5. **Extension above the model top**: `nabv = 10` extra levels from the built-in
   standard profiles `zabv/tabv/o3abv/pabv/caabv`, stretched by
   `znorm = (50 - ztop)/(50 - 20)` and linearly blended onto the model top state.
6. **Ozone column rescaling**: `o3scal(dobsnew=325., ho3=4.5, ...)` (`:2105-2133`)
   computes the column in DU (including an above-top `o3(nn)*1e5*ho3` term) and
   **rescales the whole profile so the total is exactly 325 DU**. Note `dobsi=325.`
   is hard-coded at `:1578` (the earlier `350.` is commented out just above).
7. **Regridding + temperature-dependent optics**: `subgrid` (`:2166-2545`) builds the
   internal RT grid (inserting extra layers inside clouds, so `nlevel` is *not*
   `kte+1`), computes layer-mean columns `vair/vo3/vaer/vcld/vt`, and builds the
   **level-dependent** cross sections `s(lev,kl,nr)` and quantum yields
   `qy(lev,kl,nr)` from `xs`, `xqy`, `txs`, `jpl295`, `jpl218` and the local T.
8. **Schumann-Runge**: `srband` (`:2021-2100`) — Allen & Frederick parameterization
   of the effective O2 cross section, using `sra(11,9)` / `srb(11,5)`.
9. **Spherical slant paths**: `sphers` (`:3440-3595`) -> `dsdh`, `nid`.
10. **Radiative transfer**: `optics` (`:2547-2649`) -> `ps2str` (`:2935-3394`),
    a **delta-Eddington two-stream** solve with a pseudo-spherical direct beam,
    solved by a tridiagonal system (`tridag`, `:3398-3435`), **per wavelength bin**,
    over `kl = kl0..kl1 = 30..130`, i.e. **101 bins** of the 130-bin modified-WMO grid.
    Optical properties:
    - Rayleigh: `arayl(kl) = 3.90e-28/λμ**(3.916 + 0.074λμ + 0.050/λμ)`, ω=1.0, g=0.0
    - Aerosol: `aaer(kl) = 0.379*(340/λ)` (Elterman 1968 vertical shape), ω=0.99, g=0.61, scale height 8.05 km
    - Cloud: `reff = 9.6*LWC**0.333`, `bext = (0.0275 + 1.3/reff)*1000`, ω=1.000, g=0.860
    - O3 and O2 absorption from `ao3`, `ao2`
    - **Wavelength-dependent surface albedo** (Demerjian et al. 1980), a step
      function of `wl(kl)` from 0.05 below 400 nm to 0.15 above 660 nm (`:1730-1749`)
11. **J-value integration** (`:1955-1970`):
    `df = fext(kl)*(fldir + fldn + flup)` — direct + down-diffuse + up-diffuse
    **actinic** flux — then `d(nr,lev) += df * s(lev,kl,nr) * qy(lev,kl,nr)`.
12. **Interpolation back** to the model levels by `trapez` (`:2761-2871`), then `*60`.

---

## 4. Diagnosis: everything wrong with `gaschem/fastjx/*` for this purpose

`fastjx.esm` (9059 lines) and `fastjx_interp_troposphere.esm` (37905 lines).

| # | Defect | Detail |
|---|---|---|
| **D1** | **Wrong scheme entirely** | Fast-JX is not the code that ran, and is not even a WRF-Chem option. The reference is Madronich (`module_phot_mad.F`). Nothing in `fastjx.esm` is a transcription of it. |
| **D2** | **Wrong species set; 11 of RADM2's 21 missing** | `fastjx.esm` produces 13 rates: `j_H2O2, j_H2COa, j_H2COb, j_O3, j_O31D, j_o32OH, j_NO2, j_NO3a, j_NO3b, j_N2O5, j_CH3OOH, j_ActAld, j_PAN`. Plausible correspondences to RADM2: `j_h2o2, j_ch2or, j_ch2om, j_o33p, j_o31d, j_no2, j_no3o2, j_no3o, j_ch3o2h, j_ch3cho` = **10**. **Missing: `j_hno2, j_hno3, j_hno4, j_ch3coch3, j_ch3coc2h5, j_hcocho, j_ch3cocho, j_hcochest, j_ch3coo2h, j_ch3ono2, j_hcochob` (11).** `j_N2O5` and `j_PAN` are produced but RADM2 does not use them. |
| **D3** | **0-D, but the physics is a column** | The `FastJX` model takes scalar `T, P, H2O, cos_sza` parameters. A Madronich j-value at level k depends on the *entire* overlying and underlying column through the two-stream solve: cloud in one layer changes j above (backscatter) and below (attenuation), and the surface albedo feeds `flup` at every level. A per-level 0-D lookup structurally cannot express this. The deliverable needs `lev = 59`. |
| **D4** | **Radiation is tabulated, not solved** | `fastjx_interp_troposphere.esm` gets its actinic flux from 18 registered `flux_interp_i(P, cos_sza)` bilinear lookups over a *fixed standard atmosphere*. There is **no radiative transfer**: no Rayleigh, no O3 absorption feedback, no aerosol, **no cloud**, and **no surface albedo**. The reference solves delta-Eddington per bin with all of these. |
| **D5** | **Wrong spectral resolution and binning** | 18 bands vs. Madronich's 101 active bins (`kl0=30..kl1=130`) of a 130-bin modified-WMO grid with per-bin extraterrestrial flux `fext(130)`. The band-integrated cross sections in the esm are Fast-JX's, not Madronich's `xs(130,27)`/`xqy(130,27)`. |
| **D6** | **Wrong cross sections and quantum yields** | Reference tables: `xs`, `xqy`, `txs` (T-dependence), plus `jpl295`/`jpl218` (NASA-JPL-2011 O3 at 295/218 K) and the FTUV `fo3qy` O(1D) yield (Matsumi et al. 2002). Two idiosyncrasies must be reproduced *exactly* or the rate is simply wrong: `ch3coo2h` is **not an independent spectrum** — WRF uses `0.28 * (h2o2 value)`; and `hcochest` is an **estimate with no measurement behind it**. A principled Fast-JX spectrum for either will not match. |
| **D7** | **No pseudo-spherical geometry / no Chapman function** | `fastjx.esm` uses a plain `cos_sza`. The reference switches to `1/chap(zenith)` beyond 75 deg and runs `sphers` for slant-path increments. At the large zenith angles that dominate the diurnal integral this is a first-order difference. |
| **D8** | **Wrong solar geometry** | `fastjx.esm` inlines the "NOAA Spencer-Fourier" `cos_zenith`. The reference is `calc_zenith` with its own ephemeris, and WRF calls it with **`-xlong`**. These do not agree to the tolerance the tests will demand. |
| **D9** | **No above-model extension, no Dobson rescaling** | The reference adds `nabv=10` levels above the model top from built-in standard profiles and then forces the total O3 column to **325 DU** via `o3scal`. Neither exists in the esm. The 325 DU rescaling alone sets the absolute level of `j_o31d`; the McKeen comment at the top of `module_phot_mad.F` records that the 350 -> 325 change was worth "25% increase in surf. JO1D at 17 deg SZA, ~50% at 80 deg". |
| **D10** | **No coupling to the model's own O3** | The reference reads `chem(:,:,:,p_o3)` (floored at 1e-3 ppmv) for the in-model layers, so photolysis is coupled to the chemistry's own ozone. The esm's flux lookup has ozone baked into the table. |
| **D11** | **No cloud, and the reference's cloud parameterization is specific** | `qll = 1e3*(qc+qi)*rho` g/m3, zeroed below 1e-5; `reff = 9.6*qll**0.333`; `bext = (0.0275 + 1.3/reff)*1000` km-1; ω=1.000, g=0.860; and `subgrid` **inserts extra RT layers inside clouds** so `nlevel` depends on the cloud field. |
| **D12** | **Not tested against anything** | `fastjx.esm` has **zero** inline tests. `fastjx_interp_troposphere.esm` has 3, none traceable to WRF. Per CLAUDE.md the WRF-derived tests are the whole point of stage 1. |
| **D13** | **Machine-generated bulk** | 1.4 MB / 38k lines for `fastjx_interp_troposphere.esm`, 244 variables in one flat model with inlined tables. CLAUDE.md requires compositional authoring with factored templates and equations that read as math. Even if the physics were right, this shape is not acceptable to merge. |

### Verdict

`gaschem/fastjx/*` is not a Madronich implementation with bugs in it; it is a
different scheme, of a different lineage, at a different dimensionality, missing
half the species, with no radiative transfer. **It cannot be repaired into the
reference implementation and should not be the starting skeleton.** The correct
move is a new, column-shaped `gaschem/madronich/` component tree transcribed from
`chem/module_phot_mad.F`, with `gaschem/fastjx/*` left alone (it may still be a
valid Fast-JX for GEOS-Chem-lineage mechanisms; it is simply not this).

---

## 5. What the new component has to look like

A factored tree under `components/gaschem/madronich/`, over a `lev` index set of
59 cells (+1 edge), consumer-supplied geometry, constants from `lib/wrf_constants.esm`:

* `solar_geometry.esm` — `calc_zenith` + `chap`, producing `zenith` and `zenita`.
  (Check `lib/solar.esm` and `components/atmospheric_radiation/wrf_solar` first:
  CLAUDE.md requires reuse where an implementation exists, and PLAN.md 3.1 records
  a solar-geometry component already at 63 assertions. It must be confirmed to be
  *this* ephemeris before reuse; WRF has more than one.)
* `column_setup.esm` — the `nabv=10` extension, `o3scal` 325 DU rescaling, `subgrid`
  regrid and the layer-mean columns.
* `optical_properties.esm` — Rayleigh, aerosol, cloud, O3, O2 per bin; `albedoph`.
* `schumann_runge.esm` — `srband`.
* `two_stream.esm` — `ps2str` delta-Eddington. This is the one place a discrete
  recurrence is defensible under CLAUDE.md (an up/down sweep with a tridiagonal
  solve); the description must justify it, exactly as the RRTM two-stream did.
* `j_values.esm` — the `sum_kl fext * (fldir+fldn+flup) * s * qy` integral, written
  as an integral over wavelength, and `trapez` back onto the model levels.
* `data/` — `wl`, `fext`, `xs`, `xqy`, `txs`, `jpl295`, `jpl218`, `sra`, `srb`,
  and the standard profiles, as `data_sources`/`function_tables`. CLAUDE.md permits
  a script to *translate these tables*; the equations must be authored by hand.

Output: 21 observeds `j_o31d ... j_ch3ono2` in **s^-1**, shaped `[lev]`, which is
precisely the `jv(1:21)` block of `dumps/radm2_wrf/*.json` and precisely the 21
parameters `components/gaschem/radm2/radm2.esm` already declares. That closes the
"RADM2 alone is not a column tendency" gap.

## 6. Ground truth now available

`chem/module_phot_mad.F` in the fork `ctessum-claude/WRF`
(`data/eqwefic/wrf-chem`, branch `earthsciml-instrumented`) is instrumented with
`phys/module_esm_dump.F` under scheme name **`photmad`**, dumping per column
(i = its, j = jts):

* inputs: `julday, ktau, gmt, gmtp, curr_secs, xlat, xlong, zenith, zenita,
  azimuth, zsurf, dobsi, phizz, tt, rhoa, o33, aerext, qll`
* the full static tables: `wl, fext, xs, xqy, txs, jpl295, jpl218, sra, srb,
  zabv, tabv, o3abv, pabv, caabv, zstd, tstd, airstd, o3std, aerstd`
* internals: `nlevel, nlayer, zz, t_col, air_col, o3_col, aer_col, cloud_col,
  z, zmid, vair, vo3, vaer, vcld, vt, cvo2, arayl, aaer, albedoph`, the scalar
  optical constants, and `s_sample/qy_sample/ao3_sample/ao2_sample` at 3 levels
* fluxes: `fldir, fldn, flup, endir, endn, enup` on the RT grid, and `d_rt`
* outputs: `phot1_s` (s^-1, == `jv(1:22)`), `phot1_min` (min^-1, == `ph_*`),
  `uvb_dd1, uvb_du1, uvb_dir1, uvrad`

---

## 7. Corrections and additions from building the chain (2026-09-21, later the same day)

Four things measured while implementing that change or sharpen what is written above.
Where they contradict §4, **these supersede it**.

### 7.1 D6 was right about the consequence, wrong about the mechanism

`ch3coo2h` and `hcochest` are **stored table columns**, not runtime formulas. WRF's
comments (`actually use 0.28*(h2o2 value)`, `estimate, no reliable measurement`)
describe how those columns of `xs` were *built*, decades ago, not anything the code
evaluates. Measured from the dumped `xs` table: the ratio `xs(:,20)/xs(:,10)` over the
60 bins where both are non-zero is **not constant** — it runs from **0.1414 to 0.2874**,
clustering near 0.276. So a correct implementation reproduces them by transcribing the
tables. The substantive point in D6 survives unchanged and is if anything stronger: an
implementation built from the literature — which is what `fastjx.esm` is — cannot match,
because no measured peroxyacetic-acid spectrum is 0.276 times hydrogen peroxide's.

### 7.2 The O2 Schumann-Runge parameterization is DEAD CODE in the shipped configuration

§3 step 8 lists `srband` as part of the chain. It is present, and it never runs:

* `srband`'s wavelength loop opens with `IF (wl(kl)>205.) RETURN` (`:2058`);
* the loop starts at `kl = kl0 = 30`, where `wl(30) = 254.8 nm`;
* so it **returns on the first iteration**, having applied nothing;
* and `xs(kl,1)`, the O2 cross section, is **exactly 0.0 in all 101 active bins**
  (it is non-zero only in bins 1–29, which `kl0 = 30` excludes);
* therefore `d(1,lev)` — the O2 photolysis rate — is identically zero, and
  `phot1(n) = d(n+1)` drops that row anyway.

`sra`, `srb` and the whole `cvo2` column exist to feed a routine that cannot fire.
Reactivating it means setting `kl0 = 1`, which WRF's own comment at `:1722` contemplates
("If photolysis is also desired for levels above 2, kl0 should be set equal to 1 again").
**No Schumann-Runge component was written, because there is nothing for it to contribute.**

### 7.3 A cloudy column changes the RT grid SIZE, which a static index set cannot express

`subgrid` inserts `max(int(cloud·dz·0.02), 1)` extra levels per cloudy model layer, so
`nlevel` is a function of the cloud field: **70, 72, 77 and 78** at the four dumped
daytime calls. An esm `index_sets` entry has a size fixed at load. So the cloudy case is
not a transcription problem but a design one — it needs either a fixed maximum grid with
masked levels, or an adaptive-grid discretization. The components built here are written
for, and tested on, the **clear** column (call 1410, max cloud extinction exactly 0),
where the insertion is a no-op and the regridded column *is* the interface column.

### 7.4 What is built, and what the one open edge is

Five components under `components/gaschem/madronich/`, **378 assertions, 0 fail, 0 err**
(`./esm test components/gaschem` → 1203/1203 including the 825 pre-existing RADM2):

| file | assertions | what it pins |
|---|---|---|
| `solar_geometry.esm` | 52 | `calc_zenith`, `chap`, `zenita` |
| `column_setup.esm` | 148 | `nabv` extension, `o3scal` to 325 DU, cloud extinction |
| `layer_columns.esm` | 37 | `subgrid`'s layer means, `vaer` normalization, `cvo2` |
| `spectra.esm` | 36 | T/p-dependent cross sections and the 11 quantum yields |
| `j_values.esm` | 105 | the wavelength integral → **the 21 RADM2 j-values** |

`j_values.esm` reproduces **all 21 j-values at five levels** spanning the column, against
`phot1_s` — which is bit-for-bit `jv(1:21)` of the RADM2 integrator.

**The one open edge is `optics` → `ps2str`**, the delta-Eddington two-stream solve that
produces `fldir/fldn/flup`. In `j_values.esm` those three are *parameters*, supplied from
the dump. So the chain is verified from the column state to the j-values **except the
radiative transfer itself**. What remains is a tridiagonal solve over `2·nlayer = 138`
unknowns repeated for each of the 101 bins — the one place in this scheme where a
discrete recurrence is defensible under CLAUDE.md, and the one piece not attempted here.
