# CAMELS-CL test catchments

Per-catchment daily forcing extracted from **CAMELS-CL** (Alvarez-Garreton et
al. 2018, HESS 22, 5817–5846; data: https://doi.org/10.1594/PANGAEA.894885,
CC-BY-4.0). Columns: `date`, `p` (precipitation, CR2MET, mm), `pet`
(Hargreaves, mm), `qobs` (observed streamflow, mm; `NA` = gap).
Period 1979-01-01 to 2016-12-31 (13,880 days).

| gauge_id | name | area km² | regime | swe_ratio | interv_degree | qobs coverage |
|---|---|---|---|---|---|---|
| 8123001 | Río Itata en Cholguán | 860 | pluvial | 0.00 | 0.003 | 95% |
| 7330001 | Río Perquilauquén en San Manuel | 502 | pluvial | 0.00 | 0.005 | 94% |
| 4511002 | Río Grande en Las Ramadas | 569 | nival (elev 3098 m) | 0.42 | 0.000 | 92% |
| 4703002 | Río Choapa en Cuncumén | 1132 | nival (elev 3142 m) | 0.52 | 0.025 | 99% |

Selection criteria: near-natural (no big dams, interv_degree < 0.1, glacier
< 2%), long records. The pluvial pair has snow_frac & swe_ratio < 0.05; the
snow-dominated pair (Norte Chico Andes) exercises the HBV snow routine and
carries a `tmean` column (CR2MET catchment mean).

## Split-sample results (DDS, 3000 evaluations, seed 42, warm-up 365 d)

| catchment | model | objective | cal A → val B | cal B → val A |
|---|---|---|---|---|
| 8123001 | GR4J | KGE | 0.833 → 0.794 | 0.847 → 0.761 |
| 8123001 | HBV  | KGE | 0.858 → 0.797 | 0.824 → 0.846 |
| 8123001 | GR4J | NSE | 0.691 → 0.643 | 0.709 → 0.655 |
| 8123001 | HBV  | NSE | 0.749 → 0.727 | 0.742 → 0.733 |
| 7330001 | GR4J | KGE | 0.866 → 0.819 | 0.868 → 0.823 |
| 7330001 | HBV  | KGE | 0.836 → 0.774 | 0.777 → 0.832 |
| 7330001 | GR4J | NSE | 0.740 → 0.741 | 0.744 → 0.735 |
| 7330001 | HBV  | NSE | 0.764 → 0.767 | 0.770 → 0.760 |

HBV (no snow routine; these catchments are pluvial) outperforms GR4J in NSE on
both catchments — most clearly on 8123001 (+0.07–0.08 in validation), where
its two-box response seems to handle the flashy regime better.

## Snow-dominated catchments (NSE, 3000–4000 evaluations)

| catchment | model | cal A → val B | cal B → val A |
|---|---|---|---|
| 4511002 | GR4J (no snow) | 0.041 → −0.343 | 0.158 → −0.013 |
| 4511002 | HBV w/o snow | 0.041 → 0.228 | 0.238 → 0.036 |
| 4511002 | HBV + snow | 0.473 → 0.329 | 0.679 → 0.314 |
| 4703002 | GR4J (no snow) | −0.053 → 0.023 | 0.093 → −0.163 |
| 4703002 | HBV w/o snow | 0.188 → 0.110 | 0.254 → −0.084 |
| 4703002 | HBV + snow | 0.794 → 0.616 | 0.721 → 0.344 |

The snow routine is decisive: without it both models are useless in these
Andean catchments (val NSE ≤ 0.23, often negative); with it HBV reaches
val NSE 0.31–0.62. TT and SFCF calibrate high (3–4 °C, 1.1–2.0) because the
lumped catchment-mean temperature underrepresents the cold high elevations
and CR2MET undercatches high-Andes precipitation. Absolute values are
consistent with published CAMELS-CL benchmarks for the arid Norte Chico.

## Elevation bands (semi-distributed HBV, NSE, 4000 evaluations)

Three equal-area bands per catchment with TCALT = 0.6 °C/100 m and
PCALT = 0.10/100 m (HBV-light defaults), forcing referenced to the catchment
mean elevation.

| catchment | model | cal A → val B | cal B → val A | calibrated TT |
|---|---|---|---|---|
| 4511002 | HBV lumped | 0.473 → 0.329 | 0.679 → 0.314 | 3.0 / 3.9 °C |
| 4511002 | HBV 3 bands | 0.499 → **0.513** | 0.744 → 0.327 | **0.0 / 0.5 °C** |
| 4703002 | HBV lumped | 0.794 → 0.616 | 0.721 → 0.344 | 3.8 / 3.5 °C |
| 4703002 | HBV 3 bands | 0.642 → 0.228 | 0.415 → 0.531 | 2.8 / 3.4 °C |

The headline result is **physical**, not just skill: in 4511002 the bands let
TT calibrate back to ~0 °C (the defensible rain/snow threshold) instead of the
3–4 °C the lumped model needed as a fudge for the cold high terrain — and
validation NSE improves on one fold (0.33 → 0.51). For 4703002 (an enormous
1153–5038 m relief) three equal bands with default lapse rates are too coarse.

### Calibrating the lapse rates (TCALT, PCALT)

Fixing TCALT/PCALT at the HBV-light defaults wastes the bands on the
big-relief catchment. Adding the two lapse rates to the search (14-parameter
DDS) fixes it:

| catchment | configuration | cal A → val B | cal B → val A |
|---|---|---|---|
| 4703002 | lumped HBV | 0.794 → 0.616 | 0.721 → 0.344 |
| 4703002 | 3 bands, fixed lapse | 0.642 → 0.228 | 0.415 → 0.531 |
| 4703002 | 3 bands, **fitted lapse (DDS)** | 0.867 → **0.756** | 0.799 → **0.628** |
| 4703002 | 3 bands, fitted lapse (SCE-UA) | 0.817 → 0.680 | 0.805 → 0.660 |

With fitted lapse rates the semi-distributed model finally beats the lumped
one *robustly* (validation NSE 0.63–0.76 vs 0.34–0.62, and far less spread
between folds). SCE-UA lands slightly lower peaks but more balanced across
folds — the signature of a more global optimum.

### Hypsometric (equal-area) band geometry

Hand-picking band elevations is arbitrary. `--hypsometry "min,median,max"`
builds `n` equal-area bands whose elevations are read off the catchment's
hypsometric curve (here reconstructed from the three CAMELS-CL elevation
quantiles; a curve sampled from a DEM clipped to the catchment plugs into the
same `ElevationBands::equal_area_from_hypsometry`). Five equal-area bands,
fitted lapse rates:

| catchment | geometry | cal A → val B | cal B → val A |
|---|---|---|---|
| 4511002 | 3 bands, hand-picked | 0.543 → 0.322 | 0.785 → 0.393 |
| 4511002 | 5 bands, hypsometric | 0.534 → **0.430** | 0.787 → **0.441** |
| 4703002 | 3 bands, hand-picked | 0.867 → 0.756 | 0.799 → 0.628 |
| 4703002 | 5 bands, hypsometric | 0.891 → 0.738 | 0.800 → 0.648 |

Equal-area hypsometric bands match or beat the hand-tuned geometry (clearly so
on 4511002, +0.05–0.11 validation NSE) while removing the guesswork — the band
elevations are now objective and reproducible from reported attributes. The
reconstruction from three quantiles is coarse; a DEM-derived curve is the
accuracy ceiling and uses the identical core constructor.

### DEM-derived hypsometric curves

`<id>_hypsometry.csv` holds the real hypsometric curve (101 knots) of each
snow catchment, computed from the **Copernicus DEM GLO-30** clipped to the
official CAMELS-CL catchment polygon (`scripts/hypsometry_from_dem.py`). The
clipped DEM statistics reproduce the reported CAMELS-CL attributes almost
exactly (e.g. 4703002: min/mean/max 1158/3143/5054 m vs reported
1153/3142/5038 m), confirming the clip. Feed the curve to the CLI with
`--hypsometry-file`, which builds `n` equal-area bands via the same
`ElevationBands::equal_area_from_hypsometry` core constructor.

| catchment | geometry | cal A → val B | cal B → val A |
|---|---|---|---|
| 4511002 | 5 bands, reconstructed (3 quantiles) | 0.534 → 0.430 | 0.787 → 0.441 |
| 4511002 | 5 bands, **DEM curve** | 0.528 → 0.470 | 0.789 → 0.329 |
| 4703002 | 5 bands, reconstructed (3 quantiles) | 0.891 → 0.738 | 0.800 → 0.648 |
| 4703002 | 5 bands, **DEM curve** | 0.866 → 0.747 | 0.798 → 0.667 |

The real curve is concave (most of the area sits at mid-to-high elevations, so
linear interpolation between three quantiles misplaces the outer bands by
300–400 m). With five bands and the lapse rates calibrated the *skill*
difference is small — the calibrated lapse absorbs some geometry error — but
the DEM geometry is objective and reproducible rather than reconstructed.
Pushing to 10 bands does not help (it raises the between-fold spread:
overfitting on these short records), so five equal-area bands is the
operating point. A DEM is the right ingredient; more bands is not.

Reproduce with:

```sh
rainflow split-sample --forcing data/camels-cl/7330001.csv \
    --objective kge --iterations 3000 --seed 42
```

Note: in 8123001 the two halves calibrate to quite different parameter sets
(x1 33–1128 mm) with similar KGE — equifinality plus likely non-stationarity
(post-2010 mega-drought). A good benchmark case for regularized /
gradient-based calibration later.
