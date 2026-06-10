# CAMELS-CL test catchments

Per-catchment daily forcing extracted from **CAMELS-CL** (Alvarez-Garreton et
al. 2018, HESS 22, 5817–5846; data: https://doi.org/10.1594/PANGAEA.894885,
CC-BY-4.0). Columns: `date`, `p` (precipitation, CR2MET, mm), `pet`
(Hargreaves, mm), `qobs` (observed streamflow, mm; `NA` = gap).
Period 1979-01-01 to 2016-12-31 (13,880 days).

| gauge_id | name | area km² | snow_frac | interv_degree | qobs coverage |
|---|---|---|---|---|---|
| 8123001 | Río Itata en Cholguán | 860 | 0.00 | 0.003 | 95% |
| 7330001 | Río Perquilauquén en San Manuel | 502 | 0.00 | 0.005 | 94% |

Selection criteria: near-natural (no big dams, interv_degree < 0.05),
rain-dominated (snow_frac & swe_ratio < 0.05 — GR4J has no snow module yet),
long records (n_obs > 17,000 days over the full gauge history).

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

Reproduce with:

```sh
rainflow split-sample --forcing data/camels-cl/7330001.csv \
    --objective kge --iterations 3000 --seed 42
```

Note: in 8123001 the two halves calibrate to quite different parameter sets
(x1 33–1128 mm) with similar KGE — equifinality plus likely non-stationarity
(post-2010 mega-drought). A good benchmark case for regularized /
gradient-based calibration later.
