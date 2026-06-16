# rainflow

Conceptual rainfall–runoff models in Rust — fast, operational, autodiff-first.

`rainflow` fills the gap between heavy physically-based distributed solvers and
the R ecosystem of conceptual models (airGR, TUWmodel): parsimonious
lumped/semi-distributed models (GR4J, HBV) with automatic calibration and
goodness-of-fit metrics, built for operational forecasting and massive
multi-catchment runs.

## Status (v0.1, in development)

- [x] **GR4J** (Perrin et al. 2003), generic over `Float` (autodiff-ready)
- [x] **HBV-light** (Seibert & Vis 2012): degree-day snow routine (optional
      temperature forcing), soil moisture accounting, two-box response,
      MAXBAS routing. Beats GR4J on the CAMELS-CL test catchments (val NSE
      0.73–0.77 vs 0.64–0.74)
- [x] Metrics: **NSE, KGE (+components), logNSE, PBIAS**
- [x] CSV forcing + CLI runner
- [x] **Numerical parity with airGR**: max abs diff 6e-7 mm over 10,593 daily
      steps of the L0123001 example catchment (CSV round-off level)
- [x] Calibration: **DDS** (Tolson & Shoemaker 2007) and **SCE-UA** (Duan et
      al. 1992), both deterministic per seed and selectable per run. DDS on the
      L0123001 catchment converges to the same optimum as
      `airGR::Calibration_Michel` (NSE 0.7956 vs 0.7957); DDS and SCE-UA agree
      to within rounding on the GR4J optima
- [x] **Split-sample validation on CAMELS-CL**: two near-natural pluvial
      catchments (Río Itata en Cholguán, Río Perquilauquén en San Manuel,
      1979–2016). Validation KGE 0.76–0.82 — see `data/camels-cl/README.md`
- [x] **Snow routine validated end-to-end** on two snow-dominated Andean
      catchments (Río Grande en Las Ramadas, Río Choapa en Cuncumén): the
      degree-day routine lifts validation NSE from ≤ 0.23 (no snow) to
      0.31–0.62; remaining gap is structural (lumped temperature → elevation
      bands planned for v0.2)
- [x] **Semi-distributed HBV with elevation bands** (`--model hbv-bands`):
      per-band snow + soil with temperature lapse (TCALT) and precipitation
      gradient (PCALT), shared response and routing. Single band at the
      reference elevation reproduces the lumped model exactly. With the lapse
      rates calibrated, the bands beat the lumped model robustly on Río Choapa
      en Cuncumén (validation NSE 0.63–0.76 vs 0.34–0.62) — see
      `data/camels-cl/README.md`
- [x] **Hypsometric (equal-area) band geometry** (`--hypsometry "min,median,max"`):
      band elevations read off the catchment's hypsometric curve instead of
      hand-picked. Matches or beats hand-tuned bands (+0.05–0.11 validation NSE
      on Río Grande en Las Ramadas) with no guesswork. Core constructor
      `equal_area_from_hypsometry` also accepts a dense DEM-sampled curve
- [x] **Python bindings** (PyO3 + maturin, `rainflow` on PyPI-style wheels):
      `Gr4j`, `Hbv`, the metrics, `calibrate_gr4j`/`calibrate_hbv` (DDS or
      SCE-UA) and `hypsometric_bands`. One abi3 wheel covers CPython ≥ 3.9
- [x] **CI** (GitHub Actions): fmt + clippy (`-D warnings`) + tests, plus a
      job that builds and smoke-tests the Python wheel
- [x] **DEM-derived hypsometric curves** (`--hypsometry-file`): real curve from
      Copernicus DEM GLO-30 clipped to the catchment polygon
      (`scripts/hypsometry_from_dem.py`); the clip reproduces CAMELS-CL
      elevation attributes to a few metres. Objective band geometry, same core
      constructor as the quantile reconstruction
- [x] **Spatial routing by subcatchments** (`routing` module): Muskingum
      channel routing + a drainage-tree `RiverNetwork` that accumulates and
      routes subcatchments to the outlet. On the nested Río Itata (Cholguán →
      Balsa Nueva Aldea) the 2-subcatchment + routing model beats a lumped GR4J
      by +0.06–0.08 validation NSE (`examples/route_itata.rs`)
- [x] **Parallel multi-catchment calibration** (`rainflow batch`, rayon): the
      core calibration functions are pure, so N catchments calibrate
      independently across cores — 3.86× on four catchments (16 threads), the
      operational path for the 15 BNA catchments. rayon lives only in the CLI
- [x] **Warm-start / stateful forecasting** (`run_from`): run from a given
      model state instead of the default, advancing it in place. Exact state
      continuity (a split run equals one continuous run to 1e-12). Exposed in
      Python with opaque `Gr4jState`/`HbvState` objects, so a nowcast can
      settle the state on history, snapshot it, and fan out forecast scenarios
      from the snapshot. With the `serde` feature, states (and params) serialize
      to JSON — `to_json`/`from_json` in Python — so an operational nowcast can
      persist the state to disk and resume bit-for-bit
- [x] **Gradient calibration via autodiff** — the autodiff-first design pays
      off: the `autodiff` crate's forward-mode scalar implements
      `num_traits::Float`, so it flows through `Gr4j` unchanged and yields an
      analytic loss gradient (matches finite differences to ~1e-9). Adam
      recovers the true parameters exactly (`examples/gradient_calibration.rs`).
      This is the substrate for physics+ML hybrids (δHBV-style)

## Layout

- `crates/rainflow-core` — model cores, state, metrics. No I/O, `#![no_panic]`-minded,
  everything generic over `F: num_traits::Float` so dual-number/tape scalar types
  pass through for gradient-based calibration.
- `crates/rainflow-cli` — `rainflow` binary (CSV in, simulated discharge + metrics out).
- `crates/rainflow-python` — PyO3 bindings (`rainflow` module), built with maturin.

## Quick start

```sh
cargo build --release
./target/release/rainflow run \
    --forcing data/example.csv \
    --x1 350 --x2 -1.5 --x3 90 --x4 1.7 \
    --warmup 365 --output qsim.csv

# Automatic calibration (requires a qobs column)
./target/release/rainflow calibrate \
    --forcing data/example.csv \
    --objective kge --iterations 2000 --seed 42

# Calibrate many catchments in parallel (one CSV each)
./target/release/rainflow batch \
    --forcing basin_a.csv --forcing basin_b.csv --forcing basin_c.csv \
    --model gr4j --objective kge --output summary.csv
```

The forcing CSV needs columns (flexible, case-insensitive names):
`date`, `p` (precipitation, mm), `pet` (potential evapotranspiration, mm) and
optionally `qobs` (observed discharge, mm) to compute metrics. Gaps (`NA`) are
allowed only in `qobs`.

## Python

```sh
cd crates/rainflow-python
maturin develop --release   # or: maturin build --release --out dist
```

```python
import rainflow

q = rainflow.Gr4j(350.0, -1.5, 90.0, 1.7).run(precip, pet)
print(rainflow.kge(qobs, q))

cal = rainflow.calibrate_gr4j(precip, pet, qobs, algorithm="sce", iterations=4000)
print(cal["params"], cal["value"])

# HBV with the snow routine (pass temperature)
qh = rainflow.Hbv(0.0, 3.5, 0.9, 250.0, 0.7, 2.0, 0.3, 0.1, 0.01, 20.0, 2.0, 2.5).run(
    precip, pet, temp=temperature
)

# Warm-start / chained forecasting: settle the state on the historical
# period, snapshot it, then fan out forecast scenarios from the snapshot.
g = rainflow.Gr4j(350.0, -1.5, 90.0, 1.7)
hist_q, state = g.run_from(g.initial_state(), hist_p, hist_pet)
dry_q, _ = g.run_from(state, dry_p, dry_pet)    # state is left unchanged,
wet_q, _ = g.run_from(state, wet_p, wet_pet)    # so it is reusable

# Persist the state to disk and resume in a later run (operational nowcast)
open("state.json", "w").write(state.to_json())
# ... next run ...
state = rainflow.Gr4jState.from_json(open("state.json").read())
tomorrow_q, state = g.run_from(state, today_p, today_pet)
```

## Performance

`cargo bench -p rainflow-core` (criterion). Figures below on a 12th-gen Intel
i7-1270P, `--release` with LTO:

| benchmark | time | throughput |
|---|---|---|
| GR4J, 10k days (~27 yr) | ~1.5 ms | ~6.4 M steps/s |
| HBV + snow, 10k days | ~0.68 ms | ~14.6 M steps/s |
| HBV, 5 elevation bands, 10k days | ~1.9 ms | ~5.4 M steps/s |
| Muskingum reach, 10k days | ~78 µs | ~129 M steps/s |
| GR4J calibration (DDS, 2000 evals, 6000 days) | ~5.0 s | — |
| GR4J calibration (SCE-UA, 4000 evals, 6000 days) | ~4.6 s | — |

GR4J is slower than HBV per step despite being simpler — it pays a `tanh` and
two `powf` calls per time step where HBV is mostly add/multiply.

Calibration is a few seconds per catchment, and catchments are independent, so
`rainflow batch` calibrates them in parallel (rayon). Four CAMELS-CL catchments
(DDS, 2000 evals each): **2m51s serial → 44s on 16 threads, a 3.86× speedup**,
near-linear in the four jobs. The 15 BNA catchments scale further on more cores.

## Validation

`crates/rainflow-core/tests/airgr_parity.rs` locks GR4J output to a reference
simulation generated with `airGR::RunModel_GR4J` (INRAE). Regenerate the
fixture with R + airGR if the formulation ever changes.

## License

MIT OR Apache-2.0. The airGR parity fixture derives from the GPL-2 airGR
package's example dataset and is used for testing only.
