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
- [x] Calibration: **DDS** (Tolson & Shoemaker 2007), deterministic per seed.
      On the L0123001 catchment it converges to the same optimum as
      `airGR::Calibration_Michel` (NSE 0.7956 vs 0.7957, near-identical
      parameters; 2,000 evaluations over 29 years of daily data in ~3 s)
- [x] **Split-sample validation on CAMELS-CL**: two near-natural pluvial
      catchments (Río Itata en Cholguán, Río Perquilauquén en San Manuel,
      1979–2016). Validation KGE 0.76–0.82 — see `data/camels-cl/README.md`
- [ ] SCE-UA
- [ ] Semi-distributed mode (subcatchments) + snow module (v0.2)

## Layout

- `crates/rainflow-core` — model cores, state, metrics. No I/O, `#![no_panic]`-minded,
  everything generic over `F: num_traits::Float` so dual-number/tape scalar types
  pass through for gradient-based calibration.
- `crates/rainflow-cli` — `rainflow` binary (CSV in, simulated discharge + metrics out).

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
```

The forcing CSV needs columns (flexible, case-insensitive names):
`date`, `p` (precipitation, mm), `pet` (potential evapotranspiration, mm) and
optionally `qobs` (observed discharge, mm) to compute metrics. Gaps (`NA`) are
allowed only in `qobs`.

## Validation

`crates/rainflow-core/tests/airgr_parity.rs` locks GR4J output to a reference
simulation generated with `airGR::RunModel_GR4J` (INRAE). Regenerate the
fixture with R + airGR if the formulation ever changes.

## License

MIT OR Apache-2.0. The airGR parity fixture derives from the GPL-2 airGR
package's example dataset and is used for testing only.
