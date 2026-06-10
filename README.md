# rainflow

Conceptual rainfall–runoff models in Rust — fast, operational, autodiff-first.

`rainflow` fills the gap between heavy physically-based distributed solvers and
the R ecosystem of conceptual models (airGR, TUWmodel): parsimonious
lumped/semi-distributed models (GR4J, HBV) with automatic calibration and
goodness-of-fit metrics, built for operational forecasting and massive
multi-catchment runs.

## Status (v0.1, in development)

- [x] **GR4J** (Perrin et al. 2003), generic over `Float` (autodiff-ready)
- [x] Metrics: **NSE, KGE (+components), logNSE, PBIAS**
- [x] CSV forcing + CLI runner
- [x] **Numerical parity with airGR**: max abs diff 6e-7 mm over 10,593 daily
      steps of the L0123001 example catchment (CSV round-off level)
- [ ] HBV-light core
- [ ] Calibration: DDS, SCE-UA
- [ ] Split-sample validation; CAMELS-CL cases
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
