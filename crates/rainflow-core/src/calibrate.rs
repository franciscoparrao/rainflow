//! Automatic calibration.
//!
//! Two single-objective global optimizers, selectable through [`Optimizer`]:
//!
//! - **DDS** — Dynamically Dimensioned Search (Tolson & Shoemaker 2007, WRR 43,
//!   W01413): a greedy search that scales the number of perturbed dimensions
//!   down as the budget is consumed. Parsimonious, few control parameters.
//! - **SCE-UA** — Shuffled Complex Evolution (Duan et al. 1992, WRR 28,
//!   1015–1031): a population of complexes evolved by competitive complex
//!   evolution and periodically shuffled. More robust on multimodal surfaces,
//!   at a higher evaluation cost.
//!
//! Both algorithms now live in the **forge** optimization substrate
//! (`forge-core`); rainflow consumes them rather than re-implementing. The
//! engine uses one seedable SplitMix64 RNG, so a calibration run is fully
//! reproducible for a given seed. The thin wrappers below convert between
//! rainflow's `F: Float` parameter space and forge's `f64` search space, which
//! is sound because DDS/SCE-UA are derivative-free (no autodiff needed for the
//! optimizer itself — only the model evaluation stays generic over `F`).

use forge_core::problem::{func, Maximize};
// Imported as `_`: brings the trait's `optimize` method into scope for forge's
// `Algorithm` without colliding with rainflow's own `Optimizer` enum below.
use forge_core::Optimizer as _;
use forge_core::{Algorithm, Dds, Sce, Termination};
use num_traits::Float;

use crate::error::Error;
use crate::gr4j::{Gr4j, Gr4jParams};
use crate::hbv::{ElevationBands, Hbv, HbvBands, HbvParams};
use crate::metrics;

/// DDS configuration.
#[derive(Debug, Clone, Copy)]
pub struct DdsConfig {
    /// Total objective evaluations (iteration budget), >= 2.
    pub max_iter: usize,
    /// Neighborhood perturbation size as a fraction of each parameter range.
    /// The paper's robust default is 0.2.
    pub r: f64,
    /// RNG seed; same seed + same inputs => same result.
    pub seed: u64,
}

impl Default for DdsConfig {
    fn default() -> Self {
        Self {
            max_iter: 1000,
            r: 0.2,
            seed: 42,
        }
    }
}

/// Outcome of an optimization run (shared by DDS and SCE-UA).
#[derive(Debug, Clone)]
pub struct DdsResult<F> {
    /// Best parameter vector found.
    pub params: Vec<F>,
    /// Objective value at `params` (in the maximization sense).
    pub value: F,
    /// Objective evaluations actually performed.
    pub evaluations: usize,
}

/// Converts a forge `f64` candidate into the model's `F` parameter vector.
#[inline]
fn to_params<F: Float>(x: &[f64]) -> Vec<F> {
    x.iter()
        .map(|&v| F::from(v).expect("f64 candidate must be representable in F"))
        .collect()
}

/// Builds a forge maximization problem from `bounds` and a user `objective`,
/// runs `algo`, and maps the result back into `DdsResult<F>`.
///
/// The objective is wrapped in a `RefCell` so a `FnMut` closure can be called
/// through forge's `Fn`-based [`forge_core::problem::Problem`] (the engine runs
/// it sequentially, so no aliasing occurs). Non-finite objective values are
/// passed through as `NaN`, which forge treats as rejected candidates.
fn run_forge<F: Float>(
    bounds: &[(F, F)],
    init: Option<&[F]>,
    algo: Algorithm,
    max_iter: usize,
    objective: impl FnMut(&[F]) -> F,
) -> DdsResult<F> {
    let fbounds: Vec<(f64, f64)> = bounds
        .iter()
        .map(|&(lo, hi)| {
            (
                lo.to_f64().expect("lower bound must convert to f64"),
                hi.to_f64().expect("upper bound must convert to f64"),
            )
        })
        .collect();

    let cell = std::cell::RefCell::new(objective);
    let problem = Maximize(func(fbounds, move |x: &[f64]| {
        let xf = to_params::<F>(x);
        let v = (cell.borrow_mut())(&xf);
        v.to_f64().unwrap_or(f64::NAN)
    }));

    let term = Termination::budget(max_iter);
    let report = match algo {
        Algorithm::Dds(cfg) => {
            let init_f64: Option<Vec<f64>> =
                init.map(|x0| x0.iter().map(|&v| v.to_f64().unwrap_or(f64::NAN)).collect());
            cfg.optimize_from(&problem, &term, init_f64.as_deref())
        }
        other => other.optimize(&problem, &term),
    };

    DdsResult {
        params: to_params::<F>(report.best()),
        value: F::from(report.best_value_maximized())
            .expect("objective value must be representable in F"),
        evaluations: report.evaluations,
    }
}

/// Validates box bounds: non-empty, every `lower < upper`, no NaN.
fn check_bounds<F: Float>(bounds: &[(F, F)]) -> Result<(), Error> {
    if bounds.is_empty() {
        return Err(Error::InvalidParameter {
            name: "bounds",
            reason: "at least one parameter is required".into(),
        });
    }
    for (i, &(lo, up)) in bounds.iter().enumerate() {
        // NaN bounds must be rejected too, hence not a plain `lo >= up`.
        if lo.partial_cmp(&up) != Some(std::cmp::Ordering::Less) {
            return Err(Error::InvalidParameter {
                name: "bounds",
                reason: format!("dimension {i}: lower bound must be < upper bound"),
            });
        }
    }
    Ok(())
}

/// Maximizes `objective` within `bounds` using DDS.
///
/// `init` is the starting solution; when `None`, a uniform random point is
/// drawn. Objective evaluations returning non-finite values are treated as
/// rejected candidates, so the search tolerates numerically degenerate
/// parameter combinations.
pub fn dds_maximize<F: Float>(
    bounds: &[(F, F)],
    init: Option<&[F]>,
    config: &DdsConfig,
    objective: impl FnMut(&[F]) -> F,
) -> Result<DdsResult<F>, Error> {
    check_bounds(bounds)?;
    if let Some(x0) = init
        && x0.len() != bounds.len()
    {
        return Err(Error::InvalidParameter {
            name: "init",
            reason: format!("expected {} values, got {}", bounds.len(), x0.len()),
        });
    }
    if config.max_iter < 2 {
        return Err(Error::InvalidParameter {
            name: "max_iter",
            reason: "iteration budget must be >= 2".into(),
        });
    }
    if !(config.r > 0.0 && config.r <= 1.0) {
        return Err(Error::InvalidParameter {
            name: "r",
            reason: "perturbation fraction must be in (0, 1]".into(),
        });
    }

    let algo = Algorithm::Dds(Dds {
        r: config.r,
        seed: config.seed,
    });
    Ok(run_forge(bounds, init, algo, config.max_iter, objective))
}

/// SCE-UA configuration (Shuffled Complex Evolution; Duan et al. 1992).
#[derive(Debug, Clone, Copy)]
pub struct SceConfig {
    /// Number of complexes `p`, >= 1. More complexes = more global, slower.
    pub complexes: usize,
    /// Total objective-evaluation budget (>= initial population size).
    pub max_iter: usize,
    /// RNG seed; same seed + same inputs => same result.
    pub seed: u64,
}

impl Default for SceConfig {
    fn default() -> Self {
        Self {
            complexes: 4,
            max_iter: 10_000,
            seed: 42,
        }
    }
}

/// A single optimization run, shared by both algorithms.
#[derive(Debug, Clone, Copy)]
pub enum Optimizer {
    Dds(DdsConfig),
    Sce(SceConfig),
}

impl Optimizer {
    /// Maximizes `objective` within `bounds` with the chosen algorithm.
    pub fn maximize<F: Float>(
        &self,
        bounds: &[(F, F)],
        objective: impl FnMut(&[F]) -> F,
    ) -> Result<DdsResult<F>, Error> {
        match self {
            Optimizer::Dds(c) => dds_maximize(bounds, None, c, objective),
            Optimizer::Sce(c) => sce_maximize(bounds, c, objective),
        }
    }
}

/// Maximizes `objective` within `bounds` using SCE-UA (Duan, Sorooshian &
/// Gupta 1992, Water Resources Research 28, 1015–1031): a robust global
/// optimizer that partitions a random population into complexes, evolves each
/// by competitive complex evolution (a reflection/contraction simplex), then
/// shuffles. Deterministic for a given seed.
///
/// Non-finite objective values are treated as the worst possible (rejected),
/// so degenerate parameter combinations are tolerated. The budget `max_iter`
/// counts objective evaluations, like [`dds_maximize`].
pub fn sce_maximize<F: Float>(
    bounds: &[(F, F)],
    config: &SceConfig,
    objective: impl FnMut(&[F]) -> F,
) -> Result<DdsResult<F>, Error> {
    check_bounds(bounds)?;
    if config.complexes == 0 {
        return Err(Error::InvalidParameter {
            name: "complexes",
            reason: "need at least one complex".into(),
        });
    }

    // Standard SCE-UA complex geometry (Duan et al. 1992): m = 2n+1 points per
    // complex. The initial population must fit inside the evaluation budget.
    let n = bounds.len();
    let pop_size = config.complexes * (2 * n + 1);
    if config.max_iter < pop_size {
        return Err(Error::InvalidParameter {
            name: "max_iter",
            reason: format!("budget must be >= initial population size ({pop_size})"),
        });
    }

    let algo = Algorithm::Sce(Sce {
        complexes: config.complexes,
        seed: config.seed,
    });
    Ok(run_forge(bounds, None, algo, config.max_iter, objective))
}

/// Calibration objective for GR4J.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Objective {
    Nse,
    Kge,
    LogNse,
}

impl Objective {
    fn evaluate<F: Float>(self, obs: &[F], sim: &[F]) -> Option<F> {
        match self {
            Self::Nse => metrics::nse(obs, sim).ok(),
            Self::Kge => metrics::kge(obs, sim).ok(),
            Self::LogNse => metrics::log_nse(obs, sim).ok(),
        }
    }
}

/// Default GR4J search bounds, wide enough for most catchments
/// (x1, x2, x3, x4) = ([1, 2500] mm, [-5, 5] mm/d, [1, 1000] mm, [0.5, 10] d).
pub fn gr4j_default_bounds<F: Float>() -> [(F, F); 4] {
    let lit = |v: f64| F::from(v).expect("f64 literal must be representable in F");
    [
        (lit(1.0), lit(2500.0)),
        (lit(-5.0), lit(5.0)),
        (lit(1.0), lit(1000.0)),
        (lit(0.5), lit(10.0)),
    ]
}

/// Result of a GR4J calibration.
#[derive(Debug, Clone)]
pub struct Gr4jCalibration<F> {
    pub params: Gr4jParams<F>,
    /// Objective value over the post-warm-up period.
    pub value: F,
    pub evaluations: usize,
}

/// Calibrates GR4J on (precip, pet, qobs) maximizing `objective`, skipping
/// `warmup` initial steps in the objective. `qobs` may contain NaN gaps.
pub fn calibrate_gr4j<F: Float>(
    precip: &[F],
    pet: &[F],
    qobs: &[F],
    warmup: usize,
    objective: Objective,
    bounds: &[(F, F); 4],
    optimizer: &Optimizer,
) -> Result<Gr4jCalibration<F>, Error> {
    if precip.len() != pet.len() {
        return Err(Error::ForcingLengthMismatch {
            precip: precip.len(),
            pet: pet.len(),
        });
    }
    if qobs.len() != precip.len() {
        return Err(Error::InvalidParameter {
            name: "qobs",
            reason: format!("expected {} steps, got {}", precip.len(), qobs.len()),
        });
    }
    if warmup >= precip.len() {
        return Err(Error::InvalidParameter {
            name: "warmup",
            reason: format!(
                "warm-up ({warmup}) must be shorter than the series ({})",
                precip.len()
            ),
        });
    }

    let nan = F::nan();
    let result = optimizer.maximize(bounds, |x| {
        let Ok(model) = Gr4j::new(Gr4jParams {
            x1: x[0],
            x2: x[1],
            x3: x[2],
            x4: x[3],
        }) else {
            return nan;
        };
        let Ok(qsim) = model.run(precip, pet) else {
            return nan;
        };
        objective
            .evaluate(&qobs[warmup..], &qsim[warmup..])
            .unwrap_or(nan)
    })?;

    Ok(Gr4jCalibration {
        params: Gr4jParams {
            x1: result.params[0],
            x2: result.params[1],
            x3: result.params[2],
            x4: result.params[3],
        },
        value: result.value,
        evaluations: result.evaluations,
    })
}

/// Default HBV search bounds (HBV-light conventional ranges).
///
/// Without snow (9 parameters): FC, LP, BETA, K0, K1, K2, UZL, PERC, MAXBAS.
/// With snow (12 parameters): TT, CFMAX, SFCF prepended; CFR and CWH stay at
/// their HBV-light defaults (0.05 and 0.1).
pub fn hbv_default_bounds<F: Float>(with_snow: bool) -> Vec<(F, F)> {
    let lit = |v: f64| F::from(v).expect("f64 literal must be representable in F");
    let mut bounds = Vec::new();
    if with_snow {
        // TT and SFCF ranges are wider than the HBV-light convention
        // (±2.5 °C, 0.4–1.4): with lumped catchment-mean temperature over
        // high-relief Andean basins the effective rain/snow threshold rises,
        // and SFCF must compensate gauge/product undercatch at altitude.
        bounds.extend([
            (lit(-3.0), lit(4.0)), // TT
            (lit(0.5), lit(10.0)), // CFMAX
            (lit(0.4), lit(2.0)),  // SFCF
        ]);
    }
    bounds.extend([
        (lit(50.0), lit(700.0)), // FC
        (lit(0.3), lit(1.0)),    // LP
        (lit(1.0), lit(6.0)),    // BETA
        (lit(0.05), lit(0.99)),  // K0
        (lit(0.01), lit(0.5)),   // K1
        (lit(0.001), lit(0.2)),  // K2
        (lit(0.0), lit(100.0)),  // UZL
        // Upper bound generous: wet catchments saturate the common 0-6 range.
        (lit(0.0), lit(10.0)), // PERC
        (lit(1.0), lit(7.0)),  // MAXBAS
    ]);
    bounds
}

/// Bounds for the two lapse parameters appended when calibrating
/// [`HbvBands`]: TCALT (°C/100 m) and PCALT (fraction/100 m).
fn lapse_bounds<F: Float>() -> [(F, F); 2] {
    let lit = |v: f64| F::from(v).expect("f64 literal must be representable in F");
    [
        (lit(0.0), lit(1.2)),  // TCALT — environmental lapse rate range
        (lit(0.0), lit(0.30)), // PCALT — precipitation gradient
    ]
}

fn hbv_params_from_vector<F: Float>(x: &[F], with_snow: bool) -> HbvParams<F> {
    let lit = |v: f64| F::from(v).expect("f64 literal must be representable in F");
    let (snow, rest) = if with_snow {
        ((x[0], x[1], x[2]), &x[3..])
    } else {
        // Snow routine is bypassed without temperature; placeholders only.
        ((F::zero(), lit(3.0), F::one()), x)
    };
    HbvParams {
        tt: snow.0,
        cfmax: snow.1,
        sfcf: snow.2,
        cfr: lit(0.05),
        cwh: lit(0.1),
        fc: rest[0],
        lp: rest[1],
        beta: rest[2],
        k0: rest[3],
        k1: rest[4],
        k2: rest[5],
        uzl: rest[6],
        perc: rest[7],
        maxbas: rest[8],
    }
}

/// Result of an HBV calibration.
#[derive(Debug, Clone)]
pub struct HbvCalibration<F> {
    pub params: HbvParams<F>,
    /// Objective value over the post-warm-up period.
    pub value: F,
    pub evaluations: usize,
}

/// Calibrates HBV-light on (precip, pet, temp, qobs) maximizing `objective`,
/// skipping `warmup` initial steps. When `temp` is `Some`, the snow-routine
/// parameters (TT, CFMAX, SFCF) are calibrated too. `qobs` may contain NaN.
pub fn calibrate_hbv<F: Float>(
    precip: &[F],
    pet: &[F],
    temp: Option<&[F]>,
    qobs: &[F],
    warmup: usize,
    objective: Objective,
    optimizer: &Optimizer,
) -> Result<HbvCalibration<F>, Error> {
    if precip.len() != pet.len() {
        return Err(Error::ForcingLengthMismatch {
            precip: precip.len(),
            pet: pet.len(),
        });
    }
    if qobs.len() != precip.len() {
        return Err(Error::InvalidParameter {
            name: "qobs",
            reason: format!("expected {} steps, got {}", precip.len(), qobs.len()),
        });
    }
    if warmup >= precip.len() {
        return Err(Error::InvalidParameter {
            name: "warmup",
            reason: format!(
                "warm-up ({warmup}) must be shorter than the series ({})",
                precip.len()
            ),
        });
    }

    let with_snow = temp.is_some();
    let bounds = hbv_default_bounds::<F>(with_snow);
    let nan = F::nan();
    let result = optimizer.maximize(&bounds, |x| {
        let Ok(model) = Hbv::new(hbv_params_from_vector(x, with_snow)) else {
            return nan;
        };
        let Ok(qsim) = model.run(precip, pet, temp) else {
            return nan;
        };
        objective
            .evaluate(&qobs[warmup..], &qsim[warmup..])
            .unwrap_or(nan)
    })?;

    Ok(HbvCalibration {
        params: hbv_params_from_vector(&result.params, with_snow),
        value: result.value,
        evaluations: result.evaluations,
    })
}

/// Result of a semi-distributed calibration: HBV parameters plus the fitted
/// lapse rates and the band geometry they apply to.
#[derive(Debug, Clone)]
pub struct HbvBandsCalibration<F> {
    pub params: HbvParams<F>,
    pub bands: ElevationBands<F>,
    pub value: F,
    pub evaluations: usize,
}

/// Calibrates the semi-distributed [`HbvBands`] model (temperature required).
/// Searches the full 12-parameter HBV set plus the two lapse rates TCALT and
/// PCALT (14 dimensions); the band geometry (elevations and area fractions)
/// stays fixed at `bands`, only its lapse rates are overwritten.
// One argument more than `calibrate_hbv` (the band geometry); kept as a flat
// signature to mirror the sibling calibration functions.
#[allow(clippy::too_many_arguments)]
pub fn calibrate_hbv_bands<F: Float>(
    precip: &[F],
    pet: &[F],
    temp: &[F],
    bands: &ElevationBands<F>,
    qobs: &[F],
    warmup: usize,
    objective: Objective,
    optimizer: &Optimizer,
) -> Result<HbvBandsCalibration<F>, Error> {
    if precip.len() != pet.len() {
        return Err(Error::ForcingLengthMismatch {
            precip: precip.len(),
            pet: pet.len(),
        });
    }
    if qobs.len() != precip.len() || temp.len() != precip.len() {
        return Err(Error::InvalidParameter {
            name: "qobs/temp",
            reason: format!(
                "expected {} steps, got qobs={} temp={}",
                precip.len(),
                qobs.len(),
                temp.len()
            ),
        });
    }
    if warmup >= precip.len() {
        return Err(Error::InvalidParameter {
            name: "warmup",
            reason: format!(
                "warm-up ({warmup}) must be shorter than the series ({})",
                precip.len()
            ),
        });
    }

    // 12 HBV parameters followed by [TCALT, PCALT].
    let mut bounds = hbv_default_bounds::<F>(true);
    bounds.extend_from_slice(&lapse_bounds::<F>());
    let lapse = bounds.len() - 2;

    // Reuse the supplied geometry, overriding only the lapse rates per trial.
    let with_lapse = |x: &[F]| -> ElevationBands<F> {
        let mut b = bands.clone();
        b.tcalt = x[lapse];
        b.pcalt = x[lapse + 1];
        b
    };

    let nan = F::nan();
    let result = optimizer.maximize(&bounds, |x| {
        let Ok(model) = HbvBands::new(hbv_params_from_vector(x, true), with_lapse(x)) else {
            return nan;
        };
        let Ok(qsim) = model.run(precip, pet, temp) else {
            return nan;
        };
        objective
            .evaluate(&qobs[warmup..], &qsim[warmup..])
            .unwrap_or(nan)
    })?;

    Ok(HbvBandsCalibration {
        params: hbv_params_from_vector(&result.params, true),
        bands: with_lapse(&result.params),
        value: result.value,
        evaluations: result.evaluations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dds_finds_the_sphere_optimum() {
        // Maximize -(x-3)^2 - (y+1)^2 on [-10, 10]^2; optimum at (3, -1).
        let bounds = [(-10.0, 10.0), (-10.0, 10.0)];
        let config = DdsConfig {
            max_iter: 2000,
            ..Default::default()
        };
        let res = dds_maximize(&bounds, None, &config, |x| {
            -(x[0] - 3.0).powi(2) - (x[1] + 1.0).powi(2)
        })
        .unwrap();
        assert!(res.value > -1e-2, "objective {}", res.value);
        assert!((res.params[0] - 3.0).abs() < 0.1);
        assert!((res.params[1] + 1.0).abs() < 0.1);
        assert_eq!(res.evaluations, 2000);
    }

    #[test]
    fn dds_is_reproducible_for_a_seed() {
        let bounds = [(-5.0, 5.0); 3];
        let config = DdsConfig::default();
        let f = |x: &[f64]| -x.iter().map(|v| v * v).sum::<f64>();
        let a = dds_maximize(&bounds, None, &config, f).unwrap();
        let b = dds_maximize(&bounds, None, &config, f).unwrap();
        assert_eq!(a.params, b.params);
        assert_eq!(a.value, b.value);
    }

    #[test]
    fn dds_respects_bounds() {
        let bounds = [(0.0, 1.0), (10.0, 20.0)];
        let config = DdsConfig {
            max_iter: 500,
            ..Default::default()
        };
        dds_maximize(&bounds, None, &config, |x| {
            assert!((0.0..=1.0).contains(&x[0]), "x0 out of bounds: {}", x[0]);
            assert!((10.0..=20.0).contains(&x[1]), "x1 out of bounds: {}", x[1]);
            0.0
        })
        .unwrap();
    }

    #[test]
    fn dds_rejects_bad_configuration() {
        let f = |_: &[f64]| 0.0;
        assert!(dds_maximize(&[], None, &DdsConfig::default(), f).is_err());
        assert!(dds_maximize(&[(1.0, 0.0)], None, &DdsConfig::default(), f).is_err());
        let bad_iter = DdsConfig {
            max_iter: 1,
            ..Default::default()
        };
        assert!(dds_maximize(&[(0.0, 1.0)], None, &bad_iter, f).is_err());
        assert!(dds_maximize(&[(0.0, 1.0)], Some(&[0.5, 0.5]), &DdsConfig::default(), f).is_err());
    }

    #[test]
    fn calibration_recovers_a_synthetic_truth() {
        // Generate qobs with known parameters; DDS must reach NSE ≈ 1.
        let truth = Gr4j::new(Gr4jParams {
            x1: 350.0,
            x2: -1.5,
            x3: 90.0,
            x4: 1.7,
        })
        .unwrap();
        let (p, pet) = synthetic_forcing(1500);
        let qobs = truth.run(&p, &pet).unwrap();

        let opt = Optimizer::Dds(DdsConfig {
            max_iter: 1500,
            ..Default::default()
        });
        let cal = calibrate_gr4j(
            &p,
            &pet,
            &qobs,
            365,
            Objective::Nse,
            &gr4j_default_bounds(),
            &opt,
        )
        .unwrap();
        assert!(cal.value > 0.95, "calibrated NSE too low: {}", cal.value);
    }

    #[test]
    fn sce_finds_the_sphere_optimum() {
        // Same surface as the DDS test; SCE-UA should also locate (3, -1).
        let bounds = [(-10.0, 10.0), (-10.0, 10.0)];
        let config = SceConfig {
            complexes: 3,
            max_iter: 3000,
            seed: 42,
        };
        let res = sce_maximize(&bounds, &config, |x| {
            -(x[0] - 3.0).powi(2) - (x[1] + 1.0).powi(2)
        })
        .unwrap();
        assert!(res.value > -1e-2, "objective {}", res.value);
        assert!((res.params[0] - 3.0).abs() < 0.1);
        assert!((res.params[1] + 1.0).abs() < 0.1);
    }

    #[test]
    fn sce_is_reproducible_and_respects_bounds() {
        let bounds = [(0.0, 1.0), (10.0, 20.0), (-5.0, 5.0)];
        let config = SceConfig::default();
        let f = |x: &[f64]| {
            assert!((0.0..=1.0).contains(&x[0]) && (10.0..=20.0).contains(&x[1]));
            -x.iter().map(|v| v * v).sum::<f64>()
        };
        let a = sce_maximize(&bounds, &config, f).unwrap();
        let b = sce_maximize(&bounds, &config, f).unwrap();
        assert_eq!(a.params, b.params);
        assert_eq!(a.value, b.value);
    }

    #[test]
    fn sce_rejects_bad_configuration() {
        let f = |_: &[f64]| 0.0;
        assert!(sce_maximize(&[], &SceConfig::default(), f).is_err());
        assert!(sce_maximize(&[(1.0, 0.0)], &SceConfig::default(), f).is_err());
        // Budget below the initial population must be rejected.
        let tiny = SceConfig {
            complexes: 4,
            max_iter: 5,
            seed: 1,
        };
        assert!(sce_maximize(&[(0.0, 1.0)], &tiny, f).is_err());
    }

    #[test]
    fn dds_and_sce_agree_on_gr4j_calibration() {
        let truth = Gr4j::new(Gr4jParams {
            x1: 350.0,
            x2: -1.5,
            x3: 90.0,
            x4: 1.7,
        })
        .unwrap();
        let (p, pet) = synthetic_forcing(1500);
        let qobs = truth.run(&p, &pet).unwrap();
        let sce = Optimizer::Sce(SceConfig {
            complexes: 4,
            max_iter: 4000,
            seed: 42,
        });
        let cal = calibrate_gr4j(
            &p,
            &pet,
            &qobs,
            365,
            Objective::Nse,
            &gr4j_default_bounds(),
            &sce,
        )
        .unwrap();
        assert!(
            cal.value > 0.95,
            "SCE-UA calibrated NSE too low: {}",
            cal.value
        );
    }

    #[test]
    fn hbv_calibration_recovers_a_synthetic_truth() {
        let truth = Hbv::new(HbvParams {
            tt: 0.0,
            cfmax: 3.0,
            sfcf: 1.0,
            cfr: 0.05,
            cwh: 0.1,
            fc: 250.0,
            lp: 0.7,
            beta: 2.0,
            k0: 0.3,
            k1: 0.1,
            k2: 0.01,
            uzl: 20.0,
            perc: 2.0,
            maxbas: 2.5,
        })
        .unwrap();
        let (p, pet) = synthetic_forcing(1500);
        let qobs = truth.run(&p, &pet, None).unwrap();

        let opt = Optimizer::Dds(DdsConfig {
            max_iter: 2000,
            ..Default::default()
        });
        let cal = calibrate_hbv(&p, &pet, None, &qobs, 365, Objective::Nse, &opt).unwrap();
        assert!(cal.value > 0.9, "calibrated NSE too low: {}", cal.value);
    }

    /// Same LCG-based forcing generator as the GR4J tests.
    fn synthetic_forcing(n: usize) -> (Vec<f64>, Vec<f64>) {
        let mut seed: u64 = 42;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as f64 / (1u64 << 31) as f64
        };
        let mut p = Vec::with_capacity(n);
        let mut pet = Vec::with_capacity(n);
        for i in 0..n {
            let doy = (i % 365) as f64;
            let wet = if (120.0..270.0).contains(&doy) {
                0.55
            } else {
                0.15
            };
            p.push(if next() < wet {
                -12.0 * (1.0 - next()).ln()
            } else {
                0.0
            });
            pet.push(
                (3.5 + 2.5 * (2.0 * std::f64::consts::PI * (doy - 15.0) / 365.25).sin()).max(0.1),
            );
        }
        (p, pet)
    }
}
