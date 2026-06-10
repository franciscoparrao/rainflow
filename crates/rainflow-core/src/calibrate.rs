//! Automatic calibration.
//!
//! Currently implements DDS — Dynamically Dimensioned Search (Tolson &
//! Shoemaker 2007, Water Resources Research 43, W01413): a single-objective,
//! greedy stochastic search that scales the number of perturbed dimensions
//! down as the iteration budget is consumed. Parsimonious and well suited to
//! conceptual rainfall–runoff models.
//!
//! The search is fully deterministic for a given seed (own SplitMix64 RNG, no
//! external randomness), which keeps calibration runs reproducible.

use num_traits::Float;

use crate::error::Error;
use crate::gr4j::{Gr4j, Gr4jParams};
use crate::metrics;

/// Deterministic SplitMix64 RNG — small, seedable, good enough for DDS.
struct SplitMix64(u64);

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Standard normal via Box–Muller.
    fn normal(&mut self) -> f64 {
        let u1 = (1.0 - self.uniform()).max(f64::MIN_POSITIVE); // avoid ln(0)
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

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

/// Outcome of a DDS run.
#[derive(Debug, Clone)]
pub struct DdsResult<F> {
    /// Best parameter vector found.
    pub params: Vec<F>,
    /// Objective value at `params`.
    pub value: F,
    /// Objective evaluations actually performed.
    pub evaluations: usize,
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
    mut objective: impl FnMut(&[F]) -> F,
) -> Result<DdsResult<F>, Error> {
    let dim = bounds.len();
    if dim == 0 {
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
    if let Some(x0) = init
        && x0.len() != dim
    {
        return Err(Error::InvalidParameter {
            name: "init",
            reason: format!("expected {dim} values, got {}", x0.len()),
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

    let mut rng = SplitMix64::new(config.seed);
    let lit = |v: f64| F::from(v).expect("f64 literal must be representable in F");

    let mut best: Vec<F> = match init {
        Some(x0) => x0.to_vec(),
        None => bounds
            .iter()
            .map(|&(lo, up)| lo + (up - lo) * lit(rng.uniform()))
            .collect(),
    };
    let mut best_value = objective(&best);
    let mut evaluations = 1;

    let m = config.max_iter as f64;
    let mut candidate = best.clone();
    for i in 1..config.max_iter {
        // P(perturb dim) decays from ~1 to 1/m over the budget.
        let p = 1.0 - (i as f64).ln() / m.ln();

        candidate.copy_from_slice(&best);
        let mut perturbed = 0;
        for (j, &(lo, up)) in bounds.iter().enumerate() {
            if rng.uniform() < p {
                candidate[j] = perturb(best[j], lo, up, config.r, &mut rng);
                perturbed += 1;
            }
        }
        if perturbed == 0 {
            // Always perturb at least one randomly chosen dimension.
            let j = (rng.next_u64() % dim as u64) as usize;
            let (lo, up) = bounds[j];
            candidate[j] = perturb(best[j], lo, up, config.r, &mut rng);
        }

        let value = objective(&candidate);
        evaluations += 1;
        if value.is_finite() && (!best_value.is_finite() || value > best_value) {
            best.copy_from_slice(&candidate);
            best_value = value;
        }
    }

    Ok(DdsResult {
        params: best,
        value: best_value,
        evaluations,
    })
}

/// One-dimensional DDS neighborhood move with boundary reflection.
fn perturb<F: Float>(x: F, lo: F, up: F, r: f64, rng: &mut SplitMix64) -> F {
    let lit = |v: f64| F::from(v).expect("f64 literal must be representable in F");
    let range = up - lo;
    let mut xn = x + range * lit(r * rng.normal());
    // Reflect once at each bound; if still outside, clamp to the bound
    // (Tolson & Shoemaker 2007, eq. 4).
    if xn < lo {
        xn = lo + (lo - xn);
        if xn > up {
            xn = lo;
        }
    } else if xn > up {
        xn = up - (xn - up);
        if xn < lo {
            xn = up;
        }
    }
    xn
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
    config: &DdsConfig,
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
    let result = dds_maximize(bounds, None, config, |x| {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rng_is_deterministic_and_uniform_in_range() {
        let mut a = SplitMix64::new(7);
        let mut b = SplitMix64::new(7);
        for _ in 0..1000 {
            let v = a.uniform();
            assert_eq!(v, b.uniform());
            assert!((0.0..1.0).contains(&v));
        }
    }

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

        let config = DdsConfig {
            max_iter: 1500,
            ..Default::default()
        };
        let cal = calibrate_gr4j(
            &p,
            &pet,
            &qobs,
            365,
            Objective::Nse,
            &gr4j_default_bounds(),
            &config,
        )
        .unwrap();
        assert!(cal.value > 0.95, "calibrated NSE too low: {}", cal.value);
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
