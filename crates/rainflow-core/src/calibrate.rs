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
//! Both use the same internal SplitMix64 RNG (no external randomness), so a
//! calibration run is fully reproducible for a given seed.

use num_traits::Float;

use crate::error::Error;
use crate::gr4j::{Gr4j, Gr4jParams};
use crate::hbv::{ElevationBands, Hbv, HbvBands, HbvParams};
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
    mut objective: impl FnMut(&[F]) -> F,
) -> Result<DdsResult<F>, Error> {
    let n = bounds.len();
    if n == 0 {
        return Err(Error::InvalidParameter {
            name: "bounds",
            reason: "at least one parameter is required".into(),
        });
    }
    for (i, &(lo, up)) in bounds.iter().enumerate() {
        if lo.partial_cmp(&up) != Some(std::cmp::Ordering::Less) {
            return Err(Error::InvalidParameter {
                name: "bounds",
                reason: format!("dimension {i}: lower bound must be < upper bound"),
            });
        }
    }
    if config.complexes == 0 {
        return Err(Error::InvalidParameter {
            name: "complexes",
            reason: "need at least one complex".into(),
        });
    }

    // Standard SCE-UA complex geometry (Duan et al. 1992).
    let m = 2 * n + 1; // points per complex
    let q = n + 1; // points per sub-complex (simplex)
    let beta = m; // evolution steps per complex per shuffle
    let pop_size = config.complexes * m;
    if config.max_iter < pop_size {
        return Err(Error::InvalidParameter {
            name: "max_iter",
            reason: format!("budget must be >= initial population size ({pop_size})"),
        });
    }

    let mut rng = SplitMix64::new(config.seed);
    let lit = |v: f64| F::from(v).expect("f64 literal must be representable in F");
    // Minimize -objective; rejected (non-finite) candidates score +inf.
    let inf = F::infinity();
    let mut eval = |x: &[F], objective: &mut dyn FnMut(&[F]) -> F| -> F {
        let v = objective(x);
        if v.is_finite() { -v } else { inf }
    };

    let sample = |rng: &mut SplitMix64| -> Vec<F> {
        bounds
            .iter()
            .map(|&(lo, up)| lo + (up - lo) * lit(rng.uniform()))
            .collect()
    };

    // Initial population (each entry: parameter vector + its score).
    let mut pop: Vec<(Vec<F>, F)> = Vec::with_capacity(pop_size);
    for _ in 0..pop_size {
        let x = sample(&mut rng);
        let f = eval(&x, &mut objective);
        pop.push((x, f));
    }
    let mut evaluations = pop_size;

    let cmp = |a: &(Vec<F>, F), b: &(Vec<F>, F)| {
        a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
    };

    while evaluations < config.max_iter {
        pop.sort_by(cmp);

        // Partition the rank-sorted population into p complexes round-robin.
        for k in 0..config.complexes {
            let mut complex: Vec<(Vec<F>, F)> = (0..m)
                .map(|j| pop[k + j * config.complexes].clone())
                .collect();

            // Competitive Complex Evolution: evolve the complex `beta` times.
            for _ in 0..beta {
                if evaluations >= config.max_iter {
                    break;
                }
                complex.sort_by(cmp);
                // Triangular selection: bias parents toward better points.
                let parents = triangular_subset(m, q, &mut rng);
                let worst_idx = *parents.last().unwrap();

                // Centroid of the q-1 best parents (exclude the worst).
                let mut centroid = vec![F::zero(); n];
                for &pi in &parents[..q - 1] {
                    for d in 0..n {
                        centroid[d] = centroid[d] + complex[pi].0[d];
                    }
                }
                let qm1 = lit((q - 1) as f64);
                for c in &mut centroid {
                    *c = *c / qm1;
                }

                // Reflection of the worst through the centroid.
                let worst = complex[worst_idx].0.clone();
                let mut trial: Vec<F> = (0..n)
                    .map(|d| centroid[d] + (centroid[d] - worst[d]))
                    .collect();
                let in_bounds = trial
                    .iter()
                    .zip(bounds)
                    .all(|(&v, &(lo, up))| v >= lo && v <= up);

                let mut tf;
                if in_bounds {
                    tf = eval(&trial, &mut objective);
                    evaluations += 1;
                    if tf >= complex[worst_idx].1 {
                        // Reflection failed: try contraction toward centroid.
                        let half = lit(0.5);
                        trial = (0..n).map(|d| (centroid[d] + worst[d]) * half).collect();
                        tf = eval(&trial, &mut objective);
                        evaluations += 1;
                        if tf >= complex[worst_idx].1 {
                            // Contraction failed too: random mutation in hull.
                            trial = sample_in_hull(&complex, bounds, &mut rng);
                            tf = eval(&trial, &mut objective);
                            evaluations += 1;
                        }
                    }
                } else {
                    // Out of bounds: random mutation in the complex's hull.
                    trial = sample_in_hull(&complex, bounds, &mut rng);
                    tf = eval(&trial, &mut objective);
                    evaluations += 1;
                }

                complex[worst_idx] = (trial, tf);
            }

            // Shuffle the evolved complex back into the population.
            for (j, point) in complex.into_iter().enumerate() {
                pop[k + j * config.complexes] = point;
            }
        }
    }

    pop.sort_by(cmp);
    let best = &pop[0];
    Ok(DdsResult {
        params: best.0.clone(),
        value: -best.1, // back to maximization sense
        evaluations,
    })
}

/// Picks `q` distinct indices from `0..m` with a triangular probability that
/// favors lower (better, since the complex is rank-sorted) indices, returned
/// sorted ascending (Duan et al. 1992, eq. for the trapezoidal distribution).
fn triangular_subset(m: usize, q: usize, rng: &mut SplitMix64) -> Vec<usize> {
    let mut chosen = Vec::with_capacity(q);
    let mf = m as f64;
    while chosen.len() < q {
        // P(i) = 2(m+1-i-1)/(m(m+1)); inverse-CDF sampling on rank.
        let u = rng.uniform();
        let idx = ((mf + 0.5) - ((mf + 0.5).powi(2) - mf * (mf + 1.0) * u).sqrt()) as usize;
        let idx = idx.min(m - 1);
        if !chosen.contains(&idx) {
            chosen.push(idx);
        }
    }
    chosen.sort_unstable();
    chosen
}

/// Uniform sample inside the smallest axis-aligned box containing the complex,
/// clamped to the global bounds (SCE-UA's mutation step).
fn sample_in_hull<F: Float>(
    complex: &[(Vec<F>, F)],
    bounds: &[(F, F)],
    rng: &mut SplitMix64,
) -> Vec<F> {
    let lit = |v: f64| F::from(v).expect("f64 literal must be representable in F");
    let n = bounds.len();
    (0..n)
        .map(|d| {
            let mut lo = complex[0].0[d];
            let mut hi = lo;
            for point in complex {
                lo = lo.min(point.0[d]);
                hi = hi.max(point.0[d]);
            }
            let v = lo + (hi - lo) * lit(rng.uniform());
            v.max(bounds[d].0).min(bounds[d].1)
        })
        .collect()
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
