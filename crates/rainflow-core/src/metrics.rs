//! Goodness-of-fit metrics for hydrological simulation.
//!
//! All metrics ignore time steps where either series is non-finite (missing
//! observations encoded as NaN), pairing the remaining values.

use num_traits::Float;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum MetricsError {
    #[error("series length mismatch: obs has {obs} steps, sim has {sim}")]
    LengthMismatch { obs: usize, sim: usize },

    #[error("need at least 2 valid (obs, sim) pairs, found {0}")]
    InsufficientData(usize),

    #[error("observed series has zero variance")]
    ZeroVariance,

    #[error("observed series sums to zero; relative bias is undefined")]
    ZeroObsSum,
}

/// Decomposed Kling–Gupta efficiency (Gupta et al. 2009).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KgeComponents<F> {
    /// Linear correlation between obs and sim.
    pub r: F,
    /// Variability ratio `σ_sim / σ_obs`.
    pub alpha: F,
    /// Bias ratio `μ_sim / μ_obs`.
    pub beta: F,
}

impl<F: Float> KgeComponents<F> {
    /// `KGE = 1 − √[(r−1)² + (α−1)² + (β−1)²]`.
    pub fn kge(&self) -> F {
        let one = F::one();
        let d2 = (self.r - one).powi(2) + (self.alpha - one).powi(2) + (self.beta - one).powi(2);
        one - d2.sqrt()
    }
}

/// Keeps only the pairs where both values are finite.
fn paired<F: Float>(obs: &[F], sim: &[F]) -> Result<Vec<(F, F)>, MetricsError> {
    if obs.len() != sim.len() {
        return Err(MetricsError::LengthMismatch {
            obs: obs.len(),
            sim: sim.len(),
        });
    }
    let pairs: Vec<(F, F)> = obs
        .iter()
        .zip(sim)
        .filter(|(o, s)| o.is_finite() && s.is_finite())
        .map(|(&o, &s)| (o, s))
        .collect();
    if pairs.len() < 2 {
        return Err(MetricsError::InsufficientData(pairs.len()));
    }
    Ok(pairs)
}

fn mean<F: Float>(values: impl Iterator<Item = F>, n: usize) -> F {
    let sum = values.fold(F::zero(), |acc, v| acc + v);
    sum / F::from(n).expect("pair count fits in F")
}

/// Nash–Sutcliffe efficiency. 1 is perfect, 0 matches the obs mean.
pub fn nse<F: Float>(obs: &[F], sim: &[F]) -> Result<F, MetricsError> {
    let pairs = paired(obs, sim)?;
    let om = mean(pairs.iter().map(|&(o, _)| o), pairs.len());
    let (num, den) = pairs
        .iter()
        .fold((F::zero(), F::zero()), |(n, d), &(o, s)| {
            (n + (s - o).powi(2), d + (o - om).powi(2))
        });
    if den == F::zero() {
        return Err(MetricsError::ZeroVariance);
    }
    Ok(F::one() - num / den)
}

/// NSE on log-transformed flows, `ln(q + ε)` with `ε = mean(obs)/100`
/// (airGR convention). Emphasizes low-flow fit.
pub fn log_nse<F: Float>(obs: &[F], sim: &[F]) -> Result<F, MetricsError> {
    let pairs = paired(obs, sim)?;
    let om = mean(pairs.iter().map(|&(o, _)| o), pairs.len());
    let eps = om / F::from(100).expect("literal fits in F");
    let lobs: Vec<F> = pairs.iter().map(|&(o, _)| (o + eps).ln()).collect();
    let lsim: Vec<F> = pairs.iter().map(|&(_, s)| (s + eps).ln()).collect();
    nse(&lobs, &lsim)
}

/// Kling–Gupta efficiency (2009 formulation). 1 is perfect.
pub fn kge<F: Float>(obs: &[F], sim: &[F]) -> Result<F, MetricsError> {
    kge_components(obs, sim).map(|c| c.kge())
}

/// KGE components (r, α, β) for diagnostic use.
pub fn kge_components<F: Float>(obs: &[F], sim: &[F]) -> Result<KgeComponents<F>, MetricsError> {
    let pairs = paired(obs, sim)?;
    let n = pairs.len();
    let om = mean(pairs.iter().map(|&(o, _)| o), n);
    let sm = mean(pairs.iter().map(|&(_, s)| s), n);
    let (mut cov, mut var_o, mut var_s) = (F::zero(), F::zero(), F::zero());
    for &(o, s) in &pairs {
        cov = cov + (o - om) * (s - sm);
        var_o = var_o + (o - om).powi(2);
        var_s = var_s + (s - sm).powi(2);
    }
    if var_o == F::zero() || om == F::zero() {
        return Err(MetricsError::ZeroVariance);
    }
    Ok(KgeComponents {
        r: cov / (var_o.sqrt() * var_s.sqrt()),
        alpha: var_s.sqrt() / var_o.sqrt(),
        beta: sm / om,
    })
}

/// Percent bias: `100·Σ(sim − obs)/Σ(obs)`. 0 is unbiased; positive means
/// the simulation overestimates.
pub fn pbias<F: Float>(obs: &[F], sim: &[F]) -> Result<F, MetricsError> {
    let pairs = paired(obs, sim)?;
    let (diff, total) = pairs
        .iter()
        .fold((F::zero(), F::zero()), |(d, t), &(o, s)| {
            (d + (s - o), t + o)
        });
    if total == F::zero() {
        return Err(MetricsError::ZeroObsSum);
    }
    Ok(F::from(100).expect("literal fits in F") * diff / total)
}

#[cfg(test)]
mod tests {
    use super::*;

    const OBS: [f64; 6] = [1.0, 2.5, 4.0, 3.0, 2.0, 1.5];

    #[test]
    fn perfect_simulation_scores_perfectly() {
        assert!((nse(&OBS, &OBS).unwrap() - 1.0).abs() < 1e-12);
        assert!((kge(&OBS, &OBS).unwrap() - 1.0).abs() < 1e-12);
        assert!((log_nse(&OBS, &OBS).unwrap() - 1.0).abs() < 1e-12);
        assert!(pbias(&OBS, &OBS).unwrap().abs() < 1e-12);
    }

    #[test]
    fn mean_simulation_gives_nse_zero() {
        let m = OBS.iter().sum::<f64>() / OBS.len() as f64;
        let sim = vec![m; OBS.len()];
        assert!(nse(&OBS, &sim).unwrap().abs() < 1e-12);
    }

    #[test]
    fn scaled_simulation_gives_known_pbias_and_beta() {
        let sim: Vec<f64> = OBS.iter().map(|&o| 1.1 * o).collect();
        assert!((pbias(&OBS, &sim).unwrap() - 10.0).abs() < 1e-10);
        let c = kge_components(&OBS, &sim).unwrap();
        assert!((c.beta - 1.1).abs() < 1e-12);
        assert!((c.r - 1.0).abs() < 1e-12);
        assert!((c.alpha - 1.1).abs() < 1e-12);
    }

    #[test]
    fn nan_observations_are_skipped() {
        let obs = [1.0, f64::NAN, 4.0, 3.0];
        let sim = [1.0, 99.0, 4.0, 3.0];
        assert!((nse(&obs, &sim).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn errors_on_degenerate_input() {
        assert!(matches!(
            nse(&[1.0, 2.0], &[1.0]),
            Err(MetricsError::LengthMismatch { .. })
        ));
        assert!(matches!(
            nse(&[2.0, 2.0, 2.0], &[1.0, 2.0, 3.0]),
            Err(MetricsError::ZeroVariance)
        ));
        let nans = [f64::NAN, f64::NAN];
        assert!(matches!(
            nse(&nans, &[1.0, 2.0]),
            Err(MetricsError::InsufficientData(0))
        ));
    }
}
