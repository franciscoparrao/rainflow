//! GR4J daily rainfall–runoff model.
//!
//! Reference: Perrin, C., Michel, C., Andréassian, V. (2003). *Improvement of a
//! parsimonious model for streamflow simulation*. Journal of Hydrology 279, 275–289.
//!
//! The implementation mirrors the airGR (INRAE) formulation so simulated
//! discharge can be cross-checked numerically against `airGR::RunModel_GR4J`.
//!
//! All arithmetic is generic over `F: Float` (autodiff-first): the only place a
//! concrete `f64` appears is in the conversion of literal constants, and in the
//! discrete sizing of the unit-hydrograph buffers (which is not differentiable
//! by nature; the ordinate *values* remain differentiable in `x4`).

use num_traits::Float;

use crate::error::Error;
use crate::uh::shift_front;

/// Converts an `f64` literal into the working scalar type.
///
/// Infallible for any IEEE-754-backed `Float` (f32, f64, dual numbers over
/// them); the `expect` is unreachable in practice.
#[inline]
fn lit<F: Float>(x: f64) -> F {
    F::from(x).expect("f64 literal must be representable in F")
}

/// GR4J parameters (Perrin et al. 2003 notation).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gr4jParams<F> {
    /// `x1` — production store capacity (mm), > 0.
    pub x1: F,
    /// `x2` — groundwater exchange coefficient (mm/day).
    pub x2: F,
    /// `x3` — routing store capacity (mm), > 0.
    pub x3: F,
    /// `x4` — unit hydrograph time base (days), >= 0.5.
    pub x4: F,
}

/// Mutable model state carried between time steps.
#[derive(Debug, Clone)]
pub struct Gr4jState<F> {
    /// Production store level (mm), in `[0, x1]`.
    pub s: F,
    /// Routing store level (mm), in `[0, x3]`.
    pub r: F,
    uh1: Vec<F>,
    uh2: Vec<F>,
}

/// GR4J model with pre-computed unit-hydrograph ordinates.
#[derive(Debug, Clone)]
pub struct Gr4j<F> {
    params: Gr4jParams<F>,
    ord1: Vec<F>,
    ord2: Vec<F>,
}

impl<F: Float> Gr4j<F> {
    /// Builds the model, validating parameters and pre-computing UH ordinates.
    pub fn new(params: Gr4jParams<F>) -> Result<Self, Error> {
        // NaN parameters must be rejected too, hence the explicit is_nan checks.
        if params.x1.is_nan() || params.x1 <= F::zero() {
            return Err(Error::InvalidParameter {
                name: "x1",
                reason: "production store capacity must be > 0 mm".into(),
            });
        }
        if params.x2.is_nan() {
            return Err(Error::InvalidParameter {
                name: "x2",
                reason: "groundwater exchange coefficient must be finite".into(),
            });
        }
        if params.x3.is_nan() || params.x3 <= F::zero() {
            return Err(Error::InvalidParameter {
                name: "x3",
                reason: "routing store capacity must be > 0 mm".into(),
            });
        }
        if params.x4.is_nan() || params.x4 < lit(0.5) {
            return Err(Error::InvalidParameter {
                name: "x4",
                reason: "unit hydrograph time base must be >= 0.5 days".into(),
            });
        }
        let (ord1, ord2) = ordinates(params.x4);
        Ok(Self { params, ord1, ord2 })
    }

    pub fn params(&self) -> &Gr4jParams<F> {
        &self.params
    }

    /// Default initial state: `S = 0.3·x1`, `R = 0.5·x3` (airGR convention).
    /// Use a warm-up period before evaluating goodness of fit.
    pub fn initial_state(&self) -> Gr4jState<F> {
        Gr4jState {
            s: lit::<F>(0.3) * self.params.x1,
            r: lit::<F>(0.5) * self.params.x3,
            uh1: vec![F::zero(); self.ord1.len()],
            uh2: vec![F::zero(); self.ord2.len()],
        }
    }

    /// Advances one time step with precipitation `p` and potential
    /// evapotranspiration `pet` (both mm), returning simulated discharge (mm).
    pub fn step(&self, state: &mut Gr4jState<F>, p: F, pet: F) -> F {
        let prm = &self.params;
        let zero = F::zero();
        let one = F::one();

        // Interception by net precipitation / net evapotranspiration.
        let (pn, en) = if p >= pet {
            (p - pet, zero)
        } else {
            (zero, pet - p)
        };

        // Production store: rainfall fills it, evapotranspiration drains it.
        let sr = state.s / prm.x1;
        let ps = if pn > zero {
            let w = (pn / prm.x1).tanh();
            prm.x1 * (one - sr * sr) * w / (one + sr * w)
        } else {
            zero
        };
        let es = if en > zero {
            let w = (en / prm.x1).tanh();
            state.s * (lit::<F>(2.0) - sr) * w / (one + (one - sr) * w)
        } else {
            zero
        };
        state.s = state.s + ps - es;

        // Percolation leak: Perc = S·(1 − [1 + (4S/9x1)^4]^(−1/4)).
        let c = state.s / (lit::<F>(2.25) * prm.x1);
        let perc = state.s * (one - (one + c.powi(4)).powf(lit(-0.25)));
        state.s = state.s - perc;

        // Effective rainfall, split 90% slow / 10% fast routing.
        let pr = perc + (pn - ps);
        let pr9 = lit::<F>(0.9) * pr;
        let pr1 = pr - pr9;

        for (st, o) in state.uh1.iter_mut().zip(&self.ord1) {
            *st = *st + *o * pr9;
        }
        for (st, o) in state.uh2.iter_mut().zip(&self.ord2) {
            *st = *st + *o * pr1;
        }
        let q9 = shift_front(&mut state.uh1);
        let q1 = shift_front(&mut state.uh2);

        // Groundwater exchange: F = x2·(R/x3)^(7/2).
        let f = prm.x2 * (state.r / prm.x3).powf(lit(3.5));

        // Non-linear routing store.
        state.r = (state.r + q9 + f).max(zero);
        let rr = state.r / prm.x3;
        let qr = state.r * (one - (one + rr.powi(4)).powf(lit(-0.25)));
        state.r = state.r - qr;

        // Direct branch.
        let qd = (q1 + f).max(zero);

        qr + qd
    }

    /// Runs the full series from the default initial state.
    pub fn run(&self, precip: &[F], pet: &[F]) -> Result<Vec<F>, Error> {
        if precip.len() != pet.len() {
            return Err(Error::ForcingLengthMismatch {
                precip: precip.len(),
                pet: pet.len(),
            });
        }
        let mut state = self.initial_state();
        Ok(precip
            .iter()
            .zip(pet)
            .map(|(&p, &e)| self.step(&mut state, p, e))
            .collect())
    }
}

/// S-curve of UH1: `SH1(t) = (t/x4)^(5/2)` for `0 <= t <= x4`, 1 beyond.
fn s_curve1<F: Float>(t: F, x4: F) -> F {
    if t <= F::zero() {
        F::zero()
    } else if t < x4 {
        (t / x4).powf(lit(2.5))
    } else {
        F::one()
    }
}

/// S-curve of UH2 over `[0, 2·x4]`.
fn s_curve2<F: Float>(t: F, x4: F) -> F {
    let two = lit::<F>(2.0);
    let half = lit::<F>(0.5);
    if t <= F::zero() {
        F::zero()
    } else if t <= x4 {
        half * (t / x4).powf(lit(2.5))
    } else if t < two * x4 {
        F::one() - half * (two - t / x4).powf(lit(2.5))
    } else {
        F::one()
    }
}

/// Unit-hydrograph ordinates as successive differences of the S-curves.
fn ordinates<F: Float>(x4: F) -> (Vec<F>, Vec<F>) {
    let n1 = x4.ceil().to_usize().unwrap_or(1).max(1);
    let n2 = (lit::<F>(2.0) * x4).ceil().to_usize().unwrap_or(1).max(1);
    let ord1 = (0..n1)
        .map(|j| s_curve1(lit::<F>((j + 1) as f64), x4) - s_curve1(lit::<F>(j as f64), x4))
        .collect();
    let ord2 = (0..n2)
        .map(|j| s_curve2(lit::<F>((j + 1) as f64), x4) - s_curve2(lit::<F>(j as f64), x4))
        .collect();
    (ord1, ord2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model(x1: f64, x2: f64, x3: f64, x4: f64) -> Gr4j<f64> {
        Gr4j::new(Gr4jParams { x1, x2, x3, x4 }).unwrap()
    }

    /// Deterministic synthetic forcing (LCG), winter-rain regime.
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

    #[test]
    fn uh_ordinates_sum_to_one() {
        for x4 in [0.5, 1.0, 1.7, 3.9, 10.3] {
            let (o1, o2) = ordinates::<f64>(x4);
            let s1: f64 = o1.iter().sum();
            let s2: f64 = o2.iter().sum();
            assert!((s1 - 1.0).abs() < 1e-12, "UH1 sum {s1} for x4={x4}");
            assert!((s2 - 1.0).abs() < 1e-12, "UH2 sum {s2} for x4={x4}");
            assert_eq!(o1.len(), x4.ceil() as usize);
            assert_eq!(o2.len(), (2.0 * x4).ceil() as usize);
        }
    }

    #[test]
    fn stores_stay_within_capacity() {
        let m = model(350.0, -1.5, 90.0, 1.7);
        let (p, pet) = synthetic_forcing(2000);
        let mut state = m.initial_state();
        for (&pi, &ei) in p.iter().zip(&pet) {
            let q = m.step(&mut state, pi, ei);
            assert!(q >= 0.0, "negative discharge");
            assert!(
                state.s >= 0.0 && state.s <= 350.0,
                "S out of bounds: {}",
                state.s
            );
            assert!(
                state.r >= 0.0 && state.r <= 90.0 + 1e-9,
                "R out of bounds: {}",
                state.r
            );
        }
    }

    #[test]
    fn recession_is_monotonic_without_rain() {
        let m = model(350.0, 0.0, 90.0, 1.7);
        let mut state = m.initial_state();
        // Wet spell to charge the stores, then dry recession.
        for _ in 0..30 {
            m.step(&mut state, 20.0, 1.0);
        }
        let mut prev = f64::INFINITY;
        for i in 0..100 {
            let q = m.step(&mut state, 0.0, 2.0);
            // Skip the UH tail (first 2·x4 steps still release stored pulses).
            if i > 4 {
                assert!(q <= prev + 1e-12, "recession not monotonic at step {i}");
            }
            prev = q;
        }
    }

    #[test]
    fn mass_balance_closes_without_exchange() {
        // With x2 = 0 the only fluxes are P, actual ET, Q and storage change
        // (stores + water in transit inside the UH buffers).
        let m = model(350.0, 0.0, 90.0, 1.7);
        let (p, pet) = synthetic_forcing(3000);
        let q = m.run(&p, &pet).unwrap();
        let state0 = m.initial_state();
        // Recompute final state to read storages.
        let mut state = m.initial_state();
        for (&pi, &ei) in p.iter().zip(&pet) {
            m.step(&mut state, pi, ei);
        }
        let total_p: f64 = p.iter().sum();
        let total_q: f64 = q.iter().sum();
        let in_transit: f64 = state.uh1.iter().sum::<f64>() + state.uh2.iter().sum::<f64>();
        let storage_change = (state.s - state0.s) + (state.r - state0.r) + in_transit;
        // Actual ET = P − Q − ΔS must be non-negative and bounded by PET total.
        let aet = total_p - total_q - storage_change;
        let total_pet: f64 = pet.iter().sum();
        assert!(aet >= 0.0, "negative actual ET: {aet}");
        assert!(aet <= total_pet, "actual ET {aet} exceeds PET {total_pet}");
    }

    #[test]
    fn generic_over_f32_and_f64() {
        let (p64, pet64) = synthetic_forcing(400);
        let m64 = model(350.0, -1.5, 90.0, 1.7);
        let q64 = m64.run(&p64, &pet64).unwrap();

        let m32 = Gr4j::new(Gr4jParams {
            x1: 350.0f32,
            x2: -1.5,
            x3: 90.0,
            x4: 1.7,
        })
        .unwrap();
        let p32: Vec<f32> = p64.iter().map(|&v| v as f32).collect();
        let pet32: Vec<f32> = pet64.iter().map(|&v| v as f32).collect();
        let q32 = m32.run(&p32, &pet32).unwrap();

        for (a, b) in q64.iter().zip(&q32) {
            assert!(
                (a - *b as f64).abs() < 1e-3,
                "f32/f64 divergence: {a} vs {b}"
            );
        }
    }

    #[test]
    fn rejects_invalid_parameters() {
        assert!(
            Gr4j::new(Gr4jParams {
                x1: 0.0,
                x2: 0.0,
                x3: 90.0,
                x4: 1.7
            })
            .is_err()
        );
        assert!(
            Gr4j::new(Gr4jParams {
                x1: 350.0,
                x2: 0.0,
                x3: -1.0,
                x4: 1.7
            })
            .is_err()
        );
        assert!(
            Gr4j::new(Gr4jParams {
                x1: 350.0,
                x2: 0.0,
                x3: 90.0,
                x4: 0.4
            })
            .is_err()
        );
    }
}
