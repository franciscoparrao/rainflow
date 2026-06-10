//! HBV-light daily rainfall–runoff model.
//!
//! Reference: Seibert, J., Vis, M.J.P. (2012). *Teaching hydrological modeling
//! with a user-friendly catchment-runoff-model software package*. HESS 16,
//! 3315–3325 (standard single-zone version of Bergström's HBV).
//!
//! Structure: degree-day snow routine → soil moisture accounting → two-box
//! response routine (upper zone with threshold outflow, percolation to a
//! lower zone) → triangular MAXBAS routing.
//!
//! Two deliberate, documented deviations from HBV-light:
//! - soil recharge is computed in bulk per time step (HBV-light integrates in
//!   1 mm increments); differences are small at daily resolution.
//! - upper-zone outflows are drawn sequentially (Q0 then Q1), which keeps the
//!   store non-negative and mass balance exact for any K0+K1 <= 2.
//!
//! Temperature forcing is optional: without it the snow routine is bypassed
//! and all precipitation is treated as rain (adequate for pluvial catchments).

use num_traits::Float;

use crate::error::Error;
use crate::uh::shift_front;

/// Converts an `f64` literal into the working scalar type.
/// Infallible for IEEE-754-backed `Float` types.
#[inline]
fn lit<F: Float>(x: f64) -> F {
    F::from(x).expect("f64 literal must be representable in F")
}

/// HBV-light parameters (standard notation).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HbvParams<F> {
    // Snow routine (ignored when no temperature is supplied).
    /// `TT` — rain/snow threshold temperature (°C).
    pub tt: F,
    /// `CFMAX` — degree-day melt factor (mm/°C/day), >= 0.
    pub cfmax: F,
    /// `SFCF` — snowfall correction factor (-), > 0.
    pub sfcf: F,
    /// `CFR` — refreezing coefficient (-), in [0, 1]. HBV-light default 0.05.
    pub cfr: F,
    /// `CWH` — liquid water holding capacity of the snowpack (-), >= 0.
    /// HBV-light default 0.1.
    pub cwh: F,
    // Soil moisture routine.
    /// `FC` — maximum soil moisture storage (mm), > 0.
    pub fc: F,
    /// `LP` — fraction of FC above which AET = PET (-), in (0, 1].
    pub lp: F,
    /// `BETA` — recharge non-linearity exponent (-), > 0.
    pub beta: F,
    // Response routine.
    /// `K0` — near-surface recession coefficient (1/day), in [0, 1].
    pub k0: F,
    /// `K1` — upper-zone recession coefficient (1/day), in [0, 1].
    pub k1: F,
    /// `K2` — lower-zone recession coefficient (1/day), in [0, 1].
    pub k2: F,
    /// `UZL` — upper-zone threshold for Q0 (mm), >= 0.
    pub uzl: F,
    /// `PERC` — maximum percolation to the lower zone (mm/day), >= 0.
    pub perc: F,
    // Routing.
    /// `MAXBAS` — triangular unit-hydrograph length (days), >= 1.
    pub maxbas: F,
}

impl<F: Float> HbvParams<F> {
    /// HBV-light defaults for the parameters that are usually fixed.
    pub fn with_fixed_defaults(mut self) -> Self {
        self.cfr = lit(0.05);
        self.cwh = lit(0.1);
        self
    }
}

/// Mutable model state carried between time steps.
#[derive(Debug, Clone)]
pub struct HbvState<F> {
    /// Frozen snowpack (mm water equivalent).
    pub snowpack: F,
    /// Liquid water held in the snowpack (mm).
    pub liquid_water: F,
    /// Soil moisture (mm), in `[0, FC]`.
    pub soil_moisture: F,
    /// Upper response zone (mm).
    pub suz: F,
    /// Lower response zone (mm).
    pub slz: F,
    routing: Vec<F>,
}

/// HBV-light model with pre-computed MAXBAS routing weights.
#[derive(Debug, Clone)]
pub struct Hbv<F> {
    params: HbvParams<F>,
    weights: Vec<F>,
}

impl<F: Float> Hbv<F> {
    /// Builds the model, validating parameters and pre-computing the
    /// triangular routing weights.
    pub fn new(params: HbvParams<F>) -> Result<Self, Error> {
        let zero = F::zero();
        let one = F::one();
        let check = |ok: bool, name: &'static str, reason: &str| {
            if ok {
                Ok(())
            } else {
                Err(Error::InvalidParameter {
                    name,
                    reason: reason.into(),
                })
            }
        };
        check(params.cfmax >= zero, "cfmax", "melt factor must be >= 0")?;
        check(
            params.sfcf > zero,
            "sfcf",
            "snowfall correction must be > 0",
        )?;
        check(
            params.cfr >= zero && params.cfr <= one,
            "cfr",
            "refreezing coefficient must be in [0, 1]",
        )?;
        check(params.cwh >= zero, "cwh", "holding capacity must be >= 0")?;
        check(params.fc > zero, "fc", "soil storage must be > 0 mm")?;
        check(
            params.lp > zero && params.lp <= one,
            "lp",
            "LP must be in (0, 1]",
        )?;
        check(params.beta > zero, "beta", "BETA must be > 0")?;
        for (k, name) in [(params.k0, "k0"), (params.k1, "k1"), (params.k2, "k2")] {
            check(
                k >= zero && k <= one,
                if name == "k0" {
                    "k0"
                } else if name == "k1" {
                    "k1"
                } else {
                    "k2"
                },
                "recession coefficients must be in [0, 1]",
            )?;
        }
        check(params.uzl >= zero, "uzl", "UZL must be >= 0 mm")?;
        check(params.perc >= zero, "perc", "PERC must be >= 0 mm/day")?;
        check(
            params.maxbas >= one && !params.maxbas.is_nan(),
            "maxbas",
            "MAXBAS must be >= 1 day",
        )?;
        let weights = maxbas_weights(params.maxbas);
        Ok(Self { params, weights })
    }

    pub fn params(&self) -> &HbvParams<F> {
        &self.params
    }

    /// Default initial state: soil at half capacity, empty stores.
    /// Use at least one year of warm-up so the lower zone equilibrates.
    pub fn initial_state(&self) -> HbvState<F> {
        HbvState {
            snowpack: F::zero(),
            liquid_water: F::zero(),
            soil_moisture: lit::<F>(0.5) * self.params.fc,
            suz: F::zero(),
            slz: F::zero(),
            routing: vec![F::zero(); self.weights.len()],
        }
    }

    /// Advances one time step with precipitation `p` (mm), optional mean air
    /// temperature `temp` (°C) and potential evapotranspiration `pet` (mm),
    /// returning simulated discharge (mm).
    pub fn step(&self, state: &mut HbvState<F>, p: F, temp: Option<F>, pet: F) -> F {
        let prm = &self.params;
        let zero = F::zero();
        let one = F::one();

        // --- Snow routine: split precip, melt/refreeze, release excess ---
        let to_soil = match temp {
            Some(t) => {
                if t < prm.tt {
                    state.snowpack = state.snowpack + p * prm.sfcf;
                    let refreeze = (prm.cfr * prm.cfmax * (prm.tt - t)).min(state.liquid_water);
                    state.snowpack = state.snowpack + refreeze;
                    state.liquid_water = state.liquid_water - refreeze;
                } else {
                    let melt = (prm.cfmax * (t - prm.tt)).min(state.snowpack);
                    state.snowpack = state.snowpack - melt;
                    // Rain and meltwater pass through the snowpack's liquid
                    // storage; with no snowpack the capacity is zero and all
                    // of it is released immediately.
                    state.liquid_water = state.liquid_water + melt + p;
                }
                let capacity = prm.cwh * state.snowpack;
                let release = (state.liquid_water - capacity).max(zero);
                state.liquid_water = state.liquid_water - release;
                release
            }
            None => p,
        };

        // --- Soil moisture routine (bulk recharge, then overflow, then AET) ---
        let sm_ratio = (state.soil_moisture / prm.fc).min(one);
        let mut recharge = to_soil * sm_ratio.powf(prm.beta);
        state.soil_moisture = state.soil_moisture + to_soil - recharge;
        let overflow = (state.soil_moisture - prm.fc).max(zero);
        state.soil_moisture = state.soil_moisture - overflow;
        recharge = recharge + overflow;

        let aet =
            (pet * (state.soil_moisture / (prm.fc * prm.lp)).min(one)).min(state.soil_moisture);
        state.soil_moisture = state.soil_moisture - aet;

        // --- Response routine: upper zone (Q0, Q1) and lower zone (Q2) ---
        state.suz = state.suz + recharge;
        let percolation = prm.perc.min(state.suz);
        state.suz = state.suz - percolation;
        state.slz = state.slz + percolation;

        let q0 = prm.k0 * (state.suz - prm.uzl).max(zero);
        state.suz = state.suz - q0;
        let q1 = prm.k1 * state.suz;
        state.suz = state.suz - q1;
        let q2 = prm.k2 * state.slz;
        state.slz = state.slz - q2;

        // --- Triangular routing ---
        let qgen = q0 + q1 + q2;
        for (st, w) in state.routing.iter_mut().zip(&self.weights) {
            *st = *st + *w * qgen;
        }
        shift_front(&mut state.routing)
    }

    /// Runs the full series from the default initial state. `temp` may be
    /// `None` for pluvial catchments (snow routine bypassed).
    pub fn run(&self, precip: &[F], pet: &[F], temp: Option<&[F]>) -> Result<Vec<F>, Error> {
        if precip.len() != pet.len() {
            return Err(Error::ForcingLengthMismatch {
                precip: precip.len(),
                pet: pet.len(),
            });
        }
        if let Some(t) = temp
            && t.len() != precip.len()
        {
            return Err(Error::InvalidParameter {
                name: "temp",
                reason: format!("expected {} steps, got {}", precip.len(), t.len()),
            });
        }
        let mut state = self.initial_state();
        Ok((0..precip.len())
            .map(|i| self.step(&mut state, precip[i], temp.map(|t| t[i]), pet[i]))
            .collect())
    }

    /// Total water stored in the model (stores + water in transit), for
    /// mass-balance checks.
    pub fn storage(&self, state: &HbvState<F>) -> F {
        let routed = state.routing.iter().fold(F::zero(), |acc, &v| acc + v);
        state.snowpack + state.liquid_water + state.soil_moisture + state.suz + state.slz + routed
    }
}

/// Triangular MAXBAS weights from the cumulative S-curve
/// `G(t) = 2t²/mb²` for `t <= mb/2`, `1 − 2(mb−t)²/mb²` for `t <= mb`.
fn maxbas_weights<F: Float>(maxbas: F) -> Vec<F> {
    let n = maxbas.ceil().to_usize().unwrap_or(1).max(1);
    let g = |t: F| -> F {
        let two = lit::<F>(2.0);
        if t <= F::zero() {
            F::zero()
        } else if t + t <= maxbas {
            two * t * t / (maxbas * maxbas)
        } else if t < maxbas {
            F::one() - two * (maxbas - t) * (maxbas - t) / (maxbas * maxbas)
        } else {
            F::one()
        }
    };
    (0..n)
        .map(|j| g(lit::<F>((j + 1) as f64)) - g(lit::<F>(j as f64)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> HbvParams<f64> {
        HbvParams {
            tt: 0.0,
            cfmax: 3.5,
            sfcf: 0.9,
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
        }
    }

    /// Deterministic synthetic forcing (LCG) with temperature.
    fn synthetic_forcing(n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut seed: u64 = 42;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as f64 / (1u64 << 31) as f64
        };
        let (mut p, mut pet, mut temp) = (Vec::new(), Vec::new(), Vec::new());
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
            let season = (2.0 * std::f64::consts::PI * (doy - 15.0) / 365.25).sin();
            pet.push((3.5 + 2.5 * season).max(0.1));
            // Andes-like seasonal temperature crossing 0 °C in winter.
            temp.push(8.0 + 12.0 * season + 3.0 * (next() - 0.5));
        }
        (p, pet, temp)
    }

    #[test]
    fn maxbas_weights_sum_to_one() {
        for mb in [1.0, 1.5, 2.5, 4.0, 7.3] {
            let w = maxbas_weights::<f64>(mb);
            let s: f64 = w.iter().sum();
            assert!((s - 1.0).abs() < 1e-12, "sum {s} for maxbas={mb}");
            assert_eq!(w.len(), mb.ceil() as usize);
        }
    }

    #[test]
    fn mass_balance_closes_with_snow() {
        let model = Hbv::new(params()).unwrap();
        let (p, pet, temp) = synthetic_forcing(3000);
        let mut state = model.initial_state();
        let s0 = model.storage(&state);
        let (mut total_q, mut total_in) = (0.0, 0.0);
        for i in 0..p.len() {
            // Snowfall is scaled by SFCF, so accumulate the corrected input.
            total_in += if temp[i] < model.params().tt {
                p[i] * model.params().sfcf
            } else {
                p[i]
            };
            total_q += model.step(&mut state, p[i], Some(temp[i]), pet[i]);
        }
        let aet = total_in - total_q - (model.storage(&state) - s0);
        let total_pet: f64 = pet.iter().sum();
        assert!(aet >= 0.0, "negative actual ET: {aet}");
        assert!(aet <= total_pet, "actual ET {aet} exceeds PET {total_pet}");
    }

    #[test]
    fn snow_accumulates_below_threshold_and_melts_above() {
        let model = Hbv::new(params()).unwrap();
        let mut state = model.initial_state();
        for _ in 0..30 {
            model.step(&mut state, 10.0, Some(-5.0), 0.5);
        }
        let peak = state.snowpack;
        assert!(
            (peak - 30.0 * 10.0 * 0.9).abs() < 1e-9,
            "snowpack {peak} != accumulated corrected snowfall"
        );
        for _ in 0..200 {
            model.step(&mut state, 0.0, Some(10.0), 2.0);
        }
        assert!(
            state.snowpack < 1e-9,
            "snowpack did not melt: {}",
            state.snowpack
        );
    }

    #[test]
    fn stores_stay_bounded() {
        let model = Hbv::new(params()).unwrap();
        let (p, pet, temp) = synthetic_forcing(3000);
        let mut state = model.initial_state();
        for i in 0..p.len() {
            let q = model.step(&mut state, p[i], Some(temp[i]), pet[i]);
            assert!(q >= 0.0);
            assert!(state.soil_moisture >= 0.0 && state.soil_moisture <= 250.0);
            assert!(state.suz >= 0.0 && state.slz >= 0.0);
            assert!(state.snowpack >= 0.0 && state.liquid_water >= 0.0);
        }
    }

    #[test]
    fn runs_without_temperature() {
        let model = Hbv::new(params()).unwrap();
        let (p, pet, _) = synthetic_forcing(500);
        let q = model.run(&p, &pet, None).unwrap();
        assert_eq!(q.len(), 500);
        assert!(q.iter().all(|&v| v >= 0.0 && v.is_finite()));
        // Without temperature the snow stores must stay untouched.
        let mut state = model.initial_state();
        for i in 0..p.len() {
            model.step(&mut state, p[i], None, pet[i]);
        }
        assert_eq!(state.snowpack, 0.0);
        assert_eq!(state.liquid_water, 0.0);
    }

    #[test]
    fn generic_over_f32_and_f64() {
        let (p64, pet64, temp64) = synthetic_forcing(400);
        let q64 = Hbv::new(params())
            .unwrap()
            .run(&p64, &pet64, Some(&temp64))
            .unwrap();

        let p32: Vec<f32> = p64.iter().map(|&v| v as f32).collect();
        let pet32: Vec<f32> = pet64.iter().map(|&v| v as f32).collect();
        let temp32: Vec<f32> = temp64.iter().map(|&v| v as f32).collect();
        let prm = params();
        let m32 = Hbv::new(HbvParams {
            tt: prm.tt as f32,
            cfmax: prm.cfmax as f32,
            sfcf: prm.sfcf as f32,
            cfr: prm.cfr as f32,
            cwh: prm.cwh as f32,
            fc: prm.fc as f32,
            lp: prm.lp as f32,
            beta: prm.beta as f32,
            k0: prm.k0 as f32,
            k1: prm.k1 as f32,
            k2: prm.k2 as f32,
            uzl: prm.uzl as f32,
            perc: prm.perc as f32,
            maxbas: prm.maxbas as f32,
        })
        .unwrap();
        let q32 = m32.run(&p32, &pet32, Some(&temp32)).unwrap();
        for (a, b) in q64.iter().zip(&q32) {
            assert!(
                (a - *b as f64).abs() < 5e-3,
                "f32/f64 divergence: {a} vs {b}"
            );
        }
    }

    #[test]
    fn rejects_invalid_parameters() {
        let mut bad = params();
        bad.fc = 0.0;
        assert!(Hbv::new(bad).is_err());
        let mut bad = params();
        bad.lp = 0.0;
        assert!(Hbv::new(bad).is_err());
        let mut bad = params();
        bad.k1 = 1.5;
        assert!(Hbv::new(bad).is_err());
        let mut bad = params();
        bad.maxbas = 0.5;
        assert!(Hbv::new(bad).is_err());
    }
}
