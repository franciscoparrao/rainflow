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
//!
//! [`HbvBands`] is the semi-distributed variant: the snow and soil routines
//! run per elevation band (lapsed temperature, precipitation gradient) and the
//! area-weighted recharge feeds a single shared response + routing, as in
//! HBV-light's elevation zones.

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

        let to_soil = match temp {
            Some(t) => snow_step(prm, &mut state.snowpack, &mut state.liquid_water, p, t),
            None => p,
        };
        let recharge = soil_step(prm, &mut state.soil_moisture, to_soil, pet);
        let qgen = response_step(prm, &mut state.suz, &mut state.slz, recharge);
        route(&self.weights, &mut state.routing, qgen)
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

/// Snow routine: split precipitation by TT, melt/refreeze, release the liquid
/// water exceeding the pack's holding capacity. Returns water reaching soil.
fn snow_step<F: Float>(prm: &HbvParams<F>, snowpack: &mut F, liquid: &mut F, p: F, t: F) -> F {
    let zero = F::zero();
    if t < prm.tt {
        *snowpack = *snowpack + p * prm.sfcf;
        let refreeze = (prm.cfr * prm.cfmax * (prm.tt - t)).min(*liquid);
        *snowpack = *snowpack + refreeze;
        *liquid = *liquid - refreeze;
    } else {
        let melt = (prm.cfmax * (t - prm.tt)).min(*snowpack);
        *snowpack = *snowpack - melt;
        // Rain and meltwater pass through the snowpack's liquid storage; with
        // no snowpack the capacity is zero and all of it is released at once.
        *liquid = *liquid + melt + p;
    }
    let capacity = prm.cwh * *snowpack;
    let release = (*liquid - capacity).max(zero);
    *liquid = *liquid - release;
    release
}

/// Soil moisture routine: bulk recharge, capacity overflow, then AET.
/// Returns recharge to the response routine.
fn soil_step<F: Float>(prm: &HbvParams<F>, sm: &mut F, inflow: F, pet: F) -> F {
    let zero = F::zero();
    let one = F::one();
    let sm_ratio = (*sm / prm.fc).min(one);
    let mut recharge = inflow * sm_ratio.powf(prm.beta);
    *sm = *sm + inflow - recharge;
    let overflow = (*sm - prm.fc).max(zero);
    *sm = *sm - overflow;
    recharge = recharge + overflow;

    let aet = (pet * (*sm / (prm.fc * prm.lp)).min(one)).min(*sm);
    *sm = *sm - aet;
    recharge
}

/// Response routine: upper zone (Q0 above UZL, then Q1), percolation, lower
/// zone (Q2). Returns generated runoff before routing.
fn response_step<F: Float>(prm: &HbvParams<F>, suz: &mut F, slz: &mut F, recharge: F) -> F {
    let zero = F::zero();
    *suz = *suz + recharge;
    let percolation = prm.perc.min(*suz);
    *suz = *suz - percolation;
    *slz = *slz + percolation;

    let q0 = prm.k0 * (*suz - prm.uzl).max(zero);
    *suz = *suz - q0;
    let q1 = prm.k1 * *suz;
    *suz = *suz - q1;
    let q2 = prm.k2 * *slz;
    *slz = *slz - q2;
    q0 + q1 + q2
}

/// Pushes generated runoff through the triangular routing buffer.
fn route<F: Float>(weights: &[F], buffer: &mut [F], qgen: F) -> F {
    for (st, w) in buffer.iter_mut().zip(weights) {
        *st = *st + *w * qgen;
    }
    shift_front(buffer)
}

/// One elevation band: mean elevation (m a.s.l.) and its fraction of the
/// catchment area.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElevationBand<F> {
    pub elevation: F,
    pub area_fraction: F,
}

/// Semi-distributed configuration: elevation bands plus the lapse/gradient
/// rates used to extrapolate the lumped forcing to each band.
#[derive(Debug, Clone, PartialEq)]
pub struct ElevationBands<F> {
    pub bands: Vec<ElevationBand<F>>,
    /// Elevation the forcing series represents (e.g. catchment mean for a
    /// gridded product), m a.s.l.
    pub reference_elevation: F,
    /// Temperature lapse rate (°C per 100 m, positive = cooler with height).
    /// HBV-light default (TCALT): 0.6.
    pub tcalt: F,
    /// Precipitation gradient (fraction per 100 m). HBV-light default
    /// (PCALT): 0.10.
    pub pcalt: F,
}

impl<F: Float> ElevationBands<F> {
    /// HBV-light default lapse rates.
    pub fn new(bands: Vec<ElevationBand<F>>, reference_elevation: F) -> Self {
        Self {
            bands,
            reference_elevation,
            tcalt: lit(0.6),
            pcalt: lit(0.10),
        }
    }

    /// Builds `n` equal-area bands from a hypsometric curve.
    ///
    /// `curve` is the catchment's hypsometry as `(cumulative_area_fraction,
    /// elevation)` knots, sorted by increasing fraction and spanning at least
    /// `[0, 1]` (e.g. `[(0, z_min), (0.5, z_median), (1, z_max)]`, or a dense
    /// curve sampled from a DEM). Each band covers an equal area `1/n`; its
    /// representative elevation is the curve evaluated at the band's midpoint
    /// fraction by monotone linear interpolation. Unlike hand-picked or
    /// equal-elevation bands, the elevations follow where the catchment area
    /// actually sits.
    pub fn equal_area_from_hypsometry(
        curve: &[(F, F)],
        n: usize,
        reference_elevation: F,
    ) -> Result<Self, Error> {
        if n == 0 {
            return Err(Error::InvalidParameter {
                name: "n",
                reason: "need at least one band".into(),
            });
        }
        if curve.len() < 2 {
            return Err(Error::InvalidParameter {
                name: "curve",
                reason: "hypsometric curve needs at least two knots".into(),
            });
        }
        for w in curve.windows(2) {
            // Fractions must be non-decreasing; elevations must not decrease
            // (a hypsometric curve is monotone by construction).
            if w[1].0 < w[0].0 || w[1].1 < w[0].1 || !w[0].0.is_finite() || !w[0].1.is_finite() {
                return Err(Error::InvalidParameter {
                    name: "curve",
                    reason: "knots must be sorted with non-decreasing fraction and elevation".into(),
                });
            }
        }

        let nf = lit::<F>(n as f64);
        let bands = (0..n)
            .map(|i| {
                // Midpoint of band i's area fraction: (i + 0.5)/n.
                let frac = (lit::<F>(i as f64) + lit(0.5)) / nf;
                ElevationBand {
                    elevation: interpolate_curve(curve, frac),
                    area_fraction: F::one() / nf,
                }
            })
            .collect();
        Ok(Self::new(bands, reference_elevation))
    }

    /// Convenience hypsometry from the three elevation quantiles reported by
    /// datasets such as CAMELS-CL (minimum, median, maximum). Reconstructs the
    /// curve `[(0, min), (0.5, median), (1, max)]` and delegates to
    /// [`Self::equal_area_from_hypsometry`]. The reference elevation defaults
    /// to the median.
    ///
    /// This is the data-poor fallback; a curve sampled from a DEM clipped to
    /// the catchment is strictly better and plugs into the same constructor.
    pub fn from_quantiles(min: F, median: F, max: F, n: usize) -> Result<Self, Error> {
        let curve = [
            (F::zero(), min),
            (lit(0.5), median),
            (F::one(), max),
        ];
        Self::equal_area_from_hypsometry(&curve, n, median)
    }

    fn validate(&self) -> Result<(), Error> {
        if self.bands.is_empty() {
            return Err(Error::InvalidParameter {
                name: "bands",
                reason: "at least one elevation band is required".into(),
            });
        }
        let total = self
            .bands
            .iter()
            .fold(F::zero(), |acc, b| acc + b.area_fraction);
        let tolerance = lit(1e-3);
        if (total - F::one()).abs() > tolerance {
            return Err(Error::InvalidParameter {
                name: "bands",
                reason: format!("area fractions must sum to 1 (got {:?})", total.to_f64()),
            });
        }
        for (i, b) in self.bands.iter().enumerate() {
            if b.area_fraction < F::zero() || !b.elevation.is_finite() {
                return Err(Error::InvalidParameter {
                    name: "bands",
                    reason: format!("band {i}: negative fraction or non-finite elevation"),
                });
            }
        }
        if !self.tcalt.is_finite() || !self.pcalt.is_finite() {
            return Err(Error::InvalidParameter {
                name: "bands",
                reason: "tcalt/pcalt must be finite".into(),
            });
        }
        Ok(())
    }
}

/// Per-band snow + soil state of the semi-distributed model.
#[derive(Debug, Clone)]
pub struct BandState<F> {
    pub snowpack: F,
    pub liquid_water: F,
    pub soil_moisture: F,
}

/// State of [`HbvBands`].
#[derive(Debug, Clone)]
pub struct HbvBandsState<F> {
    pub bands: Vec<BandState<F>>,
    pub suz: F,
    pub slz: F,
    routing: Vec<F>,
}

/// Semi-distributed HBV: the snow and soil routines run independently per
/// elevation band (with lapsed temperature and a precipitation gradient);
/// the area-weighted recharge feeds a single shared response routine, as in
/// HBV-light's elevation zones.
///
/// Temperature forcing is required — elevation bands without a snow routine
/// would only re-scale precipitation.
#[derive(Debug, Clone)]
pub struct HbvBands<F> {
    params: HbvParams<F>,
    config: ElevationBands<F>,
    weights: Vec<F>,
}

impl<F: Float> HbvBands<F> {
    pub fn new(params: HbvParams<F>, config: ElevationBands<F>) -> Result<Self, Error> {
        // Reuse the lumped constructor for parameter validation.
        let lumped = Hbv::new(params)?;
        config.validate()?;
        Ok(Self {
            params,
            config,
            weights: lumped.weights,
        })
    }

    pub fn params(&self) -> &HbvParams<F> {
        &self.params
    }

    pub fn config(&self) -> &ElevationBands<F> {
        &self.config
    }

    /// Default initial state: each band's soil at half capacity, empty stores.
    pub fn initial_state(&self) -> HbvBandsState<F> {
        HbvBandsState {
            bands: self
                .config
                .bands
                .iter()
                .map(|_| BandState {
                    snowpack: F::zero(),
                    liquid_water: F::zero(),
                    soil_moisture: lit::<F>(0.5) * self.params.fc,
                })
                .collect(),
            suz: F::zero(),
            slz: F::zero(),
            routing: vec![F::zero(); self.weights.len()],
        }
    }

    /// Advances one time step. `p`, `t` and `pet` are the lumped forcing at
    /// the reference elevation; each band sees lapsed temperature and a
    /// precipitation gradient (PET is kept uniform).
    pub fn step(&self, state: &mut HbvBandsState<F>, p: F, t: F, pet: F) -> F {
        let prm = &self.params;
        let cfg = &self.config;
        let per100 = lit::<F>(0.01);

        let mut recharge = F::zero();
        for (band, bs) in cfg.bands.iter().zip(state.bands.iter_mut()) {
            let dz100 = (band.elevation - cfg.reference_elevation) * per100;
            let tb = t - cfg.tcalt * dz100;
            let pb = (p * (F::one() + cfg.pcalt * dz100)).max(F::zero());
            let to_soil = snow_step(prm, &mut bs.snowpack, &mut bs.liquid_water, pb, tb);
            recharge =
                recharge + band.area_fraction * soil_step(prm, &mut bs.soil_moisture, to_soil, pet);
        }

        let qgen = response_step(prm, &mut state.suz, &mut state.slz, recharge);
        route(&self.weights, &mut state.routing, qgen)
    }

    /// Runs the full series from the default initial state.
    pub fn run(&self, precip: &[F], pet: &[F], temp: &[F]) -> Result<Vec<F>, Error> {
        if precip.len() != pet.len() {
            return Err(Error::ForcingLengthMismatch {
                precip: precip.len(),
                pet: pet.len(),
            });
        }
        if temp.len() != precip.len() {
            return Err(Error::InvalidParameter {
                name: "temp",
                reason: format!("expected {} steps, got {}", precip.len(), temp.len()),
            });
        }
        let mut state = self.initial_state();
        Ok((0..precip.len())
            .map(|i| self.step(&mut state, precip[i], temp[i], pet[i]))
            .collect())
    }
}

/// Linear interpolation of a sorted `(x, y)` curve at `x = at`, clamped to the
/// curve's endpoints outside its range.
fn interpolate_curve<F: Float>(curve: &[(F, F)], at: F) -> F {
    if at <= curve[0].0 {
        return curve[0].1;
    }
    let last = curve.len() - 1;
    if at >= curve[last].0 {
        return curve[last].1;
    }
    for w in curve.windows(2) {
        let (x0, y0) = w[0];
        let (x1, y1) = w[1];
        if at <= x1 {
            if x1 == x0 {
                return y0;
            }
            let t = (at - x0) / (x1 - x0);
            return y0 + t * (y1 - y0);
        }
    }
    curve[last].1
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
    fn single_band_at_reference_elevation_matches_lumped_model() {
        let (p, pet, temp) = synthetic_forcing(800);
        let lumped = Hbv::new(params()).unwrap();
        let q_lumped = lumped.run(&p, &pet, Some(&temp)).unwrap();

        let banded = HbvBands::new(
            params(),
            ElevationBands::new(
                vec![ElevationBand {
                    elevation: 1000.0,
                    area_fraction: 1.0,
                }],
                1000.0,
            ),
        )
        .unwrap();
        let q_banded = banded.run(&p, &pet, &temp).unwrap();
        for (a, b) in q_lumped.iter().zip(&q_banded) {
            assert!((a - b).abs() < 1e-12, "banded != lumped: {a} vs {b}");
        }
    }

    #[test]
    fn higher_bands_accumulate_more_snow() {
        let bands = ElevationBands::new(
            vec![
                ElevationBand {
                    elevation: 1500.0,
                    area_fraction: 0.5,
                },
                ElevationBand {
                    elevation: 3500.0,
                    area_fraction: 0.5,
                },
            ],
            2500.0,
        );
        let model = HbvBands::new(params(), bands).unwrap();
        let (p, pet, temp) = synthetic_forcing(2000);
        let mut state = model.initial_state();
        let mut max_low = 0.0_f64;
        let mut max_high = 0.0_f64;
        for i in 0..p.len() {
            model.step(&mut state, p[i], temp[i], pet[i]);
            max_low = max_low.max(state.bands[0].snowpack);
            max_high = max_high.max(state.bands[1].snowpack);
        }
        assert!(
            max_high > max_low,
            "high band ({max_high} mm) should out-accumulate low band ({max_low} mm)"
        );
        assert!(max_high > 0.0, "no snow accumulated at 3500 m");
    }

    #[test]
    fn equal_area_bands_from_quantiles_are_well_formed() {
        // CAMELS-CL 4703002: min 1153, median 3322, max 5038.
        let eb = ElevationBands::from_quantiles(1153.0, 3322.0, 5038.0, 3).unwrap();
        assert_eq!(eb.bands.len(), 3);
        // Equal area fractions summing to 1.
        let total: f64 = eb.bands.iter().map(|b| b.area_fraction).sum();
        assert!((total - 1.0).abs() < 1e-12);
        for b in &eb.bands {
            assert!((b.area_fraction - 1.0 / 3.0).abs() < 1e-12);
        }
        // Elevations strictly increasing and bracketed by [min, max].
        assert!(eb.bands[0].elevation < eb.bands[1].elevation);
        assert!(eb.bands[1].elevation < eb.bands[2].elevation);
        assert!(eb.bands[0].elevation > 1153.0 && eb.bands[2].elevation < 5038.0);
        // Reference elevation defaults to the median; middle band (frac 0.5)
        // sits exactly on it.
        assert!((eb.reference_elevation - 3322.0).abs() < 1e-12);
        assert!((eb.bands[1].elevation - 3322.0).abs() < 1e-9);
        // Resulting config must validate and drive a model.
        assert!(HbvBands::new(params(), eb).is_ok());
    }

    #[test]
    fn hypsometry_single_band_is_the_median() {
        let eb = ElevationBands::from_quantiles(1000.0, 2000.0, 4000.0, 1).unwrap();
        assert_eq!(eb.bands.len(), 1);
        assert!((eb.bands[0].area_fraction - 1.0).abs() < 1e-12);
        // Midpoint fraction 0.5 -> the median knot.
        assert!((eb.bands[0].elevation - 2000.0).abs() < 1e-9);
    }

    #[test]
    fn hypsometry_rejects_bad_curves() {
        assert!(ElevationBands::<f64>::from_quantiles(1000.0, 2000.0, 4000.0, 0).is_err());
        // Non-monotone elevation.
        let bad = [(0.0, 1000.0), (0.5, 800.0), (1.0, 2000.0)];
        assert!(ElevationBands::equal_area_from_hypsometry(&bad, 3, 1500.0).is_err());
        // Single knot.
        let one = [(0.0, 1000.0)];
        assert!(ElevationBands::equal_area_from_hypsometry(&one, 2, 1000.0).is_err());
    }

    #[test]
    fn dense_hypsometry_recovers_band_means() {
        // A linear hypsometry z(f) = 1000 + 2000 f over [0, 1]: equal-area
        // bands must sit at the midpoint elevations 1000 + 2000*(i+0.5)/n.
        let curve: Vec<(f64, f64)> = (0..=100)
            .map(|k| (k as f64 / 100.0, 1000.0 + 2000.0 * (k as f64 / 100.0)))
            .collect();
        let eb = ElevationBands::equal_area_from_hypsometry(&curve, 4, 2000.0).unwrap();
        for (i, b) in eb.bands.iter().enumerate() {
            let expected = 1000.0 + 2000.0 * (i as f64 + 0.5) / 4.0;
            assert!((b.elevation - expected).abs() < 1.0, "band {i}: {} vs {expected}", b.elevation);
        }
    }

    #[test]
    fn bands_validation_rejects_bad_configs() {
        let band = |e: f64, f: f64| ElevationBand {
            elevation: e,
            area_fraction: f,
        };
        // Fractions must sum to 1.
        assert!(
            HbvBands::new(
                params(),
                ElevationBands::new(vec![band(1000.0, 0.4), band(2000.0, 0.4)], 1500.0)
            )
            .is_err()
        );
        // At least one band.
        assert!(HbvBands::new(params(), ElevationBands::new(vec![], 1500.0)).is_err());
        // Negative fraction.
        assert!(
            HbvBands::new(
                params(),
                ElevationBands::new(vec![band(1000.0, 1.5), band(2000.0, -0.5)], 1500.0)
            )
            .is_err()
        );
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
