//! Channel routing and river-network accumulation.
//!
//! Turns a set of subcatchment runoff series into a single outlet hydrograph
//! by Muskingum channel routing (McCarthy 1938) along a drainage tree. This is
//! the step from "semi-distributed by elevation bands" (one outlet, vertical
//! discretization) to "semi-distributed in space" (several subcatchments
//! routed downstream).
//!
//! Everything is generic over `F: Float` (autodiff-first), like the model
//! cores. Series are unit-agnostic; route in whatever flow unit the caller
//! uses (see [`mm_per_day_to_m3s`] for the usual mm→m³/s conversion).

use num_traits::Float;

use crate::error::Error;

#[inline]
fn lit<F: Float>(x: f64) -> F {
    F::from(x).expect("f64 literal must be representable in F")
}

/// Muskingum channel routing for a single reach.
///
/// `O[t] = C0·I[t] + C1·I[t-1] + C2·O[t-1]`, with the three coefficients
/// summing to one (mass-conserving). `k` is the storage-time constant (travel
/// time, same time unit as `dt`); `x` is the weighting factor in `[0, 0.5]`
/// (0 = linear reservoir / maximum attenuation, 0.5 = pure translation).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Muskingum<F> {
    k: F,
    x: F,
}

impl<F: Float> Muskingum<F> {
    /// Validates `k > 0` and `x ∈ [0, 0.5]`.
    pub fn new(k: F, x: F) -> Result<Self, Error> {
        if k.is_nan() || k <= F::zero() {
            return Err(Error::InvalidParameter {
                name: "k",
                reason: "Muskingum K (travel time) must be > 0".into(),
            });
        }
        if x.is_nan() || x < F::zero() || x > lit(0.5) {
            return Err(Error::InvalidParameter {
                name: "x",
                reason: "Muskingum x must be in [0, 0.5]".into(),
            });
        }
        Ok(Self { k, x })
    }

    /// Routing coefficients `(C0, C1, C2)` for time step `dt`. They always sum
    /// to one; individual coefficients can go negative outside the stability
    /// window `2·K·x ≤ dt ≤ 2·K·(1−x)`, which can produce minor overshoot but
    /// never violates mass balance.
    pub fn coefficients(&self, dt: F) -> (F, F, F) {
        let half = lit::<F>(0.5);
        let kx = self.k * self.x;
        let denom = self.k - kx + half * dt;
        let c0 = (half * dt - kx) / denom;
        let c1 = (half * dt + kx) / denom;
        let c2 = (self.k - kx - half * dt) / denom;
        (c0, c1, c2)
    }

    /// Routes an inflow hydrograph through the reach with time step `dt`.
    /// The outflow is initialized to the first inflow value (steady start).
    pub fn route(&self, inflow: &[F], dt: F) -> Vec<F> {
        if inflow.is_empty() {
            return Vec::new();
        }
        let (c0, c1, c2) = self.coefficients(dt);
        let mut out = Vec::with_capacity(inflow.len());
        out.push(inflow[0]);
        for t in 1..inflow.len() {
            let o = c0 * inflow[t] + c1 * inflow[t - 1] + c2 * out[t - 1];
            // Routing is mass-conserving but can dip slightly negative under
            // overshoot; clamp to keep discharge physical.
            out.push(o.max(F::zero()));
        }
        out
    }
}

/// One subcatchment in a [`RiverNetwork`].
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Subcatchment<F> {
    /// Local runoff generated on this subcatchment (already in flow units).
    pub local_runoff: Vec<F>,
    /// Index of the downstream subcatchment this one drains into, or `None`
    /// if it is the network outlet. Must refer to a later index (the nodes
    /// are given in upstream-to-downstream topological order).
    pub downstream: Option<usize>,
    /// Muskingum reach routing this node's outflow to `downstream`. `None`
    /// passes flow on unrouted (e.g. for the outlet node, or negligible
    /// reaches).
    pub reach: Option<Muskingum<F>>,
}

/// A drainage network of subcatchments, accumulated and routed to the outlet.
#[derive(Debug, Clone)]
pub struct RiverNetwork<F> {
    nodes: Vec<Subcatchment<F>>,
    dt: F,
}

impl<F: Float> RiverNetwork<F> {
    /// Builds a network from subcatchments in upstream-to-downstream order
    /// with routing time step `dt`. Validates the topology (single outlet,
    /// strictly downstream links, equal series lengths).
    pub fn new(nodes: Vec<Subcatchment<F>>, dt: F) -> Result<Self, Error> {
        if nodes.is_empty() {
            return Err(Error::InvalidParameter {
                name: "nodes",
                reason: "network needs at least one subcatchment".into(),
            });
        }
        let len = nodes[0].local_runoff.len();
        let mut outlets = 0;
        for (i, node) in nodes.iter().enumerate() {
            if node.local_runoff.len() != len {
                return Err(Error::InvalidParameter {
                    name: "local_runoff",
                    reason: format!("node {i} series length differs from node 0"),
                });
            }
            match node.downstream {
                None => outlets += 1,
                Some(d) => {
                    if d <= i {
                        return Err(Error::InvalidParameter {
                            name: "downstream",
                            reason: format!(
                                "node {i} drains to {d}: links must point to a later (downstream) node"
                            ),
                        });
                    }
                    if d >= nodes.len() {
                        return Err(Error::InvalidParameter {
                            name: "downstream",
                            reason: format!("node {i} drains to out-of-range node {d}"),
                        });
                    }
                }
            }
        }
        if outlets != 1 {
            return Err(Error::InvalidParameter {
                name: "downstream",
                reason: format!("network must have exactly one outlet, found {outlets}"),
            });
        }
        if dt.is_nan() || dt <= F::zero() {
            return Err(Error::InvalidParameter {
                name: "dt",
                reason: "routing time step must be > 0".into(),
            });
        }
        Ok(Self { nodes, dt })
    }

    /// Accumulates and routes all subcatchments to the outlet, returning the
    /// outlet hydrograph (same length and unit as the local runoff series).
    ///
    /// Each node's total inflow is its local runoff plus the routed outflow of
    /// every upstream node draining into it; that total is then routed through
    /// the node's reach to its downstream node. Processing in topological order
    /// means every upstream contribution is ready when its target is reached.
    pub fn route_to_outlet(&self) -> Vec<F> {
        let n = self.nodes.len();
        let steps = self.nodes[0].local_runoff.len();
        // accumulated[i] = total inflow at node i (local + routed upstream).
        let mut accumulated: Vec<Vec<F>> =
            self.nodes.iter().map(|s| s.local_runoff.clone()).collect();

        let mut outlet = 0;
        for i in 0..n {
            let total = std::mem::take(&mut accumulated[i]);
            match self.nodes[i].downstream {
                Some(d) => {
                    let contribution = match &self.nodes[i].reach {
                        Some(m) => m.route(&total, self.dt),
                        None => total,
                    };
                    for (acc, c) in accumulated[d].iter_mut().zip(&contribution) {
                        *acc = *acc + *c;
                    }
                }
                None => {
                    // Outlet: its own reach (if any) is applied as a final pass.
                    outlet = i;
                    accumulated[i] = match &self.nodes[i].reach {
                        Some(m) => m.route(&total, self.dt),
                        None => total,
                    };
                }
            }
        }
        if accumulated[outlet].is_empty() {
            vec![F::zero(); steps]
        } else {
            std::mem::take(&mut accumulated[outlet])
        }
    }
}

/// Converts runoff in mm/day over `area_km2` to a volumetric flow in m³/s
/// (`q · area · 1000 / 86400`). Useful to combine subcatchments of different
/// areas before routing.
pub fn mm_per_day_to_m3s<F: Float>(q_mm_day: F, area_km2: F) -> F {
    q_mm_day * area_km2 * lit(1000.0) / lit(86400.0)
}

/// Inverse of [`mm_per_day_to_m3s`]: m³/s back to mm/day over `area_km2`.
pub fn m3s_to_mm_per_day<F: Float>(q_m3s: F, area_km2: F) -> F {
    q_m3s * lit(86400.0) / (area_km2 * lit(1000.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn triangular_hydrograph(n: usize, peak_at: usize, peak: f64) -> Vec<f64> {
        (0..n)
            .map(|t| {
                if t <= peak_at {
                    peak * t as f64 / peak_at as f64
                } else {
                    (peak * (1.0 - (t - peak_at) as f64 / (n - peak_at) as f64)).max(0.0)
                }
            })
            .collect()
    }

    #[test]
    fn coefficients_sum_to_one() {
        for (k, x) in [(2.0, 0.2), (5.0, 0.0), (1.5, 0.5), (10.0, 0.45)] {
            let (c0, c1, c2) = Muskingum::new(k, x).unwrap().coefficients(1.0);
            assert!((c0 + c1 + c2 - 1.0).abs() < 1e-12, "K={k} x={x}");
        }
    }

    #[test]
    fn k_equals_dt_x_half_is_one_step_lag() {
        // K = dt, x = 0.5 => C1 = 1, C0 = C2 = 0: output is input lagged one
        // step (pure translation by one dt — a reach always transports).
        let m = Muskingum::new(1.0, 0.5).unwrap();
        let (c0, c1, c2) = m.coefficients(1.0);
        assert!(c0.abs() < 1e-12 && (c1 - 1.0).abs() < 1e-12 && c2.abs() < 1e-12);
        let inflow = triangular_hydrograph(40, 8, 50.0);
        let out = m.route(&inflow, 1.0);
        // out[0] is the steady-start seed; out[t] == inflow[t-1] thereafter.
        for t in 1..inflow.len() {
            assert!((out[t] - inflow[t - 1]).abs() < 1e-9, "t={t}");
        }
    }

    #[test]
    fn linear_reservoir_attenuates_and_delays_peak() {
        // x = 0 is a linear reservoir: the routed peak is lower and later.
        let inflow = triangular_hydrograph(60, 10, 100.0);
        let m = Muskingum::new(3.0, 0.0).unwrap();
        let out = m.route(&inflow, 1.0);
        let in_peak = inflow.iter().cloned().fold(0.0, f64::max);
        let (out_peak_idx, out_peak) = out.iter().enumerate().fold(
            (0, 0.0),
            |(bi, bv), (i, &v)| if v > bv { (i, v) } else { (bi, bv) },
        );
        assert!(
            out_peak < in_peak,
            "peak not attenuated: {out_peak} vs {in_peak}"
        );
        assert!(out_peak_idx > 10, "peak not delayed: idx {out_peak_idx}");
    }

    #[test]
    fn routing_conserves_mass() {
        // A hydrograph returning to zero must conserve volume through routing.
        let inflow = triangular_hydrograph(200, 20, 50.0);
        for (k, x) in [(2.0, 0.2), (5.0, 0.1), (4.0, 0.0)] {
            let out = Muskingum::new(k, x).unwrap().route(&inflow, 1.0);
            let vin: f64 = inflow.iter().sum();
            let vout: f64 = out.iter().sum();
            // Small tail truncation only; volumes match closely.
            assert!(
                (vin - vout).abs() / vin < 0.02,
                "K={k} x={x}: {vin} vs {vout}"
            );
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(Muskingum::new(0.0, 0.2).is_err());
        assert!(Muskingum::new(-1.0, 0.2).is_err());
        assert!(Muskingum::new(2.0, -0.1).is_err());
        assert!(Muskingum::new(2.0, 0.6).is_err());
    }

    #[test]
    fn two_headwaters_into_outlet_sum_and_conserve() {
        // Nodes 0 and 1 drain into outlet 2. With pass-through reaches (no
        // routing) the outlet is exactly the sum of the three local series.
        let a = triangular_hydrograph(100, 10, 30.0);
        let b = triangular_hydrograph(100, 25, 20.0);
        let local_outlet = vec![1.0; 100];
        let nodes = vec![
            Subcatchment {
                local_runoff: a.clone(),
                downstream: Some(2),
                reach: None,
            },
            Subcatchment {
                local_runoff: b.clone(),
                downstream: Some(2),
                reach: None,
            },
            Subcatchment {
                local_runoff: local_outlet.clone(),
                downstream: None,
                reach: None,
            },
        ];
        let net = RiverNetwork::new(nodes, 1.0).unwrap();
        let out = net.route_to_outlet();
        let expected: Vec<f64> = (0..100).map(|t| a[t] + b[t] + local_outlet[t]).collect();
        for (o, e) in out.iter().zip(&expected) {
            assert!((o - e).abs() < 1e-9, "{o} vs {e}");
        }
    }

    #[test]
    fn routed_headwaters_conserve_volume_at_outlet() {
        // Same two-headwater network but with real Muskingum reaches: volume
        // is conserved at the outlet even though the hydrograph is reshaped.
        let a = triangular_hydrograph(300, 10, 30.0);
        let b = triangular_hydrograph(300, 25, 20.0);
        let nodes = vec![
            Subcatchment {
                local_runoff: a.clone(),
                downstream: Some(2),
                reach: Some(Muskingum::new(3.0, 0.2).unwrap()),
            },
            Subcatchment {
                local_runoff: b.clone(),
                downstream: Some(2),
                reach: Some(Muskingum::new(2.0, 0.1).unwrap()),
            },
            Subcatchment {
                local_runoff: vec![0.0; 300],
                downstream: None,
                reach: None,
            },
        ];
        let net = RiverNetwork::new(nodes, 1.0).unwrap();
        let out = net.route_to_outlet();
        let vin: f64 = a.iter().chain(&b).sum();
        let vout: f64 = out.iter().sum();
        assert!((vin - vout).abs() / vin < 0.02, "{vin} vs {vout}");
    }

    #[test]
    fn chain_routes_in_topological_order() {
        // 0 -> 1 -> 2 (outlet). Headwater pulse must reach the outlet attenuated.
        let pulse = triangular_hydrograph(120, 8, 40.0);
        let nodes = vec![
            Subcatchment {
                local_runoff: pulse.clone(),
                downstream: Some(1),
                reach: Some(Muskingum::new(2.0, 0.1).unwrap()),
            },
            Subcatchment {
                local_runoff: vec![0.0; 120],
                downstream: Some(2),
                reach: Some(Muskingum::new(2.0, 0.1).unwrap()),
            },
            Subcatchment {
                local_runoff: vec![0.0; 120],
                downstream: None,
                reach: None,
            },
        ];
        let net = RiverNetwork::new(nodes, 1.0).unwrap();
        let out = net.route_to_outlet();
        let vin: f64 = pulse.iter().sum();
        let vout: f64 = out.iter().sum();
        assert!(
            (vin - vout).abs() / vin < 0.05,
            "mass through chain: {vin} vs {vout}"
        );
        // Two reaches of attenuation: outlet peak well below the headwater peak.
        assert!(out.iter().cloned().fold(0.0, f64::max) < 40.0);
    }

    #[test]
    fn rejects_bad_topology() {
        let s = |d: Option<usize>| Subcatchment {
            local_runoff: vec![1.0; 10],
            downstream: d,
            reach: None,
        };
        // No outlet (cycle-ish: points forward but last has a target out of range).
        assert!(RiverNetwork::new(vec![s(Some(1)), s(Some(2))], 1.0).is_err());
        // Upstream link (node 1 -> 0) is not allowed.
        assert!(
            RiverNetwork::new(
                vec![
                    s(None),
                    Subcatchment {
                        local_runoff: vec![1.0; 10],
                        downstream: Some(0),
                        reach: None
                    }
                ],
                1.0
            )
            .is_err()
        );
        // Two outlets.
        assert!(RiverNetwork::new(vec![s(None), s(None)], 1.0).is_err());
    }

    #[test]
    fn unit_conversions_round_trip() {
        let q_mm = 3.5;
        let area = 860.0;
        let m3s = mm_per_day_to_m3s(q_mm, area);
        assert!((m3s_to_mm_per_day(m3s, area) - q_mm).abs() < 1e-12);
        // 1 mm/day over 86.4 km² = exactly 1 m³/s.
        assert!((mm_per_day_to_m3s(1.0, 86.4) - 1.0).abs() < 1e-12);
    }
}
