//! Spatially semi-distributed routing on the Río Itata (CAMELS-CL).
//!
//! Two subcatchments — Itata en Cholguán (headwater, 860 km²) and the
//! incremental area down to Itata en Balsa Nueva Aldea (3650 km²) — are each
//! modelled with GR4J, converted to m³/s, and the headwater is Muskingum-routed
//! to the outlet and summed with the incremental runoff. The result is compared
//! against the observed outlet discharge (Balsa Nueva Aldea, 4510 km²) and
//! against a single lumped GR4J fitted directly on the outlet forcing.
//!
//! Run with: cargo run -p rainflow-cli --example route_itata --release

use rainflow_core::routing::{
    Muskingum, RiverNetwork, Subcatchment, m3s_to_mm_per_day, mm_per_day_to_m3s,
};
use rainflow_core::{DdsConfig, Gr4j, Gr4jParams, Optimizer, calibrate, metrics};

const A_CHOL: f64 = 859.57791;
const A_INC: f64 = 3650.48161;
const A_BALSA: f64 = 4510.05952;
const WARMUP: usize = 365;

struct Forcing {
    p: Vec<f64>,
    pet: Vec<f64>,
    qobs: Vec<f64>,
}

fn read(path: &str) -> Forcing {
    let mut rdr = csv::Reader::from_path(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let headers: Vec<String> = rdr.headers().unwrap().iter().map(str::to_owned).collect();
    let qi = headers.iter().position(|h| h == "qobs");
    let (mut p, mut pet, mut qobs) = (Vec::new(), Vec::new(), Vec::new());
    for rec in rdr.records() {
        let rec = rec.unwrap();
        let num = |s: &str| -> f64 {
            let s = s.trim();
            if s.is_empty() || s.eq_ignore_ascii_case("na") {
                f64::NAN
            } else {
                s.parse().unwrap()
            }
        };
        p.push(num(&rec[1]));
        pet.push(num(&rec[2]));
        qobs.push(qi.map(|i| num(&rec[i])).unwrap_or(f64::NAN));
    }
    Forcing { p, pet, qobs }
}

/// GR4J discharge (mm/day) for a parameter slice [x1, x2, x3, x4].
fn gr4j_run(x: &[f64], p: &[f64], pet: &[f64]) -> Option<Vec<f64>> {
    let m = Gr4j::new(Gr4jParams {
        x1: x[0],
        x2: x[1],
        x3: x[2],
        x4: x[3],
    })
    .ok()?;
    m.run(p, pet).ok()
}

/// Distributed outlet discharge (mm/day over the Balsa area) for a full
/// parameter vector: GR4J_chol[0..4], GR4J_inc[4..8], Muskingum K [8] (x fixed).
fn distributed_outlet(x: &[f64], chol: &Forcing, inc: &Forcing) -> Option<Vec<f64>> {
    let q_chol = gr4j_run(&x[0..4], &chol.p, &chol.pet)?;
    let q_inc = gr4j_run(&x[4..8], &inc.p, &inc.pet)?;
    let reach = Muskingum::new(x[8], 0.2).ok()?;
    // Headwater (node 0) routes into the outlet node (node 1, incremental area).
    let chol_m3s: Vec<f64> = q_chol
        .iter()
        .map(|&q| mm_per_day_to_m3s(q, A_CHOL))
        .collect();
    let inc_m3s: Vec<f64> = q_inc.iter().map(|&q| mm_per_day_to_m3s(q, A_INC)).collect();
    let net = RiverNetwork::new(
        vec![
            Subcatchment {
                local_runoff: chol_m3s,
                downstream: Some(1),
                reach: Some(reach),
            },
            Subcatchment {
                local_runoff: inc_m3s,
                downstream: None,
                reach: None,
            },
        ],
        1.0,
    )
    .ok()?;
    Some(
        net.route_to_outlet()
            .iter()
            .map(|&q| m3s_to_mm_per_day(q, A_BALSA))
            .collect(),
    )
}

/// Masks qobs outside `range` with NaN (metrics skip it).
fn mask(qobs: &[f64], range: std::ops::Range<usize>) -> Vec<f64> {
    qobs.iter()
        .enumerate()
        .map(|(i, &v)| if range.contains(&i) { v } else { f64::NAN })
        .collect()
}

fn nse(obs: &[f64], sim: &[f64]) -> f64 {
    metrics::nse(&obs[WARMUP..], &sim[WARMUP..]).unwrap_or(f64::NAN)
}

fn main() {
    let base = "data/camels-cl/itata";
    let chol = read(&format!("{base}/cholguan.csv"));
    let inc = read(&format!("{base}/incremental.csv"));
    let outlet = read(&format!("{base}/balsa_outlet.csv"));
    let n = outlet.qobs.len();
    let mid = n / 2;
    let opt = Optimizer::Dds(DdsConfig {
        max_iter: 6000,
        seed: 42,
        ..Default::default()
    });
    let gb = calibrate::gr4j_default_bounds::<f64>();

    println!("Río Itata routing — split-sample ({n} days, split at {mid})\n");

    for (cal_name, val_name, cr, vr) in [
        ("A", "B", WARMUP..mid, mid..n),
        ("B", "A", mid..n, WARMUP..mid),
    ] {
        // --- Distributed (2 subcatchments + Muskingum), calibrated on outlet ---
        let cal_obs = mask(&outlet.qobs, cr.clone());
        let mut bounds: Vec<(f64, f64)> = gb.to_vec();
        bounds.extend_from_slice(&gb); // incremental GR4J
        bounds.push((0.1, 5.0)); // Muskingum K (days)
        let res = opt
            .maximize(&bounds, |x| match distributed_outlet(x, &chol, &inc) {
                Some(q) => metrics::nse(&cal_obs[WARMUP..], &q[WARMUP..]).unwrap_or(f64::NAN),
                None => f64::NAN,
            })
            .unwrap();
        let q_dist = distributed_outlet(&res.params, &chol, &inc).unwrap();
        let val_dist = nse(&mask(&outlet.qobs, vr.clone()), &q_dist);

        // --- Lumped baseline: one GR4J on the outlet forcing ---
        let lump = calibrate::calibrate_gr4j(
            &outlet.p,
            &outlet.pet,
            &cal_obs,
            WARMUP,
            rainflow_core::Objective::Nse,
            &gb,
            &opt,
        )
        .unwrap();
        let q_lump = gr4j_run(
            &[
                lump.params.x1,
                lump.params.x2,
                lump.params.x3,
                lump.params.x4,
            ],
            &outlet.p,
            &outlet.pet,
        )
        .unwrap();
        let val_lump = nse(&mask(&outlet.qobs, vr), &q_lump);

        println!(
            "cal {cal_name} → val {val_name}:  distributed cal NSE {:.3}, val {:.3}  (K={:.2} d)  |  lumped cal {:.3}, val {:.3}",
            res.value, val_dist, res.params[8], lump.value, val_lump
        );
    }
}
