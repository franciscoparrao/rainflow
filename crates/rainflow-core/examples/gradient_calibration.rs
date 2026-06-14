//! Gradient-based calibration of GR4J — closing the autodiff-first loop.
//!
//! The model cores are generic over `F: num_traits::Float`. The `autodiff`
//! crate's forward-mode scalar `F1` *also* implements `num_traits::Float`, so
//! it flows through `Gr4j` unchanged and yields an analytic gradient of the
//! loss with respect to the parameters — no change to the core, no finite
//! differences. This is the substrate the project was designed for: gradient
//! calibration today, physics+ML hybrids (δHBV-style) later.
//!
//! Run with: cargo run -p rainflow-core --example gradient_calibration --release

use autodiff::F1;
use num_traits::Float;
use rainflow_core::{Gr4j, Gr4jParams};

// GR4J search bounds (x1, x2, x3, x4); calibration runs in normalized [0,1]
// space so the wildly different parameter scales share one step size.
const BOUNDS: [(f64, f64); 4] = [(1.0, 2500.0), (-5.0, 5.0), (1.0, 1000.0), (0.5, 10.0)];
const WARMUP: usize = 365;

/// Maps a normalized vector in [0,1]^4 to physical GR4J parameters.
fn denormalize<F: Float>(theta: &[F; 4]) -> Gr4jParams<F> {
    let p = |i: usize| {
        let (lo, up) = BOUNDS[i];
        let (lo, up) = (F::from(lo).unwrap(), F::from(up).unwrap());
        lo + theta[i] * (up - lo)
    };
    Gr4jParams {
        x1: p(0),
        x2: p(1),
        x3: p(2),
        x4: p(3),
    }
}

/// Mean squared error of GR4J against `qobs`, generic over the scalar type.
/// With `F = f64` it is the plain loss; with `F = F1` it carries the
/// derivative through the whole simulation.
fn mse<F: Float>(theta: &[F; 4], p: &[F], pet: &[F], qobs: &[F]) -> F {
    let model = Gr4j::new(denormalize(theta)).expect("params in-range");
    let q = model.run(p, pet).expect("equal lengths");
    let mut sum = F::zero();
    let mut n = F::zero();
    for t in WARMUP..q.len() {
        let d = q[t] - qobs[t];
        sum = sum + d * d;
        n = n + F::one();
    }
    sum / n
}

/// Analytic gradient of the MSE w.r.t. the normalized parameters, by
/// forward-mode autodiff: one pass per parameter, seeding that coordinate.
fn grad(theta: &[f64; 4], p: &[F1], pet: &[F1], qobs: &[F1]) -> [f64; 4] {
    let mut g = [0.0; 4];
    for i in 0..4 {
        let mut t = [
            F1::cst(theta[0]),
            F1::cst(theta[1]),
            F1::cst(theta[2]),
            F1::cst(theta[3]),
        ];
        t[i] = F1::var(theta[i]); // seed coordinate i (dθ_i = 1)
        g[i] = mse(&t, p, pet, qobs).deriv();
    }
    g
}

fn nse(theta: &[f64; 4], p: &[f64], pet: &[f64], qobs: &[f64]) -> f64 {
    let q = Gr4j::new(denormalize(theta)).unwrap().run(p, pet).unwrap();
    rainflow_core::metrics::nse(&qobs[WARMUP..], &q[WARMUP..]).unwrap()
}

fn synthetic_forcing(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut seed: u64 = 42;
    let mut next = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 33) as f64 / (1u64 << 31) as f64
    };
    let (mut p, mut pet) = (Vec::with_capacity(n), Vec::with_capacity(n));
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
        pet.push((3.5 + 2.5 * (2.0 * std::f64::consts::PI * (doy - 15.0) / 365.25).sin()).max(0.1));
    }
    (p, pet)
}

fn main() {
    // Synthetic truth: generate qobs with known parameters, then recover them.
    let (p, pet) = synthetic_forcing(2000);
    let truth = [350.0, -1.5, 90.0, 1.7];
    let qobs = Gr4j::new(Gr4jParams {
        x1: truth[0],
        x2: truth[1],
        x3: truth[2],
        x4: truth[3],
    })
    .unwrap()
    .run(&p, &pet)
    .unwrap();

    let p1: Vec<F1> = p.iter().map(|&v| F1::cst(v)).collect();
    let pet1: Vec<F1> = pet.iter().map(|&v| F1::cst(v)).collect();
    let qobs1: Vec<F1> = qobs.iter().map(|&v| F1::cst(v)).collect();

    // --- 1. Validate the analytic gradient against central finite differences.
    let theta0 = [0.5, 0.6, 0.4, 0.3];
    let g_ad = grad(&theta0, &p1, &pet1, &qobs1);
    let h = 1e-6;
    println!("Gradient check (analytic vs finite-difference) at θ₀:");
    for i in 0..4 {
        let mut tp = theta0;
        let mut tm = theta0;
        tp[i] += h;
        tm[i] -= h;
        let fd = (mse(&tp, &p, &pet, &qobs) - mse(&tm, &p, &pet, &qobs)) / (2.0 * h);
        let rel = (g_ad[i] - fd).abs() / fd.abs().max(1e-12);
        println!(
            "  ∂L/∂θ{i}: autodiff {:+.6e}  fd {:+.6e}  rel.err {:.2e}",
            g_ad[i], fd, rel
        );
    }

    // --- 2. Gradient descent (Adam) in normalized space from a poor start.
    let mut theta = [0.1, 0.9, 0.8, 0.6]; // deliberately far from the truth
    let (mut m, mut v) = ([0.0; 4], [0.0; 4]);
    let (lr, b1, b2, eps) = (0.02, 0.9, 0.999, 1e-8);
    println!("\nGradient descent (Adam, normalized space):");
    println!("  start    NSE {:.4}", nse(&theta, &p, &pet, &qobs));
    for step in 1..=400 {
        let g = grad(&theta, &p1, &pet1, &qobs1);
        for i in 0..4 {
            m[i] = b1 * m[i] + (1.0 - b1) * g[i];
            v[i] = b2 * v[i] + (1.0 - b2) * g[i] * g[i];
            let mhat = m[i] / (1.0 - b1.powi(step));
            let vhat = v[i] / (1.0 - b2.powi(step));
            theta[i] = (theta[i] - lr * mhat / (vhat.sqrt() + eps)).clamp(0.0, 1.0);
        }
        if step % 100 == 0 {
            println!("  step {step:3}  NSE {:.4}", nse(&theta, &p, &pet, &qobs));
        }
    }

    let fitted = denormalize(&theta);
    println!(
        "\nrecovered: x1={:.1} x2={:.2} x3={:.1} x4={:.2}",
        fitted.x1, fitted.x2, fitted.x3, fitted.x4
    );
    println!(
        "truth:     x1={:.1} x2={:.2} x3={:.1} x4={:.2}",
        truth[0], truth[1], truth[2], truth[3]
    );
    println!("final NSE: {:.5}", nse(&theta, &p, &pet, &qobs));
}
