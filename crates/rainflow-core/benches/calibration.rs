//! Calibration cost: a full DDS / SCE-UA run on a realistic record.
//! This is the operational hot path (one calibration per catchment).
//!
//! Run with: cargo bench -p rainflow-core --bench calibration

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use rainflow_core::calibrate::{self, DdsConfig, Objective, Optimizer, SceConfig};
use rainflow_core::{Gr4j, Gr4jParams};

fn forcing(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut seed: u64 = 42;
    let mut next = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 33) as f64 / (1u64 << 31) as f64
    };
    let (mut p, mut pet) = (Vec::new(), Vec::new());
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

fn bench_calibration(c: &mut Criterion) {
    // ~16-year daily record; synthetic-truth observed discharge.
    let n = 6000;
    let (p, pet) = forcing(n);
    let qobs = Gr4j::new(Gr4jParams {
        x1: 350.0,
        x2: -1.5,
        x3: 90.0,
        x4: 1.7,
    })
    .unwrap()
    .run(&p, &pet)
    .unwrap();
    let bounds = calibrate::gr4j_default_bounds::<f64>();

    let mut group = c.benchmark_group("calibrate_gr4j_6000_days");
    // Calibration runs take ~seconds; keep the sample count modest.
    group.sample_size(10);

    let dds = Optimizer::Dds(DdsConfig {
        max_iter: 2000,
        seed: 42,
        ..Default::default()
    });
    group.bench_function("dds_2000_evals", |b| {
        b.iter(|| {
            calibrate::calibrate_gr4j(
                black_box(&p),
                black_box(&pet),
                black_box(&qobs),
                365,
                Objective::Nse,
                &bounds,
                &dds,
            )
            .unwrap()
        })
    });

    let sce = Optimizer::Sce(SceConfig {
        complexes: 4,
        max_iter: 4000,
        seed: 42,
    });
    group.bench_function("sce_4000_evals", |b| {
        b.iter(|| {
            calibrate::calibrate_gr4j(
                black_box(&p),
                black_box(&pet),
                black_box(&qobs),
                365,
                Objective::Nse,
                &bounds,
                &sce,
            )
            .unwrap()
        })
    });
    group.finish();
}

criterion_group!(benches, bench_calibration);
criterion_main!(benches);
