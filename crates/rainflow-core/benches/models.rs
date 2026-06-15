//! Simulation throughput: GR4J, HBV (lumped and banded) and channel routing.
//!
//! Run with: cargo bench -p rainflow-core --bench models

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use rainflow_core::routing::{Muskingum, RiverNetwork, Subcatchment};
use rainflow_core::{ElevationBands, Gr4j, Gr4jParams, Hbv, HbvBands, HbvParams};

/// Deterministic synthetic forcing (P, PET, T) — same generator as the tests.
fn forcing(n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
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
        temp.push(8.0 + 12.0 * season + 3.0 * (next() - 0.5));
    }
    (p, pet, temp)
}

fn hbv_params() -> HbvParams<f64> {
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

fn bench_models(c: &mut Criterion) {
    // 10,000 days ≈ 27 years, the scale of a CAMELS-CL record.
    let n = 10_000;
    let (p, pet, temp) = forcing(n);

    let mut group = c.benchmark_group("simulate_10k_days");
    group.throughput(criterion::Throughput::Elements(n as u64));

    let gr4j = Gr4j::new(Gr4jParams {
        x1: 350.0,
        x2: -1.5,
        x3: 90.0,
        x4: 1.7,
    })
    .unwrap();
    group.bench_function("gr4j", |b| {
        b.iter(|| gr4j.run(black_box(&p), black_box(&pet)).unwrap())
    });

    let hbv = Hbv::new(hbv_params()).unwrap();
    group.bench_function("hbv_snow", |b| {
        b.iter(|| {
            hbv.run(black_box(&p), black_box(&pet), black_box(Some(&temp)))
                .unwrap()
        })
    });

    let bands = HbvBands::new(
        hbv_params(),
        ElevationBands::from_quantiles(1153.0, 3322.0, 5038.0, 5).unwrap(),
    )
    .unwrap();
    group.bench_function("hbv_bands_5", |b| {
        b.iter(|| {
            bands
                .run(black_box(&p), black_box(&pet), black_box(&temp))
                .unwrap()
        })
    });
    group.finish();

    // Channel routing on the same series length.
    let reach = Muskingum::new(2.0, 0.2).unwrap();
    let mut rg = c.benchmark_group("route_10k_days");
    rg.throughput(criterion::Throughput::Elements(n as u64));
    rg.bench_function("muskingum_reach", |b| {
        b.iter(|| reach.route(black_box(&p), 1.0))
    });
    let net = RiverNetwork::new(
        vec![
            Subcatchment {
                local_runoff: p.clone(),
                downstream: Some(1),
                reach: Some(reach),
            },
            Subcatchment {
                local_runoff: pet.clone(),
                downstream: None,
                reach: None,
            },
        ],
        1.0,
    )
    .unwrap();
    rg.bench_function("network_2_nodes", |b| b.iter(|| net.route_to_outlet()));
    rg.finish();
}

criterion_group!(benches, bench_models);
criterion_main!(benches);
