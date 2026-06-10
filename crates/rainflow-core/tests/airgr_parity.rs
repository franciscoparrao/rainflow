//! Numerical parity against the airGR reference implementation (INRAE).
//!
//! The fixture was generated with `airGR::RunModel_GR4J` on the package's
//! L0123001 example catchment (P and E from BasinObs), parameters
//! x1=350, x2=-1.5, x3=90, x4=1.7, no warm-up, default initial store levels
//! (S = 0.3·x1, R = 0.5·x3). airGR is GPL-2 licensed; the fixture is used for
//! validation only and is not part of the published library.

use rainflow_core::{Gr4j, Gr4jParams};

#[test]
fn gr4j_matches_airgr_reference() {
    let raw = include_str!("fixtures/airgr_gr4j_l0123001.csv");
    let mut precip = Vec::new();
    let mut pet = Vec::new();
    let mut q_ref = Vec::new();
    for line in raw.lines().skip(1) {
        let mut cols = line.split(',');
        let _date = cols.next();
        let parse = |c: Option<&str>| c.and_then(|v| v.trim_matches('"').parse::<f64>().ok());
        precip.push(parse(cols.next()).expect("precip column"));
        pet.push(parse(cols.next()).expect("pet column"));
        q_ref.push(parse(cols.next()).expect("qsim_airgr column"));
    }
    assert_eq!(precip.len(), 2000, "fixture should hold 2000 steps");

    let model = Gr4j::new(Gr4jParams {
        x1: 350.0,
        x2: -1.5,
        x3: 90.0,
        x4: 1.7,
    })
    .expect("valid parameters");
    let qsim = model.run(&precip, &pet).expect("simulation runs");

    let max_diff = qsim
        .iter()
        .zip(&q_ref)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    // The residual is the airGR CSV round-off, not model divergence.
    assert!(
        max_diff < 1e-5,
        "GR4J diverges from airGR: max abs diff {max_diff:.3e} mm"
    );
}
