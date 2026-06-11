//! Python bindings for rainflow-core (PyO3).
//!
//! Mirrors the Rust API on `f64`: model classes [`Gr4j`]/[`Hbv`], the
//! goodness-of-fit metrics, and DDS/SCE-UA calibration. Sequences cross the
//! boundary as plain Python lists of floats; missing observations are encoded
//! as `nan`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use rainflow_core::calibrate::{self, DdsConfig, Objective, Optimizer, SceConfig};
use rainflow_core::{ElevationBands, Gr4j as CoreGr4j, Gr4jParams, Hbv as CoreHbv, HbvParams};

/// Maps a core error to a Python `ValueError`.
fn err(e: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn parse_objective(name: &str) -> PyResult<Objective> {
    match name.to_lowercase().as_str() {
        "nse" => Ok(Objective::Nse),
        "kge" => Ok(Objective::Kge),
        "lognse" | "log_nse" => Ok(Objective::LogNse),
        other => Err(PyValueError::new_err(format!("unknown objective {other:?}"))),
    }
}

fn build_optimizer(algorithm: &str, iterations: usize, seed: u64, complexes: usize) -> PyResult<Optimizer> {
    match algorithm.to_lowercase().as_str() {
        "dds" => Ok(Optimizer::Dds(DdsConfig {
            max_iter: iterations,
            seed,
            ..Default::default()
        })),
        "sce" => Ok(Optimizer::Sce(SceConfig {
            complexes,
            max_iter: iterations,
            seed,
        })),
        other => Err(PyValueError::new_err(format!("unknown algorithm {other:?}"))),
    }
}

/// GR4J daily rainfall–runoff model (Perrin et al. 2003).
#[pyclass]
struct Gr4j {
    inner: CoreGr4j<f64>,
}

#[pymethods]
impl Gr4j {
    #[new]
    fn new(x1: f64, x2: f64, x3: f64, x4: f64) -> PyResult<Self> {
        Ok(Self {
            inner: CoreGr4j::new(Gr4jParams { x1, x2, x3, x4 }).map_err(err)?,
        })
    }

    /// Simulates discharge (mm) from precipitation and PET series (mm).
    fn run(&self, precip: Vec<f64>, pet: Vec<f64>) -> PyResult<Vec<f64>> {
        self.inner.run(&precip, &pet).map_err(err)
    }
}

/// HBV-light daily rainfall–runoff model (Seibert & Vis 2012).
#[pyclass]
struct Hbv {
    inner: CoreHbv<f64>,
}

#[pymethods]
impl Hbv {
    #[new]
    #[pyo3(signature = (tt, cfmax, sfcf, fc, lp, beta, k0, k1, k2, uzl, perc, maxbas, cfr=0.05, cwh=0.1))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        tt: f64,
        cfmax: f64,
        sfcf: f64,
        fc: f64,
        lp: f64,
        beta: f64,
        k0: f64,
        k1: f64,
        k2: f64,
        uzl: f64,
        perc: f64,
        maxbas: f64,
        cfr: f64,
        cwh: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: CoreHbv::new(HbvParams {
                tt,
                cfmax,
                sfcf,
                cfr,
                cwh,
                fc,
                lp,
                beta,
                k0,
                k1,
                k2,
                uzl,
                perc,
                maxbas,
            })
            .map_err(err)?,
        })
    }

    /// Simulates discharge (mm). `temp` (°C) is optional; without it the snow
    /// routine is bypassed (pluvial catchments).
    #[pyo3(signature = (precip, pet, temp=None))]
    fn run(&self, precip: Vec<f64>, pet: Vec<f64>, temp: Option<Vec<f64>>) -> PyResult<Vec<f64>> {
        self.inner
            .run(&precip, &pet, temp.as_deref())
            .map_err(err)
    }
}

/// Nash–Sutcliffe efficiency.
#[pyfunction]
fn nse(obs: Vec<f64>, sim: Vec<f64>) -> PyResult<f64> {
    rainflow_core::metrics::nse(&obs, &sim).map_err(err)
}

/// Kling–Gupta efficiency (2009).
#[pyfunction]
fn kge(obs: Vec<f64>, sim: Vec<f64>) -> PyResult<f64> {
    rainflow_core::metrics::kge(&obs, &sim).map_err(err)
}

/// NSE on log-transformed flows.
#[pyfunction]
fn log_nse(obs: Vec<f64>, sim: Vec<f64>) -> PyResult<f64> {
    rainflow_core::metrics::log_nse(&obs, &sim).map_err(err)
}

/// Percent bias.
#[pyfunction]
fn pbias(obs: Vec<f64>, sim: Vec<f64>) -> PyResult<f64> {
    rainflow_core::metrics::pbias(&obs, &sim).map_err(err)
}

/// Calibrates GR4J and returns `{params: {x1,..}, value, evaluations}`.
#[pyfunction]
#[pyo3(signature = (precip, pet, qobs, warmup=365, objective="nse", algorithm="dds", iterations=2000, seed=42, complexes=4))]
#[allow(clippy::too_many_arguments)]
fn calibrate_gr4j(
    py: Python<'_>,
    precip: Vec<f64>,
    pet: Vec<f64>,
    qobs: Vec<f64>,
    warmup: usize,
    objective: &str,
    algorithm: &str,
    iterations: usize,
    seed: u64,
    complexes: usize,
) -> PyResult<PyObject> {
    let obj = parse_objective(objective)?;
    let opt = build_optimizer(algorithm, iterations, seed, complexes)?;
    let cal = calibrate::calibrate_gr4j(
        &precip,
        &pet,
        &qobs,
        warmup,
        obj,
        &calibrate::gr4j_default_bounds(),
        &opt,
    )
    .map_err(err)?;

    let params = pyo3::types::PyDict::new(py);
    params.set_item("x1", cal.params.x1)?;
    params.set_item("x2", cal.params.x2)?;
    params.set_item("x3", cal.params.x3)?;
    params.set_item("x4", cal.params.x4)?;
    let out = pyo3::types::PyDict::new(py);
    out.set_item("params", params)?;
    out.set_item("value", cal.value)?;
    out.set_item("evaluations", cal.evaluations)?;
    Ok(out.into())
}

/// Calibrates HBV-light. `temp` enables the snow routine (its parameters are
/// then calibrated too). Returns `{params: {..}, value, evaluations}`.
#[pyfunction]
#[pyo3(signature = (precip, pet, qobs, temp=None, warmup=365, objective="nse", algorithm="dds", iterations=2000, seed=42, complexes=4))]
#[allow(clippy::too_many_arguments)]
fn calibrate_hbv(
    py: Python<'_>,
    precip: Vec<f64>,
    pet: Vec<f64>,
    qobs: Vec<f64>,
    temp: Option<Vec<f64>>,
    warmup: usize,
    objective: &str,
    algorithm: &str,
    iterations: usize,
    seed: u64,
    complexes: usize,
) -> PyResult<PyObject> {
    let obj = parse_objective(objective)?;
    let opt = build_optimizer(algorithm, iterations, seed, complexes)?;
    let cal = calibrate::calibrate_hbv(&precip, &pet, temp.as_deref(), &qobs, warmup, obj, &opt)
        .map_err(err)?;

    let p = cal.params;
    let params = pyo3::types::PyDict::new(py);
    for (k, v) in [
        ("tt", p.tt), ("cfmax", p.cfmax), ("sfcf", p.sfcf), ("cfr", p.cfr), ("cwh", p.cwh),
        ("fc", p.fc), ("lp", p.lp), ("beta", p.beta), ("k0", p.k0), ("k1", p.k1),
        ("k2", p.k2), ("uzl", p.uzl), ("perc", p.perc), ("maxbas", p.maxbas),
    ] {
        params.set_item(k, v)?;
    }
    let out = pyo3::types::PyDict::new(py);
    out.set_item("params", params)?;
    out.set_item("value", cal.value)?;
    out.set_item("evaluations", cal.evaluations)?;
    Ok(out.into())
}

/// Builds `n` equal-area elevation bands from min/median/max elevations and
/// returns the band `(elevation, area_fraction)` pairs as a list of tuples.
#[pyfunction]
fn hypsometric_bands(min: f64, median: f64, max: f64, n: usize) -> PyResult<Vec<(f64, f64)>> {
    let bands = ElevationBands::from_quantiles(min, median, max, n).map_err(err)?;
    Ok(bands
        .bands
        .iter()
        .map(|b| (b.elevation, b.area_fraction))
        .collect())
}

#[pymodule]
fn rainflow(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Gr4j>()?;
    m.add_class::<Hbv>()?;
    m.add_function(wrap_pyfunction!(nse, m)?)?;
    m.add_function(wrap_pyfunction!(kge, m)?)?;
    m.add_function(wrap_pyfunction!(log_nse, m)?)?;
    m.add_function(wrap_pyfunction!(pbias, m)?)?;
    m.add_function(wrap_pyfunction!(calibrate_gr4j, m)?)?;
    m.add_function(wrap_pyfunction!(calibrate_hbv, m)?)?;
    m.add_function(wrap_pyfunction!(hypsometric_bands, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
