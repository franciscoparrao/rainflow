//! # rainflow-core
//!
//! Conceptual rainfall–runoff model cores with calibration metrics.
//!
//! Design principles:
//!
//! - **Autodiff-first**: every model and metric is generic over `F: num_traits::Float`,
//!   so dual-number or tape-based scalar types (e.g. from an AD crate) flow through
//!   the model unchanged, enabling gradient-based calibration and physics+ML hybrids.
//! - **No I/O in the core**: forcing data arrives as slices; file formats live in
//!   the CLI / bindings crates.
//! - **Numerical parity**: GR4J follows Perrin et al. (2003) and is cross-checked
//!   against the airGR reference implementation.

pub mod error;
pub mod gr4j;
pub mod metrics;

pub use error::Error;
pub use gr4j::{Gr4j, Gr4jParams, Gr4jState};
