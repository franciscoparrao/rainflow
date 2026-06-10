use thiserror::Error;

/// Errors produced by model construction and simulation.
#[derive(Debug, Error)]
pub enum Error {
    #[error("invalid parameter `{name}`: {reason}")]
    InvalidParameter { name: &'static str, reason: String },

    #[error("forcing series length mismatch: precip has {precip} steps, pet has {pet}")]
    ForcingLengthMismatch { precip: usize, pet: usize },
}
