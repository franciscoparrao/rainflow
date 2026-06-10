//! Shared unit-hydrograph convolution buffer helpers.

use num_traits::Float;

/// Pops the front of a UH convolution buffer, shifting the rest left.
pub(crate) fn shift_front<F: Float>(buf: &mut [F]) -> F {
    let out = buf[0];
    buf.copy_within(1.., 0);
    let n = buf.len();
    buf[n - 1] = F::zero();
    out
}
