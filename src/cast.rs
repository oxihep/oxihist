//! The index conversions binning needs, each justified once here: bins live
//! in a `Vec`, so every count and index is far below `isize::MAX`. They are
//! plain casts, inlined, on the fill paths.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "each function states the range that makes its cast exact or saturating"
)]

/// A bin count or index as the signed bin index (`-1` is the underflow):
/// below `isize::MAX`, as it indexes a `Vec`.
#[inline(always)]
pub(crate) fn signed(n: usize) -> isize {
    n as isize
}

/// A signed bin index the caller has checked to be ≥ 0.
#[inline(always)]
pub(crate) fn unsigned(i: isize) -> usize {
    debug_assert!(i >= 0, "a negative bin index");
    i as usize
}

/// A difference `hi − lo ≥ 0` of an integer axis's bounds, as its bin count
/// (an axis of more bins than memory holds cannot be built).
#[inline(always)]
pub(crate) fn span(d: i64) -> usize {
    debug_assert!(d >= 0, "a negative span");
    d as usize
}

/// As [`span`], as a signed bin index.
#[inline(always)]
pub(crate) fn span_signed(d: i64) -> isize {
    d as isize
}

/// `floor(x)` as an `i64`, saturating (callers handle `NaN` first).
#[inline(always)]
pub(crate) fn floor_i64(x: f64) -> i64 {
    x.floor() as i64
}

/// `floor(x)` as an `isize`, saturating (callers clamp to the bin range).
#[inline(always)]
pub(crate) fn floor_isize(x: f64) -> isize {
    x.floor() as isize
}

/// `round(x)` as an `i64`, saturating: an integer count.
#[inline(always)]
pub(crate) fn round_i64(x: f64) -> i64 {
    x.round() as i64
}
