//! The crate's error type.

use std::fmt;

/// What can go wrong building, filling or combining histograms.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// A fill value does not fit the axis it was given to (a string on a
    /// numeric axis, a float on an integer-category axis, …).
    ValueType {
        /// The axis, counted from 0.
        axis: usize,
        /// The kind of value the axis takes.
        expected: &'static str,
        /// The value that was given.
        found: String,
    },
    /// Two histograms that must share a binning do not.
    BinningMismatch(String),
    /// An axis was asked for something it cannot do (rebin a category axis, a
    /// rebin factor that does not divide the bins, an axis index out of range).
    InvalidAxis(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::ValueType {
                axis,
                expected,
                found,
            } => write!(f, "axis {axis} takes {expected}, got {found}"),
            Error::BinningMismatch(m) => write!(f, "binning mismatch: {m}"),
            Error::InvalidAxis(m) => write!(f, "invalid axis operation: {m}"),
        }
    }
}

impl std::error::Error for Error {}

/// The crate's result alias.
pub type Result<T> = std::result::Result<T, Error>;
