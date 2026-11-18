//! Histograms for high-energy physics in pure Rust, modelled on
//! [boost-histogram](https://github.com/scikit-hep/boost-histogram) and
//! scikit-hep [`hist`](https://github.com/scikit-hep/hist), with no file
//! format attached.
//!
//! ```
//! use oxihist::Hist;
//!
//! let mut h = Hist::reg(50, 0.0, 100.0).label("$p_T$ [GeV]").name("pt").weight();
//! h.fill(42.0);
//! h.fill_weight(17.0, 0.5);
//! assert_eq!(h.integral(false), 1.5);
//!
//! // A 2-D histogram with a growing string category.
//! let mut by_channel = Hist::reg(10, 0.0, 200.0).str_cat::<&str>(&[]).growth().int64();
//! by_channel.fill(91.0, "ee");
//! by_channel.fill(125.0, "mumu");
//! assert_eq!(by_channel.axis(1).size(), 2);
//! ```
//!
//! - **Axes** ([`Axis`]): regular, variable, integer, string and integer
//!   categories (optionally growing), boolean; each with optional under- and
//!   overflow bins.
//! - **Storages** ([`storage`]): [`Double`], [`Int64`], [`Weight`] (value and
//!   variance), [`Mean`] and [`WeightedMean`] (profiles).
//! - **[`Histogram<S, D>`](Histogram)** for `D` = 1, 2, 3: fills, bin access with
//!   or without flow, adding and scaling, projections, rebinning, 1-D
//!   statistics and sampling, and [`ThreadedHist`] for multithreaded fills.
//!
//! Features: `rayon` ([`fill_par`]), `serde`, `plot` (draw with
//! [oxiplot](https://github.com/oxihep/oxiplot)) and `fit` (fit with
//! [oxifit](https://github.com/oxihep/oxifit)). To read and write ROOT files,
//! [oxiroot](https://github.com/oxihep/oxiroot) converts these histograms to
//! and from ROOT's `TH1`/`TH2`/`TH3`/`TProfile`.

mod axis;
mod builder;
mod error;
mod histogram;
mod integrations;
mod stats;
pub mod storage;
mod threaded;

pub use axis::{Axis, AxisKind, Value};
pub use builder::{Build1, Build2, Build3, Hist};
pub use error::{Error, Result};
pub use histogram::Histogram;
pub use stats::{BinnedSampler, Random};
pub use storage::{
    CountStorage, Double, Int64, Mean, MeanBin, MeanStorage, Storage, Weight, WeightedMean,
    WeightedMeanBin, WeightedSum,
};
#[cfg(feature = "rayon")]
pub use threaded::fill_par;
pub use threaded::{Mergeable, ThreadedHist};
