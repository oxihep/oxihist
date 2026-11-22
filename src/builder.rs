//! The quick-construction builder, modelled on scikit-hep
//! [`hist`](https://github.com/scikit-hep/hist)'s `Hist.new`: chain the axes,
//! then pick a storage.
//!
//! ```
//! use oxihist::Hist;
//!
//! // hist:  Hist.new.Reg(50, 0, 100, label="$p_T$").Weight()
//! let h = Hist::reg(50, 0.0, 100.0).label("$p_T$ [GeV]").weight().named("pt");
//! // 2-D, a regular and a variable axis:
//! let h2 = Hist::reg(40, -4.0, 4.0).label("x").var(&[0.0, 1.0, 4.0, 10.0]).label("y").double();
//! // A string category that grows as new categories are filled:
//! let mut c = Hist::str_cat::<&str>(&[]).growth().int64();
//! c.fill("ee");
//! ```
//!
//! | `hist` | builder |
//! |---|---|
//! | `.Reg(n, lo, hi)` | [`reg`](Hist::reg) |
//! | `.Var(edges)` | [`var`](Hist::var) |
//! | `.Integer(lo, hi)` | [`int`](Hist::int) |
//! | `.StrCat(cats)` / `.IntCat(cats)` | [`str_cat`](Hist::str_cat) / [`int_cat`](Hist::int_cat) |
//! | `.Boolean()` | [`boolean`](Hist::boolean) |
//! | `.Double()` / `.Int64()` / `.Weight()` | [`double`](Build1::double) / [`int64`](Build1::int64) / [`weight`](Build1::weight) |
//! | `.Mean()` / `.WeightedMean()` | [`mean`](Build1::mean) / [`weighted_mean`](Build1::weighted_mean) |

use crate::axis::Axis;
use crate::histogram::Histogram;
use crate::storage::{Double, Int64, Mean, Storage, Weight, WeightedMean};

/// The entry point of the builder (mirrors `hist`'s `Hist.new`). Start the
/// first axis with one of its constructors.
#[derive(Debug, Clone, Copy)]
pub struct Hist;

impl Hist {
    /// A regular axis of `bins` bins over `[lo, hi)`.
    ///
    /// # Panics
    ///
    /// If the binning is invalid (see [`Axis::regular`]).
    pub fn reg(bins: usize, lo: f64, hi: f64) -> Build1 {
        Build1::start(Axis::regular(bins, lo, hi))
    }

    /// A variable axis with the given edges.
    ///
    /// # Panics
    ///
    /// If the edges are invalid (see [`Axis::variable`]).
    pub fn var(edges: &[f64]) -> Build1 {
        Build1::start(Axis::variable(edges))
    }

    /// An integer axis, one bin per integer in `lo..hi`.
    ///
    /// # Panics
    ///
    /// If `lo >= hi`.
    pub fn int(lo: i64, hi: i64) -> Build1 {
        Build1::start(Axis::integer(lo, hi))
    }

    /// A string-category axis.
    pub fn str_cat<S: AsRef<str>>(categories: &[S]) -> Build1 {
        Build1::start(Axis::str_category(categories))
    }

    /// An integer-category axis.
    pub fn int_cat(categories: &[i64]) -> Build1 {
        Build1::start(Axis::int_category(categories))
    }

    /// A boolean axis.
    pub fn boolean() -> Build1 {
        Build1::start(Axis::boolean())
    }

    /// Any axis built with [`Axis`]'s constructors.
    pub fn axis(axis: Axis) -> Build1 {
        Build1::start(axis)
    }
}

macro_rules! builder {
    ($name:ident, $d:literal, $doc:literal, next = $($next:ident)?) => {
        #[doc = $doc]
        #[derive(Debug, Clone)]
        #[must_use]
        pub struct $name {
            axes: Vec<Axis>,
        }

        impl $name {
            /// Label the last axis.
            pub fn label(mut self, label: impl Into<String>) -> Self {
                self.last().label = label.into();
                self
            }

            /// Name the last axis.
            pub fn axis_name(mut self, name: impl Into<String>) -> Self {
                self.last().name = name.into();
                self
            }

            /// Let the last axis grow (category axes only; ignored otherwise).
            pub fn growth(mut self) -> Self {
                let a = self.last().clone();
                *self.last() = a.with_growth();
                self
            }

            /// Drop the last axis's flow bins.
            pub fn no_flow(mut self) -> Self {
                let a = self.last().clone();
                *self.last() = a.with_underflow(false).with_overflow(false);
                self
            }

            fn last(&mut self) -> &mut Axis {
                self.axes.last_mut().expect("a builder always has an axis")
            }

            fn finish<S: Storage>(self) -> Histogram<S, $d> {
                Histogram::from_axes(self.axes)
            }

            /// `f64` counts ([`Double`]).
            pub fn double(self) -> Histogram<Double, $d> {
                self.finish()
            }

            /// `i64` counts ([`Int64`]).
            pub fn int64(self) -> Histogram<Int64, $d> {
                self.finish()
            }

            /// Weighted counts with their variance ([`Weight`]).
            pub fn weight(self) -> Histogram<Weight, $d> {
                self.finish()
            }

            /// A profile: per-bin mean of a sample ([`Mean`]).
            pub fn mean(self) -> Histogram<Mean, $d> {
                self.finish()
            }

            /// A weighted profile ([`WeightedMean`]).
            pub fn weighted_mean(self) -> Histogram<WeightedMean, $d> {
                self.finish()
            }

            /// Any storage, by type: `.storage::<Weight>()`.
            pub fn storage<S: Storage>(self) -> Histogram<S, $d> {
                self.finish()
            }

            $(
            /// Add a regular axis.
            ///
            /// # Panics
            ///
            /// If the binning is invalid (see [`Axis::regular`]).
            pub fn reg(self, bins: usize, lo: f64, hi: f64) -> $next {
                self.push(Axis::regular(bins, lo, hi))
            }

            /// Add a variable axis.
            ///
            /// # Panics
            ///
            /// If the edges are invalid (see [`Axis::variable`]).
            pub fn var(self, edges: &[f64]) -> $next {
                self.push(Axis::variable(edges))
            }

            /// Add an integer axis.
            ///
            /// # Panics
            ///
            /// If `lo >= hi`.
            pub fn int(self, lo: i64, hi: i64) -> $next {
                self.push(Axis::integer(lo, hi))
            }

            /// Add a string-category axis.
            pub fn str_cat<S: AsRef<str>>(self, categories: &[S]) -> $next {
                self.push(Axis::str_category(categories))
            }

            /// Add an integer-category axis.
            pub fn int_cat(self, categories: &[i64]) -> $next {
                self.push(Axis::int_category(categories))
            }

            /// Add a boolean axis.
            pub fn boolean(self) -> $next {
                self.push(Axis::boolean())
            }

            /// Add any axis.
            pub fn axis(self, axis: Axis) -> $next {
                self.push(axis)
            }

            fn push(mut self, axis: Axis) -> $next {
                self.axes.push(axis);
                $next { axes: self.axes }
            }
            )?
        }
    };
}

builder!(
    Build1,
    1,
    "A 1-D histogram under construction.",
    next = Build2
);
builder!(
    Build2,
    2,
    "A 2-D histogram under construction.",
    next = Build3
);
builder!(Build3, 3, "A 3-D histogram under construction.", next =);

impl Build1 {
    fn start(axis: Axis) -> Build1 {
        Build1 { axes: vec![axis] }
    }
}
