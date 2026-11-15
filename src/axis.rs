//! Histogram axes: how a fill value maps to a bin.
//!
//! Every axis has `size()` in-range bins, numbered `0..size`. Values below the
//! range land in the **underflow** bin (index `-1`) and values above it in the
//! **overflow** bin (index `size`), when the axis has them; otherwise they are
//! dropped. This is boost-histogram's convention.

use std::fmt;

use crate::error::{Error, Result};

/// A value to fill along one axis. Numbers, integers, strings and booleans all
/// convert into it, so `h.fill(3.2)`, `h.fill(7)`, `h.fill("ee")` and
/// `h.fill(true)` all work on an axis of the matching kind.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Value<'a> {
    /// A real number, for [`Regular`](AxisKind::Regular) and
    /// [`Variable`](AxisKind::Variable) axes (and integer axes, which floor it).
    Num(f64),
    /// An integer, for [`Integer`](AxisKind::Integer) and
    /// [`IntCategory`](AxisKind::IntCategory) axes (and numeric axes).
    Int(i64),
    /// A string, for [`StrCategory`](AxisKind::StrCategory) axes.
    Str(&'a str),
    /// A boolean, for [`Boolean`](AxisKind::Boolean) axes.
    Bool(bool),
}

impl fmt::Display for Value<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Value::Num(x) => write!(f, "{x}"),
            Value::Int(i) => write!(f, "{i}"),
            Value::Str(s) => write!(f, "{s:?}"),
            Value::Bool(b) => write!(f, "{b}"),
        }
    }
}

macro_rules! value_from {
    ($($t:ty => $v:ident as $c:ty),* $(,)?) => {$(
        impl From<$t> for Value<'_> {
            fn from(x: $t) -> Self {
                Value::$v(x as $c)
            }
        }
    )*};
}
value_from!(
    f64 => Num as f64, f32 => Num as f64,
    i64 => Int as i64, i32 => Int as i64, i16 => Int as i64, i8 => Int as i64,
    u32 => Int as i64, u16 => Int as i64, u8 => Int as i64, usize => Int as i64,
);

impl<'a> From<&'a str> for Value<'a> {
    fn from(s: &'a str) -> Self {
        Value::Str(s)
    }
}

impl<'a> From<&'a String> for Value<'a> {
    fn from(s: &'a String) -> Self {
        Value::Str(s)
    }
}

impl From<bool> for Value<'_> {
    fn from(b: bool) -> Self {
        Value::Bool(b)
    }
}

/// The binning of an [`Axis`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AxisKind {
    /// `bins` equal-width bins over `[lo, hi)`.
    Regular {
        /// The number of bins.
        bins: usize,
        /// The low edge.
        lo: f64,
        /// The high edge (exclusive).
        hi: f64,
    },
    /// Bins between the given ascending `edges` (`bins + 1` of them).
    Variable {
        /// The bin edges, strictly ascending.
        edges: Vec<f64>,
    },
    /// One bin per integer in `lo..hi`.
    Integer {
        /// The first integer.
        lo: i64,
        /// One past the last integer.
        hi: i64,
    },
    /// One bin per string category.
    StrCategory {
        /// The categories, in bin order.
        categories: Vec<String>,
        /// Whether filling an unknown category adds it as a new bin.
        growth: bool,
    },
    /// One bin per integer category.
    IntCategory {
        /// The categories, in bin order.
        categories: Vec<i64>,
        /// Whether filling an unknown category adds it as a new bin.
        growth: bool,
    },
    /// Two bins, `false` (0) and `true` (1).
    Boolean,
}

/// One histogram axis: a binning plus a name, a label and whether it keeps
/// under- and overflow bins.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Axis {
    kind: AxisKind,
    /// A short identifier, e.g. `"pt"`.
    pub name: String,
    /// A human-readable label, e.g. `"$p_T$ [GeV]"`; plotting uses it.
    pub label: String,
    underflow: bool,
    overflow: bool,
}

/// Where a value falls on an axis.
pub(crate) enum Locate {
    /// Bin index `-1..=size`.
    Index(isize),
    /// An unknown category on a growing axis: add it, then fill the new bin.
    Grow,
}

impl Axis {
    fn with_kind(kind: AxisKind, underflow: bool, overflow: bool) -> Axis {
        Axis {
            kind,
            name: String::new(),
            label: String::new(),
            underflow,
            overflow,
        }
    }

    /// `bins` equal-width bins over `[lo, hi)`, with under- and overflow.
    ///
    /// # Panics
    ///
    /// If `bins` is 0 or the range is not finite and increasing; see
    /// [`try_regular`](Axis::try_regular).
    #[must_use]
    pub fn regular(bins: usize, lo: f64, hi: f64) -> Axis {
        Self::try_regular(bins, lo, hi).unwrap_or_else(|e| panic!("{e}"))
    }

    /// Like [`regular`](Axis::regular), returning an error for an invalid
    /// binning instead of panicking.
    pub fn try_regular(bins: usize, lo: f64, hi: f64) -> Result<Axis> {
        if bins == 0 || !(lo.is_finite() && hi.is_finite() && lo < hi) {
            return Err(Error::InvalidAxis(format!(
                "a regular axis needs bins > 0 and a finite lo < hi, got {bins} bins over [{lo}, {hi})"
            )));
        }
        Ok(Self::with_kind(
            AxisKind::Regular { bins, lo, hi },
            true,
            true,
        ))
    }

    /// Bins between `edges` (at least two, strictly ascending), with under- and
    /// overflow.
    ///
    /// # Panics
    ///
    /// If the edges are invalid; see [`try_variable`](Axis::try_variable).
    #[must_use]
    pub fn variable(edges: &[f64]) -> Axis {
        Self::try_variable(edges).unwrap_or_else(|e| panic!("{e}"))
    }

    /// Like [`variable`](Axis::variable), returning an error for invalid edges
    /// instead of panicking.
    pub fn try_variable(edges: &[f64]) -> Result<Axis> {
        if edges.len() < 2
            || !edges.iter().all(|e| e.is_finite())
            || !edges.windows(2).all(|w| w[0] < w[1])
        {
            return Err(Error::InvalidAxis(
                "a variable axis needs at least two finite, strictly ascending edges".into(),
            ));
        }
        Ok(Self::with_kind(
            AxisKind::Variable {
                edges: edges.to_vec(),
            },
            true,
            true,
        ))
    }

    /// One bin per integer in `lo..hi`, with under- and overflow.
    ///
    /// # Panics
    ///
    /// If `lo >= hi`.
    #[must_use]
    pub fn integer(lo: i64, hi: i64) -> Axis {
        assert!(lo < hi, "an integer axis needs lo < hi, got {lo}..{hi}");
        Self::with_kind(AxisKind::Integer { lo, hi }, true, true)
    }

    /// One bin per string category, plus an overflow bin for any other string.
    #[must_use]
    pub fn str_category<S: AsRef<str>>(categories: &[S]) -> Axis {
        let categories = categories.iter().map(|s| s.as_ref().to_owned()).collect();
        Self::with_kind(
            AxisKind::StrCategory {
                categories,
                growth: false,
            },
            false,
            true,
        )
    }

    /// One bin per integer category, plus an overflow bin for any other value.
    #[must_use]
    pub fn int_category(categories: &[i64]) -> Axis {
        Self::with_kind(
            AxisKind::IntCategory {
                categories: categories.to_vec(),
                growth: false,
            },
            false,
            true,
        )
    }

    /// Two bins, `false` and `true`, with no flow bins.
    #[must_use]
    pub fn boolean() -> Axis {
        Self::with_kind(AxisKind::Boolean, false, false)
    }

    /// Set the axis name.
    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Axis {
        self.name = name.into();
        self
    }

    /// Set the axis label.
    #[must_use]
    pub fn with_label(mut self, label: impl Into<String>) -> Axis {
        self.label = label.into();
        self
    }

    /// Keep (or drop) the underflow bin. Category and boolean axes never have
    /// one, so this is ignored for them.
    #[must_use]
    pub fn with_underflow(mut self, on: bool) -> Axis {
        if !self.is_categorical() {
            self.underflow = on;
        }
        self
    }

    /// Keep (or drop) the overflow bin. A growing category axis never has one,
    /// and a boolean axis neither, so this is ignored for them.
    #[must_use]
    pub fn with_overflow(mut self, on: bool) -> Axis {
        if !matches!(self.kind, AxisKind::Boolean) && !self.grows() {
            self.overflow = on;
        }
        self
    }

    /// Let a category axis grow: filling an unknown category adds a bin for it
    /// (and the axis then has no overflow bin). Ignored for other axes.
    #[must_use]
    pub fn with_growth(mut self) -> Axis {
        match &mut self.kind {
            AxisKind::StrCategory { growth, .. } | AxisKind::IntCategory { growth, .. } => {
                *growth = true;
                self.overflow = false;
            }
            _ => {}
        }
        self
    }

    /// The binning.
    #[must_use]
    pub fn kind(&self) -> &AxisKind {
        &self.kind
    }

    /// The number of in-range bins.
    #[must_use]
    pub fn size(&self) -> usize {
        match &self.kind {
            AxisKind::Regular { bins, .. } => *bins,
            AxisKind::Variable { edges } => edges.len() - 1,
            AxisKind::Integer { lo, hi } => (hi - lo) as usize,
            AxisKind::StrCategory { categories, .. } => categories.len(),
            AxisKind::IntCategory { categories, .. } => categories.len(),
            AxisKind::Boolean => 2,
        }
    }

    /// Whether the axis has an underflow bin.
    #[must_use]
    pub fn underflow(&self) -> bool {
        self.underflow
    }

    /// Whether the axis has an overflow bin.
    #[must_use]
    pub fn overflow(&self) -> bool {
        self.overflow
    }

    /// The number of bins including the flow bins it has.
    #[must_use]
    pub fn extent(&self) -> usize {
        self.size() + usize::from(self.underflow) + usize::from(self.overflow)
    }

    /// Whether this is a category (string, integer) or boolean axis.
    #[must_use]
    pub fn is_categorical(&self) -> bool {
        matches!(
            self.kind,
            AxisKind::StrCategory { .. } | AxisKind::IntCategory { .. } | AxisKind::Boolean
        )
    }

    fn grows(&self) -> bool {
        matches!(
            self.kind,
            AxisKind::StrCategory { growth: true, .. } | AxisKind::IntCategory { growth: true, .. }
        )
    }

    /// The `size + 1` bin edges of a numeric axis (an integer axis's bins are
    /// `[i, i + 1)`); `None` for a category or boolean axis.
    #[must_use]
    pub fn edges(&self) -> Option<Vec<f64>> {
        match &self.kind {
            AxisKind::Regular { bins, lo, hi } => Some(
                (0..=*bins)
                    .map(|i| lo + (hi - lo) * i as f64 / *bins as f64)
                    .collect(),
            ),
            AxisKind::Variable { edges } => Some(edges.clone()),
            AxisKind::Integer { lo, hi } => Some((*lo..=*hi).map(|i| i as f64).collect()),
            _ => None,
        }
    }

    /// The bin edges for drawing: [`edges`](Axis::edges) for a numeric axis,
    /// and `0, 1, …, size` (one unit per category) for a categorical one.
    #[must_use]
    pub fn plot_edges(&self) -> Vec<f64> {
        self.edges()
            .unwrap_or_else(|| (0..=self.size()).map(|i| i as f64).collect())
    }

    /// The bin centers, from [`plot_edges`](Axis::plot_edges).
    #[must_use]
    pub fn centers(&self) -> Vec<f64> {
        self.plot_edges()
            .windows(2)
            .map(|w| 0.5 * (w[0] + w[1]))
            .collect()
    }

    /// The bin widths, from [`plot_edges`](Axis::plot_edges).
    #[must_use]
    pub fn widths(&self) -> Vec<f64> {
        self.plot_edges().windows(2).map(|w| w[1] - w[0]).collect()
    }

    /// The label of in-range bin `i`: the category for a categorical axis, the
    /// interval `[lo, hi)` for a numeric one. `None` if `i` is out of range.
    #[must_use]
    pub fn bin_label(&self, i: usize) -> Option<String> {
        if i >= self.size() {
            return None;
        }
        Some(match &self.kind {
            AxisKind::StrCategory { categories, .. } => categories[i].clone(),
            AxisKind::IntCategory { categories, .. } => categories[i].to_string(),
            AxisKind::Boolean => (i == 1).to_string(),
            AxisKind::Integer { lo, .. } => (lo + i as i64).to_string(),
            _ => {
                let e = self.plot_edges();
                format!("[{}, {})", e[i], e[i + 1])
            }
        })
    }

    /// Whether two axes bin the same way (same kind, bins and flow bins); names
    /// and labels are not compared.
    #[must_use]
    pub fn same_binning(&self, other: &Axis) -> bool {
        self.kind == other.kind
            && self.underflow == other.underflow
            && self.overflow == other.overflow
    }

    /// The bin index (`-1..=size`) a value falls in. `axis` is the axis
    /// position, for the error message.
    pub(crate) fn locate(&self, v: Value<'_>, axis: usize) -> Result<Locate> {
        let n = self.size() as isize;
        let mismatch = |expected: &'static str| Error::ValueType {
            axis,
            expected,
            found: v.to_string(),
        };
        let idx = match &self.kind {
            AxisKind::Regular { bins, lo, hi } => {
                let x = numeric(v).ok_or_else(|| mismatch("a number"))?;
                // Compare against the edges themselves, not the scaled
                // coordinate: rounding can scale a value just below `hi` to
                // exactly `bins`, which must still be the last bin.
                if x.is_nan() || x >= *hi {
                    n
                } else if x < *lo {
                    -1
                } else {
                    let z = (x - lo) / (hi - lo) * *bins as f64;
                    (z.floor() as isize).clamp(0, n - 1)
                }
            }
            AxisKind::Variable { edges } => {
                let x = numeric(v).ok_or_else(|| mismatch("a number"))?;
                if x.is_nan() || x >= edges[edges.len() - 1] {
                    n
                } else if x < edges[0] {
                    -1
                } else {
                    edges.partition_point(|&e| e <= x) as isize - 1
                }
            }
            AxisKind::Integer { lo, hi } => {
                let i = match v {
                    Value::Int(i) => i,
                    Value::Num(x) if x.is_nan() => return Ok(Locate::Index(n)),
                    Value::Num(x) => x.floor() as i64,
                    Value::Bool(b) => i64::from(b),
                    Value::Str(_) => return Err(mismatch("an integer")),
                };
                if i < *lo {
                    -1
                } else if i >= *hi {
                    n
                } else {
                    (i - lo) as isize
                }
            }
            AxisKind::StrCategory { categories, growth } => {
                let Value::Str(s) = v else {
                    return Err(mismatch("a string"));
                };
                match categories.iter().position(|c| c == s) {
                    Some(i) => i as isize,
                    None if *growth => return Ok(Locate::Grow),
                    None => n,
                }
            }
            AxisKind::IntCategory { categories, growth } => {
                let i = match v {
                    Value::Int(i) => i,
                    Value::Num(x) if x.fract() == 0.0 && x.is_finite() => x as i64,
                    _ => return Err(mismatch("an integer")),
                };
                match categories.iter().position(|&c| c == i) {
                    Some(k) => k as isize,
                    None if *growth => return Ok(Locate::Grow),
                    None => n,
                }
            }
            AxisKind::Boolean => match v {
                Value::Bool(b) => isize::from(b),
                Value::Int(0) => 0,
                Value::Int(1) => 1,
                Value::Num(0.0) => 0,
                Value::Num(1.0) => 1,
                _ => return Err(mismatch("a boolean")),
            },
        };
        Ok(Locate::Index(idx))
    }

    /// Whether two category axes can be merged by label even though their
    /// categories differ: both grow, they are the same kind of category axis,
    /// and they keep the same flow bins.
    pub(crate) fn merges_by_label(&self, other: &Axis) -> bool {
        let same_kind = matches!(
            (&self.kind, &other.kind),
            (
                AxisKind::StrCategory { growth: true, .. },
                AxisKind::StrCategory { growth: true, .. }
            ) | (
                AxisKind::IntCategory { growth: true, .. },
                AxisKind::IntCategory { growth: true, .. }
            )
        );
        same_kind && self.underflow == other.underflow && self.overflow == other.overflow
    }

    /// The fill value naming category bin `i` of a category axis.
    pub(crate) fn category_value(&self, i: usize) -> Value<'_> {
        match &self.kind {
            AxisKind::StrCategory { categories, .. } => Value::Str(&categories[i]),
            AxisKind::IntCategory { categories, .. } => Value::Int(categories[i]),
            _ => unreachable!("only category axes have category values"),
        }
    }

    /// Add the category `v` (which [`locate`](Axis::locate) reported as
    /// [`Locate::Grow`]) as a new last bin.
    pub(crate) fn grow(&mut self, v: Value<'_>) {
        match (&mut self.kind, v) {
            (AxisKind::StrCategory { categories, .. }, Value::Str(s)) => {
                categories.push(s.to_owned());
            }
            (AxisKind::IntCategory { categories, .. }, Value::Int(i)) => categories.push(i),
            (AxisKind::IntCategory { categories, .. }, Value::Num(x)) => {
                categories.push(x as i64);
            }
            _ => unreachable!("only a growing category axis grows"),
        }
    }

    /// The storage position (`0..extent`) of bin index `i` (`-1..=size`), or
    /// `None` if `i` is a flow bin the axis does not have, or out of range.
    pub(crate) fn position(&self, i: isize) -> Option<usize> {
        let n = self.size() as isize;
        if i < -1 || i > n || (i == -1 && !self.underflow) || (i == n && !self.overflow) {
            return None;
        }
        Some((i + isize::from(self.underflow)) as usize)
    }

    /// The bin index (`-1..=size`) at storage position `p`.
    pub(crate) fn index_at(&self, p: usize) -> isize {
        p as isize - isize::from(self.underflow)
    }
}

fn numeric(v: Value<'_>) -> Option<f64> {
    match v {
        Value::Num(x) => Some(x),
        Value::Int(i) => Some(i as f64),
        Value::Bool(b) => Some(f64::from(u8::from(b))),
        Value::Str(_) => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn idx<'a>(a: &Axis, v: impl Into<Value<'a>>) -> isize {
        match a.locate(v.into(), 0).unwrap() {
            Locate::Index(i) => i,
            Locate::Grow => panic!("grow"),
        }
    }

    #[test]
    fn regular_is_half_open_with_flow() {
        let a = Axis::regular(4, 0.0, 4.0);
        assert_eq!(
            [idx(&a, -0.1), idx(&a, 0.0), idx(&a, 3.999), idx(&a, 4.0)],
            [-1, 0, 3, 4]
        );
        assert_eq!(idx(&a, f64::NAN), 4);
        assert_eq!(a.edges().unwrap(), vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        assert_eq!(a.extent(), 6);
    }

    #[test]
    fn variable_and_integer() {
        let v = Axis::variable(&[0.0, 1.0, 10.0]);
        assert_eq!([idx(&v, 0.5), idx(&v, 1.0), idx(&v, 10.0)], [0, 1, 2]);
        let i = Axis::integer(-2, 3);
        assert_eq!(
            [idx(&i, -3), idx(&i, -2), idx(&i, 2), idx(&i, 3)],
            [-1, 0, 4, 5]
        );
        assert_eq!(idx(&i, -1.5), 0); // floored
    }

    #[test]
    fn categories_and_growth() {
        let c = Axis::str_category(&["a", "b"]);
        assert_eq!([idx(&c, "b"), idx(&c, "zz")], [1, 2]);
        assert!(c.locate(Value::Num(1.0), 0).is_err());
        assert_eq!(c.position(-1), None); // no underflow
        let mut g = Axis::str_category(&["a"]).with_growth();
        assert!(!g.overflow());
        assert!(matches!(g.locate("b".into(), 0).unwrap(), Locate::Grow));
        g.grow("b".into());
        assert_eq!(idx(&g, "b"), 1);
        let b = Axis::boolean();
        assert_eq!([idx(&b, false), idx(&b, true)], [0, 1]);
    }

    #[test]
    fn invalid_axes_are_rejected() {
        assert!(Axis::try_regular(0, 0.0, 1.0).is_err());
        assert!(Axis::try_regular(3, 1.0, 1.0).is_err());
        assert!(Axis::try_variable(&[0.0, 0.0]).is_err());
    }
}
