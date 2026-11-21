//! Histogram axes: how a fill value maps to a bin.
//!
//! Every axis has `size()` in-range bins, numbered `0..size`. Values below the
//! range land in the **underflow** bin (index `-1`) and values above it in the
//! **overflow** bin (index `size`), when the axis has them; otherwise they are
//! dropped. This is boost-histogram's convention.

use crate::cast;
use std::fmt;

use crate::error::{Error, Result};

/// A value to fill along one axis. Numbers, integers, strings and booleans all
/// convert into it, so `h.fill(3.2)`, `h.fill(7)`, `h.fill("ee")` and
/// `h.fill(true)` all work on an axis of the matching kind.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
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
                Value::$v(<$c>::from(x))
            }
        }
    )*};
}
value_from!(
    f64 => Num as f64, f32 => Num as f64,
    i64 => Int as i64, i32 => Int as i64, i16 => Int as i64, i8 => Int as i64,
    u32 => Int as i64, u16 => Int as i64, u8 => Int as i64,
);

impl From<usize> for Value<'_> {
    fn from(x: usize) -> Self {
        Value::Int(i64::try_from(x).unwrap_or(i64::MAX))
    }
}

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
#[non_exhaustive]
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
///
/// With the `serde` feature, a deserialised axis goes through the same checks
/// as the constructors (finite ascending edges, `lo < hi`, …), so malformed
/// input is an error rather than an axis that misbehaves later.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(try_from = "AxisRepr"))]
pub struct Axis {
    kind: AxisKind,
    /// A short identifier, e.g. `"pt"`.
    pub name: String,
    /// A human-readable label, e.g. `"$p_T$ [GeV]"`; plotting uses it.
    pub label: String,
    underflow: bool,
    overflow: bool,
    /// A hash index of the categories, kept once a category axis has more than
    /// a few: filling looks a category up on every call, and a linear scan of
    /// thousands of strings dominated the fill. Derived data, so it is not
    /// compared or serialized; deserialising rebuilds it.
    #[cfg_attr(feature = "serde", serde(skip))]
    lookup: CategoryLookup,
}

/// Whether `lo..hi` is a valid integer axis.
fn check_integer(lo: i64, hi: i64) -> Result<()> {
    // The bin count, and the count plus the two flow bins, must fit an index.
    let fits = hi
        .checked_sub(lo)
        .is_some_and(|n| n > 0 && n <= isize::MAX as i64 - 2);
    if fits {
        Ok(())
    } else {
        Err(Error::InvalidAxis(format!(
            "an integer axis needs lo < hi with hi - lo countable, got {lo}..{hi}"
        )))
    }
}

/// The serialised form of an [`Axis`]: the same fields, checked on the way in.
#[cfg(feature = "serde")]
#[derive(serde::Deserialize)]
struct AxisRepr {
    kind: AxisKind,
    name: String,
    label: String,
    underflow: bool,
    overflow: bool,
}

#[cfg(feature = "serde")]
impl TryFrom<AxisRepr> for Axis {
    type Error = Error;

    /// Accept exactly what the constructors accept.
    fn try_from(r: AxisRepr) -> Result<Axis> {
        match &r.kind {
            AxisKind::Regular { bins, lo, hi } => {
                Axis::try_regular(*bins, *lo, *hi)?;
            }
            AxisKind::Variable { edges } => {
                Axis::try_variable(edges)?;
            }
            AxisKind::Integer { lo, hi } => check_integer(*lo, *hi)?,
            AxisKind::StrCategory { .. } | AxisKind::IntCategory { .. } | AxisKind::Boolean => {}
        }
        // `with_kind` rebuilds the category index, which is not serialised.
        let mut axis = Axis::with_kind(r.kind, r.underflow, r.overflow);
        axis.name = r.name;
        axis.label = r.label;
        Ok(axis)
    }
}

/// Categories above which a category axis keeps a hash index. With foldhash
/// a lookup beats comparing strings one by one from about this many.
const LOOKUP_MIN: usize = 4;

/// The index's maps hash with foldhash: much faster than std's SipHash on
/// short keys, and still seeded per process, so crafted categories cannot
/// easily force collisions.
#[derive(Debug, Clone, Default)]
struct CategoryLookup {
    strs: foldhash::HashMap<String, usize>,
    ints: foldhash::HashMap<i64, usize>,
}

/// The index is derived from the categories, so it never makes two axes differ.
impl PartialEq for CategoryLookup {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl Axis {
    /// Rebuild the hash index from the categories, if there are enough.
    fn index_categories(&mut self) {
        self.lookup = CategoryLookup::default();
        match &self.kind {
            AxisKind::StrCategory { categories, .. } if categories.len() > LOOKUP_MIN => {
                for (i, c) in categories.iter().enumerate() {
                    self.lookup.strs.entry(c.clone()).or_insert(i);
                }
            }
            AxisKind::IntCategory { categories, .. } if categories.len() > LOOKUP_MIN => {
                for (i, &c) in categories.iter().enumerate() {
                    self.lookup.ints.entry(c).or_insert(i);
                }
            }
            _ => {}
        }
    }
}

impl Axis {
    /// The bin of string category `s`, through the index when it is built.
    #[inline]
    fn find_str(&self, categories: &[String], s: &str) -> Option<usize> {
        if self.lookup.strs.len() == categories.len() {
            self.lookup.strs.get(s).copied()
        } else {
            categories.iter().position(|c| c == s)
        }
    }

    /// The bin of integer category `i`, through the index when it is built.
    #[inline]
    fn find_int(&self, categories: &[i64], i: i64) -> Option<usize> {
        if self.lookup.ints.len() == categories.len() {
            self.lookup.ints.get(&i).copied()
        } else {
            categories.iter().position(|&c| c == i)
        }
    }
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
        let mut axis = Axis {
            kind,
            name: String::new(),
            label: String::new(),
            underflow,
            overflow,
            lookup: CategoryLookup::default(),
        };
        axis.index_categories();
        axis
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
    /// If `lo >= hi`, or the range has more integers than an index can count;
    /// see [`try_integer`](Axis::try_integer).
    #[must_use]
    pub fn integer(lo: i64, hi: i64) -> Axis {
        Self::try_integer(lo, hi).unwrap_or_else(|e| panic!("{e}"))
    }

    /// Like [`integer`](Axis::integer), returning an error instead of
    /// panicking. `hi - lo` must not overflow and must leave room for the two
    /// flow bins in an index (so `i64::MIN..i64::MAX` is refused). Each bin is
    /// stored, so a range is in practice bounded by memory long before that.
    pub fn try_integer(lo: i64, hi: i64) -> Result<Axis> {
        check_integer(lo, hi)?;
        Ok(Self::with_kind(AxisKind::Integer { lo, hi }, true, true))
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
    #[inline]
    pub fn size(&self) -> usize {
        match &self.kind {
            AxisKind::Regular { bins, .. } => *bins,
            AxisKind::Variable { edges } => edges.len() - 1,
            AxisKind::Integer { lo, hi } => cast::span(hi - lo),
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
            AxisKind::Integer { lo, .. } => (lo + i64::try_from(i).unwrap_or(i64::MAX)).to_string(),
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
    ///
    /// The common case, a number on a regular or variable axis, is handled
    /// here and inlined into the fill; every other case goes to
    /// [`locate_general`](Axis::locate_general), out of line, so the fill
    /// loop stays small.
    #[inline(always)]
    pub(crate) fn locate(&self, v: Value<'_>, axis: usize) -> Result<Locate> {
        if let Value::Num(x) = v {
            match &self.kind {
                AxisKind::Regular { bins, lo, hi } => {
                    return Ok(Locate::Index(regular_index(x, *bins, *lo, *hi)));
                }
                AxisKind::Variable { edges } => {
                    return Ok(Locate::Index(variable_index(x, edges)));
                }
                _ => {}
            }
        }
        self.locate_general(v, axis)
    }

    /// The storage position of the number `x` on a regular or variable axis
    /// (`None` for a flow bin the axis does not keep), with whether it is in
    /// range: one pass over the axis kind, so the fill does not match on it
    /// again to find out. `None` for any other axis: the general path
    /// handles those. The fill's fast path.
    #[inline(always)]
    pub(crate) fn fast_slot(&self, x: f64) -> Option<(Option<usize>, bool)> {
        let (i, n) = match &self.kind {
            AxisKind::Regular { bins, lo, hi } => (regular_index(x, *bins, *lo, *hi), *bins),
            AxisKind::Variable { edges } => (variable_index(x, edges), edges.len() - 1),
            _ => return None,
        };
        Some((self.flow_position(i, n), i >= 0 && cast::unsigned(i) < n))
    }

    /// The fast path for integer, category and boolean axes: the storage
    /// position of their native values (an integer or number on an integer
    /// axis, a string or integer on a category axis, a boolean on a boolean
    /// axis) and whether it is in range, like
    /// [`fast_slot`](Axis::fast_slot). `None` for anything the general path
    /// must handle (a new category on a growing axis, …).
    #[inline]
    pub(crate) fn other_slot(&self, v: Value<'_>) -> Option<(Option<usize>, bool)> {
        let (i, n): (isize, usize) = match (&self.kind, v) {
            (AxisKind::Integer { lo, hi }, Value::Int(k)) => {
                (integer_index(k, *lo, *hi), cast::span(hi - lo))
            }
            (AxisKind::Integer { lo, hi }, Value::Num(x)) => {
                let n = cast::span(hi - lo);
                if x.is_nan() {
                    (cast::signed(n), n)
                } else {
                    (integer_index(cast::floor_i64(x), *lo, *hi), n)
                }
            }
            (AxisKind::StrCategory { categories, growth }, Value::Str(s)) => {
                let n = categories.len();
                match self.find_str(categories, s) {
                    Some(k) => (cast::signed(k), n),
                    None if *growth => return None,
                    None => (cast::signed(n), n),
                }
            }
            (AxisKind::IntCategory { categories, growth }, Value::Int(k)) => {
                let n = categories.len();
                match self.find_int(categories, k) {
                    Some(j) => (cast::signed(j), n),
                    None if *growth => return None,
                    None => (cast::signed(n), n),
                }
            }
            (AxisKind::Boolean, Value::Bool(b)) => (isize::from(b), 2),
            _ => return None,
        };
        Some((self.flow_position(i, n), i >= 0 && cast::unsigned(i) < n))
    }

    /// The coordinate a fill of `v` landing in in-range bin `bin` adds to the
    /// histogram's [`FillStats`](crate::FillStats) on this axis: the value
    /// itself on a numeric axis, the center of the bin on a category or
    /// boolean axis, and nothing on a growing category axis, whose bins have
    /// no fixed position.
    pub(crate) fn coordinate(&self, v: Value<'_>, bin: usize) -> f64 {
        match &self.kind {
            AxisKind::Regular { .. } | AxisKind::Variable { .. } | AxisKind::Integer { .. } => {
                numeric(v).unwrap_or(0.0)
            }
            AxisKind::StrCategory { growth: true, .. }
            | AxisKind::IntCategory { growth: true, .. } => 0.0,
            _ => bin as f64 + 0.5,
        }
    }

    /// The storage position of bin `i` (`-1..=n`) on an axis of `n` bins,
    /// `None` for a flow bin the axis does not keep.
    #[inline(always)]
    fn flow_position(&self, i: isize, n: usize) -> Option<usize> {
        let under = usize::from(self.underflow);
        if i < 0 {
            self.underflow.then_some(0)
        } else if cast::unsigned(i) >= n {
            self.overflow.then_some(n + under)
        } else {
            Some(cast::unsigned(i) + under)
        }
    }

    /// [`locate`](Axis::locate) for every axis and value.
    #[inline(never)]
    fn locate_general(&self, v: Value<'_>, axis: usize) -> Result<Locate> {
        let n = cast::signed(self.size());
        let mismatch = |expected: &'static str| Error::ValueType {
            axis,
            expected,
            found: v.to_string(),
        };
        let idx = match &self.kind {
            AxisKind::Regular { bins, lo, hi } => {
                let x = numeric(v).ok_or_else(|| mismatch("a number"))?;
                regular_index(x, *bins, *lo, *hi)
            }
            AxisKind::Variable { edges } => {
                let x = numeric(v).ok_or_else(|| mismatch("a number"))?;
                variable_index(x, edges)
            }
            AxisKind::Integer { lo, hi } => {
                let i = match v {
                    Value::Int(i) => i,
                    Value::Num(x) if x.is_nan() => return Ok(Locate::Index(n)),
                    Value::Num(x) => cast::floor_i64(x),
                    Value::Bool(b) => i64::from(b),
                    Value::Str(_) => return Err(mismatch("an integer")),
                };
                if i < *lo {
                    -1
                } else if i >= *hi {
                    n
                } else {
                    cast::span_signed(i - lo)
                }
            }
            AxisKind::StrCategory { categories, growth } => {
                let Value::Str(s) = v else {
                    return Err(mismatch("a string"));
                };
                let found = self.find_str(categories, s);
                match found {
                    Some(i) => cast::signed(i),
                    None if *growth => return Ok(Locate::Grow),
                    None => n,
                }
            }
            AxisKind::IntCategory { categories, growth } => {
                let i = match v {
                    Value::Int(i) => i,
                    Value::Num(x) if x.fract() == 0.0 && x.is_finite() => cast::floor_i64(x),
                    _ => return Err(mismatch("an integer")),
                };
                let found = self.find_int(categories, i);
                match found {
                    Some(k) => cast::signed(k),
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
        let n = match (&mut self.kind, v) {
            (AxisKind::StrCategory { categories, .. }, Value::Str(s)) => {
                categories.push(s.to_owned());
                if self.lookup.strs.len() + 1 == categories.len() {
                    self.lookup
                        .strs
                        .entry(s.to_owned())
                        .or_insert(categories.len() - 1);
                }
                categories.len()
            }
            (AxisKind::IntCategory { categories, .. }, Value::Int(i)) => {
                categories.push(i);
                if self.lookup.ints.len() + 1 == categories.len() {
                    self.lookup.ints.entry(i).or_insert(categories.len() - 1);
                }
                categories.len()
            }
            (AxisKind::IntCategory { categories, .. }, Value::Num(x)) => {
                categories.push(cast::floor_i64(x));
                if self.lookup.ints.len() + 1 == categories.len() {
                    self.lookup
                        .ints
                        .entry(cast::floor_i64(x))
                        .or_insert(categories.len() - 1);
                }
                categories.len()
            }
            _ => unreachable!("only a growing category axis grows"),
        };
        // Crossing the threshold (or an axis read back without its index):
        // build the index from scratch.
        if n > LOOKUP_MIN && self.lookup.strs.len() + self.lookup.ints.len() != n {
            self.index_categories();
        }
    }

    /// The storage position (`0..extent`) of bin index `i` (`-1..=size`), or
    /// `None` if `i` is a flow bin the axis does not have, or out of range.
    #[inline]
    pub(crate) fn position(&self, i: isize) -> Option<usize> {
        let n = cast::signed(self.size());
        if i < -1 || i > n || (i == -1 && !self.underflow) || (i == n && !self.overflow) {
            return None;
        }
        Some(cast::unsigned(i + isize::from(self.underflow)))
    }
}

/// The bin (`-1..=bins`) of `x` on a regular axis over `[lo, hi)`.
#[inline(always)]
fn regular_index(x: f64, bins: usize, lo: f64, hi: f64) -> isize {
    let n = cast::signed(bins);
    // Compare against the edges themselves, not the scaled coordinate:
    // rounding can scale a value just below `hi` to exactly `bins`, which
    // must still be the last bin.
    if x.is_nan() || x >= hi {
        n
    } else if x < lo {
        -1
    } else {
        let z = (x - lo) / (hi - lo) * bins as f64;
        cast::floor_isize(z).clamp(0, n - 1)
    }
}

/// The bin (`-1..=hi - lo`) of integer `k` on an integer axis over `lo..hi`.
#[inline(always)]
fn integer_index(k: i64, lo: i64, hi: i64) -> isize {
    if k < lo {
        -1
    } else if k >= hi {
        cast::span_signed(hi - lo)
    } else {
        cast::span_signed(k - lo)
    }
}

/// The bin (`-1..=edges.len() - 1`) of `x` on a variable axis.
#[inline]
fn variable_index(x: f64, edges: &[f64]) -> isize {
    let n = cast::signed(edges.len()) - 1;
    if x.is_nan() || x >= edges[edges.len() - 1] {
        n
    } else if x < edges[0] {
        -1
    } else {
        cast::signed(edges.partition_point(|&e| e <= x)) - 1
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

    /// Wherever the fill's fast path answers, it agrees with the general path.
    #[test]
    fn fast_position_agrees_with_locate() {
        let many: Vec<String> = (0..40).map(|i| format!("c{i}")).collect(); // indexed
        let many_ints: Vec<i64> = (0..40).map(|i| i * 3).collect();
        let axes = [
            Axis::regular(4, -1.0, 3.0),
            Axis::regular(4, -1.0, 3.0)
                .with_underflow(false)
                .with_overflow(false),
            Axis::variable(&[0.0, 1.0, 1.5, 10.0]),
            Axis::integer(-2, 3),
            Axis::integer(-2, 3).with_overflow(false),
            Axis::str_category(&["a", "bb", "c"]),
            Axis::str_category(&["a", "bb", "c"]).with_overflow(false),
            Axis::str_category(&["a", "bb"]).with_growth(),
            Axis::str_category(&many),
            Axis::int_category(&[5, -1, 7]),
            Axis::int_category(&many_ints),
            Axis::int_category(&[5]).with_growth(),
            Axis::boolean(),
        ];
        let values = [
            Value::Num(-5.0),
            Value::Num(-1.0),
            Value::Num(0.0),
            Value::Num(1.25),
            Value::Num(2.9999999999999996),
            Value::Num(3.0),
            Value::Num(9.99),
            Value::Num(10.0),
            Value::Num(f64::NAN),
            Value::Num(f64::INFINITY),
            Value::Num(f64::NEG_INFINITY),
            Value::Int(-3),
            Value::Int(-2),
            Value::Int(0),
            Value::Int(2),
            Value::Int(3),
            Value::Int(5),
            Value::Int(7),
            Value::Int(57),
            Value::Int(200),
            Value::Str("a"),
            Value::Str("c"),
            Value::Str("zz"),
            Value::Str("c17"),
            Value::Str(""),
            Value::Bool(true),
            Value::Bool(false),
        ];
        let mut checked = 0;
        for axis in &axes {
            for &v in &values {
                // The fill's dispatch (see `Histogram::locate`).
                let fast = match v {
                    Value::Num(x) => axis.fast_slot(x).or_else(|| axis.other_slot(v)),
                    _ => axis.other_slot(v),
                };
                let Some((fast, inside)) = fast else {
                    continue;
                };
                checked += 1;
                match axis.locate_general(v, 0) {
                    Ok(Locate::Index(i)) => {
                        assert_eq!(fast, axis.position(i), "{:?} with {v:?}", axis.kind);
                        let in_range = i >= 0 && cast::unsigned(i) < axis.size();
                        assert_eq!(inside, in_range, "{:?} with {v:?}", axis.kind);
                    }
                    other => panic!(
                        "{:?} with {v:?}: fast path {fast:?}, general path {}",
                        axis.kind,
                        match other {
                            Ok(_) => "grows".to_owned(),
                            Err(e) => e.to_string(),
                        }
                    ),
                }
            }
        }
        assert!(checked > 60, "only {checked} fast-path cases");
    }

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
        assert!(Axis::try_integer(3, 3).is_err());
        assert!(Axis::try_integer(i64::MIN, i64::MAX).is_err());
        assert!(Axis::try_integer(i64::MIN, 0).is_err()); // 2^63 integers
        assert!(Axis::try_integer(i64::MIN, -3).is_ok()); // the widest allowed
        assert_eq!(Axis::try_integer(-2, 2).unwrap().size(), 4);
    }
}
