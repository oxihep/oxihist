//! [`Histogram`]: `D` axes and one storage accumulator per bin.

use std::marker::PhantomData;

use crate::axis::{Axis, Locate, Value};
use crate::error::{Error, Result};
use crate::storage::{CountStorage, MeanStorage, Storage};

/// A `D`-dimensional histogram (`D` = 1, 2 or 3) with storage `S`.
///
/// Bins are stored with their flow bins, the first axis varying fastest. Build
/// one with the [`Hist`](crate::Hist) builder or [`Histogram::new`].
///
/// ```
/// use oxihist::{Hist, Weight};
///
/// let mut h = Hist::reg(10, 0.0, 10.0).label("x").weight();
/// h.fill(2.5);
/// h.fill_weight(2.7, 3.0);
/// assert_eq!(h.value([2]), Some(4.0));
/// assert_eq!(h.variance([2]), Some(10.0)); // 1² + 3²
/// ```
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "serde",
    serde(bound(
        serialize = "S::Bin: serde::Serialize",
        deserialize = "S::Bin: serde::Deserialize<'de>"
    ))
)]
pub struct Histogram<S: Storage, const D: usize> {
    axes: Vec<Axis>,
    /// Storage: `cap[i]` slots along axis `i`, the first axis varying fastest.
    /// `cap[i]` is the axis's extent, except for a growing category axis,
    /// which keeps spare (default) slots so a new category rarely moves bins.
    bins: Vec<S::Bin>,
    cap: Vec<usize>,
    /// A short identifier, e.g. `"pt"`. Writing to a ROOT file uses it as the
    /// key name.
    pub name: String,
    /// A human-readable title.
    pub title: String,
    #[cfg_attr(feature = "serde", serde(skip))]
    storage: PhantomData<S>,
}

/// Two histograms are equal when their axes, name, title and every bin are;
/// spare storage kept for growing axes does not count.
impl<S: Storage, const D: usize> PartialEq for Histogram<S, D> {
    fn eq(&self, other: &Self) -> bool {
        self.axes == other.axes
            && self.name == other.name
            && self.title == other.title
            && self
                .cells()
                .zip(other.cells())
                .all(|((_, a), (_, b))| self.bins[a] == other.bins[b])
    }
}

impl<S: Storage, const D: usize> Histogram<S, D> {
    /// An empty histogram over `axes`.
    #[must_use]
    pub fn new(axes: [Axis; D]) -> Self {
        Self::from_axes(axes.into())
    }

    pub(crate) fn from_axes(axes: Vec<Axis>) -> Self {
        assert_eq!(axes.len(), D, "a {D}-D histogram needs {D} axes");
        let cap: Vec<usize> = axes.iter().map(Axis::extent).collect();
        let n = cap.iter().product();
        Histogram {
            axes,
            bins: vec![S::Bin::default(); n],
            cap,
            name: String::new(),
            title: String::new(),
            storage: PhantomData,
        }
    }

    /// Build a histogram from its axes and every bin's accumulator, flow bins
    /// included and the first axis varying fastest (the layout of
    /// [`bins`](Histogram::bins)).
    ///
    /// # Errors
    ///
    /// [`Error::BinningMismatch`] if `bins` does not have one entry per bin.
    pub fn from_bins(axes: [Axis; D], bins: Vec<S::Bin>) -> Result<Self> {
        let mut h = Self::new(axes);
        if bins.len() != h.bins.len() {
            return Err(Error::BinningMismatch(format!(
                "{} bins given for axes with {} (flow bins included)",
                bins.len(),
                h.bins.len()
            )));
        }
        h.bins = bins;
        Ok(h)
    }

    /// Set the name.
    #[must_use]
    pub fn named(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Set the title.
    #[must_use]
    pub fn titled(mut self, title: impl Into<String>) -> Self {
        self.title = title.into();
        self
    }

    /// The axes.
    #[must_use]
    pub fn axes(&self) -> &[Axis] {
        &self.axes
    }

    /// Axis `i`.
    ///
    /// # Panics
    ///
    /// If `i >= D`.
    #[must_use]
    pub fn axis(&self, i: usize) -> &Axis {
        &self.axes[i]
    }

    /// Mutable access to axis `i`'s name and label. (The binning cannot be
    /// changed in place.)
    ///
    /// # Panics
    ///
    /// If `i >= D`.
    pub fn axis_meta_mut(&mut self, i: usize) -> (&mut String, &mut String) {
        let a = &mut self.axes[i];
        (&mut a.name, &mut a.label)
    }

    /// Every bin's accumulator, flow bins included, the first axis varying
    /// fastest; [`shape`](Histogram::shape)`(true)` gives the extents.
    #[must_use]
    pub fn bins(&self) -> Vec<S::Bin> {
        self.cells().map(|(_, o)| self.bins[o].clone()).collect()
    }

    /// The number of bins along each axis, with or without the flow bins.
    #[must_use]
    pub fn shape(&self, flow: bool) -> [usize; D] {
        std::array::from_fn(|i| {
            if flow {
                self.axes[i].extent()
            } else {
                self.axes[i].size()
            }
        })
    }

    #[inline]
    fn strides(&self) -> [usize; D] {
        let mut s = [1; D];
        for i in 1..D {
            s[i] = s[i - 1] * self.cap[i - 1];
        }
        s
    }

    /// The storage offset of bin `idx` (each `-1..=size`).
    #[inline]
    fn offset(&self, idx: [isize; D]) -> Option<usize> {
        let strides = self.strides();
        let mut off = 0;
        for i in 0..D {
            off += self.axes[i].position(idx[i])? * strides[i];
        }
        Some(off)
    }

    /// Every bin (flow bins included) as its index and storage offset, the
    /// first axis varying fastest. Spare storage slots are skipped.
    fn cells(&self) -> impl Iterator<Item = ([isize; D], usize)> + '_ {
        let ext = self.shape(true);
        let strides = self.strides();
        let under: [isize; D] = std::array::from_fn(|i| isize::from(self.axes[i].underflow()));
        let total: usize = ext.iter().product();
        (0..total).map(move |mut k| {
            let mut idx = [0isize; D];
            let mut off = 0;
            for i in 0..D {
                let p = k % ext[i];
                k /= ext[i];
                idx[i] = p as isize - under[i];
                off += p * strides[i];
            }
            (idx, off)
        })
    }

    /// Bin `idx`'s accumulator (`-1` is the underflow bin, `size` the
    /// overflow), or `None` if there is no such bin.
    #[must_use]
    pub fn bin(&self, idx: [isize; D]) -> Option<&S::Bin> {
        self.offset(idx).map(|o| &self.bins[o])
    }

    /// Mutable access to bin `idx`'s accumulator.
    pub fn bin_mut(&mut self, idx: [isize; D]) -> Option<&mut S::Bin> {
        self.offset(idx).map(|o| &mut self.bins[o])
    }

    /// Bin `idx`'s value (see [`Storage::value`]).
    #[must_use]
    pub fn value(&self, idx: [isize; D]) -> Option<f64> {
        self.bin(idx).map(S::value)
    }

    /// Bin `idx`'s variance (see [`Storage::variance`]).
    #[must_use]
    pub fn variance(&self, idx: [isize; D]) -> Option<f64> {
        self.bin(idx).map(S::variance)
    }

    /// Every bin's index and accumulator, the first axis varying fastest; with
    /// `flow`, the flow bins too.
    pub fn iter(&self, flow: bool) -> impl Iterator<Item = ([isize; D], &S::Bin)> + '_ {
        self.cells().filter_map(move |(idx, o)| {
            let b = &self.bins[o];
            let in_range = idx
                .iter()
                .zip(&self.axes)
                .all(|(&i, a)| i >= 0 && (i as usize) < a.size());
            (flow || in_range).then_some((idx, b))
        })
    }

    /// Every bin's value, the first axis varying fastest; with `flow`, the
    /// flow bins too (laid out as [`shape`](Histogram::shape)`(flow)`).
    #[must_use]
    pub fn values(&self, flow: bool) -> Vec<f64> {
        self.iter(flow).map(|(_, b)| S::value(b)).collect()
    }

    /// Every bin's variance, laid out like [`values`](Histogram::values).
    #[must_use]
    pub fn variances(&self, flow: bool) -> Vec<f64> {
        self.iter(flow).map(|(_, b)| S::variance(b)).collect()
    }

    /// All bins combined into one accumulator (the total), with or without the
    /// flow bins.
    #[must_use]
    pub fn sum(&self, flow: bool) -> S::Bin {
        let mut acc = S::Bin::default();
        for (_, b) in self.iter(flow) {
            S::merge(&mut acc, b);
        }
        acc
    }

    /// Empty every bin, keeping the axes.
    pub fn reset(&mut self) {
        self.bins.fill(S::Bin::default());
    }

    /// Multiply every bin by `c` (see [`Storage::scale`]).
    pub fn scale(&mut self, c: f64) {
        self.bins.iter_mut().for_each(|b| S::scale(b, c));
    }

    /// Whether `other` has the same binning on every axis.
    #[must_use]
    pub fn same_binning(&self, other: &Self) -> bool {
        self.axes
            .iter()
            .zip(&other.axes)
            .all(|(a, b)| a.same_binning(b))
    }

    /// Add `other`'s bins to this histogram's, as if this one had also been
    /// filled with everything `other` was.
    ///
    /// Growing category axes may have met their categories in a different
    /// order, or met different ones (two threads filling copies of one
    /// template, say). They are merged by label: categories `other` has and
    /// `self` lacks are added to `self`.
    ///
    /// # Errors
    ///
    /// [`Error::BinningMismatch`] (leaving `self` unchanged) if the binnings
    /// differ in any other way.
    pub fn add(&mut self, other: &Self) -> Result<()> {
        if self.same_binning(other) {
            if self.cap == other.cap {
                for (a, b) in self.bins.iter_mut().zip(&other.bins) {
                    S::merge(a, b);
                }
            } else {
                for (idx, o) in other.cells() {
                    let off = self.offset(idx).expect("same binning");
                    S::merge(&mut self.bins[off], &other.bins[o]);
                }
            }
            return Ok(());
        }
        let compatible = self
            .axes
            .iter()
            .zip(&other.axes)
            .all(|(a, b)| a.same_binning(b) || a.merges_by_label(b));
        if !compatible {
            return Err(Error::BinningMismatch("the histograms' axes differ".into()));
        }
        // For each axis, where each of `other`'s in-range bins lands in `self`,
        // growing `self` by the categories it lacks.
        let mut maps: Vec<Option<Vec<isize>>> = vec![None; D];
        for i in 0..D {
            if self.axes[i].same_binning(&other.axes[i]) {
                continue;
            }
            let mut map = Vec::with_capacity(other.axes[i].size());
            for k in 0..other.axes[i].size() {
                let v = other.axes[i].category_value(k);
                let to = match self.axes[i].locate(v, i)? {
                    Locate::Index(j) => j,
                    Locate::Grow => {
                        self.grow(i, v);
                        self.axes[i].size() as isize - 1
                    }
                };
                map.push(to);
            }
            maps[i] = Some(map);
        }
        for (mut idx, o) in other.cells() {
            let b = &other.bins[o];
            for (i, map) in maps.iter().enumerate() {
                if let Some(map) = map {
                    // Growing axes have no flow bins, so every index is in range.
                    idx[i] = map[idx[i] as usize];
                }
            }
            let off = self
                .offset(idx)
                .expect("every bin of `other` maps into `self`");
            S::merge(&mut self.bins[off], b);
        }
        Ok(())
    }

    /// Locate a fill on every axis, growing category axes as needed, and
    /// return the bin's offset (`None` when it falls in a flow bin the axis
    /// does not keep).
    #[inline(always)]
    fn locate(&mut self, values: [Value<'_>; D]) -> Result<Option<usize>> {
        // Fast path: numbers on regular or variable axes, which neither grow
        // nor fail. Every axis is checked before a flow bin drops the fill,
        // so the general path below sees exactly what it did before.
        let strides = self.strides();
        let mut off = 0;
        let mut dropped = false;
        for i in 0..D {
            let fast = match values[i] {
                Value::Num(x) => self.axes[i].fast_position(x),
                _ => None,
            };
            match fast {
                Some(Some(p)) => off += p * strides[i],
                Some(None) => dropped = true,
                None => return self.locate_general(values),
            }
        }
        Ok((!dropped).then_some(off))
    }

    /// [`locate`](Histogram::locate) for every axis kind and value.
    #[inline(never)]
    fn locate_general(&mut self, values: [Value<'_>; D]) -> Result<Option<usize>> {
        let mut idx = [0isize; D];
        for i in 0..D {
            idx[i] = match self.axes[i].locate(values[i], i)? {
                Locate::Index(k) => k,
                Locate::Grow => {
                    self.grow(i, values[i]);
                    self.axes[i].size() as isize - 1
                }
            };
        }
        Ok(self.offset(idx))
    }

    /// Add a category to axis `i` (a growing category axis, which has no flow
    /// bins). Spare capacity is used when there is some; otherwise the axis's
    /// capacity doubles and every bin moves once, so growing to `n` categories
    /// costs `O(n)` moves per bin in all, not `O(n)` per category.
    #[cold]
    #[inline(never)]
    fn grow(&mut self, i: usize, v: Value<'_>) {
        let new_size = self.axes[i].size() + 1;
        if new_size > self.cap[i] {
            let mut cap = self.cap.clone();
            cap[i] = (self.cap[i] * 2).max(new_size).max(4);
            let n = cap.iter().product();
            let mut bins = vec![S::Bin::default(); n];
            let mut strides = [1usize; D];
            for k in 1..D {
                strides[k] = strides[k - 1] * cap[k - 1];
            }
            let moves: Vec<([isize; D], usize)> = self.cells().collect();
            let mut old = std::mem::take(&mut self.bins);
            for (idx, o) in moves {
                let mut off = 0;
                for k in 0..D {
                    off += self.axes[k].position(idx[k]).expect("a stored bin") * strides[k];
                }
                bins[off] = std::mem::take(&mut old[o]);
            }
            self.bins = bins;
            self.cap = cap;
        }
        self.axes[i].grow(v);
    }

    /// Apply `f` to the bin `values` fall in.
    #[inline(always)]
    fn apply(&mut self, values: [Value<'_>; D], f: impl FnOnce(&mut S::Bin)) -> Result<()> {
        if let Some(o) = self.locate(values)? {
            f(&mut self.bins[o]);
        }
        Ok(())
    }

    /// Sum over every axis not in `keep`, giving an `E`-dimensional histogram
    /// whose axes are `keep`'s, in that order. Flow bins are summed too.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidAxis`] if `keep` names an axis twice or one that does
    /// not exist.
    pub fn project<const E: usize>(&self, keep: [usize; E]) -> Result<Histogram<S, E>> {
        for (k, &a) in keep.iter().enumerate() {
            if a >= D || keep[..k].contains(&a) {
                return Err(Error::InvalidAxis(format!(
                    "cannot project a {D}-D histogram onto axes {keep:?}"
                )));
            }
        }
        let mut out =
            Histogram::<S, E>::from_axes(keep.iter().map(|&a| self.axes[a].clone()).collect());
        out.name.clone_from(&self.name);
        out.title.clone_from(&self.title);
        for (idx, o) in self.cells() {
            let b = &self.bins[o];
            let target: [isize; E] = std::array::from_fn(|k| idx[keep[k]]);
            let off = out.offset(target).expect("same axes, same bins");
            S::merge(&mut out.bins[off], b);
        }
        Ok(out)
    }

    /// Merge every `factor` neighbouring bins of axis `axis` into one. Flow
    /// bins are kept as they are.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidAxis`] for a categorical axis, a factor of 0, or one
    /// that does not divide the number of bins.
    pub fn rebin(&self, axis: usize, factor: usize) -> Result<Self> {
        let a = self
            .axes
            .get(axis)
            .ok_or_else(|| Error::InvalidAxis(format!("no axis {axis}")))?;
        let n = a.size();
        let edges = a.edges().ok_or_else(|| {
            Error::InvalidAxis(format!("axis {axis} is categorical and cannot be rebinned"))
        })?;
        if factor == 0 || n % factor != 0 {
            return Err(Error::InvalidAxis(format!(
                "a rebin factor of {factor} does not divide axis {axis}'s {n} bins"
            )));
        }
        let new_edges: Vec<f64> = edges.iter().step_by(factor).copied().collect();
        let mut new_axis = match a.kind() {
            crate::AxisKind::Regular { lo, hi, .. } => Axis::regular(n / factor, *lo, *hi),
            _ => Axis::variable(&new_edges),
        }
        .with_underflow(a.underflow())
        .with_overflow(a.overflow());
        new_axis.name.clone_from(&a.name);
        new_axis.label.clone_from(&a.label);

        let mut axes = self.axes.clone();
        axes[axis] = new_axis;
        let mut out = Self::from_axes(axes);
        out.name.clone_from(&self.name);
        out.title.clone_from(&self.title);
        let n = n as isize;
        for (mut idx, o) in self.cells() {
            let b = &self.bins[o];
            let i = idx[axis];
            idx[axis] = if (0..n).contains(&i) {
                i / factor as isize
            } else if i < 0 {
                -1
            } else {
                n / factor as isize
            };
            let off = out.offset(idx).expect("rebinning keeps the flow bins");
            S::merge(&mut out.bins[off], b);
        }
        Ok(out)
    }
}

impl<S: CountStorage> Histogram<S, 1> {
    /// Fill `x` with weight 1.
    ///
    /// # Panics
    ///
    /// If `x` does not fit the axis (a string on a numeric axis, …); see
    /// [`try_fill`](Histogram::try_fill).
    #[inline(always)]
    pub fn fill<'a>(&mut self, x: impl Into<Value<'a>>) {
        self.fill_weight(x, 1.0);
    }

    /// Fill `x` with weight `w`.
    ///
    /// # Panics
    ///
    /// If `x` does not fit the axis.
    #[inline(always)]
    pub fn fill_weight<'a>(&mut self, x: impl Into<Value<'a>>, w: f64) {
        self.try_fill_weight(x, w)
            .unwrap_or_else(|e| fill_failed(e));
    }

    /// Fill `x` with weight 1, or report why it does not fit the axis.
    #[inline(always)]
    pub fn try_fill<'a>(&mut self, x: impl Into<Value<'a>>) -> Result<()> {
        self.try_fill_weight(x, 1.0)
    }

    /// Fill `x` with weight `w`, or report why it does not fit the axis.
    #[inline(always)]
    pub fn try_fill_weight<'a>(&mut self, x: impl Into<Value<'a>>, w: f64) -> Result<()> {
        self.apply([x.into()], |b| S::add(b, w))
    }

    /// Fill every value of `xs` with weight 1.
    ///
    /// # Panics
    ///
    /// If a value does not fit the axis.
    pub fn fill_many<'a, V: Into<Value<'a>>>(&mut self, xs: impl IntoIterator<Item = V>) {
        for x in xs {
            self.fill(x);
        }
    }
}

impl<S: CountStorage> Histogram<S, 2> {
    /// Fill `(x, y)` with weight 1.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis; see [`try_fill`](Histogram::try_fill).
    #[inline(always)]
    pub fn fill<'a>(&mut self, x: impl Into<Value<'a>>, y: impl Into<Value<'a>>) {
        self.fill_weight(x, y, 1.0);
    }

    /// Fill `(x, y)` with weight `w`.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    #[inline(always)]
    pub fn fill_weight<'a>(&mut self, x: impl Into<Value<'a>>, y: impl Into<Value<'a>>, w: f64) {
        self.try_fill_weight(x, y, w)
            .unwrap_or_else(|e| fill_failed(e));
    }

    /// Fill `(x, y)` with weight 1, or report why a value does not fit.
    #[inline(always)]
    pub fn try_fill<'a>(&mut self, x: impl Into<Value<'a>>, y: impl Into<Value<'a>>) -> Result<()> {
        self.try_fill_weight(x, y, 1.0)
    }

    /// Fill `(x, y)` with weight `w`, or report why a value does not fit.
    #[inline(always)]
    pub fn try_fill_weight<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        w: f64,
    ) -> Result<()> {
        self.apply([x.into(), y.into()], |b| S::add(b, w))
    }
}

impl<S: CountStorage> Histogram<S, 3> {
    /// Fill `(x, y, z)` with weight 1.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis; see [`try_fill`](Histogram::try_fill).
    #[inline(always)]
    pub fn fill<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        z: impl Into<Value<'a>>,
    ) {
        self.fill_weight(x, y, z, 1.0);
    }

    /// Fill `(x, y, z)` with weight `w`.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    #[inline(always)]
    pub fn fill_weight<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        z: impl Into<Value<'a>>,
        w: f64,
    ) {
        self.try_fill_weight(x, y, z, w)
            .unwrap_or_else(|e| fill_failed(e));
    }

    /// Fill `(x, y, z)` with weight 1, or report why a value does not fit.
    #[inline(always)]
    pub fn try_fill<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        z: impl Into<Value<'a>>,
    ) -> Result<()> {
        self.try_fill_weight(x, y, z, 1.0)
    }

    /// Fill `(x, y, z)` with weight `w`, or report why a value does not fit.
    #[inline(always)]
    pub fn try_fill_weight<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        z: impl Into<Value<'a>>,
        w: f64,
    ) -> Result<()> {
        self.apply([x.into(), y.into(), z.into()], |b| S::add(b, w))
    }
}

impl<S: MeanStorage> Histogram<S, 1> {
    /// Add the sample `s` to the bin `x` falls in (a profile fill).
    ///
    /// # Panics
    ///
    /// If `x` does not fit the axis; see
    /// [`try_fill_sample`](Histogram::try_fill_sample).
    pub fn fill_sample<'a>(&mut self, x: impl Into<Value<'a>>, s: f64) {
        self.fill_sample_weight(x, s, 1.0);
    }

    /// Add the sample `s` with weight `w` to the bin `x` falls in.
    ///
    /// # Panics
    ///
    /// If `x` does not fit the axis.
    pub fn fill_sample_weight<'a>(&mut self, x: impl Into<Value<'a>>, s: f64, w: f64) {
        self.try_fill_sample_weight(x, s, w)
            .unwrap_or_else(|e| fill_failed(e));
    }

    /// Add the sample `s`, or report why `x` does not fit the axis.
    pub fn try_fill_sample<'a>(&mut self, x: impl Into<Value<'a>>, s: f64) -> Result<()> {
        self.try_fill_sample_weight(x, s, 1.0)
    }

    /// Add the sample `s` with weight `w`, or report why `x` does not fit.
    pub fn try_fill_sample_weight<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        s: f64,
        w: f64,
    ) -> Result<()> {
        self.apply([x.into()], |b| S::add(b, s, w))
    }
}

impl<S: MeanStorage> Histogram<S, 2> {
    /// Add the sample `s` to the bin `(x, y)` falls in.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    pub fn fill_sample<'a>(&mut self, x: impl Into<Value<'a>>, y: impl Into<Value<'a>>, s: f64) {
        self.fill_sample_weight(x, y, s, 1.0);
    }

    /// Add the sample `s` with weight `w` to the bin `(x, y)` falls in.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    pub fn fill_sample_weight<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        s: f64,
        w: f64,
    ) {
        self.try_fill_sample_weight(x, y, s, w)
            .unwrap_or_else(|e| fill_failed(e));
    }

    /// Add the sample `s` with weight `w`, or report why a value does not fit.
    pub fn try_fill_sample_weight<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        s: f64,
        w: f64,
    ) -> Result<()> {
        self.apply([x.into(), y.into()], |b| S::add(b, s, w))
    }
}

impl<S: MeanStorage> Histogram<S, 3> {
    /// Add the sample `s` to the bin `(x, y, z)` falls in.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    pub fn fill_sample<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        z: impl Into<Value<'a>>,
        s: f64,
    ) {
        self.fill_sample_weight(x, y, z, s, 1.0);
    }

    /// Add the sample `s` with weight `w` to the bin `(x, y, z)` falls in.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    pub fn fill_sample_weight<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        z: impl Into<Value<'a>>,
        s: f64,
        w: f64,
    ) {
        self.try_fill_sample_weight(x, y, z, s, w)
            .unwrap_or_else(|e| fill_failed(e));
    }

    /// Add the sample `s` with weight `w`, or report why a value does not fit.
    pub fn try_fill_sample_weight<'a>(
        &mut self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        z: impl Into<Value<'a>>,
        s: f64,
        w: f64,
    ) -> Result<()> {
        self.apply([x.into(), y.into(), z.into()], |b| S::add(b, s, w))
    }
}

impl<S: Storage> Histogram<S, 2> {
    /// The in-range values as a grid indexed `[x][y]`.
    #[must_use]
    pub fn values_grid(&self) -> Vec<Vec<f64>> {
        let (nx, ny) = (self.axes[0].size(), self.axes[1].size());
        (0..nx as isize)
            .map(|ix| {
                (0..ny as isize)
                    .map(|iy| self.value([ix, iy]).unwrap_or(0.0))
                    .collect()
            })
            .collect()
    }
}

/// The panic of a `fill` whose value does not fit an axis, kept out of the
/// fill so the fast path stays small.
#[cold]
#[inline(never)]
fn fill_failed(e: Error) -> ! {
    panic!("{e}")
}

#[cfg(test)]
mod tests {
    use crate::{Axis, Double, Hist, Histogram, Int64, Mean, Weight};

    #[test]
    fn layout_is_first_axis_fastest_with_flow() {
        let mut h = Hist::reg(2, 0.0, 2.0).reg(3, 0.0, 3.0).double();
        h.fill(1.5, 0.5);
        assert_eq!(h.shape(true), [4, 5]);
        assert_eq!(h.value([1, 0]), Some(1.0));
        // Offset of (1, 0) is position (2, 1): 2 + 4 * 1.
        assert_eq!(h.bins()[6], 1.0);
        assert_eq!(h.values(false), vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn flow_bins_and_dropped_values() {
        let mut h = Histogram::<Int64, 1>::new([Axis::regular(2, 0.0, 2.0).with_underflow(false)]);
        h.fill(-1.0); // no underflow bin: dropped
        h.fill(5.0);
        assert_eq!(h.value([2]), Some(1.0));
        assert_eq!(h.value([-1]), None);
        assert_eq!(h.sum(true), 1);
    }

    #[test]
    fn fast_path_matches_the_general_path() {
        // A dropped numeric value on one axis still grows a category axis
        // after it, as when every axis went through the general path.
        let mut h = Histogram::<Double, 2>::new([
            Axis::regular(2, 0.0, 2.0).with_underflow(false),
            Axis::str_category(&["a"]).with_growth(),
        ]);
        h.fill(-1.0, "b");
        assert_eq!(h.axis(1).size(), 2);
        assert_eq!(h.sum(true), 0.0);

        // Numbers on a regular axis with no flow bins drop out of range;
        // `no_flow` applies to the first axis only, so the second keeps its.
        let mut h = Hist::reg(2, 0.0, 2.0).no_flow().reg(2, 0.0, 2.0).double();
        h.fill(9.0, 0.5); // dropped: x has no overflow bin
        h.fill(0.5, 9.0); // y overflow
        h.fill(0.5, 1.5);
        assert_eq!(h.sum(true), 2.0);
        assert_eq!(h.value([0, 2]), Some(1.0));
        assert_eq!(h.value([0, 1]), Some(1.0));

        // Variable axis: edges and flow bins, like the general path.
        let mut h = Hist::var(&[0.0, 1.0, 5.0]).double();
        for x in [-1.0, 0.0, 1.0, 4.999, 5.0, f64::NAN] {
            h.fill(x);
        }
        assert_eq!(h.values(true), vec![1.0, 1.0, 2.0, 2.0]);

        // A value of the wrong type still fails instead of taking a bin.
        let mut h = Hist::reg(2, 0.0, 2.0).double();
        assert!(h.try_fill("x").is_err());
        assert_eq!(h.sum(true), 0.0);
    }

    #[test]
    fn growth_keeps_existing_bins() {
        let mut h = Histogram::<Double, 2>::new([
            Axis::str_category(&["a"]).with_growth(),
            Axis::regular(2, 0.0, 2.0),
        ]);
        h.fill("a", 0.5);
        h.fill("b", 1.5);
        h.fill("a", 1.5);
        assert_eq!(h.axis(0).size(), 2);
        assert_eq!(h.values_grid(), vec![vec![1.0, 1.0], vec![0.0, 1.0]]);
    }

    #[test]
    fn project_and_rebin() {
        let mut h = Hist::reg(4, 0.0, 4.0).reg(2, 0.0, 2.0).weight();
        for (x, y) in [(0.5, 0.5), (1.5, 1.5), (3.5, 0.5), (9.0, 0.5)] {
            h.fill(x, y);
        }
        let px = h.project([0]).unwrap();
        assert_eq!(px.values(true), vec![0.0, 1.0, 1.0, 0.0, 1.0, 1.0]);
        let py = h.project([1]).unwrap();
        assert_eq!(py.values(false), vec![3.0, 1.0]);
        let r = h.rebin(0, 2).unwrap();
        assert_eq!(r.axis(0).size(), 2);
        assert_eq!(
            r.project([0]).unwrap().values(true),
            vec![0.0, 2.0, 1.0, 1.0]
        );
        assert!(h.rebin(0, 3).is_err());
        assert!(h.project([0, 0]).is_err());
    }

    #[test]
    fn add_requires_same_binning() {
        let mut a = Hist::reg(2, 0.0, 2.0).weight();
        let b = Hist::reg(3, 0.0, 2.0).weight();
        assert!(a.add(&b).is_err());
        let mut c = Hist::reg(2, 0.0, 2.0).weight();
        c.fill(0.5);
        a.add(&c).unwrap();
        a.add(&c).unwrap();
        assert_eq!(a.value([0]), Some(2.0));
        a.scale(0.5);
        assert_eq!((a.value([0]), a.variance([0])), (Some(1.0), Some(0.5)));
        let _: &Histogram<Weight, 1> = &a;
    }

    #[test]
    fn profile_fills() {
        let mut p = Hist::reg(2, 0.0, 2.0).mean();
        p.fill_sample(0.5, 1.0);
        p.fill_sample(0.5, 3.0);
        assert_eq!(p.value([0]), Some(2.0));
        let _: &Histogram<Mean, 1> = &p;
    }

    #[test]
    #[should_panic(expected = "axis 0 takes a number")]
    fn wrong_value_type_panics_in_fill() {
        let mut h = Hist::reg(2, 0.0, 2.0).double();
        h.fill("nope");
    }
}
