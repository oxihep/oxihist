//! [`Histogram`]: `D` axes and one storage accumulator per bin.

use crate::cast;
use std::marker::PhantomData;

use crate::axis::{Axis, Locate, Value};
use crate::error::{Error, Result};
use crate::fill_stats::{FillStats, Lanes};
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
///
/// Besides its bins, a histogram keeps [`FillStats`]: running sums over its
/// fills that make [`mean`](Histogram::mean), [`std_dev`](Histogram::std_dev)
/// and [`entries`](Histogram::entries) exact. A histogram that was not built
/// by filling (from [`from_bins`](Histogram::from_bins), after
/// [`bin_mut`](Histogram::bin_mut), [`multiply`](Histogram::multiply), …) has
/// none, and those are then estimated from the bins; see
/// [`fill_stats`](Histogram::fill_stats).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "serde",
    serde(bound(
        serialize = "S::Bin: serde::Serialize",
        deserialize = "S::Bin: serde::Deserialize<'de>"
    ))
)]
// Deserialising checks what the constructors guarantee: D axes, and a storage
// laid out for them.
#[cfg_attr(feature = "serde", serde(try_from = "HistogramRepr<S::Bin, D>"))]
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
    /// The exact statistics of the fills, or `None` when the bins did not
    /// all come from fills.
    stats: Option<Lanes<D>>,
    #[cfg_attr(feature = "serde", serde(skip))]
    storage: PhantomData<S>,
}

/// The serialised form of a [`Histogram`]: the same fields, checked on the
/// way in.
#[cfg(feature = "serde")]
#[derive(serde::Deserialize)]
struct HistogramRepr<B, const D: usize> {
    axes: Vec<Axis>,
    bins: Vec<B>,
    cap: Vec<usize>,
    name: String,
    title: String,
    /// Absent in data written before 0.14, which then has no fill
    /// statistics.
    #[serde(default)]
    stats: Option<FillStats<D>>,
}

#[cfg(feature = "serde")]
impl<S: Storage, const D: usize> TryFrom<HistogramRepr<S::Bin, D>> for Histogram<S, D> {
    type Error = crate::Error;

    fn try_from(r: HistogramRepr<S::Bin, D>) -> crate::Result<Self> {
        let bad = |m: String| Err(crate::Error::InvalidData(m));
        if r.axes.len() != D || r.cap.len() != D {
            return bad(format!(
                "a {D}-D histogram needs {D} axes and storage sizes, got {} and {}",
                r.axes.len(),
                r.cap.len()
            ));
        }
        // Each axis needs at least its extent in storage; a growing axis may
        // keep spare slots beyond it.
        if let Some((i, (a, &c))) = r
            .axes
            .iter()
            .zip(&r.cap)
            .enumerate()
            .find(|(_, (a, c))| **c < a.extent())
        {
            return bad(format!(
                "axis {i} has {} bins with flow but storage for {c}",
                a.extent()
            ));
        }
        let cells = r.cap.iter().try_fold(1usize, |n, &c| n.checked_mul(c));
        if cells != Some(r.bins.len()) {
            return bad(format!(
                "storage sized {:?} holds {} bins",
                r.cap,
                r.bins.len()
            ));
        }
        Ok(Histogram {
            axes: r.axes,
            bins: r.bins,
            cap: r.cap,
            name: r.name,
            title: r.title,
            stats: r.stats.map(Lanes::from),
            storage: PhantomData,
        })
    }
}

/// Two histograms are equal when their axes, name, title and every bin are;
/// spare storage kept for growing axes does not count, and neither do the
/// [fill statistics](Histogram::fill_stats): they are sums of the fill
/// values in the order they came, so two threads' histograms merged equal a
/// serial fill bin for bin but not sum for sum, and a histogram rebuilt from
/// its bins has none. Compare `fill_stats()` as well where they matter.
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
            stats: Some(Lanes::default()),
            storage: PhantomData,
        }
    }

    /// Build a histogram from its axes and every bin's accumulator, flow bins
    /// included and the first axis varying fastest (the layout of
    /// [`bins`](Histogram::bins)).
    ///
    /// The result has no [fill statistics](Histogram::fill_stats): its mean
    /// and entry count are estimated from the bins until
    /// [`set_fill_stats`](Histogram::set_fill_stats) supplies them.
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
        h.stats = None;
        Ok(h)
    }

    /// The exact statistics of the fills (see [`FillStats`]), or `None` when
    /// the bins did not all come from fills: after
    /// [`from_bins`](Histogram::from_bins), [`bin_mut`](Histogram::bin_mut),
    /// an operation that rewrites the bins
    /// ([`multiply`](Histogram::multiply), [`divide`](Histogram::divide),
    /// `smooth`, `slice`, `cumulative`), adding a histogram that has none, or
    /// reading data serialised before 0.14.
    /// [`binned_stats`](Histogram::binned_stats) estimates them in that case.
    #[must_use]
    pub fn fill_stats(&self) -> Option<FillStats<D>> {
        self.stats.as_ref().map(Lanes::total)
    }

    /// Set or clear the fill statistics: how a reader restores the sums a
    /// file stores (ROOT's `fEntries` and `fTsumw*`) after building the
    /// histogram with [`from_bins`](Histogram::from_bins). Nothing checks
    /// them against the bins. A histogram without them does not start
    /// keeping them when filled, so clearing them on an empty histogram
    /// spares its fills the bookkeeping.
    pub fn set_fill_stats(&mut self, stats: Option<FillStats<D>>) {
        self.stats = stats.map(Lanes::from);
    }

    /// The fill statistics estimated from the bins, each in-range bin
    /// counted at its center ([`Axis::centers`]) with the sums its
    /// accumulator holds ([`Storage::sums`]). `entries` is the sum of every
    /// bin's `Σw`, flow bins included: the number of fills when they all had
    /// weight 1.
    #[must_use]
    pub fn binned_stats(&self) -> FillStats<D> {
        let centers: [Vec<f64>; D] = std::array::from_fn(|i| self.axes[i].centers());
        let mut stats = FillStats::default();
        for (idx, o) in self.cells() {
            let b = S::sums(&self.bins[o]);
            stats.entries += b.sum_w;
            if self.in_range(idx) {
                let x: [f64; D] = std::array::from_fn(|i| centers[i][cast::unsigned(idx[i])]);
                stats.accumulate(b.sum_w, b.sum_w2, &x);
                stats.sum_ws += b.sum_ws;
                stats.sum_ws2 += b.sum_ws2;
            }
        }
        stats
    }

    /// The exact fill statistics when the histogram has them, else the
    /// estimate from its bins.
    pub(crate) fn stats(&self) -> FillStats<D> {
        self.fill_stats().unwrap_or_else(|| self.binned_stats())
    }

    /// The number of fills, in range or not. Without
    /// [fill statistics](Histogram::fill_stats) it is estimated as the sum
    /// of every bin's weights, flow bins included.
    #[doc(alias = "GetEntries")]
    #[must_use]
    pub fn entries(&self) -> f64 {
        self.stats().entries
    }

    /// The effective number of entries `(Σw)² / Σw²` over the in-range
    /// fills: the number of unit-weight fills with the same relative error
    /// (0 for an empty histogram).
    #[doc(alias = "GetEffectiveEntries")]
    #[must_use]
    pub fn effective_entries(&self) -> f64 {
        self.stats().effective_entries()
    }

    /// Whether bin `idx` is an in-range bin on every axis.
    fn in_range(&self, idx: [isize; D]) -> bool {
        idx.iter()
            .zip(&self.axes)
            .all(|(&i, a)| i >= 0 && cast::unsigned(i) < a.size())
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
                idx[i] = cast::signed(p) - under[i];
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

    /// Mutable access to bin `idx`'s accumulator. The histogram loses its
    /// [fill statistics](Histogram::fill_stats): a bin edited by hand is not
    /// a fill.
    pub fn bin_mut(&mut self, idx: [isize; D]) -> Option<&mut S::Bin> {
        let o = self.offset(idx)?;
        self.stats = None;
        Some(&mut self.bins[o])
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
            (flow || self.in_range(idx)).then_some((idx, &self.bins[o]))
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

    /// Empty every bin and zero the fill statistics, keeping the axes.
    pub fn reset(&mut self) {
        self.bins.fill(S::Bin::default());
        self.stats = Some(Lanes::default());
    }

    /// Multiply every bin by `c` (see [`Storage::scale`]). The fill
    /// statistics follow: a count histogram's `Σw` and moments by `c` and
    /// its `Σw²` by `c²`, a profile's sample sums `Σw·s` by `c` and `Σw·s²`
    /// by `c²`. The number of entries is unchanged.
    pub fn scale(&mut self, c: f64) {
        self.bins.iter_mut().for_each(|b| S::scale(b, c));
        if let Some(stats) = &mut self.stats {
            stats.scale(c, S::SAMPLED);
        }
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
    /// filled with everything `other` was. The fill statistics are summed
    /// too; if either histogram has none, so does the result.
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
            self.add_stats(other);
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
                        cast::signed(self.axes[i].size()) - 1
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
                    idx[i] = map[cast::unsigned(idx[i])];
                }
            }
            let off = self
                .offset(idx)
                .expect("every bin of `other` maps into `self`");
            S::merge(&mut self.bins[off], b);
        }
        self.add_stats(other);
        Ok(())
    }

    /// Sum `other`'s fill statistics into this histogram's; without them on
    /// both sides there are none.
    fn add_stats(&mut self, other: &Self) {
        self.stats = match (self.stats, &other.stats) {
            (Some(mut a), Some(b)) => {
                a.add(b);
                Some(a)
            }
            _ => None,
        };
    }

    /// Locate a fill on every axis, growing category axes as needed: see
    /// [`Hit`].
    #[inline(always)]
    fn locate(&mut self, values: [Value<'_>; D]) -> Result<Hit<D>> {
        // Fast path, inlined: numbers on regular or variable axes, which
        // neither grow nor fail. Every axis is checked before a flow bin drops the fill,
        // so the general path below sees exactly what it did before.
        let strides = self.strides();
        let mut off = 0;
        let mut dropped = false;
        let mut inside = true;
        let mut x = [0.0; D];
        for i in 0..D {
            let fast = match values[i] {
                Value::Num(v) => {
                    x[i] = v;
                    self.axes[i].fast_slot(v)
                }
                _ => None,
            };
            let Some((pos, here)) = fast else {
                return self.locate_general(values);
            };
            match pos {
                Some(p) => off += p * strides[i],
                None => dropped = true,
            }
            inside &= here;
        }
        Ok(Hit {
            off: (!dropped).then_some(off),
            inside,
            x,
        })
    }

    /// [`locate`](Histogram::locate) for every other axis kind and value.
    /// It first tries each axis's native values (integers, known categories,
    /// booleans, numbers), which neither grow nor fail; only a new category,
    /// a converted value or a mismatch reaches the fully general code.
    #[inline(never)]
    fn locate_general(&mut self, values: [Value<'_>; D]) -> Result<Hit<D>> {
        let strides = self.strides();
        let mut off = 0;
        let mut dropped = false;
        let mut inside = true;
        let mut x = [0.0; D];
        for i in 0..D {
            let axis = &self.axes[i];
            let fast = match values[i] {
                Value::Num(v) => axis.fast_slot(v).or_else(|| axis.other_slot(values[i])),
                v => axis.other_slot(v),
            };
            let Some((pos, here)) = fast else {
                return self.locate_slow(values);
            };
            match pos {
                Some(p) => {
                    off += p * strides[i];
                    if here {
                        x[i] = axis.coordinate(values[i], p - usize::from(axis.underflow()));
                    }
                }
                None => dropped = true,
            }
            inside &= here;
        }
        Ok(Hit {
            off: (!dropped).then_some(off),
            inside,
            x,
        })
    }

    /// [`locate`](Histogram::locate) for every axis kind and value, growing
    /// category axes as needed.
    #[inline(never)]
    fn locate_slow(&mut self, values: [Value<'_>; D]) -> Result<Hit<D>> {
        let mut idx = [0isize; D];
        for i in 0..D {
            idx[i] = match self.axes[i].locate(values[i], i)? {
                Locate::Index(k) => k,
                Locate::Grow => {
                    self.grow(i, values[i]);
                    cast::signed(self.axes[i].size()) - 1
                }
            };
        }
        let inside = self.in_range(idx);
        let mut x = [0.0; D];
        if inside {
            for i in 0..D {
                x[i] = self.axes[i].coordinate(values[i], cast::unsigned(idx[i]));
            }
        }
        Ok(Hit {
            off: self.offset(idx),
            inside,
            x,
        })
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

    /// Apply `f` to the bin `values` fall in, and record the fill (weight
    /// `w`, sampled value `s` for a profile) in the fill statistics.
    #[inline(always)]
    fn apply(
        &mut self,
        values: [Value<'_>; D],
        w: f64,
        s: f64,
        f: impl FnOnce(&mut S::Bin),
    ) -> Result<()> {
        let hit = self.locate(values)?;
        if let Some(o) = hit.off {
            f(&mut self.bins[o]);
        }
        if let Some(stats) = &mut self.stats {
            // The bin picks the lane: see `Lanes`.
            let lane = stats.lane(hit.off.unwrap_or(0));
            lane.entries += 1.0;
            // A fill out of range adds zeros, not a branch: its own values
            // may be NaN or infinite and must not reach the sums.
            let (w, s, x) = if hit.inside {
                (w, s, hit.x)
            } else {
                (0.0, 0.0, [0.0; D])
            };
            lane.accumulate(w, w * w, &x);
            // Only a profile has sample sums; a count's stay 0.
            if S::SAMPLED {
                lane.sum_ws += w * s;
                lane.sum_ws2 += w * s * s;
            }
        }
        Ok(())
    }

    /// Sum over every axis not in `keep`, giving an `E`-dimensional histogram
    /// whose axes are `keep`'s, in that order. Flow bins are summed too.
    ///
    /// The fill statistics keep the kept axes' sums. They still count only
    /// the fills that were in range on every axis of this histogram, the
    /// summed-over ones included, so they can differ from the projection's
    /// own in-range content.
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
        out.stats = self.stats.map(|s| s.project(keep));
        for (idx, o) in self.cells() {
            let b = &self.bins[o];
            let target: [isize; E] = std::array::from_fn(|k| idx[keep[k]]);
            let off = out.offset(target).expect("same axes, same bins");
            S::merge(&mut out.bins[off], b);
        }
        Ok(out)
    }

    /// Merge every `factor` neighbouring bins of axis `axis` into one. Flow
    /// bins and the fill statistics are kept as they are.
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
        out.stats = self.stats;
        let n = cast::signed(n);
        for (mut idx, o) in self.cells() {
            let b = &self.bins[o];
            let i = idx[axis];
            idx[axis] = if (0..n).contains(&i) {
                i / cast::signed(factor)
            } else if i < 0 {
                -1
            } else {
                n / cast::signed(factor)
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
        self.apply([x.into()], w, 0.0, |b| S::add(b, w))
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
        self.apply([x.into(), y.into()], w, 0.0, |b| S::add(b, w))
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
        self.apply([x.into(), y.into(), z.into()], w, 0.0, |b| S::add(b, w))
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
        self.apply([x.into()], w, s, |b| S::add(b, s, w))
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
        self.apply([x.into(), y.into()], w, s, |b| S::add(b, s, w))
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
        self.apply([x.into(), y.into(), z.into()], w, s, |b| S::add(b, s, w))
    }
}

impl<S: Storage> Histogram<S, 2> {
    /// The in-range values as a grid indexed `[x][y]`.
    #[must_use]
    pub fn values_grid(&self) -> Vec<Vec<f64>> {
        let (nx, ny) = (self.axes[0].size(), self.axes[1].size());
        (0..cast::signed(nx))
            .map(|ix| {
                (0..cast::signed(ny))
                    .map(|iy| self.value([ix, iy]).unwrap_or(0.0))
                    .collect()
            })
            .collect()
    }
}

/// Where a fill lands.
struct Hit<const D: usize> {
    /// The bin's storage offset; `None` when the fill falls in a flow bin an
    /// axis does not keep.
    off: Option<usize>,
    /// Whether the fill is in range on every axis.
    inside: bool,
    /// The coordinate the fill adds to the fill statistics on each axis (see
    /// [`FillStats`]); only meaningful when `inside`.
    x: [f64; D],
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

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * (1.0 + b.abs())
    }

    #[test]
    fn fill_stats_match_root_2d() {
        // ROOT 6.40: TH2D(4, 0, 4, 3, 0, 3) with these fills; GetEntries and
        // GetStats.
        let mut h = Hist::reg(4, 0.0, 4.0).reg(3, 0.0, 3.0).weight();
        h.fill_weight(0.3, 0.2, 2.0);
        h.fill_weight(1.7, 2.9, 0.5);
        h.fill(3.9, 1.1);
        h.fill_weight(9.0, 1.0, 3.0); // x overflow
        h.fill(2.2, -1.0); // y underflow
        let s = h.fill_stats().unwrap();
        assert_eq!(s.entries, 5.0);
        for (got, want) in [
            (s.sum_w, 3.5),
            (s.sum_w2, 5.25),
            (s.sum_wx[0], 5.35),
            (s.sum_wxy[0][0], 16.835),
            (s.sum_wx[1], 2.95),
            (s.sum_wxy[1][1], 5.495),
            (s.sum_wxy[0][1], 6.875),
        ] {
            assert!(close(got, want), "got {got}, want {want}");
        }
        assert_eq!(s.sum_wxy[0][1], s.sum_wxy[1][0]);
        assert_eq!((s.sum_ws, s.sum_ws2), (0.0, 0.0));
        assert_eq!(h.entries(), 5.0);
        assert!(close(h.effective_entries(), 3.5 * 3.5 / 5.25));
    }

    #[test]
    fn fill_stats_match_root_profile() {
        // ROOT 6.40: TProfile(4, 0, 4) with Fill(0.5, 10, 2), Fill(1.25, 5),
        // Fill(7, 3); GetStats.
        let mut p = Hist::reg(4, 0.0, 4.0).weighted_mean();
        p.fill_sample_weight(0.5, 10.0, 2.0);
        p.fill_sample(1.25, 5.0);
        p.fill_sample(7.0, 3.0); // overflow
        let s = p.fill_stats().unwrap();
        assert_eq!(
            (s.entries, s.sum_w, s.sum_w2, s.sum_wx, s.sum_wxy),
            (3.0, 3.0, 5.0, [2.25], [[2.0625]])
        );
        assert_eq!((s.sum_ws, s.sum_ws2), (25.0, 225.0));
        assert_eq!(p.mean(), 0.75);
        // Scaling a profile scales the sampled values, not the weights.
        p.scale(2.0);
        let t = p.fill_stats().unwrap();
        assert_eq!(
            (t.sum_w, t.sum_wx, t.sum_ws, t.sum_ws2),
            (3.0, [2.25], 50.0, 900.0)
        );
    }

    #[test]
    fn three_dimensional_cross_terms() {
        let mut h = Hist::reg(2, 0.0, 2.0)
            .reg(2, 0.0, 2.0)
            .reg(2, 0.0, 2.0)
            .double();
        h.fill_weight(0.5, 1.0, 1.5, 2.0);
        h.fill(0.5, 1.0, 9.0); // z overflow
        let s = h.fill_stats().unwrap();
        assert_eq!((s.entries, s.sum_w), (2.0, 2.0));
        assert_eq!(s.sum_wx, [1.0, 2.0, 3.0]);
        assert_eq!(
            s.sum_wxy,
            [[0.5, 1.0, 1.5], [1.0, 2.0, 3.0], [1.5, 3.0, 4.5]]
        );
        let mut p = Hist::reg(2, 0.0, 2.0)
            .reg(2, 0.0, 2.0)
            .reg(2, 0.0, 2.0)
            .mean();
        p.fill_sample(0.5, 1.0, 1.5, 4.0);
        assert_eq!(p.fill_stats().unwrap().sum_ws, 4.0);
        let mut q = Hist::reg(2, 0.0, 2.0).reg(2, 0.0, 2.0).mean();
        q.fill_sample_weight(0.5, 1.0, 4.0, 0.5);
        let s = q.fill_stats().unwrap();
        assert_eq!(
            (s.sum_w, s.sum_ws, s.sum_ws2, s.sum_wxy[0][1]),
            (0.5, 2.0, 8.0, 0.25)
        );
    }

    #[test]
    fn every_fill_counts_but_only_in_range_ones_weigh() {
        let mut h = Histogram::<Int64, 1>::new([Axis::regular(2, 0.0, 2.0).with_underflow(false)]);
        h.fill(-1.0); // dropped: no underflow bin
        h.fill(5.0); // overflow
        h.fill(f64::NAN); // overflow
        h.fill(0.25);
        assert!(h.try_fill("x").is_err()); // not a fill at all
        let s = h.fill_stats().unwrap();
        assert_eq!((s.entries, s.sum_w, s.sum_wx), (4.0, 1.0, [0.25]));
        // What an out-of-range fill carries never reaches the sums.
        let mut d = Hist::reg(2, 0.0, 2.0).double();
        d.fill(0.25);
        d.fill_weight(5.0, f64::INFINITY);
        d.fill_weight(f64::NAN, f64::NAN);
        let s = d.fill_stats().unwrap();
        assert_eq!(
            (s.entries, s.sum_w, s.sum_w2, s.sum_wx),
            (3.0, 1.0, 1.0, [0.25])
        );
        let mut p = Hist::reg(2, 0.0, 2.0).mean();
        p.fill_sample(0.5, 3.0);
        p.fill_sample(9.0, f64::NAN);
        let t = p.fill_stats().unwrap();
        assert_eq!(
            (t.entries, t.sum_w, t.sum_ws, t.sum_ws2),
            (2.0, 1.0, 3.0, 9.0)
        );
        // The statistics take the weight as given; an integer bin rounds it.
        h.fill_weight(1.5, 0.25);
        assert_eq!(h.fill_stats().unwrap().sum_w, 1.25);
        assert_eq!(h.value([1]), Some(0.0));
    }

    #[test]
    fn many_fills_in_many_bins_total_exactly() {
        // The sums are kept in lanes picked by the bin; the record is their
        // total, whichever bins the fills went to.
        let mut h = Hist::reg(16, 0.0, 16.0).reg(3, 0.0, 3.0).double();
        let (mut n, mut sx, mut sy, mut sxy) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..200u32 {
            let (x, y) = (f64::from(i % 17), f64::from(i % 4) * 0.75);
            h.fill(x, y);
            if x < 16.0 {
                n += 1.0;
                sx += x;
                sy += y;
                sxy += x * y;
            }
        }
        let s = h.fill_stats().unwrap();
        assert_eq!((s.entries, s.sum_w, s.sum_w2), (200.0, n, n));
        assert_eq!(
            (s.sum_wx, s.sum_wxy[0][1], s.sum_wxy[1][0]),
            ([sx, sy], sxy, sxy)
        );
        // Setting them back is the identity, however they were spread.
        let mut again = h.clone();
        again.set_fill_stats(h.fill_stats());
        assert_eq!(again.fill_stats(), h.fill_stats());
        assert!(format!("{h:?}").contains("stats: Some(FillStats {"));
    }

    #[test]
    fn coordinates_of_non_numeric_axes() {
        // An integer axis takes the value, also a fractional one.
        let mut i = Hist::int(-2, 3).double();
        i.fill(-2);
        i.fill(1.75);
        assert_eq!(i.fill_stats().unwrap().sum_wx, [-0.25]);
        // A category or boolean axis takes the bin's center.
        let mut c = Hist::str_cat(&["a", "b", "c"]).boolean().double();
        c.fill("c", true);
        c.fill("zz", true); // overflow
        let s = c.fill_stats().unwrap();
        assert_eq!((s.entries, s.sum_w, s.sum_wx), (2.0, 1.0, [2.5, 1.5]));
        let mut k = Hist::int_cat(&[7, 9]).no_flow().double();
        k.fill(9);
        assert_eq!(k.fill_stats().unwrap().sum_wx, [1.5]);
        // A growing category axis has no positions to sum.
        let mut g = Hist::str_cat::<&str>(&[])
            .growth()
            .reg(2, 0.0, 2.0)
            .double();
        g.fill("a", 0.5);
        g.fill("b", 1.5);
        let s = g.fill_stats().unwrap();
        assert_eq!((s.sum_w, s.sum_wx), (2.0, [0.0, 2.0]));
        assert_eq!(s.sum_wxy, [[0.0, 0.0], [0.0, 2.5]]);
    }

    #[test]
    fn statistics_are_exact_or_estimated_from_bins() {
        let mut h = Hist::reg(4, 0.0, 4.0).weight();
        h.fill(0.25);
        h.fill_weight(1.0, 3.0);
        h.fill(9.0);
        assert_eq!(h.mean(), 0.8125);
        assert_eq!(h.entries(), 3.0);

        // The same bins without the fills: the estimate from bin centers.
        let mut rebuilt = Histogram::<Weight, 1>::from_bins([h.axis(0).clone()], h.bins()).unwrap();
        assert_eq!(rebuilt, h);
        assert!(rebuilt.fill_stats().is_none());
        assert_eq!(rebuilt.mean(), 1.25);
        assert_eq!(rebuilt.entries(), 5.0); // Σw with the overflow
        let binned = h.binned_stats();
        assert_eq!(
            (binned.sum_w, binned.sum_w2, binned.sum_wx),
            (4.0, 10.0, [5.0])
        );
        assert_eq!(rebuilt.effective_entries(), 1.6);

        // A reader restores what the file stored.
        rebuilt.set_fill_stats(h.fill_stats());
        assert_eq!(rebuilt.fill_stats(), h.fill_stats());
        assert_eq!(rebuilt.mean(), 0.8125);

        // Editing a bin by hand ends the exact statistics; asking for a bin
        // that does not exist does not.
        assert!(h.bin_mut([99]).is_none());
        assert!(h.fill_stats().is_some());
        h.bin_mut([0]).unwrap().value += 1.0;
        assert!(h.fill_stats().is_none());
        // Emptying the histogram makes them exact again.
        h.reset();
        assert_eq!(h.fill_stats(), Some(crate::FillStats::default()));
        assert!(h.mean().is_nan() && h.std_dev().is_nan());
    }

    #[test]
    fn a_profile_without_statistics_weights_centers_by_entries() {
        let mut p = Hist::reg(4, 0.0, 4.0).mean();
        for (x, s) in [(0.5, 10.0), (0.5, 20.0), (1.5, 5.0), (2.5, 30.0)] {
            p.fill_sample(x, s);
        }
        // x values 0.5, 0.5, 1.5, 2.5: all at bin centers.
        let var = (0.25 + 0.25 + 2.25 + 6.25) / 4.0 - 1.25f64.powi(2);
        let exact = p.fill_stats().unwrap();
        assert_eq!(p.mean(), 1.25);
        assert!(close(p.std_dev(), var.sqrt()));
        assert_eq!((exact.sum_ws, exact.sum_ws2), (65.0, 1425.0));
        let binned = p.binned_stats();
        assert_eq!(
            (binned.entries, binned.sum_w, binned.sum_wx),
            (4.0, 4.0, [5.0])
        );
        assert!(close(binned.sum_ws, 65.0) && close(binned.sum_ws2, 1425.0));
        p.set_fill_stats(None);
        assert_eq!(p.mean(), 1.25);
    }

    #[test]
    fn operations_keep_the_statistics_consistent() {
        let fill = |xs: &[(f64, f64, f64)]| {
            let mut h = Hist::reg(4, 0.0, 4.0).reg(2, 0.0, 2.0).weight();
            for &(x, y, w) in xs {
                h.fill_weight(x, y, w);
            }
            h
        };
        let (first, second) = (
            [(0.5, 0.5, 1.0), (1.5, 1.5, 2.0)],
            [(3.5, 0.25, 4.0), (9.0, 0.5, 1.0)],
        );
        let whole = fill(&[first, second].concat());

        // add: the sums of both.
        let mut a = fill(&first);
        a.add(&fill(&second)).unwrap();
        assert_eq!(a.fill_stats(), whole.fill_stats());
        // ...unless one side has none.
        let mut b = fill(&first);
        let mut bare = fill(&second);
        bare.set_fill_stats(None);
        b.add(&bare).unwrap();
        assert!(b.fill_stats().is_none());
        let mut c = bare.clone();
        c.add(&fill(&first)).unwrap();
        assert!(c.fill_stats().is_none());

        // scale: Σw and the moments by c, Σw² by c², the entries unchanged.
        let s = whole.fill_stats().unwrap();
        let mut scaled = whole.clone();
        scaled.scale(-2.0);
        let t = scaled.fill_stats().unwrap();
        assert_eq!(
            (t.entries, t.sum_w, t.sum_w2),
            (s.entries, -2.0 * s.sum_w, 4.0 * s.sum_w2)
        );
        assert_eq!(t.sum_wx, s.sum_wx.map(|v| -2.0 * v));
        assert_eq!(t.sum_wxy, s.sum_wxy.map(|r| r.map(|v| -2.0 * v)));
        assert_eq!(
            scaled.project([0]).unwrap().mean(),
            whole.project([0]).unwrap().mean()
        );

        // project: the kept axes' terms, in the order asked for.
        let px = whole.project([0]).unwrap();
        let sx = px.fill_stats().unwrap();
        assert_eq!(
            (sx.entries, sx.sum_w, sx.sum_w2),
            (s.entries, s.sum_w, s.sum_w2)
        );
        assert_eq!(
            (sx.sum_wx, sx.sum_wxy),
            ([s.sum_wx[0]], [[s.sum_wxy[0][0]]])
        );
        let swapped = whole.project([1, 0]).unwrap();
        let sw = swapped.fill_stats().unwrap();
        assert_eq!(sw.sum_wx, [s.sum_wx[1], s.sum_wx[0]]);
        assert_eq!(sw.sum_wxy[0], [s.sum_wxy[1][1], s.sum_wxy[1][0]]);
        // The x overflow fill is not in the sums, so the exact mean of y
        // differs from the estimate over the projection's bins.
        let py = whole.project([1]).unwrap();
        assert!(close(py.mean(), 4.5 / 7.0));
        assert!(close(py.binned_stats().mean(0), 6.0 / 8.0));

        // rebin: unchanged.
        assert_eq!(whole.rebin(0, 2).unwrap().fill_stats(), whole.fill_stats());

        // reset: zero.
        let mut r = whole.clone();
        r.reset();
        assert_eq!(r.fill_stats(), Some(crate::FillStats::default()));
        assert_eq!(r.entries(), 0.0);
    }

    #[test]
    fn growing_axes_merge_their_statistics_by_label() {
        let template = Hist::str_cat::<&str>(&[])
            .growth()
            .reg(2, 0.0, 2.0)
            .double();
        let (mut a, mut b, mut serial) = (template.clone(), template.clone(), template);
        a.fill("x", 0.5);
        b.fill("y", 1.5);
        b.fill("x", 1.5);
        for (c, v) in [("x", 0.5), ("y", 1.5), ("x", 1.5)] {
            serial.fill(c, v);
        }
        b.add(&a).unwrap(); // b met "y" first
        assert_eq!(b.fill_stats(), serial.fill_stats());
    }

    #[test]
    #[should_panic(expected = "axis 0 takes a number")]
    fn wrong_value_type_panics_in_fill() {
        let mut h = Hist::reg(2, 0.0, 2.0).double();
        h.fill("nope");
    }
}
