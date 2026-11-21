//! [`FillStats`]: the running sums a histogram keeps of what it was filled
//! with, so its mean, standard deviation and entry count are exact instead of
//! estimated from the bin centers.

/// What a [`Histogram`](crate::Histogram) was filled with, as running sums
/// over its fills: the statistics ROOT keeps as `fEntries` and `fTsumw*`.
///
/// Every fill counts in [`entries`](FillStats::entries). The weighted sums
/// count only fills that land in range on every axis (a fill in an under- or
/// overflow bin, or dropped by an axis without one, adds nothing to them).
///
/// The coordinate `xᵢ` a fill contributes on axis `i` is:
///
/// - on a regular, variable or integer axis, the value filled;
/// - on a category or boolean axis, the center of the bin it lands in, in
///   the axis's plot coordinates (`0.5`, `1.5`, …);
/// - on a growing category axis, nothing (0): its bins have no fixed
///   position, so two copies that met the categories in different orders can
///   still be added.
///
/// The weight is the one given to the fill, also for an
/// [`Int64`](crate::Int64) histogram, whose bins round it.
///
/// ```
/// use oxihist::Hist;
///
/// let mut h = Hist::reg(4, 0.0, 4.0).double();
/// h.fill(0.25);
/// h.fill_weight(1.0, 3.0);
/// h.fill(9.0); // overflow: an entry, but not in the sums
/// let stats = h.fill_stats().unwrap();
/// assert_eq!((stats.entries, stats.sum_w, stats.sum_w2), (3.0, 4.0, 10.0));
/// assert_eq!(stats.sum_wx, [3.25]);
/// assert_eq!(stats.mean(0), 0.8125); // not the 1.25 of the bin centers
/// ```
#[doc(alias = "fTsumw", alias = "fEntries", alias = "GetStats")]
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "serde",
    serde(into = "FillStatsRepr", try_from = "FillStatsRepr")
)]
pub struct FillStats<const D: usize> {
    /// The number of fills, in range or not.
    pub entries: f64,
    /// `Σw` over the in-range fills.
    pub sum_w: f64,
    /// `Σw²` over the in-range fills.
    pub sum_w2: f64,
    /// `Σw·xᵢ` for each axis `i`.
    pub sum_wx: [f64; D],
    /// `Σw·xᵢ·xⱼ` for each pair of axes: the diagonal is `Σw·xᵢ²`, and the
    /// matrix is symmetric (`[i][j]` equals `[j][i]`).
    pub sum_wxy: [[f64; D]; D],
    /// `Σw·s` of the sampled value `s`, for a profile
    /// ([`Mean`](crate::Mean), [`WeightedMean`](crate::WeightedMean)); 0 for
    /// a count histogram.
    pub sum_ws: f64,
    /// `Σw·s²` of the sampled value; 0 for a count histogram.
    pub sum_ws2: f64,
}

impl<const D: usize> Default for FillStats<D> {
    fn default() -> Self {
        FillStats {
            entries: 0.0,
            sum_w: 0.0,
            sum_w2: 0.0,
            sum_wx: [0.0; D],
            sum_wxy: [[0.0; D]; D],
            sum_ws: 0.0,
            sum_ws2: 0.0,
        }
    }
}

impl<const D: usize> FillStats<D> {
    /// The weighted mean `Σw·x / Σw` along `axis` (`NaN` when nothing was
    /// filled in range).
    ///
    /// # Panics
    ///
    /// If `axis >= D`.
    #[must_use]
    pub fn mean(&self, axis: usize) -> f64 {
        self.sum_wx[axis] / self.sum_w
    }

    /// The weighted standard deviation `√(Σw·x²/Σw − mean²)` along `axis`
    /// (population form; `NaN` when nothing was filled in range).
    ///
    /// # Panics
    ///
    /// If `axis >= D`.
    #[must_use]
    pub fn std_dev(&self, axis: usize) -> f64 {
        let mean = self.mean(axis);
        let var = self.sum_wxy[axis][axis] / self.sum_w - mean * mean;
        if var.is_nan() {
            return f64::NAN;
        }
        var.max(0.0).sqrt()
    }

    /// The effective number of entries `(Σw)² / Σw²`: the number of
    /// unit-weight fills with the same relative error (0 when nothing was
    /// filled in range).
    #[must_use]
    pub fn effective_entries(&self) -> f64 {
        if self.sum_w2 > 0.0 {
            self.sum_w * self.sum_w / self.sum_w2
        } else {
            0.0
        }
    }

    /// Add an in-range contribution at `x`: `w` and `w2` to the weight sums,
    /// and `w` times the coordinates to the moments.
    #[inline(always)]
    pub(crate) fn accumulate(&mut self, w: f64, w2: f64, x: &[f64; D]) {
        self.sum_w += w;
        self.sum_w2 += w2;
        for i in 0..D {
            let wx = w * x[i];
            self.sum_wx[i] += wx;
            // The upper triangle, mirrored, so the matrix stays symmetric to
            // the last bit.
            for j in i..D {
                let wxy = wx * x[j];
                self.sum_wxy[i][j] += wxy;
                if j != i {
                    self.sum_wxy[j][i] += wxy;
                }
            }
        }
    }

    /// Add `other`'s sums, as if this record had also seen every fill
    /// `other` did.
    pub(crate) fn add(&mut self, other: &Self) {
        self.entries += other.entries;
        self.sum_w += other.sum_w;
        self.sum_w2 += other.sum_w2;
        for i in 0..D {
            self.sum_wx[i] += other.sum_wx[i];
            for j in 0..D {
                self.sum_wxy[i][j] += other.sum_wxy[i][j];
            }
        }
        self.sum_ws += other.sum_ws;
        self.sum_ws2 += other.sum_ws2;
    }

    /// The sums after every bin was multiplied by `c`. A count histogram's
    /// weights scale (`Σw` and the moments by `c`, `Σw²` by `c²`); a
    /// profile's sampled values do (`Σw·s` by `c`, `Σw·s²` by `c²`). The
    /// number of fills is unchanged.
    pub(crate) fn scale(&mut self, c: f64, sampled: bool) {
        if sampled {
            self.sum_ws *= c;
            self.sum_ws2 *= c * c;
            return;
        }
        self.sum_w *= c;
        self.sum_w2 *= c * c;
        self.sum_wx.iter_mut().for_each(|v| *v *= c);
        self.sum_wxy
            .as_flattened_mut()
            .iter_mut()
            .for_each(|v| *v *= c);
    }

    /// The sums of the axes `keep`, in that order (each `< D`).
    pub(crate) fn project<const E: usize>(&self, keep: [usize; E]) -> FillStats<E> {
        FillStats {
            entries: self.entries,
            sum_w: self.sum_w,
            sum_w2: self.sum_w2,
            sum_wx: std::array::from_fn(|i| self.sum_wx[keep[i]]),
            sum_wxy: std::array::from_fn(|i| {
                std::array::from_fn(|j| self.sum_wxy[keep[i]][keep[j]])
            }),
            sum_ws: self.sum_ws,
            sum_ws2: self.sum_ws2,
        }
    }
}

/// How many lanes a histogram spreads its running sums over.
const LANES: usize = 4;

/// A histogram's running sums while it is being filled: [`LANES`] partial
/// [`FillStats`], whose total is the record.
///
/// One set of sums would make every fill wait for the one before it: each
/// adds to the same few numbers, and an addition cannot start until the last
/// one's result is stored. A fill therefore adds to the lane its bin picks,
/// so consecutive fills in different bins usually touch different sums, as
/// they touch different bins.
#[derive(Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[cfg_attr(feature = "serde", serde(into = "FillStats<D>"))]
pub(crate) struct Lanes<const D: usize>([Lane<D>; LANES]);

/// One lane, on its own cache line so that its address is a shift away.
#[derive(Clone, Copy, Default)]
#[repr(align(64))]
struct Lane<const D: usize>(FillStats<D>);

impl<const D: usize> std::fmt::Debug for Lanes<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.total().fmt(f)
    }
}

impl<const D: usize> Default for Lanes<D> {
    fn default() -> Self {
        Lanes([Lane::default(); LANES])
    }
}

impl<const D: usize> From<FillStats<D>> for Lanes<D> {
    fn from(total: FillStats<D>) -> Self {
        let mut lanes = Lanes::default();
        lanes.0[0].0 = total;
        lanes
    }
}

impl<const D: usize> From<Lanes<D>> for FillStats<D> {
    fn from(lanes: Lanes<D>) -> Self {
        lanes.total()
    }
}

impl<const D: usize> Lanes<D> {
    /// The record: the lanes summed, in order.
    pub(crate) fn total(&self) -> FillStats<D> {
        let mut total = self.0[0].0;
        for lane in &self.0[1..] {
            total.add(&lane.0);
        }
        total
    }

    /// The partial sums a fill landing at storage offset `off` adds to.
    #[inline(always)]
    pub(crate) fn lane(&mut self, off: usize) -> &mut FillStats<D> {
        &mut self.0[off % LANES].0
    }

    /// Add `other`'s sums, lane by lane.
    pub(crate) fn add(&mut self, other: &Self) {
        for (a, b) in self.0.iter_mut().zip(&other.0) {
            a.0.add(&b.0);
        }
    }

    /// See [`FillStats::scale`].
    pub(crate) fn scale(&mut self, c: f64, sampled: bool) {
        self.0.iter_mut().for_each(|lane| lane.0.scale(c, sampled));
    }

    /// See [`FillStats::project`].
    pub(crate) fn project<const E: usize>(&self, keep: [usize; E]) -> Lanes<E> {
        Lanes(self.0.map(|lane| Lane(lane.0.project(keep))))
    }
}

/// The serialised form of a [`FillStats`]: serde has no arrays of a generic
/// length, so the per-axis sums travel as vectors and are checked on the way
/// in.
#[cfg(feature = "serde")]
#[derive(serde::Serialize, serde::Deserialize)]
struct FillStatsRepr {
    entries: f64,
    sum_w: f64,
    sum_w2: f64,
    sum_wx: Vec<f64>,
    sum_wxy: Vec<Vec<f64>>,
    sum_ws: f64,
    sum_ws2: f64,
}

#[cfg(feature = "serde")]
impl<const D: usize> From<FillStats<D>> for FillStatsRepr {
    fn from(s: FillStats<D>) -> Self {
        FillStatsRepr {
            entries: s.entries,
            sum_w: s.sum_w,
            sum_w2: s.sum_w2,
            sum_wx: s.sum_wx.to_vec(),
            sum_wxy: s.sum_wxy.iter().map(|row| row.to_vec()).collect(),
            sum_ws: s.sum_ws,
            sum_ws2: s.sum_ws2,
        }
    }
}

#[cfg(feature = "serde")]
impl<const D: usize> TryFrom<FillStatsRepr> for FillStats<D> {
    type Error = crate::Error;

    fn try_from(r: FillStatsRepr) -> crate::Result<Self> {
        let bad = || {
            crate::Error::InvalidData(format!(
                "the fill statistics of a {D}-D histogram need {D} sums per axis"
            ))
        };
        let sum_wx: [f64; D] = r.sum_wx.try_into().map_err(|_| bad())?;
        let rows: Vec<[f64; D]> = r
            .sum_wxy
            .into_iter()
            .map(|row| row.try_into().map_err(|_| bad()))
            .collect::<crate::Result<_>>()?;
        let sum_wxy: [[f64; D]; D] = rows.try_into().map_err(|_| bad())?;
        Ok(FillStats {
            entries: r.entries,
            sum_w: r.sum_w,
            sum_w2: r.sum_w2,
            sum_wx,
            sum_wxy,
            sum_ws: r.sum_ws,
            sum_ws2: r.sum_ws2,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::FillStats;

    #[test]
    fn moments_and_cross_terms_accumulate() {
        let mut s = FillStats::<2>::default();
        s.accumulate(2.0, 4.0, &[1.0, 3.0]);
        s.accumulate(1.0, 1.0, &[4.0, 0.5]);
        assert_eq!((s.sum_w, s.sum_w2), (3.0, 5.0));
        assert_eq!(s.sum_wx, [6.0, 6.5]);
        assert_eq!(s.sum_wxy, [[18.0, 8.0], [8.0, 18.25]]);
        assert_eq!(s.mean(0), 2.0);
        assert!((s.std_dev(0) - 2.0f64.sqrt()).abs() < 1e-12);
        assert!((s.effective_entries() - 1.8).abs() < 1e-12);

        let swapped = s.project([1, 0]);
        assert_eq!(swapped.sum_wx, [6.5, 6.0]);
        assert_eq!(swapped.sum_wxy, [[18.25, 8.0], [8.0, 18.0]]);
        assert_eq!(s.project([1]).sum_wxy, [[18.25]]);
    }

    #[test]
    fn nothing_in_range_has_no_mean() {
        let s = FillStats::<1>::default();
        assert!(s.mean(0).is_nan() && s.std_dev(0).is_nan());
        assert_eq!(s.effective_entries(), 0.0);
    }

    #[test]
    fn scaling_weights_or_samples() {
        let mut s = FillStats::<1>::default();
        s.accumulate(2.0, 4.0, &[3.0]);
        (s.sum_ws, s.sum_ws2) = (10.0, 50.0);
        let mut counts = s;
        counts.scale(-2.0, false);
        assert_eq!((counts.sum_w, counts.sum_w2), (-4.0, 16.0));
        assert_eq!((counts.sum_wx, counts.sum_wxy), ([-12.0], [[-36.0]]));
        assert_eq!(counts.mean(0), s.mean(0));
        let mut profile = s;
        profile.scale(-2.0, true);
        assert_eq!((profile.sum_w, profile.sum_wx), (2.0, [6.0]));
        assert_eq!((profile.sum_ws, profile.sum_ws2), (-20.0, 200.0));
    }
}
