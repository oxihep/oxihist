//! Analysis of filled histograms: bin-by-bin products and ratios, extrema,
//! compatibility tests, smoothing, interpolation, slicing and profiles.

use std::ops::Range;

use oxistats::{chi_square_prob, kolmogorov_prob};

use crate::axis::{Axis, AxisKind};
use crate::cast;
use crate::error::{Error, Result};
use crate::fill_stats::FillStats;
use crate::histogram::Histogram;
use crate::storage::{CountStorage, MeanStorage as _, Storage, WeightedMean, WeightedMeanBin};

/// What a chi-square test assumes about the two histograms it compares
/// (the options of ROOT's `Chi2Test`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum Chi2TestKind {
    /// Both histograms hold unweighted counts (`"UU"`).
    #[default]
    UnweightedUnweighted,
    /// The first holds unweighted counts, the second weighted ones (`"UW"`).
    UnweightedWeighted,
    /// Both hold weighted counts (`"WW"`).
    WeightedWeighted,
}

/// The result of [`Histogram::chi2_test`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Chi2TestResult {
    /// The probability of a chi-square at least this large if the two
    /// histograms follow the same distribution.
    pub p_value: f64,
    /// The chi-square statistic.
    pub chi2: f64,
    /// The degrees of freedom: the bins that entered the sum, minus 1.
    pub ndf: usize,
}

/// The result of [`Histogram::kolmogorov_test`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KsTestResult {
    /// The Kolmogorov probability that the two histograms follow the same
    /// distribution (near 1 when indistinguishable, near 0 when clearly
    /// different).
    pub prob: f64,
    /// The largest distance between the two normalised cumulatives.
    pub distance: f64,
}

/// [`Error::BinningMismatch`] unless `a` and `b` bin the same way on every
/// axis.
fn check_binning<S: Storage, T: Storage, const D: usize>(
    a: &Histogram<S, D>,
    b: &Histogram<T, D>,
    what: &str,
) -> Result<()> {
    if a.axes()
        .iter()
        .zip(b.axes())
        .all(|(x, y)| x.same_binning(y))
    {
        Ok(())
    } else {
        Err(Error::BinningMismatch(format!(
            "{what}: the histograms' axes differ"
        )))
    }
}

impl<S: Storage, const D: usize> Histogram<S, D> {
    /// The in-range bin with the largest value, as its index and that value.
    /// Among equal values the first wins (the first axis varying fastest);
    /// `NaN` values are skipped. `None` if there is no such bin.
    #[doc(alias = "GetMaximum", alias = "GetMaximumBin")]
    #[must_use]
    pub fn maximum(&self) -> Option<([isize; D], f64)> {
        self.extremum(|v, best| v > best)
    }

    /// The in-range bin with the smallest value, as its index and that
    /// value. Among equal values the first wins (the first axis varying
    /// fastest); `NaN` values are skipped. `None` if there is no such bin.
    #[doc(alias = "GetMinimum", alias = "GetMinimumBin")]
    #[must_use]
    pub fn minimum(&self) -> Option<([isize; D], f64)> {
        self.extremum(|v, best| v < best)
    }

    fn extremum(&self, better: impl Fn(f64, f64) -> bool) -> Option<([isize; D], f64)> {
        let mut best: Option<([isize; D], f64)> = None;
        for (idx, b) in self.iter(false) {
            let v = S::value(b);
            if !v.is_nan() && best.is_none_or(|(_, w)| better(v, w)) {
                best = Some((idx, v));
            }
        }
        best
    }
}

impl<S: CountStorage, const D: usize> Histogram<S, D> {
    /// Multiply bin by bin by `other`, flow bins included, propagating the
    /// variances as for independent values:
    /// `σ² = σ₁²·v₂² + σ₂²·v₁²`.
    ///
    /// The result has no [fill statistics](Histogram::fill_stats). A storage
    /// that does not track variances keeps the product alone (see
    /// [`CountStorage::bin`]).
    ///
    /// # Errors
    ///
    /// [`Error::BinningMismatch`] (leaving `self` unchanged) if the binnings
    /// differ.
    pub fn multiply<T: Storage>(&mut self, other: &Histogram<T, D>) -> Result<()> {
        self.combine(other, "multiply", |v1, e1, v2, e2| {
            (v1 * v2, e1 * v2 * v2 + e2 * v1 * v1)
        })
    }

    /// Divide bin by bin by `other`, flow bins included, propagating the
    /// variances as for independent values:
    /// `σ² = (σ₁²·v₂² + σ₂²·v₁²) / v₂⁴`. A bin whose denominator is 0
    /// becomes 0, with variance 0.
    ///
    /// The result has no [fill statistics](Histogram::fill_stats). A storage
    /// that does not track variances keeps the ratio alone (see
    /// [`CountStorage::bin`]); [`Int64`](crate::Int64) rounds it.
    ///
    /// # Errors
    ///
    /// [`Error::BinningMismatch`] (leaving `self` unchanged) if the binnings
    /// differ.
    pub fn divide<T: Storage>(&mut self, other: &Histogram<T, D>) -> Result<()> {
        self.combine(other, "divide", |v1, e1, v2, e2| {
            if v2 == 0.0 {
                return (0.0, 0.0);
            }
            let v2sq = v2 * v2;
            (v1 / v2, (e1 * v2sq + e2 * v1 * v1) / (v2sq * v2sq))
        })
    }

    /// Replace each bin by `f(value, variance, other's value, other's
    /// variance)`, a new value and variance.
    fn combine<T: Storage>(
        &mut self,
        other: &Histogram<T, D>,
        what: &str,
        f: impl Fn(f64, f64, f64, f64) -> (f64, f64),
    ) -> Result<()> {
        check_binning(self, other, what)?;
        for (idx, b) in other.iter(true) {
            let bin = self.bin_mut(idx).expect("the same binning");
            let (value, variance) = f(S::value(bin), S::variance(bin), T::value(b), T::variance(b));
            *bin = S::bin(value, variance);
        }
        Ok(())
    }
}

impl<S: CountStorage> Histogram<S, 1> {
    /// A chi-square test of whether this histogram and `other` follow the
    /// same distribution, over the in-range bins (ROOT's `Chi2Test`).
    /// `kind` says which of the two hold weighted counts; the weighted
    /// forms use the bins' variances.
    ///
    /// # Errors
    ///
    /// [`Error::BinningMismatch`] if the binnings differ.
    #[doc(alias = "Chi2Test")]
    pub fn chi2_test<T: CountStorage>(
        &self,
        other: &Histogram<T, 1>,
        kind: Chi2TestKind,
    ) -> Result<Chi2TestResult> {
        check_binning(self, other, "chi2_test")?;
        let (c1, c2) = (self.values(false), other.values(false));
        let positive = |v: Vec<f64>| v.into_iter().map(|e| e.max(0.0)).collect::<Vec<f64>>();
        let (e1, e2) = (
            positive(self.variances(false)),
            positive(other.variances(false)),
        );
        let n = c1.len();
        let sum1: f64 = c1.iter().sum();
        let sum2: f64 = c2.iter().sum();

        let mut chi2 = 0.0;
        let mut used = 0usize;
        match kind {
            Chi2TestKind::UnweightedUnweighted => {
                for i in 0..n {
                    let (a, b) = (c1[i], c2[i]);
                    if a == 0.0 && b == 0.0 {
                        continue;
                    }
                    used += 1;
                    let delta = sum2 * a - sum1 * b;
                    chi2 += delta * delta / (a + b);
                }
                if sum1 * sum2 > 0.0 {
                    chi2 /= sum1 * sum2;
                }
            }
            Chi2TestKind::WeightedWeighted => {
                for i in 0..n {
                    if e1[i] == 0.0 && e2[i] == 0.0 {
                        continue;
                    }
                    used += 1;
                    let delta = sum1 * c2[i] - sum2 * c1[i];
                    let denom = sum1 * sum1 * e2[i] + sum2 * sum2 * e1[i];
                    if denom > 0.0 {
                        chi2 += delta * delta / denom;
                    }
                }
            }
            Chi2TestKind::UnweightedWeighted => {
                let sumw2: f64 = e2.iter().sum();
                for i in 0..n {
                    let (a, b) = (c1[i], c2[i]);
                    let mut s2 = e2[i];
                    if a == 0.0 && b == 0.0 {
                        continue;
                    }
                    if b == 0.0 && s2 == 0.0 {
                        // An empty weighted bin: take the average weight as
                        // its variance.
                        if sumw2 > 0.0 && sum2 > 0.0 {
                            s2 = sumw2 / sum2;
                        } else {
                            continue;
                        }
                    }
                    used += 1;
                    // The probability of the bin that maximises the
                    // likelihood of both observations: the root of a
                    // quadratic.
                    let var1 = sum2 * b - sum1 * s2;
                    let var2 = (var1 * var1 + 4.0 * sum2 * sum2 * a * s2).sqrt();
                    let prob = (var1 + var2) / (2.0 * sum2 * sum2);
                    let (expected1, expected2) = (prob * sum1, prob * sum2);
                    if expected1 > 0.0 {
                        chi2 += (a - expected1).powi(2) / expected1;
                    }
                    if s2 > 0.0 {
                        chi2 += (b - expected2).powi(2) / s2;
                    }
                }
            }
        }
        let ndf = used.saturating_sub(1);
        Ok(Chi2TestResult {
            p_value: chi_square_prob(chi2, ndf),
            chi2,
            ndf,
        })
    }

    /// A Kolmogorov–Smirnov test of whether this histogram and `other`
    /// follow the same distribution, from the largest distance between the
    /// cumulatives of their in-range bins (ROOT's `KolmogorovTest`). The
    /// probability is the one for unweighted counts; it is 0, with distance
    /// 0, when either histogram has no positive total.
    ///
    /// # Errors
    ///
    /// [`Error::BinningMismatch`] if the binnings differ.
    #[doc(alias = "KolmogorovTest")]
    pub fn kolmogorov_test<T: CountStorage>(
        &self,
        other: &Histogram<T, 1>,
    ) -> Result<KsTestResult> {
        check_binning(self, other, "kolmogorov_test")?;
        let (c1, c2) = (self.values(false), other.values(false));
        let sum1: f64 = c1.iter().sum();
        let sum2: f64 = c2.iter().sum();
        if sum1 <= 0.0 || sum2 <= 0.0 {
            return Ok(KsTestResult {
                prob: 0.0,
                distance: 0.0,
            });
        }
        let (mut cum1, mut cum2, mut distance) = (0.0, 0.0, 0.0f64);
        for (a, b) in c1.iter().zip(&c2) {
            cum1 += a / sum1;
            cum2 += b / sum2;
            distance = distance.max((cum1 - cum2).abs());
        }
        let z = distance * (sum1 * sum2 / (sum1 + sum2)).sqrt();
        Ok(KsTestResult {
            prob: kolmogorov_prob(z),
            distance,
        })
    }

    /// Smooth the in-range values `times` times with the 353QH-twice
    /// running-median smoother (ROOT's `TH1::Smooth`). The variances are
    /// kept where the storage tracks them, the flow bins are untouched, and
    /// the result has no [fill statistics](Histogram::fill_stats). A
    /// histogram of fewer than 3 bins is left as it is, as is any histogram
    /// when `times` is 0.
    #[doc(alias = "Smooth")]
    pub fn smooth(&mut self, times: usize) {
        let mut values = self.values(false);
        if values.len() < 3 || times == 0 {
            return;
        }
        smooth_353qh(&mut values, times);
        for (i, v) in values.into_iter().enumerate() {
            let bin = self.bin_mut([cast::signed(i)]).expect("in range");
            *bin = S::bin(v, S::variance(bin));
        }
    }
}

impl<S: Storage> Histogram<S, 1> {
    /// The value at `x`, interpolated linearly between the centers of the
    /// two bins around it (ROOT's `TH1::Interpolate`). Between an end of the
    /// axis and the center of its first or last bin, that bin's value.
    ///
    /// `None` for an `x` outside the axis (its upper edge counts as inside),
    /// a `NaN`, or a categorical axis.
    #[doc(alias = "Interpolate")]
    #[must_use]
    pub fn interpolate(&self, x: f64) -> Option<f64> {
        let edges = self.axis(0).edges()?;
        let n = edges.len() - 1;
        // Written so that a NaN fails it.
        if !(x >= edges[0] && x <= edges[n]) {
            return None;
        }
        let centers = self.axis(0).centers();
        let values = self.values(false);
        if x <= centers[0] {
            return Some(values[0]);
        }
        if x >= centers[n - 1] {
            return Some(values[n - 1]);
        }
        // centers[k - 1] <= x < centers[k]
        let k = centers.partition_point(|&c| c <= x);
        let (x0, x1) = (centers[k - 1], centers[k]);
        let (y0, y1) = (values[k - 1], values[k]);
        Some(y0 + (x - x0) * (y1 - y0) / (x1 - x0))
    }

    /// The part of the histogram between two values: the bins overlapping
    /// `range` (from the one holding `range.start` to the one holding
    /// `range.end`, or the one before when `range.end` is its low edge), on
    /// an axis cut down to them. A range reaching past an end of the axis
    /// stops there.
    ///
    /// The bins cut off are merged into the flow bins, so the total is
    /// kept; on an axis without a flow bin on that side they are dropped.
    /// The result has no [fill statistics](Histogram::fill_stats): its mean
    /// and standard deviation are those of the kept bins.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidAxis`] for a categorical axis, a range that is empty,
    /// reversed or has a `NaN` bound, or one that lies wholly outside the
    /// axis.
    pub fn slice(&self, range: Range<f64>) -> Result<Self> {
        let axis = self.axis(0);
        let edges = axis.edges().ok_or_else(|| {
            Error::InvalidAxis("a categorical axis cannot be sliced by value".into())
        })?;
        let n = edges.len() - 1;
        let (lo, hi) = (range.start, range.end);
        if lo.is_nan() || hi.is_nan() || lo >= hi {
            return Err(Error::InvalidAxis(format!(
                "cannot slice the empty range {lo}..{hi}"
            )));
        }
        if hi <= edges[0] || lo >= edges[n] {
            return Err(Error::InvalidAxis(format!(
                "the range {lo}..{hi} lies outside the axis [{}, {})",
                edges[0], edges[n]
            )));
        }
        // The bin holding `lo`, and the last bin starting below `hi`.
        let first = edges.partition_point(|&e| e <= lo).saturating_sub(1);
        let last = edges[..n].partition_point(|&e| e < hi) - 1;

        let mut cut = match axis.kind() {
            AxisKind::Regular { .. } => {
                Axis::regular(last - first + 1, edges[first], edges[last + 1])
            }
            AxisKind::Integer { lo, .. } => {
                let at = |i: usize| lo + i64::try_from(i).expect("a bin index");
                Axis::integer(at(first), at(last + 1))
            }
            _ => Axis::variable(&edges[first..=last + 1]),
        }
        .with_underflow(axis.underflow())
        .with_overflow(axis.overflow());
        cut.name.clone_from(&axis.name);
        cut.label.clone_from(&axis.label);

        let mut out = Histogram::<S, 1>::new([cut]);
        out.name.clone_from(&self.name);
        out.title.clone_from(&self.title);
        let (first, last) = (cast::signed(first), cast::signed(last));
        for (idx, b) in self.iter(true) {
            let to = if idx[0] < first {
                -1
            } else if idx[0] > last {
                last - first + 1
            } else {
                idx[0] - first
            };
            // `None`: cut off towards a flow bin the axis does not keep.
            if let Some(bin) = out.bin_mut([to]) {
                S::merge(bin, b);
            }
        }
        out.set_fill_stats(None);
        Ok(out)
    }
}

impl<S: CountStorage> Histogram<S, 2> {
    /// The profile along `axis` (0 or 1): for each of its bins, flow bins
    /// included, the mean and spread of the other axis's coordinate over
    /// that bin's entries (ROOT's `ProfileX` and `ProfileY`). Each in-range
    /// bin of the other axis counts at its center, with its value as the
    /// weight and its variance as `Σw²`.
    ///
    /// The profile inherits the [fill statistics](Histogram::fill_stats):
    /// the sums of the other axis's coordinate become those of the sampled
    /// value.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidAxis`] if `axis` is not 0 or 1.
    #[doc(alias = "ProfileX", alias = "ProfileY")]
    pub fn profile(&self, axis: usize) -> Result<Histogram<WeightedMean, 1>> {
        if axis > 1 {
            return Err(Error::InvalidAxis(format!(
                "a 2-D histogram has no axis {axis} to profile along"
            )));
        }
        let other = 1 - axis;
        let kept = self.axis(axis);
        let centers = self.axis(other).centers();
        let mut bins = Vec::with_capacity(kept.extent());
        let n = cast::signed(kept.size());
        for k in -1..=n {
            if kept.position(k).is_none() {
                continue; // a flow bin the axis does not keep
            }
            let mut acc = WeightedMeanBin::default();
            let mut sum_w2 = 0.0;
            for (o, &center) in centers.iter().enumerate() {
                let mut idx = [0isize; 2];
                idx[axis] = k;
                idx[other] = cast::signed(o);
                let bin = self.bin(idx).expect("an existing bin");
                WeightedMean::add(&mut acc, center, S::value(bin));
                sum_w2 += S::variance(bin);
            }
            acc.sum_of_weights_squared = sum_w2;
            bins.push(acc);
        }
        let mut out = Histogram::from_bins([kept.clone()], bins)?;
        out.name.clone_from(&self.name);
        out.title.clone_from(&self.title);
        out.set_fill_stats(self.fill_stats().map(|s| FillStats {
            entries: s.entries,
            sum_w: s.sum_w,
            sum_w2: s.sum_w2,
            sum_wx: [s.sum_wx[axis]],
            sum_wxy: [[s.sum_wxy[axis][axis]]],
            sum_ws: s.sum_wx[other],
            sum_ws2: s.sum_wxy[other][other],
        }));
        Ok(out)
    }
}

/// The median of `values` (3 or 5 of them). A total order, so a `NaN` bin
/// does not panic the smoother.
fn median(values: &[f64]) -> f64 {
    let mut v = [0.0f64; 5];
    let n = values.len();
    v[..n].copy_from_slice(values);
    v[..n].sort_by(f64::total_cmp);
    v[n / 2]
}

/// The 353QH-twice smoother of ROOT's `TH1::SmoothArray`, applied `times`
/// times to `xx` (at least 3 values): running medians of 3, 5 and 3,
/// quadratic interpolation over flat stretches, a Hanning running mean, then
/// the same again on the residuals.
fn smooth_353qh(xx: &mut [f64], times: usize) {
    let nn = xx.len();
    let mut hh = [0.0f64; 5];
    let mut yy = vec![0.0f64; nn];
    let mut zz = vec![0.0f64; nn];
    let mut rr = vec![0.0f64; nn];

    for _ in 0..times {
        yy.copy_from_slice(xx);

        for pass in 0..2 {
            // 353: running medians of window 3, 5, 3.
            for kk in 0..3 {
                zz.copy_from_slice(&yy);
                let window = if kk == 1 { 5 } else { 3 };
                let first = if kk == 1 { 2 } else { 1 };
                let last = if kk == 1 { nn - 2 } else { nn - 1 };
                for ii in first..last {
                    zz[ii] = median(&yy[ii - first..ii - first + window]);
                }
                if kk == 0 {
                    // The end points of the first median of 3.
                    hh[0] = zz[1];
                    hh[1] = zz[0];
                    hh[2] = 3.0 * zz[1] - 2.0 * zz[2];
                    zz[0] = median(&hh[..3]);
                    hh[0] = zz[nn - 2];
                    hh[1] = zz[nn - 1];
                    hh[2] = 3.0 * zz[nn - 2] - 2.0 * zz[nn - 3];
                    zz[nn - 1] = median(&hh[..3]);
                }
                if kk == 1 {
                    // The points next to the ends of the median of 5.
                    zz[1] = median(&yy[..3]);
                    zz[nn - 2] = median(&yy[nn - 3..]);
                }
                yy.copy_from_slice(&zz);
            }

            // Quadratic interpolation over flat stretches of 3 values.
            zz.copy_from_slice(&yy);
            for ii in 2..nn.saturating_sub(2) {
                if zz[ii - 1] != zz[ii] || zz[ii] != zz[ii + 1] {
                    continue;
                }
                let (left, right) = (zz[ii - 2] - zz[ii], zz[ii + 2] - zz[ii]);
                if left * right <= 0.0 {
                    continue;
                }
                // `near` is the neighbour on the side of the smaller step.
                let (far, back, near) = if right.abs() > left.abs() {
                    (ii - 2, ii + 2, ii - 1)
                } else {
                    (ii + 2, ii - 2, ii + 1)
                };
                yy[ii] = -0.5 * zz[back] + zz[ii] / 0.75 + zz[far] / 6.0;
                yy[near] = 0.5 * (zz[far] - zz[back]) + zz[ii];
            }

            // Hanning running mean.
            for ii in 1..nn - 1 {
                zz[ii] = 0.25 * yy[ii - 1] + 0.5 * yy[ii] + 0.25 * yy[ii + 1];
            }
            zz[0] = yy[0];
            zz[nn - 1] = yy[nn - 1];

            if pass == 0 {
                // Keep the smoothed values and smooth the residuals next.
                rr.copy_from_slice(&zz);
                for ii in 0..nn {
                    yy[ii] = xx[ii] - zz[ii];
                }
            }
        }

        // Counts stay non-negative unless some were negative to begin with.
        let lowest = xx.iter().copied().fold(f64::INFINITY, f64::min);
        for ii in 0..nn {
            let v = rr[ii] + zz[ii];
            xx[ii] = if lowest < 0.0 { v } else { v.max(0.0) };
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{Axis, Chi2TestKind, Hist, Histogram, Weight, WeightedSum};

    /// The parabola `h` and its perturbed copy `g` of oxiroot's
    /// `analysis.root` fixture: 20 bins over `[0, 20)`. The expected numbers
    /// below are ROOT 6.40's for these contents.
    const H: [f64; 20] = [
        20.0, 38.0, 54.0, 68.0, 80.0, 90.0, 98.0, 104.0, 108.0, 110.0, 110.0, 108.0, 104.0, 98.0,
        90.0, 80.0, 68.0, 54.0, 38.0, 20.0,
    ];
    const G: [f64; 20] = [
        23.0, 41.0, 52.0, 71.0, 83.0, 88.0, 101.0, 107.0, 106.0, 113.0, 113.0, 106.0, 107.0, 101.0,
        88.0, 83.0, 71.0, 52.0, 41.0, 23.0,
    ];

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9 * (1.0 + b.abs())
    }

    fn counts(values: &[f64]) -> Histogram<crate::Double, 1> {
        let mut bins = vec![0.0];
        bins.extend_from_slice(values);
        bins.push(0.0);
        let n = values.len();
        Histogram::from_bins([Axis::regular(n, 0.0, n as f64)], bins).unwrap()
    }

    /// Weighted bins with the given variance for each value.
    fn weighted(values: &[f64], variance: impl Fn(f64) -> f64) -> Histogram<Weight, 1> {
        let mut bins = vec![WeightedSum::default()];
        bins.extend(values.iter().map(|&value| WeightedSum {
            value,
            variance: variance(value),
        }));
        bins.push(WeightedSum::default());
        let n = values.len();
        Histogram::from_bins([Axis::regular(n, 0.0, n as f64)], bins).unwrap()
    }

    #[test]
    fn chi2_test_matches_root() {
        let (h, g) = (counts(&H), counts(&G));
        // Chi2TestX("UU"): chi2 = 1.01136014673, ndf = 19, p = 0.999999999141.
        let uu = h.chi2_test(&g, Chi2TestKind::UnweightedUnweighted).unwrap();
        assert_eq!(uu.ndf, 19);
        assert!(close(uu.chi2, 1.01136014673), "{}", uu.chi2);
        assert!(close(uu.p_value, 0.999999999141), "{}", uu.p_value);

        // The fixture's weighted pair: `hw1` has Poisson errors, `hw2` the
        // errors √content + 1.
        let hw1 = weighted(&H, |c| c);
        let hw2 = weighted(&G, |c| (c.sqrt() + 1.0).powi(2));
        let uw = hw1
            .chi2_test(&hw2, Chi2TestKind::UnweightedWeighted)
            .unwrap();
        let ww = hw1.chi2_test(&hw2, Chi2TestKind::WeightedWeighted).unwrap();
        assert_eq!((uw.ndf, ww.ndf), (19, 19));
        assert!(close(uw.chi2, 0.870192920943), "{}", uw.chi2);
        assert!(close(uw.p_value, 0.99999999978), "{}", uw.p_value);
        assert!(close(ww.chi2, 0.872240891087), "{}", ww.chi2);
        assert!(close(ww.p_value, 0.999999999776), "{}", ww.p_value);
        // The kind decides the statistic, not the storages: counts against
        // weighted bins give the same answer.
        let mixed = h.chi2_test(&hw2, Chi2TestKind::UnweightedWeighted).unwrap();
        assert_eq!(mixed, uw);
    }

    #[test]
    fn kolmogorov_test_matches_root() {
        let (h, g) = (counts(&H), counts(&G));
        // KolmogorovTest: distance 0.00314335346182, probability 1.
        let ks = h.kolmogorov_test(&g).unwrap();
        assert!(close(ks.distance, 0.00314335346182), "{}", ks.distance);
        assert!(close(ks.prob, 1.0), "{}", ks.prob);
        // A histogram against itself, and against nothing.
        let same = h.kolmogorov_test(&h).unwrap();
        assert_eq!((same.distance, same.prob), (0.0, 1.0));
        let empty = h.kolmogorov_test(&counts(&[0.0; 20])).unwrap();
        assert_eq!((empty.distance, empty.prob), (0.0, 0.0));
    }

    #[test]
    fn comparisons_need_the_same_binning() {
        let h = counts(&H);
        let other = Hist::reg(10, 0.0, 10.0).double();
        assert!(h.chi2_test(&other, Chi2TestKind::default()).is_err());
        assert!(h.kolmogorov_test(&other).is_err());
    }

    #[test]
    fn interpolate_matches_root() {
        let h = counts(&H);
        // TH1::Interpolate.
        for (x, want) in [(5.0, 85.0), (5.5, 90.0), (10.3, 110.0), (0.7, 23.6)] {
            assert!(close(h.interpolate(x).unwrap(), want), "at {x}");
        }
        // Between an end of the axis and the outermost bin center.
        assert_eq!(h.interpolate(0.0), Some(20.0));
        assert_eq!(h.interpolate(19.8), Some(20.0));
        assert_eq!(h.interpolate(20.0), Some(20.0));
        // Outside the axis there is nothing to interpolate.
        assert_eq!(h.interpolate(-0.1), None);
        assert_eq!(h.interpolate(100.0), None);
        assert_eq!(h.interpolate(f64::NAN), None);
        assert_eq!(Hist::str_cat(&["a", "b"]).double().interpolate(0.5), None);
        // Variable bins interpolate between their own centers.
        let mut v = Hist::var(&[0.0, 2.0, 3.0]).double();
        v.fill_weight(1.0, 4.0);
        v.fill_weight(2.5, 1.0);
        assert!(close(v.interpolate(2.0).unwrap(), 2.0)); // 2/3 of the way from 1 to 2.5
        // A single bin is flat.
        let mut one = Hist::reg(1, 0.0, 1.0).double();
        one.fill_weight(0.2, 7.0);
        assert_eq!(one.interpolate(0.9), Some(7.0));
    }

    #[test]
    fn smooth_matches_root() {
        // TH1D(10, 0, 10) with these contents; Smooth(1).
        let mut h = counts(&[2.0, 8.0, 4.0, 9.0, 3.0, 7.0, 1.0, 6.0, 5.0, 2.0]);
        h.smooth(1);
        let want = [2.0, 3.5, 4.0, 4.25, 4.75, 5.0, 5.0, 5.0, 5.0, 5.0];
        for (got, want) in h.values(false).iter().zip(want) {
            assert!(close(*got, want), "got {got}, want {want}");
        }
        assert_eq!((h.value([-1]), h.value([10])), (Some(0.0), Some(0.0)));
    }

    #[test]
    fn smooth_edge_cases() {
        // Fewer than 3 bins, or no passes: unchanged, statistics included.
        let mut two = Hist::reg(2, 0.0, 2.0).double();
        two.fill(0.5);
        two.smooth(3);
        assert_eq!(two.values(false), vec![1.0, 0.0]);
        assert!(two.fill_stats().is_some());
        let mut h = Hist::reg(5, 0.0, 5.0).weight();
        h.fill_many([0.5, 1.5, 1.5, 4.5]);
        let before = h.clone();
        h.smooth(0);
        assert_eq!(h, before);
        assert!(h.fill_stats().is_some());
        // Smoothing keeps the variances and drops the statistics.
        h.smooth(1);
        assert_eq!(h.variances(true), before.variances(true));
        assert!(h.fill_stats().is_none());
        // Exactly 3 bins, and a NaN bin, do not panic.
        let mut three = counts(&[1.0, 5.0, 2.0]);
        three.smooth(2);
        assert_eq!(three.values(false).len(), 3);
        let mut nan = counts(&[0.0, 1.0, 2.0, f64::NAN, 4.0, 5.0]);
        nan.smooth(2);
        assert_eq!(nan.values(false).len(), 6);
    }

    #[test]
    fn multiply_and_divide_propagate_variances() {
        let mut num = Hist::reg(2, 0.0, 2.0).weight();
        num.fill_many([0.5, 0.5, 1.5]);
        let mut den = Hist::reg(2, 0.0, 2.0).weight();
        den.fill_many([0.5; 4]);
        den.fill_many([1.5; 4]);

        let mut eff = num.clone();
        eff.divide(&den).unwrap();
        assert_eq!(eff.values(false), vec![0.5, 0.25]);
        // (σ₁²·v₂² + σ₂²·v₁²) / v₂⁴ = (2·16 + 4·4) / 256.
        assert!(close(eff.variance([0]).unwrap(), 48.0 / 256.0));
        assert!(eff.fill_stats().is_none());

        eff.multiply(&den).unwrap();
        assert!(close(eff.value([0]).unwrap(), 2.0));
        assert!(close(eff.value([1]).unwrap(), 1.0));
        // σ₁²·v₂² + σ₂²·v₁² = (48/256)·16 + 4·0.25.
        assert!(close(eff.variance([0]).unwrap(), 4.0));

        // A zero denominator gives an empty bin, flow bins included.
        let mut r = num.clone();
        let mut zero = Hist::reg(2, 0.0, 2.0).weight();
        zero.fill(0.5);
        r.fill(-1.0);
        r.divide(&zero).unwrap();
        assert_eq!(r.values(true), vec![0.0, 2.0, 0.0, 0.0]);
        assert_eq!(r.variances(true), vec![0.0, 6.0, 0.0, 0.0]);

        // Storages may differ; a storage without variances keeps the value.
        let mut d = Hist::reg(2, 0.0, 2.0).double();
        d.fill_many([0.5, 0.5, 1.5]);
        d.divide(&den).unwrap();
        assert_eq!(d.values(false), vec![0.5, 0.25]);
        let mut i = Hist::reg(2, 0.0, 2.0).int64();
        i.fill_many([0.5, 0.5, 1.5]);
        i.multiply(&den).unwrap();
        assert_eq!(i.values(false), vec![8.0, 4.0]);
    }

    #[test]
    fn multiply_and_divide_need_the_same_binning() {
        let mut a = Hist::reg(2, 0.0, 2.0).weight();
        a.fill(0.5);
        let before = a.clone();
        let other = Hist::reg(3, 0.0, 2.0).weight();
        assert!(a.multiply(&other).is_err());
        assert!(a.divide(&other).is_err());
        assert_eq!(a, before);
        assert_eq!(a.fill_stats(), before.fill_stats());
    }

    #[test]
    fn two_dimensional_ratio() {
        let mut num = Hist::reg(2, 0.0, 2.0).reg(2, 0.0, 2.0).weight();
        let mut den = num.clone();
        num.fill(0.5, 1.5);
        den.fill(0.5, 1.5);
        den.fill(0.5, 1.5);
        num.divide(&den).unwrap();
        assert_eq!(num.value([0, 1]), Some(0.5));
        assert_eq!(num.value([1, 1]), Some(0.0));
    }

    #[test]
    fn extrema_are_in_range_bins() {
        let mut h = Hist::reg(4, 0.0, 4.0).double();
        h.fill_many([0.5, 0.5, 1.5, 2.5, 2.5, 2.5]);
        for _ in 0..9 {
            h.fill(7.0); // the overflow is not a candidate
        }
        assert_eq!(h.maximum(), Some(([2], 3.0)));
        assert_eq!(h.minimum(), Some(([3], 0.0)));

        let mut h2 = Hist::reg(2, 0.0, 2.0).reg(2, 0.0, 2.0).double();
        h2.fill(0.5, 0.5);
        h2.fill(1.5, 1.5);
        h2.fill(1.5, 1.5);
        assert_eq!(h2.maximum(), Some(([1, 1], 2.0)));
        // Ties go to the first bin, the first axis varying fastest.
        assert_eq!(h2.minimum(), Some(([1, 0], 0.0)));

        // NaN bins are skipped; with no bin at all there is no extremum.
        let nan = counts(&[f64::NAN, 1.0, f64::NAN]);
        assert_eq!(nan.maximum(), Some(([1], 1.0)));
        assert_eq!(nan.minimum(), Some(([1], 1.0)));
        assert_eq!(counts(&[f64::NAN]).maximum(), None);
        let empty = Hist::str_cat::<&str>(&[]).growth().double();
        assert_eq!((empty.maximum(), empty.minimum()), (None, None));
    }

    #[test]
    fn slice_keeps_the_total_in_the_flow_bins() {
        let mut g = Hist::reg(10, 0.0, 10.0).double().named("g").titled("t");
        g.axis_meta_mut(0).1.push_str("x [cm]");
        for i in 0..10 {
            g.fill(f64::from(i) + 0.5);
        }
        let s = g.slice(2.5..4.5).unwrap();
        assert_eq!(s.values(true), vec![2.0, 1.0, 1.0, 1.0, 5.0]);
        assert_eq!(s.axis(0), &Axis::regular(3, 2.0, 5.0).with_label("x [cm]"));
        assert_eq!((s.name.as_str(), s.title.as_str()), ("g", "t"));
        assert_eq!(s.integral(true), g.integral(true));
        // The slice describes its own bins.
        assert!(s.fill_stats().is_none());
        assert!(close(s.mean(), 3.5));
        // An upper bound on a bin edge leaves that bin out.
        assert_eq!(g.slice(2.0..5.0).unwrap(), s);
        // A range past the ends stops at them.
        assert_eq!(
            g.slice(-5.0..2.0).unwrap().values(true),
            [0.0, 1.0, 1.0, 8.0]
        );
        assert_eq!(g.slice(-1.0..99.0).unwrap().values(true), g.values(true));
    }

    #[test]
    fn slice_merges_variances_and_follows_the_axis_kind() {
        let mut h = Hist::reg(4, 0.0, 4.0).weight();
        for (x, w) in [(0.5, 2.0), (1.5, 2.0), (2.5, 3.0), (3.5, 3.0)] {
            h.fill_weight(x, w);
        }
        let s = h.slice(1.5..2.5).unwrap();
        assert_eq!(s.values(true), vec![2.0, 2.0, 3.0, 3.0]);
        assert_eq!(s.variances(true), vec![4.0, 4.0, 9.0, 9.0]);

        let mut v = Hist::var(&[0.0, 1.0, 5.0, 6.0, 10.0]).double();
        v.fill_many([0.5, 2.0, 5.5, 7.0]);
        let s = v.slice(1.0..6.0).unwrap();
        assert_eq!(s.axis(0), &Axis::variable(&[1.0, 5.0, 6.0]));
        assert_eq!(s.values(true), vec![1.0, 1.0, 1.0, 1.0]);

        let mut i = Hist::int(-2, 3).int64();
        i.fill_many([-2, -1, 0, 1, 2]);
        let s = i.slice(-1.0..1.5).unwrap();
        assert_eq!(s.axis(0), &Axis::integer(-1, 2));
        assert_eq!(s.values(true), vec![1.0, 1.0, 1.0, 1.0, 1.0]);

        // Without flow bins, what is cut off is gone.
        let mut n = Hist::reg(4, 0.0, 4.0).no_flow().double();
        n.fill_many([0.5, 1.5, 2.5, 3.5]);
        let s = n.slice(1.0..3.0).unwrap();
        assert_eq!(s.values(true), vec![1.0, 1.0]);
    }

    #[test]
    fn slice_refuses_what_it_cannot_cut() {
        let h = Hist::reg(4, 0.0, 4.0).double();
        for bad in [2.0..2.0, 3.0..1.0, f64::NAN..1.0, 0.0..f64::NAN] {
            assert!(h.slice(bad).is_err());
        }
        assert!(h.slice(4.0..9.0).is_err());
        assert!(h.slice(-3.0..0.0).is_err());
        assert!(Hist::str_cat(&["a"]).double().slice(0.0..1.0).is_err());
        assert!(Hist::boolean().double().slice(0.0..1.0).is_err());
    }

    #[test]
    fn profile_is_the_mean_of_the_other_axis() {
        let mut h = Hist::reg(2, 0.0, 2.0).reg(2, 0.0, 2.0).double().named("h");
        h.fill(0.5, 0.5);
        h.fill(0.5, 1.5);
        for _ in 0..3 {
            h.fill(1.5, 0.5);
        }
        let px = h.profile(0).unwrap();
        assert_eq!(px.values(false), vec![1.0, 0.5]);
        assert_eq!(px.name, "h");
        let py = h.profile(1).unwrap();
        // y bin 0: x centers 0.5 once and 1.5 three times.
        assert_eq!(py.values(false), vec![1.25, 0.5]);
        assert!(h.profile(2).is_err());

        // The statistics carry over: x as the axis, y as the sampled value.
        let stats = px.fill_stats().unwrap();
        let from = h.fill_stats().unwrap();
        assert_eq!(stats.entries, 5.0);
        assert_eq!(stats.sum_wx, [from.sum_wx[0]]);
        assert_eq!(stats.sum_ws, from.sum_wx[1]);
        assert_eq!(stats.sum_ws2, from.sum_wxy[1][1]);
        assert_eq!(px.mean(), 1.1);
    }

    #[test]
    fn profile_of_weighted_bins_equals_a_weighted_profile() {
        let mut h2 = Hist::reg(2, 0.0, 2.0).reg(2, 0.0, 2.0).weight();
        h2.fill_weight(0.5, 0.5, 2.0);
        h2.fill_weight(0.5, 1.5, 2.0);
        h2.fill_weight(0.5, 7.0, 5.0); // the y overflow is left out
        let mut direct = Hist::reg(2, 0.0, 2.0).weighted_mean();
        direct.fill_sample_weight(0.5, 0.5, 2.0);
        direct.fill_sample_weight(0.5, 1.5, 2.0);
        let p = h2.profile(0).unwrap();
        assert_eq!(p, direct);
        // Σw = 4 and Σw² = 8: 2 effective entries, not 4.
        assert_eq!(p.bin([0]).unwrap().effective_count(), 2.0);

        // An axis without flow bins profiles to one without.
        let mut n = Hist::reg(2, 0.0, 2.0).no_flow().reg(2, 0.0, 2.0).double();
        n.fill(1.5, 1.5);
        assert_eq!(n.profile(0).unwrap().values(true), vec![0.0, 1.5]);
    }
}
