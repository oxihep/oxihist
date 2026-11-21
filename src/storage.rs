//! What each bin accumulates.
//!
//! | Storage | Bin | Use |
//! |---|---|---|
//! | [`Double`] | `f64` | unweighted or weighted counts, no variance tracking |
//! | [`Int64`] | `i64` | exact integer counts |
//! | [`Weight`] | [`WeightedSum`] | weighted counts with their variance (`Σw`, `Σw²`) |
//! | [`Mean`] | [`MeanBin`] | a profile: the mean and spread of a sample per bin |
//! | [`WeightedMean`] | [`WeightedMeanBin`] | a weighted profile |
//!
//! The first three are [`CountStorage`]s (filled with a weight); the last two
//! are [`MeanStorage`]s (filled with a sample value and a weight).
//!
//! A bin knows nothing of where inside it an entry fell: the histogram keeps
//! that, for all bins together, in its [`FillStats`](crate::FillStats).

use crate::cast;
use std::fmt::Debug;

/// What a histogram's bins hold. The storage types themselves are empty
/// markers; [`Storage::Bin`] is the per-bin accumulator.
pub trait Storage: Clone + Debug + Default + Send + Sync + 'static {
    /// The per-bin accumulator.
    type Bin: Clone + Debug + Default + PartialEq + Send + Sync;

    /// The bin's value: the count (or sum of weights), or the mean for a
    /// profile.
    fn value(bin: &Self::Bin) -> f64;

    /// The variance of the bin's value. Counts without tracked weights use the
    /// Poisson estimate (the count itself); a profile reports the variance of
    /// its mean (the sample variance over the effective count).
    fn variance(bin: &Self::Bin) -> f64;

    /// Combine `b` into `a`, as if `a` had also been filled with everything
    /// `b` was.
    fn merge(a: &mut Self::Bin, b: &Self::Bin);

    /// Multiply the bin by `c` (a count's value by `c` and its variance by
    /// `c²`; a profile's mean by `c` and its spread by `c²`).
    fn scale(bin: &mut Self::Bin, c: f64);

    /// Whether a bin holds samples of a value (a profile) instead of a
    /// count. Scaling a profile scales the sampled values, not the weights,
    /// which decides what happens to a histogram's
    /// [`FillStats`](crate::FillStats).
    const SAMPLED: bool = false;

    /// What the bin adds to its histogram's statistics when they are
    /// estimated from the bins
    /// ([`Histogram::binned_stats`](crate::Histogram::binned_stats)). The
    /// default is a count's: the value as `Σw` and the variance as `Σw²`.
    fn sums(bin: &Self::Bin) -> BinSums {
        BinSums {
            sum_w: Self::value(bin),
            sum_w2: Self::variance(bin),
            sum_ws: 0.0,
            sum_ws2: 0.0,
        }
    }
}

/// The sums one bin holds, as far as its accumulator knows them: see
/// [`Storage::sums`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BinSums {
    /// `Σw` of the bin's entries.
    pub sum_w: f64,
    /// `Σw²` of the bin's entries.
    pub sum_w2: f64,
    /// `Σw·s` of a profile bin's samples (0 for a count).
    pub sum_ws: f64,
    /// `Σw·s²` of a profile bin's samples (0 for a count).
    pub sum_ws2: f64,
}

/// A storage filled with a weight per entry: [`Double`], [`Int64`], [`Weight`].
pub trait CountStorage: Storage {
    /// Add an entry of weight `w`.
    fn add(bin: &mut Self::Bin, w: f64);

    /// A bin with the given value and variance, as close as the storage can
    /// hold them: [`Double`] and [`Int64`] keep the value only (their
    /// variance is the value), and [`Int64`] rounds it.
    fn bin(value: f64, variance: f64) -> Self::Bin;
}

/// A storage filled with a sample value and a weight per entry: [`Mean`],
/// [`WeightedMean`].
pub trait MeanStorage: Storage {
    /// Add the sample `x` with weight `w`.
    fn add(bin: &mut Self::Bin, x: f64, w: f64);
    /// The (effective) number of entries in the bin: `Σw` for a weighted mean.
    fn count(bin: &Self::Bin) -> f64;
}

/// Plain `f64` counts, weighted or not; the variance is taken to be the count.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Double;

impl Storage for Double {
    type Bin = f64;
    fn value(bin: &f64) -> f64 {
        *bin
    }
    fn variance(bin: &f64) -> f64 {
        *bin
    }
    fn merge(a: &mut f64, b: &f64) {
        *a += b;
    }
    fn scale(bin: &mut f64, c: f64) {
        *bin *= c;
    }
}

impl CountStorage for Double {
    fn add(bin: &mut f64, w: f64) {
        *bin += w;
    }
    fn bin(value: f64, _variance: f64) -> f64 {
        value
    }
}

/// Exact `i64` counts. A weight is rounded to the nearest integer, and scaling
/// rounds the result.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Int64;

impl Storage for Int64 {
    type Bin = i64;
    fn value(bin: &i64) -> f64 {
        *bin as f64
    }
    fn variance(bin: &i64) -> f64 {
        *bin as f64
    }
    fn merge(a: &mut i64, b: &i64) {
        *a += b;
    }
    fn scale(bin: &mut i64, c: f64) {
        *bin = cast::round_i64(*bin as f64 * c);
    }
}

impl CountStorage for Int64 {
    fn add(bin: &mut i64, w: f64) {
        *bin += cast::round_i64(w);
    }
    fn bin(value: f64, _variance: f64) -> i64 {
        cast::round_i64(value)
    }
}

/// Weighted counts with their variance: each bin holds `Σw` and `Σw²`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Weight;

/// A [`Weight`] bin: the sum of weights and the sum of squared weights.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WeightedSum {
    /// `Σw`.
    pub value: f64,
    /// `Σw²`, the variance of `value`.
    pub variance: f64,
}

impl Storage for Weight {
    type Bin = WeightedSum;
    fn value(bin: &WeightedSum) -> f64 {
        bin.value
    }
    fn variance(bin: &WeightedSum) -> f64 {
        bin.variance
    }
    fn merge(a: &mut WeightedSum, b: &WeightedSum) {
        a.value += b.value;
        a.variance += b.variance;
    }
    fn scale(bin: &mut WeightedSum, c: f64) {
        bin.value *= c;
        bin.variance *= c * c;
    }
}

impl CountStorage for Weight {
    fn add(bin: &mut WeightedSum, w: f64) {
        bin.value += w;
        bin.variance += w * w;
    }
    fn bin(value: f64, variance: f64) -> WeightedSum {
        WeightedSum { value, variance }
    }
}

/// A profile: the count, mean and spread of the samples in each bin. A weight
/// counts as a frequency (an entry of weight 2 is two identical samples).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Mean;

/// A [`Mean`] bin, accumulated with Welford's algorithm.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MeanBin {
    /// The number of samples.
    pub count: f64,
    /// Their mean.
    pub value: f64,
    /// `Σ(x − mean)²`.
    pub sum_of_deltas_squared: f64,
}

impl MeanBin {
    /// The sample variance `Σ(x − mean)² / (n − 1)` (0 with fewer than two
    /// samples).
    #[must_use]
    pub fn sample_variance(&self) -> f64 {
        if self.count > 1.0 {
            self.sum_of_deltas_squared / (self.count - 1.0)
        } else {
            0.0
        }
    }
}

impl Storage for Mean {
    type Bin = MeanBin;
    fn value(bin: &MeanBin) -> f64 {
        bin.value
    }
    fn variance(bin: &MeanBin) -> f64 {
        if bin.count > 0.0 {
            bin.sample_variance() / bin.count
        } else {
            0.0
        }
    }
    fn merge(a: &mut MeanBin, b: &MeanBin) {
        let n = a.count + b.count;
        if n == 0.0 {
            return;
        }
        let delta = b.value - a.value;
        a.value += delta * b.count / n;
        a.sum_of_deltas_squared += b.sum_of_deltas_squared + delta * delta * a.count * b.count / n;
        a.count = n;
    }
    fn scale(bin: &mut MeanBin, c: f64) {
        bin.value *= c;
        bin.sum_of_deltas_squared *= c * c;
    }
    const SAMPLED: bool = true;
    fn sums(bin: &MeanBin) -> BinSums {
        // Frequency weights: without a record of Σw², take it to be Σw.
        BinSums {
            sum_w: bin.count,
            sum_w2: bin.count,
            sum_ws: bin.count * bin.value,
            sum_ws2: bin.sum_of_deltas_squared + bin.count * bin.value * bin.value,
        }
    }
}

impl MeanStorage for Mean {
    fn add(bin: &mut MeanBin, x: f64, w: f64) {
        bin.count += w;
        let delta = x - bin.value;
        bin.value += w * delta / bin.count;
        bin.sum_of_deltas_squared += w * delta * (x - bin.value);
    }
    fn count(bin: &MeanBin) -> f64 {
        bin.count
    }
}

/// A weighted profile: the weighted mean and spread of the samples in each
/// bin, with weights as reliability weights.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WeightedMean;

/// A [`WeightedMean`] bin.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WeightedMeanBin {
    /// `Σw`.
    pub sum_of_weights: f64,
    /// `Σw²`.
    pub sum_of_weights_squared: f64,
    /// The weighted mean.
    pub value: f64,
    /// `Σw(x − mean)²`.
    pub sum_of_weighted_deltas_squared: f64,
}

impl WeightedMeanBin {
    /// The effective number of entries `(Σw)² / Σw²`.
    #[must_use]
    pub fn effective_count(&self) -> f64 {
        if self.sum_of_weights_squared > 0.0 {
            self.sum_of_weights * self.sum_of_weights / self.sum_of_weights_squared
        } else {
            0.0
        }
    }

    /// The unbiased weighted sample variance
    /// `Σw(x − mean)² / (Σw − Σw²/Σw)` (0 when it is undefined).
    #[must_use]
    pub fn sample_variance(&self) -> f64 {
        if self.sum_of_weights <= 0.0 {
            return 0.0;
        }
        let denom = self.sum_of_weights - self.sum_of_weights_squared / self.sum_of_weights;
        if denom > 0.0 {
            self.sum_of_weighted_deltas_squared / denom
        } else {
            0.0
        }
    }
}

impl Storage for WeightedMean {
    type Bin = WeightedMeanBin;
    fn value(bin: &WeightedMeanBin) -> f64 {
        bin.value
    }
    fn variance(bin: &WeightedMeanBin) -> f64 {
        let n = bin.effective_count();
        if n > 0.0 {
            bin.sample_variance() / n
        } else {
            0.0
        }
    }
    fn merge(a: &mut WeightedMeanBin, b: &WeightedMeanBin) {
        let sw = a.sum_of_weights + b.sum_of_weights;
        if sw == 0.0 {
            return;
        }
        let delta = b.value - a.value;
        a.value += delta * b.sum_of_weights / sw;
        a.sum_of_weighted_deltas_squared += b.sum_of_weighted_deltas_squared
            + delta * delta * a.sum_of_weights * b.sum_of_weights / sw;
        a.sum_of_weights = sw;
        a.sum_of_weights_squared += b.sum_of_weights_squared;
    }
    fn scale(bin: &mut WeightedMeanBin, c: f64) {
        bin.value *= c;
        bin.sum_of_weighted_deltas_squared *= c * c;
    }
    const SAMPLED: bool = true;
    fn sums(bin: &WeightedMeanBin) -> BinSums {
        BinSums {
            sum_w: bin.sum_of_weights,
            sum_w2: bin.sum_of_weights_squared,
            sum_ws: bin.sum_of_weights * bin.value,
            sum_ws2: bin.sum_of_weighted_deltas_squared
                + bin.sum_of_weights * bin.value * bin.value,
        }
    }
}

impl MeanStorage for WeightedMean {
    fn add(bin: &mut WeightedMeanBin, x: f64, w: f64) {
        bin.sum_of_weights += w;
        bin.sum_of_weights_squared += w * w;
        if bin.sum_of_weights == 0.0 {
            return;
        }
        let delta = x - bin.value;
        bin.value += w * delta / bin.sum_of_weights;
        bin.sum_of_weighted_deltas_squared += w * delta * (x - bin.value);
    }
    fn count(bin: &WeightedMeanBin) -> f64 {
        bin.sum_of_weights
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12 * (1.0 + a.abs().max(b.abs()))
    }

    #[test]
    fn mean_matches_two_pass_and_merges() {
        let xs = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let mut all = MeanBin::default();
        for &x in &xs {
            Mean::add(&mut all, x, 1.0);
        }
        assert!(close(all.value, 5.0));
        assert!(close(all.sample_variance(), 32.0 / 7.0));
        let (mut a, mut b) = (MeanBin::default(), MeanBin::default());
        for &x in &xs[..3] {
            Mean::add(&mut a, x, 1.0);
        }
        for &x in &xs[3..] {
            Mean::add(&mut b, x, 1.0);
        }
        Mean::merge(&mut a, &b);
        assert!(
            close(a.value, all.value) && close(a.sum_of_deltas_squared, all.sum_of_deltas_squared)
        );
    }

    #[test]
    fn weighted_mean_with_unit_weights_is_the_mean() {
        let mut w = WeightedMeanBin::default();
        let mut m = MeanBin::default();
        for x in [1.0, 2.0, 6.0] {
            WeightedMean::add(&mut w, x, 1.0);
            Mean::add(&mut m, x, 1.0);
        }
        assert!(close(w.value, m.value));
        assert!(close(w.sample_variance(), m.sample_variance()));
        assert!(close(WeightedMean::variance(&w), Mean::variance(&m)));
    }

    #[test]
    fn weight_tracks_sum_of_squares() {
        let mut b = WeightedSum::default();
        Weight::add(&mut b, 2.0);
        Weight::add(&mut b, 3.0);
        assert_eq!((b.value, b.variance), (5.0, 13.0));
        Weight::scale(&mut b, 2.0);
        assert_eq!((b.value, b.variance), (10.0, 52.0));
    }
}
