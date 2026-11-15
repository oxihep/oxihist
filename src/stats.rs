//! Statistics of 1-D histograms, and drawing random values from them.

use crate::axis::Axis;
use crate::histogram::Histogram;
use crate::storage::{CountStorage, Storage};

/// A small, seedable pseudo-random generator (SplitMix64) for reproducible
/// sampling without a `rand` dependency.
///
/// ```
/// use oxihist::Random;
/// let mut rng = Random::seed(5);
/// let x = rng.sample_binned(&[0.0, 1.0, 0.0], &[0.0, 1.0, 2.0, 3.0]).unwrap();
/// assert!((1.0..2.0).contains(&x)); // only the middle bin has weight
/// assert_eq!(rng.sample_binned(&[1.0], &[0.0]), None);
/// ```
#[doc(alias = "Rng", alias = "SplitMix64")]
#[derive(Debug, Clone)]
pub struct Random {
    state: u64,
}

impl Random {
    /// A generator seeded with `seed` (any value, including 0).
    #[must_use]
    pub fn seed(seed: u64) -> Random {
        Random { state: seed }
    }

    /// The next uniform `f64` in `[0, 1)`.
    pub fn uniform(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        // Top 53 bits → a double in [0, 1).
        (z >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Draw one value from a binned density: bin `k` spans
    /// `[edges[k], edges[k + 1])` and is picked with probability proportional to
    /// `weights[k]`, and the value is interpolated within it. Makes exactly one
    /// [`uniform`](Random::uniform) draw.
    ///
    /// Returns `None`, without drawing, when `edges.len() != weights.len() + 1`
    /// or the weights do not sum to a positive, finite total.
    pub fn sample_binned(&mut self, weights: &[f64], edges: &[f64]) -> Option<f64> {
        if edges.len() != weights.len() + 1 {
            return None;
        }
        let cdf = build_cdf(weights)?;
        Some(sample_cdf(&cdf, edges, self.uniform()))
    }
}

impl Default for Random {
    fn default() -> Random {
        Random::seed(0)
    }
}

/// The normalized cumulative of `weights`: `cdf[0] = 0`,
/// `cdf[k] = Σ weights[..k] / total`. `None` if the total is not positive and
/// finite.
fn build_cdf(weights: &[f64]) -> Option<Vec<f64>> {
    let mut cdf = vec![0.0; weights.len() + 1];
    let mut acc = 0.0;
    for (i, &w) in weights.iter().enumerate() {
        acc += w;
        cdf[i + 1] = acc;
    }
    if !(acc.is_finite() && acc > 0.0) {
        return None;
    }
    for c in &mut cdf {
        *c /= acc;
    }
    Some(cdf)
}

/// Inverse-transform sample: map `r ∈ [0, 1)` through the normalized `cdf`,
/// interpolating linearly within the chosen bin.
fn sample_cdf(cdf: &[f64], edges: &[f64], r: f64) -> f64 {
    let nbins = cdf.len() - 1;
    // The first k in 1..=nbins with cdf[k] > r.
    let k = 1 + cdf[1..nbins].partition_point(|&c| c <= r);
    let (c0, c1) = (cdf[k - 1], cdf[k]);
    let low = edges[k - 1];
    let width = edges[k] - edges[k - 1];
    if c1 > c0 {
        low + (r - c0) / (c1 - c0) * width
    } else {
        low + 0.5 * width
    }
}

/// Statistics over the in-range bins of a 1-D histogram, treating each bin as
/// concentrated at its center. On a categorical axis the centers are the bin
/// positions `0.5, 1.5, …`.
impl<S: Storage> Histogram<S, 1> {
    fn centers_values(&self) -> (Vec<f64>, Vec<f64>) {
        (self.axis(0).centers(), self.values(false))
    }

    /// The sum of the in-range values (`flow = true` adds the flow bins).
    #[must_use]
    pub fn integral(&self, flow: bool) -> f64 {
        self.values(flow).iter().sum()
    }

    /// The value-weighted mean of the bin centers (`NaN` for an empty
    /// histogram).
    #[must_use]
    pub fn mean(&self) -> f64 {
        let (c, v) = self.centers_values();
        let sw: f64 = v.iter().sum();
        c.iter().zip(&v).map(|(c, w)| c * w).sum::<f64>() / sw
    }

    /// The value-weighted standard deviation of the bin centers (population
    /// form; `NaN` for an empty histogram).
    #[must_use]
    pub fn std(&self) -> f64 {
        let (c, v) = self.centers_values();
        let sw: f64 = v.iter().sum();
        let m = c.iter().zip(&v).map(|(c, w)| c * w).sum::<f64>() / sw;
        let var = c
            .iter()
            .zip(&v)
            .map(|(c, w)| w * (c - m) * (c - m))
            .sum::<f64>()
            / sw;
        var.max(0.0).sqrt()
    }

    /// The `x` where the cumulative distribution of the in-range values reaches
    /// each probability in `probs`, interpolating linearly within a bin (as
    /// ROOT's `GetQuantiles` does). A probability landing exactly on a bin
    /// boundary returns that bin's center. An empty histogram gives its low
    /// edge for every probability.
    #[must_use]
    pub fn quantiles(&self, probs: &[f64]) -> Vec<f64> {
        let edges = self.axis(0).plot_edges();
        let v = self.values(false);
        let n = v.len();
        let mut cum = vec![0.0; n + 1];
        for i in 0..n {
            cum[i + 1] = cum[i] + v[i];
        }
        let total = cum[n];
        if total <= 0.0 {
            return vec![edges[0]; probs.len()];
        }
        cum.iter_mut().for_each(|c| *c /= total);
        probs
            .iter()
            .map(|&p| {
                let i = cum[..n].iter().rposition(|&c| c <= p).unwrap_or(0);
                if i > 0 && cum[i] == p {
                    return 0.5 * (edges[i - 1] + edges[i]);
                }
                let d = cum[i + 1] - cum[i];
                if d > 0.0 {
                    edges[i] + (p - cum[i]) / d * (edges[i + 1] - edges[i])
                } else {
                    edges[i]
                }
            })
            .collect()
    }

    /// The running sum of the in-range values, from the left (`forward`) or the
    /// right, as a histogram on the same axis; flow bins are left empty.
    #[must_use]
    pub fn cumulative(&self, forward: bool) -> Self {
        let mut out = Histogram::new([self.axis(0).clone()]);
        out.name.clone_from(&self.name);
        out.title.clone_from(&self.title);
        let n = self.axis(0).size() as isize;
        let order: Vec<isize> = if forward {
            (0..n).collect()
        } else {
            (0..n).rev().collect()
        };
        let mut acc = S::Bin::default();
        for i in order {
            S::merge(&mut acc, self.bin([i]).expect("in range"));
            *out.bin_mut([i]).expect("in range") = acc.clone();
        }
        out
    }

    /// Draw one value distributed like the in-range values (each bin uniform
    /// across its width). `None` for a histogram with no positive total.
    pub fn sample(&self, rng: &mut Random) -> Option<f64> {
        rng.sample_binned(&self.values(false), &self.axis(0).plot_edges())
    }
}

impl<S: CountStorage> Histogram<S, 1> {
    /// Fill `n` values drawn from `source`'s distribution (ROOT's
    /// `FillRandom`). Does nothing if `source` has no positive total.
    pub fn fill_random<T: Storage>(
        &mut self,
        source: &Histogram<T, 1>,
        n: usize,
        rng: &mut Random,
    ) {
        let edges = source.axis(0).plot_edges();
        let Some(cdf) = build_cdf(&source.values(false)) else {
            return;
        };
        for _ in 0..n {
            self.fill(sample_cdf(&cdf, &edges, rng.uniform()));
        }
    }

    /// Fill `n` values drawn from the density `f` over this histogram's own
    /// (numeric) axis: `f` is evaluated at each bin center, and a value is
    /// drawn from the resulting binned density. Does nothing for a categorical
    /// axis or an `f` with no positive total there.
    pub fn fill_random_fn(&mut self, f: impl Fn(f64) -> f64, n: usize, rng: &mut Random) {
        let axis: &Axis = self.axis(0);
        let Some(edges) = axis.edges() else {
            return;
        };
        let weights: Vec<f64> = axis.centers().into_iter().map(|c| f(c).max(0.0)).collect();
        let Some(cdf) = build_cdf(&weights) else {
            return;
        };
        for _ in 0..n {
            self.fill(sample_cdf(&cdf, &edges, rng.uniform()));
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{Hist, Random};

    #[test]
    fn mean_std_integral() {
        let mut h = Hist::reg(4, 0.0, 4.0).double();
        h.fill_many([0.5, 1.5, 1.5, 2.5]);
        h.fill(10.0);
        assert_eq!(h.integral(false), 4.0);
        assert_eq!(h.integral(true), 5.0);
        assert!((h.mean() - 1.5).abs() < 1e-12);
        assert!((h.std() - 0.5f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn quantiles_interpolate_and_hit_centers() {
        let mut h = Hist::reg(4, 0.0, 4.0).double();
        h.fill_many([0.5, 1.5, 2.5, 3.5]);
        let q = h.quantiles(&[0.0, 0.25, 0.5, 0.625, 1.0]);
        assert_eq!(q, vec![0.0, 0.5, 1.5, 2.5, 4.0]);
    }

    #[test]
    fn cumulative_both_ways() {
        let mut h = Hist::reg(3, 0.0, 3.0).weight();
        h.fill_many([0.5, 1.5, 1.5, 2.5]);
        assert_eq!(h.cumulative(true).values(false), vec![1.0, 3.0, 4.0]);
        assert_eq!(h.cumulative(false).values(false), vec![4.0, 3.0, 1.0]);
    }

    #[test]
    fn sampling_follows_the_source_and_is_seeded() {
        let mut src = Hist::reg(2, 0.0, 2.0).double();
        src.fill_weight(1.5, 3.0);
        src.fill_weight(0.5, 1.0);
        let mut a = Hist::reg(2, 0.0, 2.0).int64();
        a.fill_random(&src, 4000, &mut Random::seed(1));
        let frac = a.value([1]).unwrap() / 4000.0;
        assert!((frac - 0.75).abs() < 0.03, "{frac}");
        let mut b = Hist::reg(2, 0.0, 2.0).int64();
        b.fill_random(&src, 4000, &mut Random::seed(1));
        assert_eq!(a, b);
        let mut c = Hist::reg(10, 0.0, 1.0).double();
        c.fill_random_fn(|x| x, 1000, &mut Random::seed(3));
        assert_eq!(c.integral(false), 1000.0);
        assert!(c.value([9]) > c.value([0]));
    }
}
