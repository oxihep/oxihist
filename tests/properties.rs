//! Invariants that must hold for any data.

use oxihist::{Chi2TestKind, FillStats, Hist, Histogram, Weight};
use proptest::prelude::*;

fn filled(xs: &[(f64, f64, f64)]) -> Histogram<Weight, 2> {
    let mut h = Hist::reg(6, -3.0, 3.0).var(&[-2.0, 0.0, 0.5, 2.0]).weight();
    for &(x, y, w) in xs {
        h.fill_weight(x, y, w);
    }
    h
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

fn entries() -> impl Strategy<Value = Vec<(f64, f64, f64)>> {
    prop::collection::vec((-4.0..4.0f64, -3.0..3.0f64, 0.1..5.0f64), 0..200)
}

proptest! {
    #[test]
    fn splitting_and_merging_equals_one_fill(xs in entries(), cut in 0usize..200) {
        let cut = cut.min(xs.len());
        let mut a = filled(&xs[..cut]);
        a.add(&filled(&xs[cut..])).unwrap();
        let whole = filled(&xs);
        for (p, q) in a.bins().iter().zip(whole.bins()) {
            prop_assert!(close(p.value, q.value) && close(p.variance, q.variance));
        }
    }

    #[test]
    fn projections_and_rebinning_keep_the_total(xs in entries()) {
        let h = filled(&xs);
        let total: f64 = xs.iter().map(|e| e.2).sum();
        prop_assert!(close(h.sum(true).value, total));
        prop_assert!(close(h.project([0]).unwrap().sum(true).value, total));
        prop_assert!(close(h.project([1]).unwrap().sum(true).value, total));
        prop_assert!(close(h.rebin(0, 3).unwrap().sum(true).value, total));
        let swapped = h.project([1, 0]).unwrap();
        prop_assert_eq!(swapped.value([1, 2]), h.value([2, 1]));
    }

    #[test]
    fn scaling_scales_value_and_variance(xs in entries(), c in -3.0..3.0f64) {
        let mut h = filled(&xs);
        let before = h.sum(true);
        h.scale(c);
        let after = h.sum(true);
        prop_assert!(close(after.value, c * before.value));
        prop_assert!(close(after.variance, c * c * before.variance));
    }
}

/// The fill statistics of `xs` by definition: every fill an entry, the
/// in-range ones in the sums.
fn stats_of(xs: &[(f64, f64, f64)]) -> FillStats<2> {
    let mut s = FillStats::<2>::default();
    for &(x, y, w) in xs {
        s.entries += 1.0;
        if (-3.0..3.0).contains(&x) && (-2.0..2.0).contains(&y) {
            s.sum_w += w;
            s.sum_w2 += w * w;
            s.sum_wx[0] += w * x;
            s.sum_wx[1] += w * y;
            s.sum_wxy[0][0] += w * x * x;
            s.sum_wxy[1][1] += w * y * y;
            s.sum_wxy[0][1] += w * x * y;
            s.sum_wxy[1][0] += w * x * y;
        }
    }
    s
}

fn stats_close(a: &FillStats<2>, b: &FillStats<2>) -> bool {
    a.entries == b.entries
        && close(a.sum_w, b.sum_w)
        && close(a.sum_w2, b.sum_w2)
        && (0..2).all(|i| {
            close(a.sum_wx[i], b.sum_wx[i])
                && (0..2).all(|j| close(a.sum_wxy[i][j], b.sum_wxy[i][j]))
        })
}

fn counts() -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(0.0..50.0f64, 3..30)
}

fn from_counts(values: &[f64]) -> Histogram<Weight, 1> {
    let mut h = Hist::reg(values.len(), 0.0, values.len() as f64).weight();
    for (i, &v) in values.iter().enumerate() {
        h.fill_weight(i as f64 + 0.5, v);
    }
    h
}

proptest! {
    #[test]
    fn fill_stats_are_the_sums_of_the_fills(xs in entries()) {
        let h = filled(&xs);
        let want = stats_of(&xs);
        let got = h.fill_stats().unwrap();
        prop_assert!(stats_close(got, &want), "{got:?} != {want:?}");
        prop_assert_eq!(got.sum_wxy[0][1], got.sum_wxy[1][0]);
        prop_assert_eq!(h.entries(), xs.len() as f64);
        // The bins agree with the fills on the weights, not on where they
        // fell inside a bin.
        let binned = h.binned_stats();
        prop_assert!(close(binned.sum_w, want.sum_w) && close(binned.sum_w2, want.sum_w2));
    }

    #[test]
    fn fill_stats_follow_add_scale_project_and_rebin(
        xs in entries(),
        cut in 0usize..200,
        c in -3.0..3.0f64,
    ) {
        let cut = cut.min(xs.len());
        let whole = filled(&xs);
        let stats = *whole.fill_stats().unwrap();

        let mut a = filled(&xs[..cut]);
        a.add(&filled(&xs[cut..])).unwrap();
        prop_assert!(stats_close(a.fill_stats().unwrap(), &stats));

        let mut scaled = whole.clone();
        scaled.scale(c);
        let s = scaled.fill_stats().unwrap();
        prop_assert_eq!(s.entries, stats.entries);
        prop_assert!(close(s.sum_w, c * stats.sum_w) && close(s.sum_w2, c * c * stats.sum_w2));
        prop_assert!(close(s.sum_wx[1], c * stats.sum_wx[1]));
        prop_assert!(close(s.sum_wxy[0][1], c * stats.sum_wxy[0][1]));

        let py = whole.project([1]).unwrap();
        let p = py.fill_stats().unwrap();
        prop_assert_eq!((p.entries, p.sum_w, p.sum_w2), (stats.entries, stats.sum_w, stats.sum_w2));
        prop_assert_eq!((p.sum_wx, p.sum_wxy), ([stats.sum_wx[1]], [[stats.sum_wxy[1][1]]]));
        if stats.sum_w > 0.0 {
            prop_assert!(close(py.mean(), stats.sum_wx[1] / stats.sum_w));
            // An exact mean lies inside the axis; so does the spread.
            prop_assert!((-2.0..2.0).contains(&py.mean()));
            prop_assert!(py.std_dev() <= 4.0);
        }

        let rebinned = whole.rebin(0, 3).unwrap();
        prop_assert_eq!(rebinned.fill_stats(), whole.fill_stats());

        let mut emptied = whole.clone();
        emptied.reset();
        prop_assert_eq!(emptied.fill_stats(), Some(&FillStats::default()));
    }

    /// Without fill statistics the mean is that of the bin centers, which is
    /// within half a bin of the exact one.
    #[test]
    fn the_binned_mean_is_within_half_a_bin(xs in prop::collection::vec(0.0..10.0f64, 1..200)) {
        let mut h = Hist::reg(20, 0.0, 10.0).double();
        h.fill_many(xs.iter().copied());
        let exact = xs.iter().sum::<f64>() / xs.len() as f64;
        prop_assert!(close(h.mean(), exact));
        let mut rebuilt = Histogram::<oxihist::Double, 1>::from_bins([h.axis(0).clone()], h.bins()).unwrap();
        prop_assert_eq!(&rebuilt, &h);
        prop_assert!(rebuilt.fill_stats().is_none());
        prop_assert!((rebuilt.mean() - exact).abs() <= 0.25 + 1e-9);
        prop_assert_eq!(rebuilt.entries(), h.entries());
        rebuilt.set_fill_stats(h.fill_stats().copied());
        prop_assert_eq!(rebuilt.mean(), h.mean());
    }

    #[test]
    fn dividing_then_multiplying_gives_the_values_back(a in counts(), b in counts()) {
        let n = a.len().min(b.len());
        let (num, den) = (from_counts(&a[..n]), from_counts(&b[..n]));
        let mut r = num.clone();
        r.divide(&den).unwrap();
        prop_assert!(r.fill_stats().is_none());
        for i in 0..n {
            let (v, d) = (r.value([i.cast_signed()]).unwrap(), b[i]);
            if d == 0.0 {
                prop_assert_eq!(v, 0.0);
            } else {
                prop_assert!(close(v, a[i] / d));
            }
            prop_assert!(r.variance([i.cast_signed()]).unwrap() >= 0.0);
        }
        r.multiply(&den).unwrap();
        for i in 0..n {
            let want = if b[i] == 0.0 { 0.0 } else { a[i] };
            prop_assert!(close(r.value([i.cast_signed()]).unwrap(), want));
        }
    }

    #[test]
    fn a_histogram_is_compatible_with_itself_and_its_multiples(a in counts(), c in 0.1..10.0f64) {
        let h = from_counts(&a);
        let mut scaled = h.clone();
        scaled.scale(c);
        for kind in [
            Chi2TestKind::UnweightedUnweighted,
            Chi2TestKind::UnweightedWeighted,
            Chi2TestKind::WeightedWeighted,
        ] {
            let r = h.chi2_test(&h, kind).unwrap();
            prop_assert!(r.chi2.abs() < 1e-6, "{kind:?}: {}", r.chi2);
            prop_assert!(r.ndf < a.len());
        }
        // The weighted test compares shapes: a multiple is the same shape.
        let ww = h.chi2_test(&scaled, Chi2TestKind::WeightedWeighted).unwrap();
        prop_assert!(ww.chi2.abs() < 1e-6, "{}", ww.chi2);
        let ks = h.kolmogorov_test(&scaled).unwrap();
        prop_assert!(ks.distance < 1e-9);
        prop_assert!((0.0..=1.0).contains(&ks.prob));
    }

    #[test]
    fn smoothing_keeps_counts_non_negative_and_bounded(a in counts(), times in 1usize..4) {
        let mut h = from_counts(&a);
        let top = a.iter().copied().fold(0.0, f64::max);
        h.smooth(times);
        for v in h.values(false) {
            prop_assert!(v >= 0.0 && v.is_finite());
            // Medians and running means stay within the data; the twice
            // pass adds back at most the residuals.
            prop_assert!(v <= 3.0 * top + 1e-9, "{v} > 3 × {top}");
        }
        prop_assert_eq!(h.variances(true), from_counts(&a).variances(true));
        // A flat histogram is a fixed point.
        let mut flat = from_counts(&vec![7.0; a.len()]);
        flat.smooth(times);
        prop_assert!(flat.values(false).iter().all(|&v| close(v, 7.0)));
    }

    #[test]
    fn interpolation_stays_between_its_neighbours(a in counts(), t in 0.0..1.0f64) {
        let h = from_counts(&a);
        let n = a.len() as f64;
        let x = t * n;
        let y = h.interpolate(x).unwrap();
        let (lo, hi) = a.iter().fold((f64::INFINITY, 0.0f64), |(l, h), &v| (l.min(v), h.max(v)));
        prop_assert!(y >= lo - 1e-9 && y <= hi + 1e-9);
        // At a bin center it is the bin's value.
        for (i, &v) in a.iter().enumerate() {
            prop_assert!(close(h.interpolate(i as f64 + 0.5).unwrap(), v));
        }
        prop_assert_eq!(h.interpolate(-0.5 - t), None);
        prop_assert_eq!(h.interpolate(n + 0.5 + t), None);
    }

    #[test]
    fn a_slice_keeps_the_total_and_its_bins(a in counts(), lo in 0.0..30.0f64, len in 0.01..30.0f64) {
        let h = from_counts(&a);
        let n = a.len() as f64;
        let hi = lo + len;
        let Ok(s) = h.slice(lo..hi) else {
            prop_assert!(lo >= n);
            return Ok(());
        };
        prop_assert!(close(s.sum(true).value, h.sum(true).value));
        prop_assert!(close(s.sum(true).variance, h.sum(true).variance));
        // Every kept bin overlaps the range and holds what it held.
        let edges = s.axis(0).edges().unwrap();
        prop_assert!(edges[0] <= lo && lo < edges[1]);
        let last = edges.len() - 1;
        prop_assert!(edges[last - 1] < hi && (hi <= edges[last] || edges[last] == n));
        for (k, e) in edges[..last].iter().enumerate() {
            let i = h.axis(0).edges().unwrap().iter().position(|x| x == e).unwrap();
            prop_assert_eq!(s.bin([k.cast_signed()]), h.bin([i.cast_signed()]));
        }
    }

    #[test]
    fn the_extrema_bound_every_bin(xs in entries()) {
        let h = filled(&xs);
        let ((imax, max), (imin, min)) = (h.maximum().unwrap(), h.minimum().unwrap());
        prop_assert_eq!(h.value(imax), Some(max));
        prop_assert_eq!(h.value(imin), Some(min));
        prop_assert!(h.values(false).iter().all(|&v| min <= v && v <= max));
    }

    /// A profile of a 2-D histogram is the profile filled at the y bin
    /// centers.
    #[test]
    fn a_profile_matches_filling_at_bin_centers(xs in entries()) {
        let h = filled(&xs);
        let centers = h.axis(1).centers();
        let edges = h.axis(1).edges().unwrap();
        let mut direct = Hist::reg(6, -3.0, 3.0).weighted_mean();
        for &(x, y, w) in &xs {
            if let Some(k) = (0..3).find(|&k| edges[k] <= y && y < edges[k + 1]) {
                direct.fill_sample_weight(x, centers[k], w);
            }
        }
        let p = h.profile(0).unwrap();
        for (a, b) in p.bins().iter().zip(direct.bins()) {
            prop_assert!(close(a.sum_of_weights, b.sum_of_weights));
            prop_assert!(close(a.sum_of_weights_squared, b.sum_of_weights_squared));
            prop_assert!(close(a.value, b.value));
            prop_assert!((a.sum_of_weighted_deltas_squared - b.sum_of_weighted_deltas_squared).abs() < 1e-7);
        }
        let (s, from) = (p.fill_stats().unwrap(), h.fill_stats().unwrap());
        prop_assert_eq!((s.sum_wx[0], s.sum_ws, s.sum_ws2), (from.sum_wx[0], from.sum_wx[1], from.sum_wxy[1][1]));
    }
}

#[cfg(feature = "serde")]
#[test]
fn serde_round_trip() {
    let mut h = Hist::reg(3, 0.0, 3.0)
        .str_cat::<&str>(&[])
        .growth()
        .label("channel")
        .weight();
    h.fill_weight(1.5, "ee", 2.0);
    h.fill(0.5, "mumu");
    let h = h.named("h").titled("a title");
    let json = serde_json::to_string(&h).unwrap();
    let back: Histogram<Weight, 2> = serde_json::from_str(&json).unwrap();
    assert_eq!(back, h);
    assert_eq!(back.fill_stats(), h.fill_stats());
}

proptest! {
    /// A value below `hi` is never put in the overflow bin, and the last bin's
    /// upper edge belongs to the overflow.
    #[test]
    fn regular_axis_upper_edge(lo in -1e6..1e6f64, width in 1e-3..1e6f64, bins in 1usize..500) {
        let hi = lo + width;
        let mut h = Hist::reg(bins, lo, hi).int64();
        let below = hi.next_down(); // the largest f64 < hi
        if below >= lo {
            h.fill(below);
            prop_assert_eq!(h.value([bins.cast_signed() - 1]), Some(1.0));
        }
        h.fill(hi);
        prop_assert_eq!(h.value([bins.cast_signed()]), Some(1.0));
    }

    /// Filling copies of a growing-category template on different "threads",
    /// in any order, then merging, equals one serial fill.
    #[test]
    fn growing_categories_merge_like_a_serial_fill(
        events in prop::collection::vec((0usize..6, 0.0..4.0f64), 0..120),
        cut in 0usize..120,
    ) {
        let names = ["ee", "mumu", "emu", "tautau", "gg", "qq"];
        let template = Hist::str_cat::<&str>(&[]).growth().reg(4, 0.0, 4.0).weight();
        let fill = |evs: &[(usize, f64)]| {
            let mut h = template.clone();
            for &(c, x) in evs {
                h.fill(names[c], x);
            }
            h
        };
        let cut = cut.min(events.len());
        let (mut a, b) = (fill(&events[cut..]), fill(&events[..cut]));
        a.add(&b).unwrap();
        let serial = fill(&events);
        prop_assert_eq!(a.axis(0).size(), serial.axis(0).size());
        for c in 0..serial.axis(0).size() {
            let name = serial.axis(0).bin_label(c).unwrap();
            let ca = (0..a.axis(0).size()).find(|&k| a.axis(0).bin_label(k).as_deref() == Some(&name)).unwrap();
            for x in -1..=4 {
                prop_assert_eq!(a.value([ca.cast_signed(), x]), serial.value([c.cast_signed(), x]));
            }
        }
    }
}

#[test]
fn threaded_fill_with_growing_categories() {
    use oxihist::ThreadedHist;
    let acc = ThreadedHist::new(Hist::str_cat::<&str>(&[]).growth().int64());
    std::thread::scope(|s| {
        for t in 0..4 {
            let acc = &acc;
            s.spawn(move || {
                for i in 0..50 {
                    acc.fill(["a", "b", "c", "d"][(i + t) % 4]);
                }
            });
        }
    });
    let h = acc.merge().unwrap();
    assert_eq!(h.axis(0).size(), 4);
    assert_eq!(h.sum(true), 200);
}

#[cfg(feature = "rayon")]
#[test]
fn fill_par_with_growing_categories() {
    let data: Vec<usize> = (0..10_000).map(|i| i % 7).collect();
    let names = ["a", "b", "c", "d", "e", "f", "g"];
    let h = oxihist::fill_par(
        &Hist::str_cat::<&str>(&[]).growth().int64(),
        &data,
        |h, &i| {
            h.fill(names[i]);
        },
    );
    assert_eq!(h.axis(0).size(), 7);
    assert_eq!(h.sum(true), 10_000);
}

proptest! {
    /// Growing categories while filling gives the same histogram as declaring
    /// every category up front (in first-seen order), whatever the number of
    /// categories: below and above the hash-index threshold, and across
    /// several capacity doublings.
    #[test]
    fn growth_equals_categories_declared_up_front(
        events in prop::collection::vec((0usize..60, 0.0..3.0f64), 0..400),
    ) {
        let mut grown = Hist::str_cat::<&str>(&[]).growth().reg(3, 0.0, 3.0).weight();
        let names: Vec<String> = (0..60).map(|i| format!("c{i}")).collect();
        let mut seen: Vec<&str> = Vec::new();
        for &(c, x) in &events {
            grown.fill(names[c].as_str(), x);
            if !seen.contains(&names[c].as_str()) {
                seen.push(&names[c]);
            }
        }
        let mut declared = Hist::str_cat(&seen).growth().reg(3, 0.0, 3.0).weight();
        for &(c, x) in &events {
            declared.fill(names[c].as_str(), x);
        }
        prop_assert_eq!(&grown, &declared);
        prop_assert_eq!(grown.bins(), declared.bins());
        prop_assert_eq!(grown.values(true), declared.values(true));
    }
}
