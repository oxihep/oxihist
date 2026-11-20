//! Invariants that must hold for any data.

use oxihist::{Hist, Histogram, Weight};
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
