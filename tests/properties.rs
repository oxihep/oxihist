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
