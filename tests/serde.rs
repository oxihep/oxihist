//! Deserialising goes through the constructors' checks: malformed data is an
//! error, never a histogram that misbehaves later.
#![cfg(feature = "serde")]

use oxihist::{Axis, Double, Hist, Histogram};
use serde_json::{Value, json};

/// A serialised 1-D `Double` histogram over `axis`, as JSON to tamper with.
fn histogram_json(axis: Axis) -> Value {
    let h: Histogram<Double, 1> = Histogram::new([axis]);
    serde_json::to_value(&h).unwrap()
}

fn load(v: Value) -> Result<Histogram<Double, 1>, String> {
    serde_json::from_value(v).map_err(|e| e.to_string())
}

fn load_axis(v: Value) -> Result<Axis, String> {
    serde_json::from_value(v).map_err(|e| e.to_string())
}

#[test]
fn a_valid_histogram_round_trips() {
    for axis in [
        Axis::regular(4, 0.0, 1.0),
        Axis::variable(&[0.0, 1.0, 5.0]),
        Axis::integer(-3, 3),
        Axis::str_category(&["a", "b"]),
        Axis::int_category(&[7, 9]),
        Axis::boolean(),
    ] {
        let v = histogram_json(axis.clone());
        let back = load(v).unwrap();
        assert_eq!(back.axes()[0], axis);
    }
}

#[test]
fn malformed_axes_are_refused() {
    let axis = |kind: Value| json!({ "kind": kind, "name": "", "label": "", "underflow": true, "overflow": true });
    for (kind, why) in [
        (
            json!({ "Regular": { "bins": 0, "lo": 0.0, "hi": 1.0 } }),
            "no bins",
        ),
        (
            json!({ "Regular": { "bins": 3, "lo": 1.0, "hi": 1.0 } }),
            "empty range",
        ),
        (
            json!({ "Regular": { "bins": 3, "lo": 2.0, "hi": 1.0 } }),
            "reversed",
        ),
        (
            json!({ "Variable": { "edges": [0.0, 2.0, 1.0] } }),
            "descending edges",
        ),
        (json!({ "Variable": { "edges": [0.0] } }), "one edge"),
        (
            json!({ "Integer": { "lo": 5, "hi": 5 } }),
            "empty integer range",
        ),
        (
            json!({ "Integer": { "lo": i64::MIN, "hi": i64::MAX } }),
            "a span that overflows",
        ),
    ] {
        assert!(load_axis(axis(kind)).is_err(), "accepted: {why}");
    }
    // The valid shapes of the same fields load.
    assert!(
        load_axis(axis(
            json!({ "Regular": { "bins": 3, "lo": 0.0, "hi": 1.0 } })
        ))
        .is_ok()
    );
    assert!(load_axis(axis(json!({ "Integer": { "lo": -2, "hi": 2 } }))).is_ok());
}

#[test]
fn storage_that_does_not_match_its_axes_is_refused() {
    let good = histogram_json(Axis::regular(3, 0.0, 1.0)); // 3 bins + 2 flow
    assert!(load(good.clone()).is_ok());

    let mut short = good.clone();
    short["bins"].as_array_mut().unwrap().pop();
    assert!(load(short).unwrap_err().contains("holds 4 bins"));

    let mut small_cap = good.clone();
    small_cap["cap"] = json!([4]);
    small_cap["bins"] = json!([0.0, 0.0, 0.0, 0.0]);
    assert!(load(small_cap).unwrap_err().contains("storage for 4"));

    let mut two_axes = good.clone();
    let a = two_axes["axes"][0].clone();
    two_axes["axes"].as_array_mut().unwrap().push(a);
    assert!(
        load(two_axes)
            .unwrap_err()
            .contains("1-D histogram needs 1 axes")
    );

    let mut huge = good;
    huge["cap"] = json!([usize::MAX]);
    assert!(load(huge).is_err(), "an overflowing storage size");
}

#[test]
fn a_loaded_category_axis_finds_its_categories() {
    // Past a few categories the axis keeps a hash index, which is not
    // serialised: loading rebuilds it, so fills land in the right bins.
    let cats: Vec<String> = (0..20).map(|i| format!("c{i}")).collect();
    let mut h = Hist::str_cat(&cats).double();
    h.fill("c13");
    let json = serde_json::to_string(&h).unwrap();
    let mut back: Histogram<Double, 1> = serde_json::from_str(&json).unwrap();
    back.fill("c17");
    assert_eq!(back.value([13]), Some(1.0));
    assert_eq!(back.value([17]), Some(1.0));
}
