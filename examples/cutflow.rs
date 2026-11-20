//! A weighted cutflow, by channel, with category axes.
//!
//! An analysis applies its selection step by step and reports how many
//! (weighted) events survive each step, per channel. A 2-D histogram over a
//! string category for the step and a growing string category for the
//! channel holds the whole table: no bookkeeping by hand, and it merges,
//! scales and converts like any other histogram. Run with:
//!
//! ```sh
//! cargo run --example cutflow
//! ```

use oxihist::Hist;

struct Toy(u64);

impl Toy {
    fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
}

/// A simulated event: its channel and the quantities the cuts look at.
struct Event {
    channel: &'static str,
    weight: f64,
    triggered: bool,
    n_leptons: u32,
    opposite_sign: bool,
    mass: f64,
    met: f64,
}

fn main() {
    let steps = [
        "all",
        "trigger",
        "2 leptons",
        "opposite sign",
        "Z window",
        "MET < 40",
    ];
    let mut cutflow = Hist::str_cat(&steps)
        .label("selection step")
        .str_cat::<&str>(&[])
        .growth()
        .label("channel")
        .weight();

    let mut toy = Toy(5);
    for _ in 0..100_000 {
        let r = toy.uniform();
        let channel = if r < 0.45 {
            "ee"
        } else if r < 0.9 {
            "μμ"
        } else {
            "eμ"
        };
        let e = Event {
            channel,
            // Cross-section weight per event: the eμ sample is a different
            // (top-like) process, simulated with more events per fb⁻¹.
            weight: if channel == "eμ" { 0.2 } else { 1.0 },
            triggered: toy.uniform() < if channel == "ee" { 0.92 } else { 0.95 },
            n_leptons: if toy.uniform() < 0.8 { 2 } else { 1 },
            // Z → ℓℓ is opposite-sign except for charge mis-measurement; the
            // eμ (top-like) pairs are uncorrelated.
            opposite_sign: toy.uniform() < if channel == "eμ" { 0.5 } else { 0.995 },
            // A Breit–Wigner Z peak (with tails outside the window), and a
            // flat continuum for eμ.
            mass: if channel == "eμ" {
                40.0 + 120.0 * toy.uniform()
            } else {
                91.19 + 1.25 * (std::f64::consts::PI * (toy.uniform() - 0.5)).tan()
            },
            // Missing ET: small for Z events (resolution only), large with
            // the neutrinos of a top-like decay.
            met: -(if channel == "eμ" { 60.0 } else { 12.0 }) * toy.uniform().ln(),
        };
        // Each step is filled only if every earlier one passed.
        let passed = [
            true,
            e.triggered,
            e.n_leptons == 2,
            e.opposite_sign,
            (76.0..106.0).contains(&e.mass),
            e.met < 40.0,
        ];
        for (step, _) in steps.iter().zip(passed).take_while(|(_, ok)| *ok) {
            cutflow.fill_weight(*step, e.channel, e.weight);
        }
    }

    let channels: Vec<String> = (0..cutflow.axis(1).size())
        .filter_map(|i| cutflow.axis(1).bin_label(i))
        .collect();
    print!("{:<15}", "step");
    for c in &channels {
        print!(" {:>10} {:>7}", c, "eff");
    }
    println!();
    for (i, step) in steps.iter().enumerate() {
        print!("{step:<15}");
        for j in 0..channels.len() {
            let n = cutflow.value([i.cast_signed(), j.cast_signed()]).unwrap();
            let first = cutflow.value([0, j.cast_signed()]).unwrap();
            print!(" {n:>10.1} {:>6.1}%", 100.0 * n / first);
        }
        println!();
    }
    println!("\nThe Z window removes most of the eμ (top-like) sample and keeps the Z channels.");
}
