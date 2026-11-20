//! Jet energy response and resolution versus jet pT, with a profile.
//!
//! For each simulated jet, the detector measures `pT_reco` for a true
//! `pT_true`. The response `R = pT_reco / pT_true` is below 1 (the
//! calorimeter is non-compensating) and its spread shrinks with energy as
//! `σ/E = a/√E ⊕ b`. A profile (`Mean` storage) over `pT_true` gives the mean
//! response per bin — the jet energy scale correction is its inverse — and
//! each bin's spread is the resolution. Run with:
//!
//! ```sh
//! cargo run --example jet_response
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
    fn gauss(&mut self, mean: f64, sigma: f64) -> f64 {
        let (u1, u2) = (self.uniform(), self.uniform());
        mean + sigma * (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

fn main() {
    let edges = [
        20.0, 30.0, 45.0, 70.0, 100.0, 150.0, 220.0, 320.0, 500.0, 800.0,
    ];
    let mut profile = Hist::var(&edges).label("jet $p_T^{true}$ [GeV]").mean();
    let mut toy = Toy(11);

    for _ in 0..300_000 {
        // A falling jet spectrum, 1/pT² above 20 GeV.
        let pt_true = 20.0 / (1.0 - toy.uniform() * (1.0 - 20.0 / 800.0));
        // Response: a scale rising from ~0.85 to ~0.97, and a resolution of
        // 100 %/√E ⊕ 5 %.
        let scale = 0.98 - 0.4 / pt_true.sqrt();
        let resolution = (1.0 / pt_true + 0.05 * 0.05).sqrt();
        let pt_reco = pt_true * toy.gauss(scale, resolution);
        profile.fill_sample(pt_true, pt_reco / pt_true);
    }

    println!(
        "{:>14} {:>8} {:>16} {:>14} {:>14}",
        "pT_true [GeV]", "jets", "response ⟨R⟩", "resolution σ_R", "JES correction"
    );
    for i in 0..edges.len() - 1 {
        let bin = profile.bin([i.cast_signed()]).expect("in range");
        let err = profile.variance([i.cast_signed()]).unwrap().sqrt();
        let spread = bin.sample_variance().sqrt();
        println!(
            "{:>6.0} – {:<5.0} {:>8.0} {:>8.4} ± {:.4} {:>14.3} {:>14.3}",
            edges[i],
            edges[i + 1],
            bin.count,
            bin.value,
            err,
            spread / bin.value,
            1.0 / bin.value
        );
    }
    println!("\nThe relative resolution falls with pT towards the 5 % constant term;");
    println!("the correction factor 1/⟨R⟩ is what a jet energy scale calibration applies.");
}
