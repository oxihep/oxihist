//! The dimuon invariant mass of Z → μμ events, from muon (pT, η, φ).
//!
//! A toy generator produces Z bosons with a spread of rapidity, decays them
//! isotropically, applies a CMS-like acceptance (`pT > 20 GeV`, `|η| < 2.4`)
//! and a 2 % pT resolution, and computes the pair mass the way an analysis
//! does from the reconstructed muons:
//!
//! ```text
//! m² = 2 pT₁ pT₂ [cosh(η₁ − η₂) − cos(φ₁ − φ₂)]      (massless muons)
//! ```
//!
//! Like an NLO generator, 15 % of events carry weight −1: a `Weight` storage
//! keeps `Σw` and `Σw²`, so the errors account for them. Run with:
//!
//! ```sh
//! cargo run --example dimuon_mass
//! ```

use oxihist::Hist;

/// A small reproducible generator.
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
    /// A Breit–Wigner (Cauchy) draw: the Z lineshape.
    fn breit_wigner(&mut self, mass: f64, width: f64) -> f64 {
        mass + 0.5 * width * (std::f64::consts::PI * (self.uniform() - 0.5)).tan()
    }
}

/// A muon as (pT, η, φ).
type Muon = (f64, f64, f64);

/// Decay a Z of mass `m` and rapidity `y` (no transverse momentum) into two
/// massless muons, isotropically in its rest frame.
fn decay(toy: &mut Toy, m: f64, y: f64) -> [Muon; 2] {
    let cos_theta = 2.0 * toy.uniform() - 1.0;
    let phi = std::f64::consts::TAU * toy.uniform();
    let p = m / 2.0;
    let pt = p * (1.0 - cos_theta * cos_theta).sqrt();
    // In the rest frame each muon has pseudorapidity atanh(cos θ); a boost
    // along z by rapidity y adds y to it (exact for massless particles).
    let eta = cos_theta.atanh();
    let wrap = |a: f64| {
        (a + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
    };
    [
        (pt, eta + y, wrap(phi)),
        (pt, -eta + y, wrap(phi + std::f64::consts::PI)),
    ]
}

fn mass(a: Muon, b: Muon) -> f64 {
    (2.0 * a.0 * b.0 * ((a.1 - b.1).cosh() - (a.2 - b.2).cos())).sqrt()
}

fn main() {
    let mut toy = Toy(91);
    let mut m_mumu = Hist::reg(60, 60.0, 120.0)
        .label("$m_{\\mu\\mu}$ [GeV]")
        .name("m_mumu")
        .weight();
    let mut leading_pt = Hist::reg(40, 0.0, 80.0)
        .label("leading muon $p_T$ [GeV]")
        .weight();

    let (mut generated, mut accepted) = (0u32, 0u32);
    for _ in 0..200_000 {
        generated += 1;
        let m_true = toy.breit_wigner(91.1876, 2.4952);
        let y = toy.gauss(0.0, 1.8);
        let mut muons = decay(&mut toy, m_true, y);
        // Reconstruction: smear pT by 2 %, then the acceptance cuts.
        for mu in &mut muons {
            mu.0 *= toy.gauss(1.0, 0.02);
        }
        if muons
            .iter()
            .any(|&(pt, eta, _)| pt < 20.0 || eta.abs() > 2.4)
        {
            continue;
        }
        accepted += 1;
        let weight = if toy.uniform() < 0.15 { -1.0 } else { 1.0 };
        m_mumu.fill_weight(mass(muons[0], muons[1]), weight);
        leading_pt.fill_weight(muons[0].0.max(muons[1].0), weight);
    }

    println!(
        "generated {generated} Z → μμ, accepted {accepted} ({:.1} %)",
        100.0 * f64::from(accepted) / f64::from(generated)
    );
    let total = m_mumu.sum(false);
    // Negative weights cost statistics: the effective number of events
    // (Σw)²/Σw² is what the yield is as precise as.
    let n_eff = total.value * total.value / total.variance;
    println!(
        "weighted yield in 60–120 GeV: {:.0} ± {:.0} ({:.2} %)",
        total.value,
        total.variance.sqrt(),
        100.0 * total.variance.sqrt() / total.value,
    );
    println!(
        "  the 15 % of negative weights leave {n_eff:.0} effective events ({:.0} % of the entries).",
        100.0 * n_eff / f64::from(accepted)
    );
    let q = m_mumu.quantiles(&[0.16, 0.5, 0.84]);
    println!(
        "mass: mean {:.2} GeV, median {:.2}, 68 % of events in [{:.2}, {:.2}] GeV",
        m_mumu.mean(),
        q[1],
        q[0],
        q[2]
    );
    // The Jacobian peak: a muon from a Z at rest carries up to m_Z/2.
    let peak_bin = leading_pt
        .values(false)
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .expect("bins");
    println!(
        "leading-muon pT peaks in bin [{}, {}) GeV: the Jacobian peak below m_Z/2 = 45.6 GeV",
        2 * peak_bin,
        2 * peak_bin + 2
    );

    println!("\n{:>10} {:>10} {:>8}", "m [GeV]", "Σw", "√Σw²");
    for i in (20..40).step_by(2) {
        println!(
            "{:>10.1} {:>10.0} {:>8.1}",
            60.0 + f64::from(i) + 0.5,
            m_mumu.value([i as isize]).unwrap(),
            m_mumu.variance([i as isize]).unwrap().sqrt()
        );
    }
}
