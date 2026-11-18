//! A detector occupancy map in (η, φ), to spot dead and hot regions.
//!
//! Tracks hit a barrel detector uniformly in φ and roughly uniformly in η.
//! Fill a 2-D η × φ histogram, then compare each cell with the expectation
//! from its η row: a cell far below is dead (here, a module that is off), a
//! cell far above is hot (here, a noisy channel). Projections give the η and φ
//! distributions on their own. Run with:
//!
//! ```sh
//! cargo run --example occupancy
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

fn main() {
    let pi = std::f64::consts::PI;
    let (n_eta, n_phi) = (24, 32);
    let mut occupancy = Hist::reg(n_eta, -2.4, 2.4)
        .label("$\\eta$")
        .reg(n_phi, -pi, pi)
        .label("$\\phi$")
        .int64();
    let mut toy = Toy(3);

    for _ in 0..500_000 {
        let eta = -2.4 + 4.8 * toy.uniform();
        let phi = -pi + 2.0 * pi * toy.uniform();
        // A dead module: no hits in 0.4 < η < 0.8, 1.0 < φ < 1.4.
        if (0.4..0.8).contains(&eta) && (1.0..1.4).contains(&phi) {
            continue;
        }
        occupancy.fill(eta, phi);
    }
    // A noisy channel fires on its own.
    for _ in 0..3000 {
        occupancy.fill(-1.5, -2.0);
    }

    // Expectation per cell: the median of its η row, which a single hot or
    // dead cell cannot drag (a mean would). Flag cells beyond 5σ.
    println!("cells more than 5σ from their row's median:");
    let grid = occupancy.values_grid(); // [η][φ]
    for (ie, row) in grid.iter().enumerate() {
        let mut sorted = row.clone();
        sorted.sort_by(f64::total_cmp);
        let mean = sorted[sorted.len() / 2];
        for (ip, &n) in row.iter().enumerate() {
            let pull = (n - mean) / mean.sqrt();
            if pull.abs() > 5.0 {
                let eta = -2.4 + (ie as f64 + 0.5) * 0.2;
                let phi = -pi + (ip as f64 + 0.5) * 2.0 * pi / n_phi as f64;
                let kind = if pull < 0.0 { "dead" } else { "hot" };
                println!(
                    "  η = {eta:+.1}, φ = {phi:+.2}: {n:>6.0} hits (row median {mean:.0}, {pull:+.0}σ) → {kind}"
                );
            }
        }
    }

    // Project onto η and rebin by 4: the dead module shows as a dip there too.
    let eta = occupancy
        .project([0])
        .expect("axis 0 exists")
        .rebin(0, 4)
        .expect("24 bins divide by 4");
    println!(
        "\nhits per η slice (rebinned to {} bins):",
        eta.axis(0).size()
    );
    for (i, n) in eta.values(false).iter().enumerate() {
        let lo = -2.4 + i as f64 * 0.8;
        println!("  [{lo:+.1}, {:+.1})  {n:>7.0}", lo + 0.8);
    }
}
