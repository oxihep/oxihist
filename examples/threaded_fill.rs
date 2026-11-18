//! Fill one histogram from many threads, then merge.
//!
//! Event loops parallelise well, but a shared histogram would serialise them
//! on a lock. `ThreadedHist` gives each thread its own copy (found without a
//! shared lock) and adds them up at the end, ROOT's `TThreadedObject`
//! pattern. The result equals a single-threaded fill. Run with:
//!
//! ```sh
//! cargo run --release --example threaded_fill
//! ```

use std::time::Instant;

use oxihist::{Hist, ThreadedHist};

/// A per-event "measurement": a deterministic function of the event number,
/// with enough arithmetic to make the loop the bottleneck, as reconstruction
/// would be.
fn reconstruct(event: u64) -> f64 {
    let mut s = event.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let mut x = 0.0;
    for _ in 0..12 {
        s ^= s >> 31;
        s = s.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        x += (s >> 11) as f64 / (1u64 << 53) as f64;
    }
    91.0 + 4.0 * (x - 6.0) // ~Gaussian around 91
}

fn main() {
    let n_events: u64 = 4_000_000;
    let template = Hist::reg(100, 70.0, 110.0).label("$m$ [GeV]").double();

    let start = Instant::now();
    let mut serial = template.clone();
    for e in 0..n_events {
        serial.fill(reconstruct(e));
    }
    let t_serial = start.elapsed();

    let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
    let start = Instant::now();
    let shared = ThreadedHist::new(template);
    std::thread::scope(|scope| {
        for t in 0..threads as u64 {
            let shared = &shared;
            scope.spawn(move || {
                // Each thread takes every `threads`-th event.
                let mut e = t;
                while e < n_events {
                    shared.fill(reconstruct(e));
                    e += threads as u64;
                }
            });
        }
    });
    let merged = shared.merge().expect("copies of one template merge");
    let t_threaded = start.elapsed();

    println!("{n_events} events");
    println!("  1 thread : {t_serial:>10.2?}");
    println!(
        "  {threads} threads: {t_threaded:>10.2?}  ({:.1}× faster)",
        t_serial.as_secs_f64() / t_threaded.as_secs_f64()
    );
    println!(
        "same histogram either way: {} (entries {} vs {}, mean {:.3} vs {:.3})",
        merged == serial,
        merged.integral(true),
        serial.integral(true),
        merged.mean(),
        serial.mean()
    );
}
