# oxihist

[![CI](https://github.com/oxihep/oxihist/actions/workflows/ci.yml/badge.svg)](https://github.com/oxihep/oxihist/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Rust 1.100+](https://img.shields.io/badge/rust-1.100%2B-orange.svg)](https://www.rust-lang.org)

Histograms for high-energy physics in pure Rust, modelled on
[boost-histogram](https://github.com/scikit-hep/boost-histogram) and
scikit-hep [`hist`](https://github.com/scikit-hep/hist). No file format is
attached: [oxiroot](https://github.com/oxihep/oxiroot) reads and writes these
histograms as ROOT `TH1`/`TH2`/`TH3`/`TProfile` objects, the way uproot does
for `hist`.

[oxiplot](https://github.com/oxihep/oxiplot) draws them (the `plot` feature)
and [oxifit](https://github.com/oxihep/oxifit) fits them (the `fit` feature).

```toml
[dependencies]
oxihist = { git = "https://github.com/oxihep/oxihist", tag = "v0.4.3" }
```

```rust
use oxihist::Hist;

// hist:  Hist.new.Reg(50, 0, 100, label="$p_T$").Weight()
let mut h = Hist::reg(50, 0.0, 100.0).label("$p_T$ [GeV]").name("pt").weight();
h.fill(42.0);
h.fill_weight(17.0, 0.5);
assert_eq!(h.value([8]), Some(0.5));       // bins are numbered from 0
assert_eq!(h.variance([8]), Some(0.25));   // Σw²

// 2-D, with a string category that grows as new categories arrive.
let mut m = Hist::reg(40, 0.0, 200.0).str_cat::<&str>(&[]).growth().int64();
m.fill(91.2, "ee");
m.fill(125.0, "mumu");
let mass = m.project([0]).unwrap();       // sum over the channels
```

## What is in it

| | |
| --- | --- |
| **Axes** | `Regular`, `Variable`, `Integer`, `StrCategory`, `IntCategory` (both optionally growing), `Boolean`; under- and overflow bins on or off per axis |
| **Storages** | `Double`, `Int64`, `Weight` (`Σw`, `Σw²`), `Mean` and `WeightedMean` (profiles) |
| **Histograms** | `Histogram<S, D>` for 1, 2 and 3 dimensions: `fill`/`fill_weight`, `fill_sample` for profiles, bin access with index `-1` for underflow and `size` for overflow |
| **Operations** | `add`, `scale`, `project`, `rebin`, `sum`, `reset`; 1-D `mean`, `std`, `quantiles`, `integral`, `cumulative` |
| **Sampling** | a seedable `Random`, `sample`, `fill_random`, `fill_random_fn` |
| **Threads** | `ThreadedHist` (one copy per thread, found without a shared lock, merged at the end) and `fill_par` with the `rayon` feature; both work for any `Mergeable` type, so oxiroot uses them for ROOT histograms too |
| **Speed** | a number on a regular axis fills at ~770 M/s, as fast as oxiroot's `Hist1D`; integers, known categories (a foldhash index) and booleans skip the general growing/converting path |

Bins are numbered boost-histogram's way: `0..size` in range, `-1` for the
underflow and `size` for the overflow bin.

## Features

| Feature | Adds |
| --- | --- |
| `rayon` | `fill_par`, a parallel fill from a slice |
| `serde` | `Serialize`/`Deserialize` for histograms, axes and bins |
| `plot` | oxiplot's data traits: `Axes::hist(&h)`, `Axes::profile(&p)`, `Axes::hist2d(&h2)` |
| `fit` | oxifit's `FitData`: `h.fit(&model)` on a 1-D histogram |

## The oxihep family

Rust crates for high-energy physics, split the way Python's
[Scikit-HEP](https://scikit-hep.org) is. Each works on its own; together they
cover reading data, histogramming, fitting and plotting.

| Crate | What it does | Python counterpart |
| --- | --- | --- |
| [oxiroot](https://github.com/oxihep/oxiroot) | Read and write ROOT files: `TTree`, RNTuple, ROOT histograms and graphs | uproot |
| **oxihist** (this crate) | Histograms: typed axes, weighted and profile storages | boost-histogram, hist |
| [oxiplot](https://github.com/oxihep/oxiplot) | Plotting to SVG/PNG/PDF with an mplhep look | matplotlib, mplhep |
| [oxifit](https://github.com/oxihep/oxifit) | Curve fitting with Minuit2, and a `TFormula` engine | iminuit |
| [oxistats](https://github.com/oxihep/oxistats) | Statistics, HEP lineshapes and intervals | scipy.stats, hepstats |
| [oxiparticle](https://github.com/oxihep/oxiparticle) | PDG particle data and ID decoding | particle |

The crates are not on crates.io yet: depend on them from GitHub, pinned to a
release tag, as shown above.

## Not yet

N-dimensional (more than 3 axes) and sparse histograms, slicing by value, and
axis growth on numeric axes.

## License

MIT — see [LICENSE](LICENSE).
