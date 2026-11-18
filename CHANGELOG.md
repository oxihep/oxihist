# Changelog

Notable changes to oxihist, by release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.5.0]

### Changed

- **Breaking: the `plot` feature uses oxiplot 0.3.0.** oxiplot 0.3 adds log
  axes and legend placement, lays out multi-panel figures around their
  labels, and fixes axis limits, non-finite data and the PDF page size (see
  its changelog). The version moves because oxiplot is a public dependency
  (its traits are implemented for `Histogram`).
- **Rust 1.100 and edition 2024.** The minimum supported Rust version
  (`rust-version`) moves from 1.95 to 1.100, the current stable release that
  CI already builds with, and the crates move to the 2024 edition (migrated
  with `cargo fix --edition`; formatting follows the 2024 style).

## [0.4.3]

### Changed

- **Faster fills on category and integer axes.** A fill of an integer, a
  known category or a boolean now skips the general locate (which also
  handles growth, value conversion and errors), as numbers on regular and
  variable axes already did. Per fill, string categories 12–14 → 6–7 ns
  (8, 20 and 200 categories; the `fill` bench: 69 → 125 M/s), integer
  categories 5.5–12 → 3.9 ns, integer axis 4.8 → 3.4 ns, boolean axis
  7.0 → 3.5 ns. Regular and variable axes are unchanged.
- The category index hashes with foldhash instead of std's SipHash (still
  seeded per process), and is kept from 5 categories instead of 17: with a
  fast hash a lookup beats comparing strings one by one.

## [0.4.2]

### Changed

- **Filling numbers into regular and variable axes is 2–3× faster**, on par
  with oxiroot's `Hist1D` (1-D regular, `Double`: 185 → 762 M fills/s;
  `Hist1D`: 751 M/s). The fill used to call an out-of-line `locate` that
  handled every axis kind, return its result through memory, re-match the
  axis kind to find the storage position, and carry the category-growth
  allocation and the panic formatting in its stack frame. Now a small,
  always-inlined path maps a number on a regular or variable axis straight
  to its storage position; category axes, other value types, growth and the
  panic stay out of line. Results are unchanged: every axis is still
  located (and grown) before a fill is dropped into a missing flow bin.

## [0.4.1]

### Changed

- The `plot` feature uses oxiplot 0.2.1, which places math subscripts where
  the font asks (it read the wrong MATH constant before).

## [0.4.0]

### Changed

- **Breaking: the `plot` feature uses oxiplot 0.2.0.** oxiplot now reads
  fonts with skrifa instead of the unmaintained ttf-parser
  (RUSTSEC-2026-0192), whose advisory ignore is removed here. Plots render
  the same; the version moves because oxiplot is a public dependency (its
  traits are implemented for `Histogram`).

## [0.3.0]

### Changed

- **Breaking: `Mergeable` has an associated `Error` type.** `ThreadedHist`
  and `fill_par` now work for any type that implements it, with that type's
  own error, so oxiroot's ROOT histograms use this accumulator instead of a
  copy of it. `ThreadedHist::merge` returns `Result<H, H::Error>`; for
  oxihist histograms that is still `oxihist::Error`.
- `ThreadedHist::merge` folds the per-thread copies only. Each copy starts as
  a clone of the template, so a non-empty template is no longer counted once
  more than the copies that hold it.

### Added

- `Mergeable::merge_all` folds an iterator into one value (`Ok(None)` when it
  is empty), an in-memory `hadd`.

## [0.2.0]

### Changed

- **Breaking:** `Histogram::bins` returns a `Vec` (a compact copy, first axis
  fastest, flow bins included) instead of a slice: storage now keeps spare
  slots for growing category axes.
- **Breaking:** `AxisKind` and `Value` are `#[non_exhaustive]`, so new axis
  kinds and value types can be added without a breaking release. Match them
  with a `_` arm.

### Fixed

- **Growing a category axis was quadratic.** Every new category moved every
  bin, and every fill scanned the categories linearly. Growing axes now keep
  doubling spare capacity, and category axes with more than 16 categories keep
  a hash index. 16 000 categories × 50 bins, 96 000 fills: 36 s before, 14 ms
  now.

## [0.1.1]

### Changed

- `ThreadedHist::with_local` panics if its closure re-enters the same
  accumulator (it used to deadlock).

### Fixed

- **`ThreadedHist` got slower as threads were added.** Every fill took a
  shared read lock to find the thread's copy. The copies now live in
  thread-local storage (the `thread_local` crate), so a fill touches no shared
  state. 24 M fills with 12 threads: 15 s before, 32 ms now.
- **`fill_par` cloned the template once per rayon split** (about 1 300 copies
  of a large histogram). It now fills a few chunks per thread. With an 8 MB
  template: 533 ms before, 20 ms now.
- A value just below a regular axis's upper edge could land in the overflow
  bin (in about a quarter of cases, depending on rounding). The axis now
  compares against its edges directly.
- Merging histograms whose growing category axes met different categories, or
  met them in a different order, failed with `BinningMismatch`; so did
  `ThreadedHist::merge`, and `fill_par` panicked. `add` now merges such axes by
  label.

## [0.1.0]

First release: 1-, 2- and 3-D histograms with regular, variable, integer,
category and boolean axes; `Double`, `Int64`, `Weight`, `Mean` and
`WeightedMean` storages; filling, adding, scaling, projecting, rebinning, 1-D
statistics and sampling; `ThreadedHist` and `fill_par`; the `serde`, `plot`
and `fit` features. The builder and the sampling generator come from
[oxiroot](https://github.com/oxihep/oxiroot)'s histogram crate, which now
converts these histograms to and from ROOT's.
