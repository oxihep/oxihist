# Changelog

Notable changes to oxihist, by release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

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

### Fixed

- **`ThreadedHist` got slower as threads were added.** Every fill took a
  shared read lock to find the thread's copy. The copies now live in
  thread-local storage (the `thread_local` crate), so a fill touches no shared
  state. 24 M fills with 12 threads: 15 s before, 32 ms now.
- **`fill_par` cloned the template once per rayon split** (about 1 300 copies
  of a large histogram). It now fills a few chunks per thread. With an 8 MB
  template: 533 ms before, 20 ms now.

### Changed

- `ThreadedHist::with_local` panics if its closure re-enters the same
  accumulator (it used to deadlock).

### Fixed

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
