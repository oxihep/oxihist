# Changelog

Notable changes to oxihist, by release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.14.0]

Exact fill statistics and the analysis methods oxiroot's ROOT histograms
had and these lacked, so that `Histogram<S, D>` can be the one histogram
type and oxiroot's `Hist1D`/`Hist2D`/`Hist3D`/`Profile*` plain I/O records.
Breaking; see **Changed**.

### Added

- **`FillStats<D>`**: what a histogram was filled with, as running sums
  (ROOT's `fEntries` and `fTsumw*`). `entries` counts every fill; `sum_w`,
  `sum_w2`, `sum_wx[i]` and `sum_wxy[i][j]` (`Σw·xᵢ·xⱼ`: the diagonal is
  `Σw·xᵢ²`, the rest the cross terms) count the fills in range on every
  axis; `sum_ws` and `sum_ws2` are a profile's sums of the sampled value.
  `FillStats::mean`, `std_dev` and `effective_entries` read them. Every fill
  path accumulates them: `fill*`, `try_fill*`, `fill_many`, `fill_sample*`,
  `fill_random*`, `ThreadedHist` and `fill_par`.
  - The coordinate of a fill is the value filled on a regular, variable or
    integer axis, and the bin's center (`0.5`, `1.5`, …) on a category or
    boolean axis. A growing category axis adds nothing: its bins have no
    fixed position, which is what lets two copies that met the categories in
    different orders be added.
  - `add` sums them, `scale` scales them (`Σw` and the moments by `c`, `Σw²`
    by `c²`; a profile's sample sums instead), `reset` zeroes them,
    `project` keeps the kept axes' terms and `rebin` leaves them alone.
  - A histogram whose bins did not all come from fills has none:
    `Histogram::fill_stats` returns `None` after `from_bins`, `bin_mut`,
    `multiply`, `divide`, `smooth`, `slice` and `cumulative`, after adding a
    histogram that has none, and for data serialised by 0.13.
    `Histogram::set_fill_stats` restores them (a ROOT reader, from the
    file's sums), and `Histogram::binned_stats` is the estimate from the
    bin centers that takes their place.
- **`Histogram::entries`** and **`effective_entries`** (`(Σw)² / Σw²`), exact
  with fill statistics and estimated from the bins without.
- **`Histogram::multiply`** and **`divide`**: bin by bin, flow bins
  included, by a histogram of any storage, with the variances propagated as
  for independent values. A zero denominator gives an empty bin.
- **`Histogram::maximum`** and **`minimum`**: the in-range bin with the
  largest or smallest value, as `Some((index, value))`.
- **`Histogram::chi2_test(other, kind)`** and **`kolmogorov_test(other)`**
  for 1-D count histograms, returning `Chi2TestResult` and `KsTestResult`.
  `Chi2TestKind` chooses ROOT's `"UU"`, `"UW"` or `"WW"` test. The numbers
  are ROOT's `Chi2Test` and `KolmogorovTest`.
- **`Histogram::smooth(times)`**: ROOT's 353QH-twice smoother
  (`TH1::Smooth`) for 1-D counts.
- **`Histogram::interpolate(x)`**: linear between bin centers
  (`TH1::Interpolate`); `None` outside the axis.
- **`Histogram::slice(lo..hi)`**: the bins of a 1-D histogram overlapping a
  value range, on an axis cut down to them; what is cut off goes to the flow
  bins.
- **`Histogram::profile(axis)`**: a 2-D count histogram's profile along one
  axis (`ProfileX`, `ProfileY`) as a 1-D `WeightedMean` histogram.
- **`Storage::sums`** and **`BinSums`** (what a bin adds to the statistics
  estimated from bins), **`Storage::SAMPLED`** (whether a bin holds samples
  of a value) and **`CountStorage::bin`** (a bin from a value and a
  variance).

### Changed

- **`mean` and `std_dev` are those of the values filled**, not of the bin
  centers, whenever the histogram has fill statistics. A histogram without
  them is estimated from its bin centers as before, with two differences:
  a profile's bins are weighted by their number of samples instead of their
  mean, and `std_dev` of an empty histogram is `NaN`, as documented, not 0.
- **Equality does not compare the fill statistics.** Two histograms are
  equal when their axes, name, title and bins are, as before. The sums
  depend on the order of the fills in the last bit, so two threads'
  histograms merged would no longer equal a serial fill, which they do bin
  for bin with integer counts; and a histogram rebuilt from its bins, which
  has no statistics, would not equal the one it was rebuilt from. Compare
  `fill_stats()` as well where the sums matter; `set_fill_stats` lets a
  round trip through a file restore them exactly.
- **`bin_mut` drops the fill statistics**: a bin edited by hand is not a
  fill.
- **serde**: a histogram serialises a `stats` field (`null` without fill
  statistics). Data written by 0.13 has no such field and still loads, with
  no fill statistics. The bins' layout is unchanged.
- **`CountStorage` has a new required method, `bin`**, for storages
  implemented outside this crate.
- **Fills cost more**: each one also updates the running sums, as in ROOT.
  Measured with oxiroot's `fill` benchmark (a million values, Apple
  silicon), against 0.13.1:

  | Fill | 0.13.1 | 0.14.0 |
  | --- | --- | --- |
  | 1-D regular, `Double` | 750–800 M/s | 475 M/s |
  | 1-D regular, `Weight` | 680–705 M/s | 460 M/s |
  | 2-D regular, `Double` | 255–265 M/s | 223 M/s |
  | 1-D profile, `Mean` | 443–453 M/s | 318 M/s |
  | string category, 20 categories | 124–132 M/s | 106 M/s |

  Two things keep it to that. The fast path finds the bin, whether the
  fill is in range and its coordinate in one pass over the axis, instead
  of matching on the axis kind again for each. And the sums are kept in
  four partial sets, the bin choosing which, so that consecutive fills in
  different bins do not wait on each other's additions; `fill_stats`
  returns their total. This costs 256 bytes per 1-D histogram, 512 per
  2-D and 768 per 3-D. `set_fill_stats(None)` on the empty histogram
  turns the bookkeeping off for a fill that does not need it.
- **oxistats 0.2.0** is a dependency (the chi-square and Kolmogorov
  probabilities), at the tag oxifit and oxiroot name.

## [0.13.1]

### Changed

- **oxiplot 0.11.1** behind the `plot` feature: the legend's automatic
  placement keeps off text annotations.

## [0.13.0]

API review toward 1.0 (breaking).

### Changed

- **`Histogram::std` is now `std_dev`**, the name the oxihep crates use for
  a standard deviation everywhere.
- **oxifit 0.9.0 and oxiplot 0.11.0** behind the `fit` and `plot` features
  (renamed toggles, `XyPoints`, `BinomialData`; see their changelogs).

## [0.12.0]

### Changed

- **oxifit 0.8.0 and oxiplot 0.10.0** behind the `fit` and `plot` features
  (oxifit: NUTS sampling and faster gradient fits; oxiplot: result plots,
  and `#[non_exhaustive]` on its public enums).
- **`ThreadedHist` implements `Debug`**.
- **Lints**: the cast lints are on. The binning code's `usize` ↔ `isize`
  bin-index conversions go through a few inlined helpers that state why
  each is exact (bins live in a `Vec`), so the fill paths are unchanged;
  `From<usize> for Value` saturates instead of wrapping.

## [0.11.0]

### Added

- **Bin-integrated fits use the real bin edges.** A 1-D histogram's
  `FitData` now reports each bin's `(lo, hi)`, so oxifit's
  `FitOptions::with_bin_integration(true)` averages the model over the
  actual bins of a variable-width axis instead of guessing them halfway
  between centers.

### Changed

- **oxifit 0.7.0 and oxiplot 0.9.0** behind the `fit` and `plot` features:
  oxifit gains a `Cost` trait and `Fitter`, CLs limits, template fits,
  sPlot, efficiency fits, bin-integrated fits and an ensemble sampler. No
  oxihist code changes.

## [0.10.0]

### Added

- **2-D and 3-D histograms can be fitted** (the `fit` feature): they
  implement oxifit's `FitDataNd`, each in-range bin a point at its center
  with its value and `√variance`, so `h.fit(&model)` works with a model of
  as many variables (`Model::new2`, `Model::product`, …).

### Changed

- **oxifit 0.5.0 and oxiplot 0.8.0** behind the `fit` and `plot` features.

## [0.9.1]

### Changed

- **oxifit 0.4.1 and oxiplot 0.7.1**: oxifit's numerical `integrate` no
  longer hangs on a narrow peak or a `NaN` integrand.

## [0.9.0]

### Changed

- **oxifit 0.4.0 and oxiplot 0.7.0** behind the `fit` and `plot` features:
  unbinned and extended likelihood fits (on `oxifit::Events`), the SumW2
  correction, toy studies and simultaneous fits — histograms take part in a
  `Simultaneous` fit as channels.

## [0.8.0]

### Changed

- **oxifit 0.3.0 and oxiplot 0.6.0** behind the `fit` and `plot` features:
  CMA-ES and Levenberg–Marquardt minimizers, sum models (`sig + bkg`),
  Gaussian parameter constraints, profile likelihoods and 2-D contours
  (`h.profile(…)`, `h.contour(…)` on a histogram through `FitExt`), and
  oxiplot's `fit_contour`/`fit_profile` to draw them.

## [0.7.0]

### Changed

- **oxifit 0.2.0 and oxiplot 0.5.0** behind the `fit` and `plot` features.
  oxifit 0.2 adds global minimizers (multistart, differential evolution),
  runs HESSE after MIGRAD by default, and makes `FitResult`
  `#[non_exhaustive]`; `h.fit(&model)` is unchanged. A minor version because
  those crates' types appear in this crate's API.

## [0.6.3]

### Changed

- The `plot` feature uses oxiplot 0.4.3 (its `Loc::Best` legend placement
  keeps off filled histograms).

### Added

- Examples: `dimuon_mass`, `cutflow`, `jet_response`, `occupancy` and
  `threaded_fill` (`cargo run --example …`), listed in the README.

## [0.6.2]

### Added

- `BinnedSampler`: a binned density prepared once and sampled many times
  (one uniform draw and a binary search per value), what `fill_random` and
  `Random::sample_binned` now use. Downstream code (oxiroot's `Hist1D`)
  shares it instead of keeping its own copy. The draws are unchanged.

### Changed

- The `plot` feature uses oxiplot 0.4.2 (its `Error::Io` reports its cause
  through `source()`).

## [0.6.1]

### Changed

- The `fit` feature uses oxifit 0.1.2, the version oxiplot 0.4.0 uses, so a
  build with both resolves one oxifit (oxifit 0.1.2's Nelder–Mead fallback
  is valid only when it converged).

## [0.6.0]

### Changed

- **Breaking: the `plot` feature uses oxiplot 0.4.0** (oxiplot is a public
  dependency: its traits are implemented for `Histogram`). Its plots gain
  stacks, labelled category axes, SVGs sized in points and large heatmaps
  as one image (see its changelog).

### Added

- A 1-D histogram over a string, integer-category or boolean axis names its
  bins in plots (`Hist1dData::x_labels`): each bin is labelled with its
  category instead of a numeric tick.
- `Axis::try_integer`, and `Error::InvalidData` for deserialised data that
  breaks an invariant.

### Fixed

- **Deserialising skipped the constructors' checks.** With the `serde`
  feature, `Axis` and `Histogram` derived `Deserialize` directly, so a
  hand-edited or hostile file could build an axis with descending edges or
  no bins, or a histogram whose storage does not match its axes (which then
  panicked or misindexed on fill). Loading now checks exactly what the
  constructors do, and the storage size against the axes, and returns an
  error. The serialised format is unchanged. A loaded category axis also
  rebuilds its lookup index at once rather than scanning until it grows.
- **`Axis::integer` over a range wider than an index.** `hi - lo` was not
  checked, so `Axis::integer(i64::MIN, i64::MAX)` overflowed. It now panics
  with a clear message (`try_integer` returns an error).

## [0.5.1]

### Changed

- The `plot` feature uses oxiplot 0.3.1, whose histograms no longer drop
  their y axis below zero for error bars reaching below an empty bin.

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
