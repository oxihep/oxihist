# Changelog

Notable changes to oxihist, by release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0]

First release: 1-, 2- and 3-D histograms with regular, variable, integer,
category and boolean axes; `Double`, `Int64`, `Weight`, `Mean` and
`WeightedMean` storages; filling, adding, scaling, projecting, rebinning, 1-D
statistics and sampling; `ThreadedHist` and `fill_par`; the `serde`, `plot`
and `fit` features. The builder and the sampling generator come from
[oxiroot](https://github.com/oxihep/oxiroot)'s histogram crate, which now
converts these histograms to and from ROOT's.
