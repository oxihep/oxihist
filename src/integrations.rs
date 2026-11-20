//! Optional integrations with the other oxihep crates.

/// [oxiplot](https://github.com/oxihep/oxiplot) draws 1-D and 2-D histograms
/// (the `plot` feature): a profile ([`Mean`](crate::Mean) storage) draws its
/// per-bin means with their standard errors, a categorical axis one unit per
/// category, labelled with the category (a 1-D histogram's string, integer
/// or boolean categories name its bins on the plot).
#[cfg(feature = "plot")]
mod plot {
    use oxiplot::{Hist1dData, Hist2dData};

    use crate::{Histogram, Storage};

    impl<S: Storage> Hist1dData for Histogram<S, 1> {
        fn edges(&self) -> Vec<f64> {
            self.axis(0).plot_edges()
        }
        fn values(&self) -> Vec<f64> {
            Histogram::values(self, false)
        }
        fn error(&self, i: usize) -> f64 {
            self.variance([crate::cast::signed(i)])
                .unwrap_or(0.0)
                .max(0.0)
                .sqrt()
        }
        fn x_labels(&self) -> Option<Vec<String>> {
            let axis = self.axis(0);
            axis.is_categorical()
                .then(|| (0..axis.size()).filter_map(|i| axis.bin_label(i)).collect())
        }
    }

    impl<S: Storage> Hist2dData for Histogram<S, 2> {
        fn x_edges(&self) -> Vec<f64> {
            self.axis(0).plot_edges()
        }
        fn y_edges(&self) -> Vec<f64> {
            self.axis(1).plot_edges()
        }
        fn grid(&self) -> Vec<Vec<f64>> {
            self.values_grid()
        }
    }

    #[cfg(test)]
    mod tests {
        use oxiplot::Axes;

        use crate::Hist;

        #[test]
        fn histograms_render() {
            let mut h = Hist::reg(10, 0.0, 10.0).weight();
            h.fill(3.0);
            let mut p = Hist::reg(10, 0.0, 10.0).mean();
            p.fill_sample(3.0, 1.0);
            let mut h2 = Hist::reg(4, 0.0, 1.0).str_cat(&["a", "b"]).double();
            h2.fill(0.2, "b");
            let mut ax = Axes::new();
            ax.hist(&h).profile(&p);
            let mut ax2 = Axes::new();
            ax2.hist2d(&h2);
            for svg in [ax.to_svg_string(), ax2.to_svg_string()] {
                assert!(svg.starts_with("<svg") && svg.ends_with("</svg>"));
            }
        }

        #[test]
        fn category_axes_name_their_bins() {
            use oxiplot::Hist1dData;
            let mut h = Hist::str_cat(&["ee", "μμ"]).growth().double();
            h.fill("eμ");
            assert_eq!(
                h.x_labels(),
                Some(vec!["ee".into(), "μμ".into(), "eμ".into()])
            );
            let b = Hist::boolean().double();
            assert_eq!(b.x_labels(), Some(vec!["false".into(), "true".into()]));
            let i = Hist::int_cat(&[11, 13]).double();
            assert_eq!(i.x_labels(), Some(vec!["11".into(), "13".into()]));
            // Numeric axes, integer ones included, keep numeric ticks.
            assert_eq!(Hist::reg(4, 0.0, 1.0).double().x_labels(), None);
            assert_eq!(Hist::int(0, 5).double().x_labels(), None);
            // And the plot draws.
            let mut ax = Axes::new();
            ax.hist(&h);
            assert!(ax.to_svg_string().ends_with("</svg>"));
        }
    }
}

/// [oxifit](https://github.com/oxihep/oxifit) fits 1-, 2- and 3-D histograms
/// (the `fit` feature): each in-range bin is a point at its center with the
/// bin's value and `√variance`. An empty count bin has σ = 0, so a Neyman χ²
/// fit skips it, as ROOT does. A 2-D or 3-D histogram is fitted with a model
/// of as many variables (`Model::new2`, `Model::product`, …).
#[cfg(feature = "fit")]
mod fit {
    use oxifit::{FitData, FitDataNd, Point, PointNd};

    use crate::{Histogram, Storage};

    impl<S: Storage> FitData for Histogram<S, 1> {
        fn points(&self) -> Vec<Point> {
            let centers = self.axis(0).centers();
            (0..self.axis(0).size())
                .map(|i| {
                    let b = self.bin([crate::cast::signed(i)]).expect("in range");
                    Point::new(centers[i], S::value(b), S::variance(b).max(0.0).sqrt())
                })
                .collect()
        }

        /// The bins' `(lo, hi)` for bin-integrated fits
        /// (`FitOptions::with_bin_integration`); `None` for a categorical axis.
        fn bin_edges(&self) -> Option<Vec<(f64, f64)>> {
            let edges = self.axis(0).edges()?;
            Some(edges.windows(2).map(|w| (w[0], w[1])).collect())
        }
    }

    impl<S: Storage> FitDataNd for Histogram<S, 2> {
        fn dim(&self) -> usize {
            2
        }
        fn points_nd(&self) -> Vec<PointNd> {
            let (cx, cy) = (self.axis(0).centers(), self.axis(1).centers());
            let mut out = Vec::with_capacity(cx.len() * cy.len());
            for (i, &x) in cx.iter().enumerate() {
                for (j, &y) in cy.iter().enumerate() {
                    let b = self
                        .bin([crate::cast::signed(i), crate::cast::signed(j)])
                        .expect("in range");
                    out.push(PointNd::new(
                        [x, y, 0.0],
                        S::value(b),
                        S::variance(b).max(0.0).sqrt(),
                    ));
                }
            }
            out
        }
    }

    impl<S: Storage> FitDataNd for Histogram<S, 3> {
        fn dim(&self) -> usize {
            3
        }
        fn points_nd(&self) -> Vec<PointNd> {
            let c = [
                self.axis(0).centers(),
                self.axis(1).centers(),
                self.axis(2).centers(),
            ];
            let mut out = Vec::with_capacity(c[0].len() * c[1].len() * c[2].len());
            for (i, &x) in c[0].iter().enumerate() {
                for (j, &y) in c[1].iter().enumerate() {
                    for (k, &z) in c[2].iter().enumerate() {
                        let b = self
                            .bin([
                                crate::cast::signed(i),
                                crate::cast::signed(j),
                                crate::cast::signed(k),
                            ])
                            .expect("in range");
                        out.push(PointNd::new(
                            [x, y, z],
                            S::value(b),
                            S::variance(b).max(0.0).sqrt(),
                        ));
                    }
                }
            }
            out
        }
    }

    #[cfg(test)]
    mod tests {
        use oxifit::{FitData, FitExt, FitOptions, Model};

        use crate::Hist;

        #[test]
        fn a_line_fits() {
            let mut h = Hist::reg(10, 0.0, 10.0).weight();
            for i in 0..10 {
                let x = f64::from(i) + 0.5;
                h.fill_weight(x, 2.0 * x + 1.0);
            }
            let r = h.fit(&Model::polynomial("line", 1).with_params(vec![0.0, 1.0]));
            assert!(r.valid);
            assert!((r.params[1] - 2.0).abs() < 1e-6, "{:?}", r.params);
        }

        #[test]
        fn variable_bins_are_integrated_over_their_edges() {
            // Contents are the exact bin averages of 100·e^(−x/2) over very
            // uneven bins: only a fit averaging the model over each bin's real
            // edges recovers the slope.
            let edges = [0.0, 0.5, 1.5, 3.0, 5.0, 8.0];
            let mut h = Hist::var(&edges).weight();
            for w in edges.windows(2) {
                let integral = 200.0 * ((-w[0] / 2.0).exp() - (-w[1] / 2.0).exp());
                // Many small fills: errors of 1 % of the content.
                for _ in 0..10_000 {
                    h.fill_weight(0.5 * (w[0] + w[1]), integral / (w[1] - w[0]) / 1e4);
                }
            }
            let got = h.bin_edges().expect("a numeric axis");
            assert_eq!(
                got,
                vec![(0.0, 0.5), (0.5, 1.5), (1.5, 3.0), (3.0, 5.0), (5.0, 8.0)]
            );
            let model = Model::new("e", &["n", "tau"], vec![80.0, 1.5], |x, p| {
                p[0] * (-x / p[1]).exp()
            });
            let opts = FitOptions::new().bin_integration(true);
            let r = h.fit_opts(&model, &opts);
            assert!(
                r.valid && (r.params[1] - 2.0).abs() < 1e-4,
                "{:?}",
                r.params
            );
            // At the bin centers the same data pull the slope away.
            let r = h.fit(&model);
            assert!((r.params[1] - 2.0).abs() > 1e-2, "{:?}", r.params);
            assert!(Hist::str_cat(&["a", "b"]).double().bin_edges().is_none());
        }

        #[test]
        fn a_2d_and_a_3d_histogram_fit() {
            // A plane a + b·x + c·y filled as weights at the bin centers.
            let mut h = Hist::reg(8, 0.0, 8.0).reg(6, 0.0, 6.0).weight();
            for i in 0..8 {
                for j in 0..6 {
                    let (x, y) = (f64::from(i) + 0.5, f64::from(j) + 0.5);
                    h.fill_weight(x, y, 1.0 + 2.0 * x - 0.5 * y);
                }
            }
            let plane = Model::new2("plane", &["a", "b", "c"], vec![0.0, 1.0, 0.0], |x, y, p| {
                p[0] + p[1] * x + p[2] * y
            });
            let r = h.fit(&plane);
            assert!(r.valid && (r.params[2] + 0.5).abs() < 1e-6, "{r:?}");
            assert_eq!(r.ndf, 8 * 6 - 3);

            let mut h3 = Hist::reg(4, 0.0, 4.0)
                .reg(4, 0.0, 4.0)
                .reg(4, 0.0, 4.0)
                .weight();
            for i in 0..4 {
                for j in 0..4 {
                    for k in 0..4 {
                        let (x, y, z) =
                            (f64::from(i) + 0.5, f64::from(j) + 0.5, f64::from(k) + 0.5);
                        h3.fill_weight(x, y, z, 3.0 + x + y + z);
                    }
                }
            }
            let flat = Model::new3("f", &["a", "b"], vec![0.0, 0.5], |x, y, z, p| {
                p[0] + p[1] * (x + y + z)
            });
            let r = h3.fit(&flat);
            assert!(r.valid && (r.params[1] - 1.0).abs() < 1e-6, "{r:?}");
        }
    }
}
