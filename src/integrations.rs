//! Optional integrations with the other oxihep crates.

/// [oxiplot](https://github.com/oxihep/oxiplot) draws 1-D and 2-D histograms
/// (the `plot` feature): a profile ([`Mean`](crate::Mean) storage) draws its
/// per-bin means with their standard errors, a categorical axis one unit per
/// category.
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
            self.variance([i as isize]).unwrap_or(0.0).max(0.0).sqrt()
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
    }
}

/// [oxifit](https://github.com/oxihep/oxifit) fits 1-D histograms (the `fit`
/// feature): each in-range bin is a point `(center, value, √variance)`. An
/// empty count bin has σ = 0, so a Neyman χ² fit skips it, as ROOT does.
#[cfg(feature = "fit")]
mod fit {
    use oxifit::{FitData, Point};

    use crate::{Histogram, Storage};

    impl<S: Storage> FitData for Histogram<S, 1> {
        fn points(&self) -> Vec<Point> {
            let centers = self.axis(0).centers();
            (0..self.axis(0).size())
                .map(|i| {
                    let b = self.bin([i as isize]).expect("in range");
                    Point::new(centers[i], S::value(b), S::variance(b).max(0.0).sqrt())
                })
                .collect()
        }
    }

    #[cfg(test)]
    mod tests {
        use oxifit::{FitExt, Model};

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
    }
}
