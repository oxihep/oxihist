//! Filling one histogram from many threads: [`ThreadedHist`] and, with the
//! `rayon` feature, [`fill_par`].

use std::collections::HashMap;
use std::sync::{Arc, Mutex, PoisonError, RwLock};
use std::thread::ThreadId;

use crate::axis::Value;
use crate::error::Result;
use crate::histogram::Histogram;
use crate::storage::{CountStorage, MeanStorage, Storage};

/// Histograms that combine into one: the reduction behind multithreaded
/// fills.
pub trait Mergeable: Clone + Send + Sized {
    /// Combine `other` into `self`, as if `self` had also been filled with
    /// everything `other` was.
    ///
    /// # Errors
    ///
    /// If the two cannot be combined (different binnings), leaving `self`
    /// unchanged.
    fn merge(&mut self, other: &Self) -> Result<()>;
}

impl<S: Storage, const D: usize> Mergeable for Histogram<S, D> {
    fn merge(&mut self, other: &Self) -> Result<()> {
        self.add(other)
    }
}

/// A multithreaded fill accumulator, the analog of ROOT's `TThreadedObject`.
///
/// Share `&ThreadedHist` across threads and fill from any of them: each thread
/// fills its own copy of the template (made on first use), so threads never
/// contend. [`merge`](ThreadedHist::merge) combines the copies at the end.
///
/// ```
/// use oxihist::{Hist, ThreadedHist};
///
/// let data: Vec<f64> = (0..1000).map(|i| f64::from(i % 100)).collect();
/// let hist = ThreadedHist::new(Hist::reg(100, 0.0, 100.0).int64());
/// std::thread::scope(|s| {
///     for chunk in data.chunks(250) {
///         let hist = &hist;
///         s.spawn(move || chunk.iter().for_each(|&x| hist.fill(x)));
///     }
/// });
/// assert_eq!(hist.merge().unwrap().sum(true), 1000);
/// ```
pub struct ThreadedHist<H: Mergeable> {
    template: H,
    slots: RwLock<HashMap<ThreadId, Arc<Mutex<H>>>>,
}

impl<H: Mergeable> ThreadedHist<H> {
    /// An accumulator whose per-thread copies start as clones of `template`
    /// (normally empty).
    pub fn new(template: H) -> Self {
        ThreadedHist {
            template,
            slots: RwLock::new(HashMap::new()),
        }
    }

    /// Run `f` on the calling thread's own copy.
    pub fn with_local<R>(&self, f: impl FnOnce(&mut H) -> R) -> R {
        let id = std::thread::current().id();
        let slot = {
            let read = self.slots.read().unwrap_or_else(PoisonError::into_inner);
            read.get(&id).cloned()
        };
        let slot = slot.unwrap_or_else(|| {
            let mut write = self.slots.write().unwrap_or_else(PoisonError::into_inner);
            write
                .entry(id)
                .or_insert_with(|| Arc::new(Mutex::new(self.template.clone())))
                .clone()
        });
        let mut h = slot.lock().unwrap_or_else(PoisonError::into_inner);
        f(&mut h)
    }

    /// The number of per-thread copies made so far.
    pub fn num_slots(&self) -> usize {
        self.slots
            .read()
            .unwrap_or_else(PoisonError::into_inner)
            .len()
    }

    /// Combine every thread's copy into one (the template, if no thread
    /// filled).
    ///
    /// # Errors
    ///
    /// Never for copies of one template; the error of [`Mergeable::merge`]
    /// otherwise.
    pub fn merge(self) -> Result<H> {
        let mut acc = self.template;
        let slots = self
            .slots
            .into_inner()
            .unwrap_or_else(PoisonError::into_inner);
        for (_, slot) in slots {
            let h = Arc::try_unwrap(slot).map_or_else(
                |shared| {
                    shared
                        .lock()
                        .unwrap_or_else(PoisonError::into_inner)
                        .clone()
                },
                |m| m.into_inner().unwrap_or_else(PoisonError::into_inner),
            );
            acc.merge(&h)?;
        }
        Ok(acc)
    }
}

impl<S: CountStorage> ThreadedHist<Histogram<S, 1>> {
    /// Fill `x` into this thread's copy.
    ///
    /// # Panics
    ///
    /// If `x` does not fit the axis.
    pub fn fill<'a>(&self, x: impl Into<Value<'a>>) {
        let x = x.into();
        self.with_local(|h| h.fill(x));
    }

    /// Fill `x` with weight `w` into this thread's copy.
    ///
    /// # Panics
    ///
    /// If `x` does not fit the axis.
    pub fn fill_weight<'a>(&self, x: impl Into<Value<'a>>, w: f64) {
        let x = x.into();
        self.with_local(|h| h.fill_weight(x, w));
    }
}

impl<S: CountStorage> ThreadedHist<Histogram<S, 2>> {
    /// Fill `(x, y)` into this thread's copy.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    pub fn fill<'a>(&self, x: impl Into<Value<'a>>, y: impl Into<Value<'a>>) {
        let (x, y) = (x.into(), y.into());
        self.with_local(|h| h.fill(x, y));
    }

    /// Fill `(x, y)` with weight `w` into this thread's copy.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    pub fn fill_weight<'a>(&self, x: impl Into<Value<'a>>, y: impl Into<Value<'a>>, w: f64) {
        let (x, y) = (x.into(), y.into());
        self.with_local(|h| h.fill_weight(x, y, w));
    }
}

impl<S: CountStorage> ThreadedHist<Histogram<S, 3>> {
    /// Fill `(x, y, z)` into this thread's copy.
    ///
    /// # Panics
    ///
    /// If a value does not fit its axis.
    pub fn fill<'a>(
        &self,
        x: impl Into<Value<'a>>,
        y: impl Into<Value<'a>>,
        z: impl Into<Value<'a>>,
    ) {
        let (x, y, z) = (x.into(), y.into(), z.into());
        self.with_local(|h| h.fill(x, y, z));
    }
}

impl<S: MeanStorage> ThreadedHist<Histogram<S, 1>> {
    /// Add the sample `s` at `x` in this thread's copy.
    ///
    /// # Panics
    ///
    /// If `x` does not fit the axis.
    pub fn fill_sample<'a>(&self, x: impl Into<Value<'a>>, s: f64) {
        let x = x.into();
        self.with_local(|h| h.fill_sample(x, s));
    }
}

/// Fill a histogram from a slice in parallel (the `rayon` feature): each task
/// folds into a clone of `template`, and the partial histograms are merged.
///
/// ```
/// use oxihist::{fill_par, Hist};
/// let data: Vec<f64> = (0..1000).map(|i| f64::from(i % 100)).collect();
/// let h = fill_par(&Hist::reg(100, 0.0, 100.0).int64(), &data, |h, &x| h.fill(x));
/// assert_eq!(h.sum(true), 1000);
/// ```
#[cfg(feature = "rayon")]
pub fn fill_par<H, T, F>(template: &H, data: &[T], fill: F) -> H
where
    H: Mergeable + Sync,
    T: Sync,
    F: Fn(&mut H, &T) + Sync,
{
    use rayon::prelude::*;
    data.par_iter()
        .fold(
            || template.clone(),
            |mut h, item| {
                fill(&mut h, item);
                h
            },
        )
        .reduce(
            || template.clone(),
            |mut a, b| {
                a.merge(&b)
                    .expect("fill_par: every partial is a clone of one template");
                a
            },
        )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Hist;

    #[test]
    fn threads_merge_to_a_serial_fill() {
        let data: Vec<f64> = (0..400).map(|i| f64::from(i % 40) / 4.0).collect();
        let acc = ThreadedHist::new(Hist::reg(10, 0.0, 10.0).weight());
        std::thread::scope(|s| {
            for chunk in data.chunks(100) {
                let acc = &acc;
                s.spawn(move || chunk.iter().for_each(|&x| acc.fill(x)));
            }
        });
        assert!(acc.num_slots() >= 1);
        let merged = acc.merge().unwrap();
        let mut serial = Hist::reg(10, 0.0, 10.0).weight();
        serial.fill_many(data.iter().copied());
        assert_eq!(merged, serial);
    }

    #[test]
    fn mismatched_binnings_do_not_merge() {
        let mut a = Hist::reg(2, 0.0, 1.0).double();
        assert!(a.merge(&Hist::reg(3, 0.0, 1.0).double()).is_err());
    }
}
