//! Filling one histogram from many threads: [`ThreadedHist`] and, with the
//! `rayon` feature, [`fill_par`].

use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering};

use thread_local::ThreadLocal;

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
/// fills its own copy of the template (made on first use), and finding that
/// copy takes no shared lock, so fills scale with the number of threads.
/// [`merge`](ThreadedHist::merge) combines the copies at the end.
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
    /// One copy per thread. `ThreadLocal` finds the calling thread's slot
    /// without taking a shared lock, so filling threads never touch a common
    /// cache line (a shared `RwLock` here made fills slower with more threads).
    slots: ThreadLocal<RefCell<H>>,
    count: AtomicUsize,
}

impl<H: Mergeable> ThreadedHist<H> {
    /// An accumulator whose per-thread copies start as clones of `template`
    /// (normally empty).
    pub fn new(template: H) -> Self {
        ThreadedHist {
            template,
            slots: ThreadLocal::new(),
            count: AtomicUsize::new(0),
        }
    }

    /// Run `f` on the calling thread's own copy.
    ///
    /// # Panics
    ///
    /// If `f` calls `with_local` (or a `fill`) on the same accumulator again:
    /// the thread's copy is already borrowed.
    pub fn with_local<R>(&self, f: impl FnOnce(&mut H) -> R) -> R {
        let slot = self.slots.get_or(|| {
            self.count.fetch_add(1, Ordering::Relaxed);
            RefCell::new(self.template.clone())
        });
        f(&mut slot.borrow_mut())
    }

    /// The number of per-thread copies made so far.
    pub fn num_slots(&self) -> usize {
        self.count.load(Ordering::Relaxed)
    }

    /// Combine every thread's copy into one (the template, if no thread
    /// filled).
    ///
    /// # Errors
    ///
    /// The error of [`Mergeable::merge`] if two copies cannot be combined;
    /// copies of one template always can, including when growing category
    /// axes met different categories on different threads.
    pub fn merge(self) -> Result<H> {
        let mut acc = self.template;
        for slot in self.slots {
            acc.merge(&slot.into_inner())?;
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

/// Fill a histogram from a slice in parallel (the `rayon` feature).
///
/// The slice is cut into a few chunks per worker thread; each chunk is filled
/// into one clone of `template` and the clones are merged. The cost is a
/// handful of clones per thread, so even a large template stays cheap.
///
/// ```
/// use oxihist::{fill_par, Hist};
/// let data: Vec<f64> = (0..1000).map(|i| f64::from(i % 100)).collect();
/// let h = fill_par(&Hist::reg(100, 0.0, 100.0).int64(), &data, |h, &x| h.fill(x));
/// assert_eq!(h.sum(true), 1000);
/// ```
///
/// # Panics
///
/// If two partial histograms cannot be merged, which cannot happen for copies
/// of one template.
#[cfg(feature = "rayon")]
pub fn fill_par<H, T, F>(template: &H, data: &[T], fill: F) -> H
where
    H: Mergeable + Sync,
    T: Sync,
    F: Fn(&mut H, &T) + Sync,
{
    use rayon::prelude::*;
    // A few chunks per thread balance the load without many clones.
    let chunks = rayon::current_num_threads().max(1) * 4;
    let chunk_len = data.len().div_ceil(chunks).max(1);
    data.par_chunks(chunk_len)
        .map(|chunk| {
            let mut h = template.clone();
            for item in chunk {
                fill(&mut h, item);
            }
            h
        })
        .reduce_with(|mut a, b| {
            a.merge(&b)
                .expect("fill_par: partial histograms are copies of one template");
            a
        })
        .unwrap_or_else(|| template.clone())
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
