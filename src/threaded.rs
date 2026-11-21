//! Filling one histogram from many threads: [`ThreadedHist`] and, with the
//! `rayon` feature, [`fill_par`].
//!
//! Both work for any [`Mergeable`] type, not just oxihist's own histograms:
//! a crate implements `Mergeable` (with its own error type) for its histogram
//! types and gets the accumulator and the parallel fill for free. oxiroot does
//! this for ROOT's `TH1`/`TH2`/`TH3`, profiles and graphs.

use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering};

use thread_local::ThreadLocal;

use crate::axis::Value;
use crate::error::Result;
use crate::histogram::Histogram;
use crate::storage::{CountStorage, MeanStorage, Storage};

/// Histograms that combine into one: the reduction behind multithreaded
/// fills and in-memory merges.
///
/// Implement it for your own histogram type, with your own error type, to use
/// it with [`ThreadedHist`] and [`fill_par`].
#[doc(alias = "Merge")]
pub trait Mergeable: Clone + Send + Sized {
    /// Why two values cannot be combined.
    type Error: std::error::Error + Send + Sync + 'static;

    /// Combine `other` into `self`, as if `self` had also been filled with
    /// everything `other` was.
    ///
    /// # Errors
    ///
    /// If the two cannot be combined (different binnings), leaving `self`
    /// unchanged.
    fn merge(&mut self, other: &Self) -> std::result::Result<(), Self::Error>;

    /// Fold an iterator into one value, starting from the first item: an
    /// in-memory `hadd`. `Ok(None)` for an empty iterator.
    ///
    /// # Errors
    ///
    /// The error of the first pair that cannot be combined.
    fn merge_all<I: IntoIterator<Item = Self>>(
        items: I,
    ) -> std::result::Result<Option<Self>, Self::Error> {
        let mut it = items.into_iter();
        let Some(mut acc) = it.next() else {
            return Ok(None);
        };
        for item in it {
            acc.merge(&item)?;
        }
        Ok(Some(acc))
    }
}

impl<S: Storage, const D: usize> Mergeable for Histogram<S, D> {
    type Error = crate::Error;

    fn merge(&mut self, other: &Self) -> Result<()> {
        self.add(other)
    }
}

/// A multithreaded fill accumulator, the analog of ROOT's `TThreadedObject`.
///
/// Share `&ThreadedHist` across threads and fill from any of them: each thread
/// fills its own copy of the template (made on first use), and finding that
/// copy takes no shared lock, so fills scale with the number of threads.
/// [`merge`](ThreadedHist::merge) combines the copies at the end, their
/// [fill statistics](Histogram::fill_stats) included.
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

impl<H: Mergeable> std::fmt::Debug for ThreadedHist<H> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ThreadedHist")
            .field(
                "copies",
                &self.count.load(std::sync::atomic::Ordering::Relaxed),
            )
            .finish_non_exhaustive()
    }
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

    /// The number of per-thread copies made so far: at most one per thread
    /// that filled. A thread started after another has exited may reuse the
    /// exited thread's copy, which keeps that thread's fills, so the merged
    /// result is unaffected.
    pub fn num_slots(&self) -> usize {
        self.count.load(Ordering::Relaxed)
    }

    /// Combine every thread's copy into one (the template, if no thread
    /// filled). Each copy started as a clone of the template, so the template
    /// itself is not added again.
    ///
    /// # Errors
    ///
    /// The error of [`Mergeable::merge`] if two copies cannot be combined;
    /// copies of one template always can, including when growing category
    /// axes met different categories on different threads.
    pub fn merge(self) -> std::result::Result<H, H::Error> {
        let copies = self.slots.into_iter().map(RefCell::into_inner);
        Ok(H::merge_all(copies)?.unwrap_or(self.template))
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
        // The data are multiples of 1/4, so the sums are exact in any order.
        assert_eq!(merged.fill_stats(), serial.fill_stats());
        assert_eq!(merged.entries(), 400.0);
    }

    #[test]
    fn profiles_merge_their_statistics() {
        let acc = ThreadedHist::new(Hist::reg(4, 0.0, 4.0).mean());
        std::thread::scope(|s| {
            for t in 0..4u8 {
                let acc = &acc;
                s.spawn(move || {
                    for i in 0..25u8 {
                        acc.fill_sample(f64::from(i % 4) + 0.5, f64::from(t));
                    }
                });
            }
        });
        let merged = acc.merge().unwrap();
        let stats = merged.fill_stats().unwrap();
        assert_eq!((stats.entries, stats.sum_w), (100.0, 100.0));
        assert_eq!((stats.sum_ws, stats.sum_ws2), (150.0, 350.0));
    }

    #[cfg(feature = "rayon")]
    #[test]
    fn fill_par_keeps_the_statistics() {
        let data: Vec<f64> = (0..4000).map(|i| f64::from(i % 40) / 4.0).collect();
        let template = Hist::reg(8, 0.0, 8.0).reg(8, 0.0, 8.0).weight();
        let h = fill_par(&template, &data, |h, &x| h.fill_weight(x, 9.75 - x, 0.5));
        let mut serial = template.clone();
        for &x in &data {
            serial.fill_weight(x, 9.75 - x, 0.5);
        }
        assert_eq!(h, serial);
        assert_eq!(h.fill_stats(), serial.fill_stats());
        assert_eq!(h.entries(), 4000.0);
    }

    /// A type from another crate, with its own error, plugs in through
    /// `Mergeable` alone.
    #[derive(Clone, Debug, PartialEq)]
    struct Counter(u64);

    #[derive(Debug)]
    struct Never;
    impl std::fmt::Display for Never {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str("never")
        }
    }
    impl std::error::Error for Never {}

    impl Mergeable for Counter {
        type Error = Never;
        fn merge(&mut self, other: &Self) -> std::result::Result<(), Never> {
            self.0 += other.0;
            Ok(())
        }
    }

    #[test]
    fn foreign_types_use_the_accumulator() {
        let acc = ThreadedHist::new(Counter(0));
        std::thread::scope(|s| {
            for _ in 0..4 {
                let acc = &acc;
                s.spawn(move || (0..250).for_each(|_| acc.with_local(|c| c.0 += 1)));
            }
        });
        assert_eq!(acc.merge().unwrap(), Counter(1000));
        assert_eq!(Counter::merge_all(Vec::new()).unwrap(), None);
        assert_eq!(
            Counter::merge_all([Counter(1), Counter(2)]).unwrap(),
            Some(Counter(3))
        );
        // The template is the start of every copy, not an extra addend.
        let acc = ThreadedHist::new(Counter(5));
        acc.with_local(|c| c.0 += 1);
        assert_eq!(acc.merge().unwrap(), Counter(6));
    }

    #[test]
    fn mismatched_binnings_do_not_merge() {
        let mut a = Hist::reg(2, 0.0, 1.0).double();
        assert!(a.merge(&Hist::reg(3, 0.0, 1.0).double()).is_err());
    }
}
