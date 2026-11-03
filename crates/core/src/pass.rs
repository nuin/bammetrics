//! The reader loop every metric shares. A BAM is read once; BGZF blocks are
//! inflated on worker threads (`noodles` `MultithreadedReader`), and records
//! are decoded and handed to an [`Accumulator`] on the calling thread, in
//! file order. Results therefore do not depend on the worker count.
//!
//! Several tools can share one read ([`run_many`]): each one is a boxed
//! [`Collector`] fed from batches on a thread of its own, the way Picard's
//! `SinglePassSamProgram.makeItSo` hands every record to every program of
//! `CollectMultipleMetrics`.

use std::fs::File;
use std::num::NonZero;
use std::path::Path;
use std::sync::Arc;
use std::sync::mpsc::{self, SyncSender};
use std::thread::{self, JoinHandle};

use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::{bam, bgzf, sam};

use crate::error::{Error, Result};
use crate::metric::ToolOutput;

/// Cap on the default worker count. BGZF inflation stops being the
/// bottleneck well before this on the corpus benchmark (docs/benchmarks.md).
pub const MAX_DEFAULT_THREADS: usize = 8;

/// A metric computed from one pass over the records of a BAM, the role of
/// Picard's `SinglePassSamProgram.acceptRead` and `finish`.
pub trait Accumulator {
    type Output;

    /// Take one record, in file order.
    fn accept(&mut self, record: &bam::Record) -> Result<()>;

    /// True once the accumulator wants no more records (`STOP_AFTER`); the
    /// pass then stops reading.
    fn done(&self) -> bool {
        false
    }

    /// The metric after the last record.
    fn finish(self) -> Result<Self::Output>;
}

/// An [`Accumulator`] behind a trait object, so tools of different types can
/// share one pass. Its result is the tool's files.
pub trait Collector: Send {
    fn accept(&mut self, record: &bam::Record) -> Result<()>;
    fn finish(self: Box<Self>) -> Result<ToolOutput>;
}

impl Accumulator for Box<dyn Collector> {
    type Output = ToolOutput;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        (**self).accept(record)
    }

    fn finish(self) -> Result<ToolOutput> {
        Collector::finish(self)
    }
}

/// An accumulator plus the step from its output to the files: the way a
/// tool becomes a [`Collector`].
struct Adapter<A, F> {
    accumulator: A,
    finish: Option<F>,
}

impl<A, F> Collector for Adapter<A, F>
where
    A: Accumulator + Send,
    F: FnOnce(A::Output) -> Result<ToolOutput> + Send,
{
    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        self.accumulator.accept(record)
    }

    fn finish(mut self: Box<Self>) -> Result<ToolOutput> {
        let finish = self.finish.take().expect("finish runs once");
        finish(self.accumulator.finish()?)
    }
}

/// `accumulator` as a [`Collector`] whose files `finish` builds from its
/// output.
pub fn boxed<A, F>(accumulator: A, finish: F) -> Box<dyn Collector>
where
    A: Accumulator + Send + 'static,
    F: FnOnce(A::Output) -> Result<ToolOutput> + Send + 'static,
{
    Box::new(Adapter {
        accumulator,
        finish: Some(finish),
    })
}

/// The available cores, at most [`MAX_DEFAULT_THREADS`].
pub fn default_threads() -> usize {
    thread::available_parallelism()
        .map_or(1, NonZero::get)
        .min(MAX_DEFAULT_THREADS)
}

/// Read `path` once with `threads` BGZF workers (0 counts as 1), build the
/// accumulator from the header, feed it every record and finish it.
pub fn run<A, F>(path: &Path, threads: usize, build: F) -> Result<A::Output>
where
    A: Accumulator,
    F: FnOnce(&sam::Header) -> Result<A>,
{
    run_limited(path, threads, Limits::default(), build)
}

/// The `SinglePassSamProgram` options every pass takes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Limits {
    /// `STOP_AFTER`: stop after this many records, every flag counted; 0
    /// reads them all.
    pub stop_after: u64,
    /// `ASSUME_SORTED`: when false, a header whose sort order is not
    /// `coordinate` (or that has none) is an error, as Picard throws.
    pub assume_sorted: bool,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            stop_after: 0,
            assume_sorted: true,
        }
    }
}

/// [`run`] with `STOP_AFTER` and `ASSUME_SORTED` applied.
pub fn run_limited<A, F>(path: &Path, threads: usize, limits: Limits, build: F) -> Result<A::Output>
where
    A: Accumulator,
    F: FnOnce(&sam::Header) -> Result<A>,
{
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = File::open(path).map_err(input)?;
    let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
    let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
        workers, file,
    ));
    let header = reader.read_header().map_err(input)?;
    if !limits.assume_sorted {
        check_coordinate_sorted(path, &header)?;
    }
    let mut accumulator = build(&header)?;
    let mut record = bam::Record::default();
    let mut records = 0;
    while !accumulator.done() && reader.read_record(&mut record)? != 0 {
        accumulator.accept(&record)?;
        records += 1;
        if limits.stop_after > 0 && records >= limits.stop_after {
            break;
        }
    }
    accumulator.finish()
}

/// One pass over `path` shared by every collector `build` returns, each on a
/// thread of its own. The outputs come back in the same order. The first
/// collector to fail fails the run, and no file of any of them is written,
/// as one exception ends a Picard `CollectMultipleMetrics` run.
pub fn run_many<F>(path: &Path, threads: usize, limits: Limits, build: F) -> Result<Vec<ToolOutput>>
where
    F: FnOnce(&sam::Header) -> Result<Vec<Box<dyn Collector>>>,
{
    let results = run_limited(path, threads, limits, |header| {
        Ok(Fanout::new(build(header)?))
    })?;
    results.into_iter().collect()
}

/// Records per batch handed to the collector threads.
pub(crate) const BATCH_RECORDS: usize = 4096;
/// Batches a collector thread may have queued before the reader waits.
const QUEUED_BATCHES: usize = 8;

pub(crate) type Batch = Arc<Vec<bam::Record>>;

/// An accumulator, or the error that stopped it. Once failed it takes no
/// more records.
pub(crate) enum Slot<A> {
    Running(A),
    Failed(Error),
}

impl<A: Accumulator> Slot<A> {
    pub(crate) fn accept(&mut self, record: &bam::Record) {
        if let Slot::Running(accumulator) = self
            && let Err(e) = accumulator.accept(record)
        {
            *self = Slot::Failed(e);
        }
    }

    pub(crate) fn finish(self) -> Result<A::Output> {
        match self {
            Slot::Running(accumulator) => accumulator.finish(),
            Slot::Failed(e) => Err(e),
        }
    }
}

/// One accumulator on its own thread. It takes every batch in file order,
/// so it sees the records exactly as a single-metric run does, then
/// finishes there too.
pub(crate) struct Worker<A: Accumulator> {
    sender: SyncSender<Batch>,
    handle: JoinHandle<Result<A::Output>>,
}

impl<A> Worker<A>
where
    A: Accumulator + Send + 'static,
    A::Output: Send + 'static,
{
    pub(crate) fn spawn(accumulator: A) -> Self {
        let (sender, batches) = mpsc::sync_channel::<Batch>(QUEUED_BATCHES);
        let handle = thread::spawn(move || {
            let mut slot = Slot::Running(accumulator);
            for batch in batches {
                for record in batch.iter() {
                    slot.accept(record);
                }
            }
            slot.finish()
        });
        Self { sender, handle }
    }

    pub(crate) fn send(&self, batch: &Batch) {
        // The receiver is gone only when the thread panicked; `finish`
        // re-raises that panic, so the batch can be dropped here.
        let _ = self.sender.send(Arc::clone(batch));
    }

    /// Close the queue and wait for the metric. A panic on its thread
    /// carries on here.
    pub(crate) fn finish(self) -> Result<A::Output> {
        drop(self.sender);
        self.handle
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
    }
}

/// Collectors of one pass, each on its own thread, fed in batches.
struct Fanout {
    batch: Vec<bam::Record>,
    workers: Vec<Worker<Box<dyn Collector>>>,
}

impl Fanout {
    fn new(collectors: Vec<Box<dyn Collector>>) -> Self {
        Self {
            batch: Vec::with_capacity(BATCH_RECORDS),
            workers: collectors.into_iter().map(Worker::spawn).collect(),
        }
    }

    fn flush(&mut self) {
        let batch: Batch = Arc::new(std::mem::replace(
            &mut self.batch,
            Vec::with_capacity(BATCH_RECORDS),
        ));
        for worker in &self.workers {
            worker.send(&batch);
        }
    }
}

impl Accumulator for Fanout {
    type Output = Vec<Result<ToolOutput>>;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        self.batch.push(record.clone());
        if self.batch.len() == BATCH_RECORDS {
            self.flush();
        }
        Ok(())
    }

    fn finish(mut self) -> Result<Vec<Result<ToolOutput>>> {
        if !self.batch.is_empty() {
            self.flush();
        }
        Ok(self.workers.into_iter().map(Worker::finish).collect())
    }
}

/// `SinglePassSamProgram.makeItSo` with `ASSUME_SORTED=false`: the header
/// must say `SO:coordinate`.
fn check_coordinate_sorted(path: &Path, header: &sam::Header) -> Result<()> {
    let order = header
        .header()
        .and_then(|h| h.other_fields().get(&SORT_ORDER))
        .map_or_else(|| "unsorted".to_owned(), ToString::to_string);
    if order == "coordinate" {
        return Ok(());
    }
    Err(Error::Usage(format!(
        "{} should be coordinate sorted but the header says the sort order is {order}; \
         pass ASSUME_SORTED=true to read it anyway",
        path.display()
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fails on its second record and counts the records it took.
    struct FailsSecond(usize);

    impl Accumulator for FailsSecond {
        type Output = usize;

        fn accept(&mut self, _: &bam::Record) -> Result<()> {
            self.0 += 1;
            if self.0 == 2 {
                return Err(Error::Metric("second record".to_owned()));
            }
            Ok(())
        }

        fn finish(self) -> Result<usize> {
            Ok(self.0)
        }
    }

    #[test]
    fn a_failed_slot_keeps_its_error_and_takes_no_more_records() {
        let record = bam::Record::default();
        let mut slot = Slot::Running(FailsSecond(0));
        for _ in 0..3 {
            slot.accept(&record);
        }
        assert!(matches!(slot.finish(), Err(Error::Metric(m)) if m == "second record"));
        let mut healthy = Slot::Running(FailsSecond(5));
        healthy.accept(&record);
        assert_eq!(healthy.finish().expect("ok"), 6);
    }

    #[test]
    fn a_boxed_collector_runs_its_finish_once_on_the_output() {
        let record = bam::Record::default();
        let mut collector = boxed(FailsSecond(0), |n| {
            assert_eq!(n, 1);
            Ok(ToolOutput::NONE)
        });
        collector.accept(&record).expect("first");
        assert!(Collector::finish(collector).expect("finish").main.is_none());
        let mut fanout = Fanout::new(vec![
            boxed(FailsSecond(0), |_| Ok(ToolOutput::NONE)),
            boxed(FailsSecond(10), |n| {
                assert_eq!(n, 12);
                Ok(ToolOutput::NONE)
            }),
        ]);
        fanout.accept(&record).expect("one");
        fanout.accept(&record).expect("two");
        let results = fanout.finish().expect("results");
        assert!(results[0].is_err());
        assert!(results[1].is_ok());
    }
}
