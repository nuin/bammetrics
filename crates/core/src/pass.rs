//! The reader loop every metric shares. A BAM is read once; BGZF blocks are
//! inflated on worker threads (`noodles` `MultithreadedReader`), and records
//! are decoded and handed to an [`Accumulator`] on the calling thread, in
//! file order. Results therefore do not depend on the worker count.
//!
//! Several tools can share one read ([`run_many`]): each one is a boxed
//! [`Collector`] fed from batches on a thread of its own, the way Picard's
//! `SinglePassSamProgram.makeItSo` hands every record to every program of
//! `CollectMultipleMetrics`. [`run_isolated`] is the same read for `qc`:
//! a collector that fails or panics fails alone, and the others finish.

use std::any::Any;
use std::fs::File;
use std::io::{self, Read};
use std::num::NonZero;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
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

/// The header of the BAM at `path`, no record read.
pub fn read_header(path: &Path) -> Result<sam::Header> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = File::open(path).map_err(input)?;
    bam::io::Reader::new(file).read_header().map_err(input)
}

/// [`run`] with `STOP_AFTER` and `ASSUME_SORTED` applied.
pub fn run_limited<A, F>(path: &Path, threads: usize, limits: Limits, build: F) -> Result<A::Output>
where
    A: Accumulator,
    F: FnOnce(&sam::Header) -> Result<A>,
{
    read(path, threads, limits, None, build)
}

/// How far a pass has got: records handed to the accumulator and
/// compressed bytes read from the file.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Progress {
    pub records: u64,
    pub bytes: u64,
}

/// A progress callback: called every [`BATCH_RECORDS`] records and once
/// after the last.
pub type OnProgress<'a> = &'a mut dyn FnMut(Progress);

/// The file, counting the bytes the BGZF reader takes from it.
struct Counted {
    file: File,
    bytes: Arc<AtomicU64>,
}

impl Read for Counted {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.file.read(buf)?;
        self.bytes.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

/// The read loop behind [`run_limited`] and [`run_isolated`].
fn read<A, F>(
    path: &Path,
    threads: usize,
    limits: Limits,
    mut progress: Option<OnProgress>,
    build: F,
) -> Result<A::Output>
where
    A: Accumulator,
    F: FnOnce(&sam::Header) -> Result<A>,
{
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let bytes = Arc::new(AtomicU64::new(0));
    let file = Counted {
        file: File::open(path).map_err(input)?,
        bytes: Arc::clone(&bytes),
    };
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
    let mut report = |records: u64| {
        if let Some(callback) = progress.as_mut() {
            callback(Progress {
                records,
                bytes: bytes.load(Ordering::Relaxed),
            });
        }
    };
    while !accumulator.done() && reader.read_record(&mut record)? != 0 {
        accumulator.accept(&record)?;
        records += 1;
        if records % BATCH_RECORDS as u64 == 0 {
            report(records);
        }
        if limits.stop_after > 0 && records >= limits.stop_after {
            break;
        }
    }
    report(records);
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
        Ok(Fanout::new(
            build(header)?.into_iter().map(Ok).collect(),
            false,
        ))
    })?;
    results.into_iter().collect()
}

/// One pass over `path` shared by the collectors `build` returns, each
/// isolated from the others: a collector that fails to build, fails on a
/// record, fails to finish or panics gives its own `Err`, and the others
/// finish. Only a failure to open or read the BAM fails the pass. The
/// results come back in the order of `build`.
pub fn run_isolated<F>(
    path: &Path,
    threads: usize,
    limits: Limits,
    progress: Option<OnProgress>,
    build: F,
) -> Result<Vec<Result<ToolOutput>>>
where
    F: FnOnce(&sam::Header) -> Vec<Result<Box<dyn Collector>>>,
{
    read(path, threads, limits, progress, |header| {
        Ok(Fanout::new(build(header), true))
    })
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
        self.join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
    }

    /// Close the queue and wait for the thread: the metric, or the panic
    /// that ended the thread.
    fn join(self) -> thread::Result<Result<A::Output>> {
        drop(self.sender);
        self.handle.join()
    }
}

/// The text of a panic payload, for the error that replaces it.
fn panic_text(panic: &(dyn Any + Send)) -> String {
    panic
        .downcast_ref::<&str>()
        .map(|s| (*s).to_owned())
        .or_else(|| panic.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "no message".to_owned())
}

/// Collectors of one pass, each on its own thread, fed in batches. A
/// collector that failed to build holds its error and gets no batch.
struct Fanout {
    batch: Vec<bam::Record>,
    workers: Vec<Result<Worker<Box<dyn Collector>>>>,
    /// True: a panic on a collector thread becomes that collector's error.
    /// False: it carries on to the caller, as a Picard exception would.
    isolate: bool,
}

impl Fanout {
    fn new(collectors: Vec<Result<Box<dyn Collector>>>, isolate: bool) -> Self {
        Self {
            batch: Vec::with_capacity(BATCH_RECORDS),
            workers: collectors
                .into_iter()
                .map(|c| c.map(Worker::spawn))
                .collect(),
            isolate,
        }
    }

    fn flush(&mut self) {
        let batch: Batch = Arc::new(std::mem::replace(
            &mut self.batch,
            Vec::with_capacity(BATCH_RECORDS),
        ));
        for worker in self.workers.iter().flatten() {
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
        let isolate = self.isolate;
        Ok(self
            .workers
            .into_iter()
            .map(|worker| {
                let worker = worker?;
                if !isolate {
                    return worker.finish();
                }
                worker.join().unwrap_or_else(|panic| {
                    Err(Error::Metric(format!(
                        "the collector stopped: {}",
                        panic_text(panic.as_ref())
                    )))
                })
            })
            .collect())
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
        let mut fanout = Fanout::new(
            vec![
                Ok(boxed(FailsSecond(0), |_| Ok(ToolOutput::NONE))),
                Ok(boxed(FailsSecond(10), |n| {
                    assert_eq!(n, 12);
                    Ok(ToolOutput::NONE)
                })),
            ],
            false,
        );
        fanout.accept(&record).expect("one");
        fanout.accept(&record).expect("two");
        let results = fanout.finish().expect("results");
        assert!(results[0].is_err());
        assert!(results[1].is_ok());
    }

    /// Counts records and panics on the second when asked to.
    struct Counter {
        seen: usize,
        panic_on_second: bool,
    }

    impl Accumulator for Counter {
        type Output = usize;

        fn accept(&mut self, _: &bam::Record) -> Result<()> {
            self.seen += 1;
            assert!(!(self.panic_on_second && self.seen == 2), "planted panic");
            Ok(())
        }

        fn finish(self) -> Result<usize> {
            Ok(self.seen)
        }
    }

    fn fixture() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/insert_size.bam")
    }

    /// A collector whose finish reports its record count as a note.
    fn counting(panic_on_second: bool) -> Box<dyn Collector> {
        boxed(
            Counter {
                seen: 0,
                panic_on_second,
            },
            |n| Err(Error::Metric(format!("{n} records"))),
        )
    }

    #[test]
    fn an_isolated_pass_keeps_each_failure_to_its_own_collector() {
        let mut updates = Vec::new();
        let mut on_progress = |p: Progress| updates.push(p);
        let results = run_isolated(
            &fixture(),
            2,
            Limits::default(),
            Some(&mut on_progress),
            |_| {
                vec![
                    Ok(boxed(FailsSecond(0), |_| Ok(ToolOutput::NONE))),
                    Err(Error::Usage("no such interval file".to_owned())),
                    Ok(counting(true)),
                    Ok(counting(false)),
                ]
            },
        )
        .expect("the BAM reads");
        assert_eq!(results.len(), 4);
        assert!(matches!(&results[0], Err(Error::Metric(m)) if m == "second record"));
        assert!(matches!(&results[1], Err(Error::Usage(m)) if m == "no such interval file"));
        assert!(
            matches!(&results[2], Err(Error::Metric(m)) if m == "the collector stopped: planted panic")
        );
        assert!(matches!(&results[3], Err(Error::Metric(m)) if m == "10 records"));
        let last = updates.last().expect("a final update");
        assert_eq!(last.records, 10);
        assert!(last.bytes > 0);
    }

    #[test]
    fn an_unreadable_bam_still_fails_the_isolated_pass() {
        let missing = fixture().with_file_name("no_such.bam");
        let err = run_isolated(&missing, 1, Limits::default(), None, |_| {
            vec![Ok(counting(false))]
        })
        .expect_err("no file");
        assert!(matches!(err, Error::Input { .. }));
    }
}
