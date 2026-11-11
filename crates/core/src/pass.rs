//! The reader loop every metric shares. A BAM is read once; BGZF blocks are
//! inflated on worker threads (`noodles` `MultithreadedReader`), and records
//! are decoded and handed to an [`Accumulator`] on the calling thread, in
//! file order. Results therefore do not depend on the worker count. A CRAM
//! (told by its magic bytes) is decoded by [`crate::cram`] into the same
//! record stream, its containers on worker threads.
//!
//! Several tools can share one read ([`run_many`]): each one is a boxed
//! [`Collector`] fed from batches on a thread of its own, the way Picard's
//! `SinglePassSamProgram.makeItSo` hands every record to every program of
//! `CollectMultipleMetrics`. [`run_isolated`] is the same read for `qc`:
//! a collector that fails or panics fails alone, and the others finish.

use std::any::Any;
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::num::NonZero;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{self, SyncSender};
use std::thread::{self, JoinHandle};

use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::{bam, bgzf, sam};

use crate::cram;
use crate::error::{Error, Result};
use crate::metric::ToolOutput;

/// Cap on the default worker count. BGZF inflation stops being the
/// bottleneck well before this on the corpus benchmark (docs/benchmarks.md).
pub const MAX_DEFAULT_THREADS: usize = 8;

/// The most threads a run may ask for. Each read, each `write` and each
/// CRAM decode starts that many workers at once. A count the system cannot
/// start panics inside `noodles-bgzf` (exit 101), so it is refused here.
pub const MAX_THREADS: usize = 256;

/// `threads` when a run may ask for it, else what is wrong with it, for a
/// usage error that names the value in front (`--threads 300 is ...`).
/// Every place that parses a thread count (the verbs, `chain`, a plan file,
/// `qc`, `collect`, `threads=` on the Picard face, Python) asks here.
pub fn check_threads(threads: usize) -> std::result::Result<usize, String> {
    if threads > MAX_THREADS {
        return Err(format!("more than the maximum of {MAX_THREADS} threads"));
    }
    Ok(threads)
}

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
/// accumulator from the header, feed it every record and finish it. A
/// CRAM needs [`run_limited`] with its reference.
pub fn run<A, F>(path: &Path, threads: usize, build: F) -> Result<A::Output>
where
    A: Accumulator,
    F: FnOnce(&sam::Header) -> Result<A>,
{
    run_limited(path, None, threads, Limits::default(), build)
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

/// The header of the BAM or CRAM at `path`, no record read. A CRAM header
/// needs no reference.
pub fn read_header(path: &Path) -> Result<sam::Header> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let sniffed = cram::sniff(File::open(path).map_err(input)?).map_err(input)?;
    if sniffed.cram {
        noodles::cram::io::Reader::new(BufReader::new(sniffed.reader))
            .read_header()
            .map_err(|e| input(cram::truncated(e, true)))
    } else {
        bam::io::Reader::new(sniffed.reader)
            .read_header()
            .map_err(|e| input(cram::bam_header_error(e)))
    }
}

/// [`run`] with `STOP_AFTER` and `ASSUME_SORTED` applied, and the
/// `REFERENCE_SEQUENCE` a CRAM is decoded with (a BAM ignores it).
pub fn run_limited<A, F>(
    path: &Path,
    reference: Option<&Path>,
    threads: usize,
    limits: Limits,
    build: F,
) -> Result<A::Output>
where
    A: Accumulator,
    F: FnOnce(&sam::Header) -> Result<A>,
{
    read(path, reference, threads, limits, None, build)
}

/// Open the BAM or CRAM at `path` and hand `consume` its header and the
/// inflated record stream: each record as `block_size` and its bytes
/// ([`crate::bam_raw::RawRecord::read_from`]). A BAM is inflated by
/// `threads` BGZF workers, a CRAM decoded by `threads` container workers
/// with `reference`. The second flag is true for a CRAM.
pub fn with_record_stream<T>(
    path: &Path,
    reference: Option<&Path>,
    threads: usize,
    consume: impl FnOnce(&sam::Header, &mut (dyn Read + 'static), bool) -> Result<T>,
) -> Result<T> {
    let bytes = Arc::new(AtomicU64::new(0));
    with_counted_stream(path, reference, threads, &bytes, consume)
}

/// [`with_record_stream`], adding the compressed bytes read from the file
/// to `bytes` as the reader takes them.
pub(crate) fn with_counted_stream<T>(
    path: &Path,
    reference: Option<&Path>,
    threads: usize,
    bytes: &Arc<AtomicU64>,
    consume: impl FnOnce(&sam::Header, &mut (dyn Read + 'static), bool) -> Result<T>,
) -> Result<T> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = Counted {
        file: File::open(path).map_err(input)?,
        bytes: Arc::clone(bytes),
    };
    let sniffed = cram::sniff(file).map_err(input)?;
    let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
    if sniffed.cram {
        let mut reader =
            noodles::cram::io::Reader::new(cram::Peek::new(BufReader::new(sniffed.reader)));
        let header = reader
            .read_header()
            .map_err(|e| input(cram::truncated(e, true)))?;
        cram::decode(reader, &header, reference, workers, |records| {
            consume(&header, records.get_mut(), true)
        })
    } else {
        let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
            workers,
            sniffed.reader,
        ));
        let header = reader
            .read_header()
            .map_err(|e| input(cram::bam_header_error(e)))?;
        consume(&header, reader.get_mut(), false)
    }
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

/// The file, counting the bytes the BGZF or CRAM reader takes from it.
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

/// The read loop behind [`run_limited`] and [`run_isolated`]. A BAM is
/// inflated by `threads` BGZF workers; a CRAM ([`crate::cram`]) is decoded
/// by `threads` container workers into the same BAM record stream.
fn read<A, F>(
    path: &Path,
    reference: Option<&Path>,
    threads: usize,
    limits: Limits,
    progress: Option<OnProgress>,
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
    let sniffed = cram::sniff(file).map_err(input)?;
    let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
    let pass = Pass {
        path,
        cram: sniffed.cram,
        limits,
        progress,
        bytes: &bytes,
    };
    if sniffed.cram {
        let mut reader =
            noodles::cram::io::Reader::new(cram::Peek::new(BufReader::new(sniffed.reader)));
        let header = reader
            .read_header()
            .map_err(|e| input(cram::truncated(e, true)))?;
        cram::decode(reader, &header, reference, workers, |records| {
            pass.drive(records, &header, build)
        })
    } else {
        let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
            workers,
            sniffed.reader,
        ));
        let header = reader
            .read_header()
            .map_err(|e| input(cram::bam_header_error(e)))?;
        pass.drive(&mut reader, &header, build)
    }
}

/// One read's settings: the record loop over a BAM record stream.
struct Pass<'a, 'p> {
    path: &'a Path,
    cram: bool,
    limits: Limits,
    progress: Option<OnProgress<'p>>,
    bytes: &'a AtomicU64,
}

impl Pass<'_, '_> {
    /// Build the accumulator from `header`, feed it the records of
    /// `reader` in file order and finish it. A CRAM decoding error names
    /// the input file; a BAM one stays a bad record.
    fn drive<R, A, F>(
        mut self,
        reader: &mut bam::io::Reader<R>,
        header: &sam::Header,
        build: F,
    ) -> Result<A::Output>
    where
        R: Read,
        A: Accumulator,
        F: FnOnce(&sam::Header) -> Result<A>,
    {
        if !self.limits.assume_sorted {
            check_coordinate_sorted(self.path, header)?;
        }
        let mut accumulator = build(header)?;
        let mut record = bam::Record::default();
        let mut records = 0;
        let mut report = |records: u64| {
            if let Some(callback) = self.progress.as_mut() {
                callback(Progress {
                    records,
                    bytes: self.bytes.load(Ordering::Relaxed),
                });
            }
        };
        let read_error = |source: io::Error| read_error(self.path, self.cram, source);
        while !accumulator.done() && reader.read_record(&mut record).map_err(read_error)? != 0 {
            accumulator.accept(&record)?;
            records += 1;
            if records % BATCH_RECORDS as u64 == 0 {
                report(records);
            }
            if self.limits.stop_after > 0 && records >= self.limits.stop_after {
                break;
            }
        }
        report(records);
        accumulator.finish()
    }
}

/// The error of a failed record read. A CRAM decoding error names the
/// input file; a BAM that ends inside a record or a block was cut; any
/// other BAM error stays a bad record.
pub(crate) fn read_error(path: &Path, cram: bool, source: io::Error) -> Error {
    if cram {
        Error::Input {
            path: path.to_path_buf(),
            source,
        }
    } else if source.kind() == io::ErrorKind::UnexpectedEof {
        // A BAM ends with an empty BGZF block, so a stream that ends
        // inside a record or a block was cut.
        Error::Input {
            path: path.to_path_buf(),
            source: io::Error::new(
                source.kind(),
                "truncated BAM: the file ends inside a record or a block",
            ),
        }
    } else {
        Error::Record(source)
    }
}

/// One pass over `path` shared by every collector `build` returns, each on a
/// thread of its own. The outputs come back in the same order. The first
/// collector to fail fails the run, and no file of any of them is written,
/// as one exception ends a Picard `CollectMultipleMetrics` run.
pub fn run_many<F>(
    path: &Path,
    reference: Option<&Path>,
    threads: usize,
    limits: Limits,
    build: F,
) -> Result<Vec<ToolOutput>>
where
    F: FnOnce(&sam::Header) -> Result<Vec<Box<dyn Collector>>>,
{
    let results = run_limited(path, reference, threads, limits, |header| {
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
    reference: Option<&Path>,
    threads: usize,
    limits: Limits,
    progress: Option<OnProgress>,
    build: F,
) -> Result<Vec<Result<ToolOutput>>>
where
    F: FnOnce(&sam::Header) -> Vec<Result<Box<dyn Collector>>>,
{
    read(path, reference, threads, limits, progress, |header| {
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
pub(crate) fn panic_text(panic: &(dyn Any + Send)) -> String {
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

/// The header's `@HD SO` value, `unsorted` without one (htsjdk's
/// `getSortOrder` on a header without the tag).
pub fn sort_order(header: &sam::Header) -> String {
    header
        .header()
        .and_then(|h| h.other_fields().get(&SORT_ORDER))
        .map_or_else(|| "unsorted".to_owned(), ToString::to_string)
}

/// `SinglePassSamProgram.makeItSo` with `ASSUME_SORTED=false`: the header
/// must say `SO:coordinate`.
pub(crate) fn check_coordinate_sorted(path: &Path, header: &sam::Header) -> Result<()> {
    let order = sort_order(header);
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
            None,
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
        let err = run_isolated(&missing, None, 1, Limits::default(), None, |_| {
            vec![Ok(counting(false))]
        })
        .expect_err("no file");
        assert!(matches!(err, Error::Input { .. }));
    }

    /// A file that is not a BAM (a SAM text, an empty file) and a BAM cut
    /// inside its records name the failure, in place of `failed to fill
    /// whole buffer`.
    #[test]
    fn a_file_without_a_bam_header_and_a_cut_bam_name_the_failure() {
        use std::io::Write as _;
        use std::os::unix::fs::{DirBuilderExt, OpenOptionsExt};

        let dir = std::env::temp_dir().join(format!("bammetrics-pass-{}", std::process::id()));
        std::fs::DirBuilder::new()
            .mode(0o700)
            .create(&dir)
            .expect("temp dir");
        let write = |name: &str, bytes: &[u8]| {
            let path = dir.join(name);
            std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .mode(0o600)
                .open(&path)
                .and_then(|mut f| f.write_all(bytes))
                .expect(name);
            path
        };
        let sam = write(
            "hdr_only.sam",
            b"@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n",
        );
        let empty = write("empty.bam", b"");
        let whole = std::fs::read(fixture()).expect("fixture");
        let cut = write("cut.bam", &whole[..whole.len() - 40]);
        for (path, text) in [
            (&sam, "not a BAM or CRAM file"),
            (&empty, "not a BAM or CRAM file"),
            (&cut, "truncated BAM"),
        ] {
            let err = run(path, 1, |_| {
                Ok(Counter {
                    seen: 0,
                    panic_on_second: false,
                })
            })
            .expect_err("fails");
            assert!(
                matches!(err, Error::Input { .. }),
                "{}: {err}",
                path.display()
            );
            assert!(err.to_string().contains(text), "{}: {err}", path.display());
        }
        std::fs::remove_dir_all(&dir).expect("cleanup");
    }
}
