//! The reader loop every metric shares. A BAM is read once; BGZF blocks are
//! inflated on worker threads (`noodles` `MultithreadedReader`), and records
//! are decoded and handed to an [`Accumulator`] on the calling thread, in
//! file order. Results therefore do not depend on the worker count.

use std::fs::File;
use std::num::NonZero;
use std::path::Path;
use std::thread;

use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::{bam, bgzf, sam};

use crate::error::{Error, Result};

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
