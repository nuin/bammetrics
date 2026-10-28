//! The reader loop every metric shares. A BAM is read once; BGZF blocks are
//! inflated on worker threads (`noodles` `MultithreadedReader`), and records
//! are decoded and handed to an [`Accumulator`] on the calling thread, in
//! file order. Results therefore do not depend on the worker count.

use std::fs::File;
use std::num::NonZero;
use std::path::Path;
use std::thread;

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
    let mut accumulator = build(&header)?;
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        accumulator.accept(&record)?;
    }
    accumulator.finish()
}
