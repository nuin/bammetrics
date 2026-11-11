//! A sort that keeps at most a set number of records in memory and spills
//! sorted runs to files under a temporary directory, the role of htsjdk's
//! `SortingCollection` and `SortingLongCollection` in `MarkDuplicates`.
//!
//! Records are sorted stably, and runs merge with ties going to the earlier
//! run, so equal records come back in the order they were added, as
//! htsjdk's collections give them. A run is sorted by [`Spill::sort_run`]:
//! the standard library's stable sort for a total order, the port of Java's
//! `TimSort` ([`crate::timsort`]) for an order that may not be one. A run
//! file holds each record as [`Spill::encode`] writes it, plain or in BGZF
//! blocks ([`Sorter::compressed`]). The buffer holds at most a set number
//! of records and, with [`Sorter::byte_limit`], a set number of bytes. The
//! directory and every run file are created here and removed when the
//! [`TempDir`] is dropped; nothing else is deleted.

use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::sync::{Mutex, OnceLock};

use noodles::bgzf;

use crate::error::{Error, Result};

/// A record that can be written to and read back from a run file.
pub trait Spill: Sized {
    fn encode(&self, out: &mut Vec<u8>);
    /// The next record, `None` at the end of the file.
    fn decode(input: &mut impl BufRead) -> io::Result<Option<Self>>;

    /// Bytes the record holds in memory, counted against
    /// [`Sorter::byte_limit`].
    fn memory(&self) -> usize {
        std::mem::size_of::<Self>()
    }

    /// Sort one run in memory, stably. The default suits a total order. A
    /// type whose order may not be one sorts as Java does
    /// ([`crate::timsort`]), so it never panics.
    fn sort_run(run: &mut Vec<Self>) -> Result<()>
    where
        Self: Ord,
    {
        run.sort();
        Ok(())
    }
}

/// A directory of run files, removed with everything this process put in it.
/// The directory is created by the first [`TempDir::file`], so a run that
/// spills nothing never touches `TMP_DIR`. Its name holds the process ID
/// and a process-wide counter; a name already taken is skipped.
#[derive(Debug)]
pub struct TempDir {
    base: PathBuf,
    /// The directory, once created.
    path: OnceLock<PathBuf>,
    files: Mutex<Vec<PathBuf>>,
    next: AtomicU64,
}

/// Directories made by this process so far, for unique names.
static DIRECTORIES: AtomicU64 = AtomicU64::new(0);

impl TempDir {
    /// A directory to be made under `base` (the first `TMP_DIR`, or the
    /// system temporary directory).
    pub fn new(base: Option<&Path>) -> Self {
        Self {
            base: base.map_or_else(std::env::temp_dir, Path::to_path_buf),
            path: OnceLock::new(),
            files: Mutex::new(Vec::new()),
            next: AtomicU64::new(0),
        }
    }

    /// The directory, `None` before the first file.
    pub fn path(&self) -> Option<&Path> {
        self.path.get().map(PathBuf::as_path)
    }

    fn create(&self) -> Result<&Path> {
        if let Some(path) = self.path.get() {
            return Ok(path);
        }
        let error = |path: &Path| {
            let path = path.to_path_buf();
            move |source| Error::Output { path, source }
        };
        fs::create_dir_all(&self.base).map_err(error(&self.base))?;
        loop {
            let n = DIRECTORIES.fetch_add(1, AtomicOrdering::Relaxed);
            let path = self
                .base
                .join(format!("bammetrics-markdup.{}.{n}", std::process::id()));
            match create_private_dir(&path) {
                Ok(()) => {
                    crate::cleanup::register(&path);
                    return Ok(self.path.get_or_init(|| path));
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {}
                Err(e) => return Err(error(&path)(e)),
            }
        }
    }

    /// A new file name in the directory, recorded for removal. The first
    /// call creates the directory.
    pub fn file(&self, stem: &str) -> Result<PathBuf> {
        let dir = self.create()?;
        let n = self.next.fetch_add(1, AtomicOrdering::Relaxed);
        let path = dir.join(format!("{stem}.{n}.tmp"));
        crate::cleanup::register(&path);
        if let Ok(mut files) = self.files.lock() {
            files.push(path.clone());
        }
        Ok(path)
    }

    /// A new file in the directory, created for writing and readable by
    /// this user only (mode 0600): spilled read ends and spooled inputs hold
    /// read data, and `TMP_DIR` is often shared.
    pub fn create_file(&self, stem: &str) -> Result<(PathBuf, File)> {
        let path = self.file(stem)?;
        let file = create_private_file(&path).map_err(write_error(&path))?;
        Ok((path, file))
    }

    /// Remove one file made by [`TempDir::file`] before the directory goes.
    pub fn remove(&self, path: &Path) {
        if let Ok(mut files) = self.files.lock()
            && let Some(i) = files.iter().position(|f| f == path)
        {
            files.swap_remove(i);
            let _ = fs::remove_file(path);
            crate::cleanup::unregister(path);
        }
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        if let Ok(files) = self.files.lock() {
            for file in files.iter() {
                let _ = fs::remove_file(file);
                crate::cleanup::unregister(file);
            }
        }
        if let Some(path) = self.path.get() {
            let _ = fs::remove_dir(path);
            crate::cleanup::unregister(path);
        }
    }
}

/// A directory this user alone can enter (mode 0700). Fails when the path
/// exists.
fn create_private_dir(path: &Path) -> io::Result<()> {
    let mut builder = fs::DirBuilder::new();
    #[cfg(unix)]
    {
        use std::os::unix::fs::DirBuilderExt;
        builder.mode(0o700);
    }
    builder.create(path)
}

/// A new file this user alone can read (mode 0600). Fails when the path
/// exists, so a name planted in `TMP_DIR` is never written through.
fn create_private_file(path: &Path) -> io::Result<File> {
    let mut options = File::options();
    options.write(true).create_new(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.mode(0o600);
    }
    options.open(path)
}

fn write_error(path: &Path) -> impl Fn(io::Error) -> Error + '_ {
    move |source| Error::Output {
        path: path.to_path_buf(),
        source,
    }
}

/// How run files are written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RunFormat {
    Plain,
    /// BGZF blocks at this level, compressed on `workers` threads.
    Bgzf {
        level: u8,
        workers: usize,
    },
}

/// Records added in any order, given back sorted.
pub struct Sorter<'d, T> {
    dir: &'d TempDir,
    stem: &'static str,
    limit: usize,
    /// Most bytes ([`Spill::memory`]) the buffer holds; `usize::MAX` for no
    /// bound.
    byte_limit: usize,
    /// Bytes the buffer holds now.
    held: usize,
    format: RunFormat,
    buffer: Vec<T>,
    runs: Vec<PathBuf>,
    /// How many merges made each run: 0 for a spilled buffer.
    levels: Vec<u32>,
    count: u64,
    /// Bytes written to run files so far, merged runs included.
    spilled: u64,
    peak_records: usize,
    peak_bytes: usize,
    /// The deepest level of merges made while spilling.
    merge_levels: u32,
}

impl<'d, T: Spill + Ord> Sorter<'d, T> {
    /// A sorter holding at most `limit` records (at least 1) in memory.
    pub fn new(dir: &'d TempDir, stem: &'static str, limit: usize) -> Self {
        Self {
            dir,
            stem,
            limit: limit.max(1),
            byte_limit: usize::MAX,
            held: 0,
            format: RunFormat::Plain,
            buffer: Vec::new(),
            runs: Vec::new(),
            levels: Vec::new(),
            count: 0,
            spilled: 0,
            peak_records: 0,
            peak_bytes: 0,
            merge_levels: 0,
        }
    }

    /// Also spill before the buffer would hold more than `bytes` bytes of
    /// records ([`Spill::memory`]). A record larger than `bytes` is held
    /// alone.
    #[must_use]
    pub fn byte_limit(mut self, bytes: usize) -> Self {
        self.byte_limit = bytes.max(1);
        self
    }

    /// Write run files as BGZF blocks at `level` (0 to 9), compressed on
    /// `workers` threads.
    #[must_use]
    pub fn compressed(mut self, level: u8, workers: usize) -> Self {
        self.format = RunFormat::Bgzf {
            level: level.min(9),
            workers: workers.max(1),
        };
        self
    }

    pub fn push(&mut self, item: T) -> Result<()> {
        let size = item.memory();
        let over_bytes =
            !self.buffer.is_empty() && self.held.saturating_add(size) > self.byte_limit;
        if self.buffer.len() >= self.limit || over_bytes {
            self.spill()?;
        }
        self.buffer.push(item);
        self.held += size;
        self.count += 1;
        self.peak_records = self.peak_records.max(self.buffer.len());
        self.peak_bytes = self.peak_bytes.max(self.held);
        Ok(())
    }

    /// What the sort has done so far.
    pub fn stats(&self) -> SortStats {
        SortStats {
            records: self.count,
            runs: self.runs.len(),
            spilled_bytes: self.spilled,
            peak_records: self.peak_records,
            peak_bytes: self.peak_bytes,
            merge_levels: self.merge_levels,
        }
    }

    /// Records added so far.
    pub fn len(&self) -> u64 {
        self.count
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    fn spill(&mut self) -> Result<()> {
        T::sort_run(&mut self.buffer)?;
        let (path, file) = self.dir.create_file(self.stem)?;
        let mut out = RunWriter::new(file, self.format);
        let mut bytes = Vec::with_capacity(64);
        for item in self.buffer.drain(..) {
            bytes.clear();
            item.encode(&mut bytes);
            out.write_all(&bytes).map_err(write_error(&path))?;
        }
        self.held = 0;
        self.spilled += out.finish().map_err(write_error(&path))?;
        self.runs.push(path);
        self.levels.push(0);
        // The last MAX_OPEN_RUNS runs of one level merge into one run of
        // the next, so a small buffer cannot pile up thousands of files;
        // consecutive runs merge, so ties keep their order.
        while self.levels.len() >= MAX_OPEN_RUNS {
            let from = self.levels.len() - MAX_OPEN_RUNS;
            let level = self.levels[from];
            if self.levels[from..].iter().any(|l| *l != level) {
                break;
            }
            let group = self.runs[from..].to_vec();
            let merged = self.merge_into_run(&group)?;
            self.runs.truncate(from);
            self.levels.truncate(from);
            self.runs.push(merged);
            self.levels.push(level + 1);
            self.merge_levels = self.merge_levels.max(level + 1);
        }
        Ok(())
    }

    /// Every record in order.
    pub fn finish(self) -> Result<Sorted<'d, T>> {
        self.finish_with_stats().map(|(sorted, _)| sorted)
    }

    /// [`Sorter::finish`], and what the sort did with the last merge
    /// counted.
    pub fn finish_with_stats(mut self) -> Result<(Sorted<'d, T>, SortStats)> {
        if self.runs.is_empty() {
            T::sort_run(&mut self.buffer)?;
            let stats = self.stats();
            let items = std::mem::take(&mut self.buffer).into_iter();
            return Ok((
                Sorted {
                    source: Source::Memory(items),
                },
                stats,
            ));
        }
        if !self.buffer.is_empty() {
            self.spill()?;
        }
        // Consecutive runs merge into one until few enough files stay open;
        // a merged run takes its place in the list, so ties keep their order.
        let mut level = self.merge_levels;
        while self.runs.len() > MAX_OPEN_RUNS {
            let mut merged = Vec::with_capacity(self.runs.len() / MAX_OPEN_RUNS + 1);
            for group in std::mem::take(&mut self.runs).chunks(MAX_OPEN_RUNS) {
                merged.push(self.merge_into_run(group)?);
            }
            self.runs = merged;
            level += 1;
        }
        self.merge_levels = level + 1;
        self.levels.clear();
        let stats = self.stats();
        let runs = std::mem::take(&mut self.runs);
        Ok((
            Sorted {
                source: open_runs(self.dir, runs, self.format)?,
            },
            stats,
        ))
    }

    /// Merge `group` into one new run file; the files of `group` are removed.
    fn merge_into_run(&mut self, group: &[PathBuf]) -> Result<PathBuf> {
        if let [only] = group {
            return Ok(only.clone());
        }
        let mut source: Sorted<'_, T> = Sorted {
            source: open_runs(self.dir, group.to_vec(), self.format)?,
        };
        let (path, file) = self.dir.create_file(self.stem)?;
        let mut out = RunWriter::new(file, self.format);
        let mut bytes = Vec::with_capacity(64);
        while let Some(item) = source.next_item()? {
            bytes.clear();
            item.encode(&mut bytes);
            out.write_all(&bytes).map_err(write_error(&path))?;
        }
        self.spilled += out.finish().map_err(write_error(&path))?;
        Ok(path)
    }
}

/// What a [`Sorter`] did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SortStats {
    pub records: u64,
    /// Run files read by the last merge.
    pub runs: usize,
    /// Bytes written to run files, merged runs included.
    pub spilled_bytes: u64,
    /// Most records the buffer held at once.
    pub peak_records: usize,
    /// Most bytes ([`Spill::memory`]) the buffer held at once.
    pub peak_bytes: usize,
    /// 0 for a sort in memory, 1 for one merge of the spilled runs, one
    /// more per level of merges made before that.
    pub merge_levels: u32,
}

/// Most run files read at once; more are merged in groups first.
const MAX_OPEN_RUNS: usize = 64;

/// A run file being written: plain, or BGZF blocks.
enum RunWriter {
    Plain(BufWriter<File>),
    Bgzf(bgzf::io::MultithreadedWriter<BufWriter<File>>),
}

impl RunWriter {
    fn new(file: File, format: RunFormat) -> Self {
        let out = BufWriter::with_capacity(1 << 20, file);
        match format {
            RunFormat::Plain => Self::Plain(out),
            RunFormat::Bgzf { level, workers } => {
                let level = bgzf::io::writer::CompressionLevel::new(level)
                    .unwrap_or(bgzf::io::writer::CompressionLevel::FAST);
                let workers = NonZero::new(workers).unwrap_or(NonZero::<usize>::MIN);
                Self::Bgzf(
                    bgzf::io::multithreaded_writer::Builder::default()
                        .set_compression_level(level)
                        .set_worker_count(workers)
                        .build_from_writer(out),
                )
            }
        }
    }

    fn write_all(&mut self, bytes: &[u8]) -> io::Result<()> {
        match self {
            Self::Plain(out) => out.write_all(bytes),
            Self::Bgzf(out) => out.write_all(bytes),
        }
    }

    /// Flush the file; its size in bytes.
    fn finish(self) -> io::Result<u64> {
        let mut file = match self {
            Self::Plain(out) => out,
            Self::Bgzf(mut out) => out.finish()?,
        };
        file.flush()?;
        file.get_ref().metadata().map(|m| m.len())
    }
}

/// A run file being read.
enum RunReader {
    Plain(BufReader<File>),
    Bgzf(bgzf::io::Reader<BufReader<File>>),
}

impl Read for RunReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Plain(r) => r.read(buf),
            Self::Bgzf(r) => r.read(buf),
        }
    }
}

impl BufRead for RunReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Self::Plain(r) => r.fill_buf(),
            Self::Bgzf(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amount: usize) {
        match self {
            Self::Plain(r) => r.consume(amount),
            Self::Bgzf(r) => r.consume(amount),
        }
    }
}

/// The k-way merge over `paths`, ties to the earlier run.
fn open_runs<T: Spill + Ord>(
    dir: &TempDir,
    paths: Vec<PathBuf>,
    format: RunFormat,
) -> Result<Source<'_, T>> {
    let mut readers = Vec::with_capacity(paths.len());
    let mut heap = BinaryHeap::with_capacity(paths.len());
    for (run, path) in paths.iter().enumerate() {
        let file = File::open(path).map_err(|source| Error::Input {
            path: path.clone(),
            source,
        })?;
        let buffered = BufReader::with_capacity(1 << 18, file);
        let mut reader = match format {
            RunFormat::Plain => RunReader::Plain(buffered),
            RunFormat::Bgzf { .. } => RunReader::Bgzf(bgzf::io::Reader::new(buffered)),
        };
        if let Some(item) = T::decode(&mut reader).map_err(Error::Record)? {
            heap.push(Reverse(Head { item, run }));
        }
        readers.push(reader);
    }
    Ok(Source::Runs {
        dir,
        paths,
        readers,
        heap,
    })
}

/// The smallest unread record of one run; ties go to the earlier run.
struct Head<T> {
    item: T,
    run: usize,
}

impl<T: Ord> PartialEq for Head<T> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<T: Ord> Eq for Head<T> {}

impl<T: Ord> PartialOrd for Head<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: Ord> Ord for Head<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.item.cmp(&other.item).then(self.run.cmp(&other.run))
    }
}

enum Source<'d, T> {
    Memory(std::vec::IntoIter<T>),
    Runs {
        dir: &'d TempDir,
        paths: Vec<PathBuf>,
        readers: Vec<RunReader>,
        heap: BinaryHeap<Reverse<Head<T>>>,
    },
}

/// The sorted records; run files are removed once read.
pub struct Sorted<'d, T> {
    source: Source<'d, T>,
}

impl<T: Spill + Ord> Sorted<'_, T> {
    /// The next record in order.
    pub fn next_item(&mut self) -> Result<Option<T>> {
        match &mut self.source {
            Source::Memory(items) => Ok(items.next()),
            Source::Runs { readers, heap, .. } => {
                let Some(Reverse(Head { item, run })) = heap.pop() else {
                    return Ok(None);
                };
                if let Some(next) = T::decode(&mut readers[run]).map_err(Error::Record)? {
                    heap.push(Reverse(Head { item: next, run }));
                }
                Ok(Some(item))
            }
        }
    }
}

impl<T> Drop for Sorted<'_, T> {
    fn drop(&mut self) {
        if let Source::Runs { dir, paths, .. } = &self.source {
            for path in paths {
                dir.remove(path);
            }
        }
    }
}

/// LEB128 of an unsigned value.
pub fn put_u64(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

/// Zigzag LEB128 of a signed value.
pub fn put_i64(out: &mut Vec<u8>, v: i64) {
    put_u64(out, ((v << 1) ^ (v >> 63)) as u64);
}

/// The next LEB128 value; `None` at a clean end of input.
pub fn get_u64(input: &mut impl BufRead) -> io::Result<Option<u64>> {
    let mut value = 0u64;
    let mut shift = 0;
    loop {
        let mut byte = [0u8];
        if input.read(&mut byte)? == 0 {
            return if shift == 0 {
                Ok(None)
            } else {
                Err(io::ErrorKind::UnexpectedEof.into())
            };
        }
        value |= u64::from(byte[0] & 0x7f) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(Some(value));
        }
        shift += 7;
        if shift > 63 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "varint too long",
            ));
        }
    }
}

/// [`get_u64`] that must find a value.
pub fn need_u64(input: &mut impl BufRead) -> io::Result<u64> {
    get_u64(input)?.ok_or_else(|| io::ErrorKind::UnexpectedEof.into())
}

/// The zigzag value [`put_i64`] wrote.
pub fn need_i64(input: &mut impl BufRead) -> io::Result<i64> {
    let v = need_u64(input)?;
    Ok(((v >> 1) as i64) ^ -((v & 1) as i64))
}

impl Spill for u64 {
    fn encode(&self, out: &mut Vec<u8>) {
        put_u64(out, *self);
    }

    fn decode(input: &mut impl BufRead) -> io::Result<Option<Self>> {
        get_u64(input)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runs_merge_into_one_sorted_stream_and_are_removed() {
        let dir = TempDir::new(None);
        let mut sorter = Sorter::new(&dir, "n", 3);
        let values = [9u64, 1, 7, 3, 3, 0, 1 << 40, 5, 2, 8];
        for v in values {
            sorter.push(v).expect("push");
        }
        assert_eq!(sorter.len(), 10);
        let mut sorted = sorter.finish().expect("finish");
        let mut out = Vec::new();
        while let Some(v) = sorted.next_item().expect("next") {
            out.push(v);
        }
        let mut expected = values.to_vec();
        expected.sort_unstable();
        assert_eq!(out, expected);
        let runs: Vec<_> = fs::read_dir(dir.path().expect("created"))
            .expect("dir")
            .collect();
        assert_eq!(runs.len(), 4);
        drop(sorted);
        assert_eq!(
            fs::read_dir(dir.path().expect("created"))
                .expect("dir")
                .count(),
            0
        );
        let path = dir.path().expect("created").to_path_buf();
        drop(dir);
        assert!(!path.exists());
    }

    #[cfg(unix)]
    #[test]
    fn the_directory_and_its_files_belong_to_this_user_only() {
        use std::os::unix::fs::PermissionsExt;
        let dir = TempDir::new(None);
        let (path, file) = dir.create_file("n").expect("create");
        drop(file);
        let mode = |p: &Path| fs::metadata(p).expect("metadata").permissions().mode() & 0o777;
        assert_eq!(mode(dir.path().expect("created")), 0o700);
        assert_eq!(mode(&path), 0o600);
        // A name that already exists is never written through.
        assert!(create_private_file(&path).is_err());
    }

    #[test]
    fn many_runs_merge_in_groups() {
        let dir = TempDir::new(None);
        let mut sorter = Sorter::new(&dir, "many", 2);
        let values: Vec<u64> = (0..300u64).map(|i| (i * 7919) % 301).collect();
        for v in &values {
            sorter.push(*v).expect("push");
        }
        let mut sorted = sorter.finish().expect("finish");
        let mut out = Vec::new();
        while let Some(v) = sorted.next_item().expect("next") {
            out.push(v);
        }
        let mut expected = values;
        expected.sort_unstable();
        assert_eq!(out, expected);
        assert!(
            fs::read_dir(dir.path().expect("created"))
                .expect("dir")
                .count()
                <= MAX_OPEN_RUNS
        );
    }

    #[test]
    fn a_one_record_buffer_keeps_few_files() {
        let dir = TempDir::new(None);
        let mut sorter = Sorter::new(&dir, "ones", 1);
        let values: Vec<u64> = (0..5000u64).map(|i| ((i % 3) << 32) | i).collect();
        for v in &values {
            sorter.push(*v).expect("push");
            let files = dir
                .path()
                .map_or(0, |p| fs::read_dir(p).expect("dir").count());
            assert!(files < 3 * MAX_OPEN_RUNS, "{files} run files");
        }
        let mut sorted = sorter.finish().expect("finish");
        let mut out = Vec::new();
        while let Some(v) = sorted.next_item().expect("next") {
            out.push(v);
        }
        let mut expected = values;
        expected.sort_unstable();
        assert_eq!(out, expected);
    }

    #[test]
    fn nothing_is_created_until_a_run_spills() {
        let base = std::env::temp_dir().join(format!("bammetrics-lazy.{}", std::process::id()));
        let dir = TempDir::new(Some(&base));
        let mut sorter = Sorter::new(&dir, "lazy", 10);
        sorter.push(1u64).expect("push");
        let mut sorted = sorter.finish().expect("finish");
        assert_eq!(sorted.next_item().expect("next"), Some(1));
        assert!(dir.path().is_none() && !base.exists());
    }

    #[test]
    fn signed_varints_round_trip() {
        let mut bytes = Vec::new();
        for v in [0i64, -1, 1, i64::MIN, i64::MAX, -300, 300] {
            put_i64(&mut bytes, v);
        }
        let mut input = bytes.as_slice();
        for v in [0i64, -1, 1, i64::MIN, i64::MAX, -300, 300] {
            assert_eq!(need_i64(&mut input).expect("value"), v);
        }
        assert_eq!(get_u64(&mut input).expect("end"), None);
    }
}
