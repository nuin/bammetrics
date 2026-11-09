//! A sort that keeps at most a set number of records in memory and spills
//! sorted runs to files under a temporary directory, the role of htsjdk's
//! `SortingCollection` and `SortingLongCollection` in `MarkDuplicates`.
//!
//! Records are sorted stably, and runs merge with ties going to the earlier
//! run, so equal records come back in the order they were added, as
//! htsjdk's collections give them. A run file holds each record as LEB128
//! varints ([`Spill`]). The directory and every run file are created here
//! and removed when the [`TempDir`] is dropped; nothing else is deleted.

use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

use crate::error::{Error, Result};

/// A record that can be written to and read back from a run file.
pub trait Spill: Sized {
    fn encode(&self, out: &mut Vec<u8>);
    /// The next record, `None` at the end of the file.
    fn decode(input: &mut impl BufRead) -> io::Result<Option<Self>>;
}

/// A directory of run files, removed with everything this process put in it.
#[derive(Debug)]
pub struct TempDir {
    path: PathBuf,
    files: Mutex<Vec<PathBuf>>,
    next: AtomicU64,
}

impl TempDir {
    /// A new directory under `base` (the first `TMP_DIR`, or the system
    /// temporary directory).
    pub fn new(base: Option<&Path>) -> Result<Self> {
        let base = base.map_or_else(std::env::temp_dir, Path::to_path_buf);
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let path = base.join(format!("bammetrics-markdup.{}.{stamp}", std::process::id()));
        fs::create_dir_all(&base)
            .and_then(|()| fs::create_dir(&path))
            .map_err(|source| Error::Output {
                path: path.clone(),
                source,
            })?;
        Ok(Self {
            path,
            files: Mutex::new(Vec::new()),
            next: AtomicU64::new(0),
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// A new file name in the directory, recorded for removal.
    pub fn file(&self, stem: &str) -> PathBuf {
        let n = self.next.fetch_add(1, AtomicOrdering::Relaxed);
        let path = self.path.join(format!("{stem}.{n}.tmp"));
        if let Ok(mut files) = self.files.lock() {
            files.push(path.clone());
        }
        path
    }

    /// Remove one file made by [`TempDir::file`] before the directory goes.
    pub fn remove(&self, path: &Path) {
        if let Ok(mut files) = self.files.lock()
            && let Some(i) = files.iter().position(|f| f == path)
        {
            files.swap_remove(i);
            let _ = fs::remove_file(path);
        }
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        if let Ok(files) = self.files.lock() {
            for file in files.iter() {
                let _ = fs::remove_file(file);
            }
        }
        let _ = fs::remove_dir(&self.path);
    }
}

fn write_error(path: &Path) -> impl Fn(io::Error) -> Error + '_ {
    move |source| Error::Output {
        path: path.to_path_buf(),
        source,
    }
}

/// Records added in any order, given back sorted.
pub struct Sorter<'d, T> {
    dir: &'d TempDir,
    stem: &'static str,
    limit: usize,
    buffer: Vec<T>,
    runs: Vec<PathBuf>,
    count: u64,
}

impl<'d, T: Spill + Ord> Sorter<'d, T> {
    /// A sorter holding at most `limit` records (at least 1) in memory.
    pub fn new(dir: &'d TempDir, stem: &'static str, limit: usize) -> Self {
        Self {
            dir,
            stem,
            limit: limit.max(1),
            buffer: Vec::new(),
            runs: Vec::new(),
            count: 0,
        }
    }

    pub fn push(&mut self, item: T) -> Result<()> {
        if self.buffer.len() >= self.limit {
            self.spill()?;
        }
        self.buffer.push(item);
        self.count += 1;
        Ok(())
    }

    /// Records added so far.
    pub fn len(&self) -> u64 {
        self.count
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    fn spill(&mut self) -> Result<()> {
        self.buffer.sort();
        let path = self.dir.file(self.stem);
        let file = File::create(&path).map_err(write_error(&path))?;
        let mut out = BufWriter::with_capacity(1 << 20, file);
        let mut bytes = Vec::with_capacity(64);
        for item in self.buffer.drain(..) {
            bytes.clear();
            item.encode(&mut bytes);
            out.write_all(&bytes).map_err(write_error(&path))?;
        }
        out.flush().map_err(write_error(&path))?;
        self.runs.push(path);
        Ok(())
    }

    /// Every record in order.
    pub fn finish(mut self) -> Result<Sorted<'d, T>> {
        if self.runs.is_empty() {
            self.buffer.sort();
            return Ok(Sorted {
                source: Source::Memory(self.buffer.into_iter()),
            });
        }
        if !self.buffer.is_empty() {
            self.spill()?;
        }
        let mut readers = Vec::with_capacity(self.runs.len());
        let mut heap = BinaryHeap::with_capacity(self.runs.len());
        for (run, path) in self.runs.iter().enumerate() {
            let file = File::open(path).map_err(|source| Error::Input {
                path: path.clone(),
                source,
            })?;
            let mut reader = BufReader::with_capacity(1 << 18, file);
            if let Some(item) = T::decode(&mut reader).map_err(Error::Record)? {
                heap.push(Reverse(Head { item, run }));
            }
            readers.push(reader);
        }
        Ok(Sorted {
            source: Source::Runs {
                dir: self.dir,
                paths: std::mem::take(&mut self.runs),
                readers,
                heap,
            },
        })
    }
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
        readers: Vec<BufReader<File>>,
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
        let dir = TempDir::new(None).expect("temp dir");
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
        let runs: Vec<_> = fs::read_dir(dir.path()).expect("dir").collect();
        assert_eq!(runs.len(), 4);
        drop(sorted);
        assert_eq!(fs::read_dir(dir.path()).expect("dir").count(), 0);
        let path = dir.path().to_path_buf();
        drop(dir);
        assert!(!path.exists());
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
