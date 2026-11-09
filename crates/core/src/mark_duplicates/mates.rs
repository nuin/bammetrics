//! Mates waiting for their pair during the first read, keyed by the `RG`
//! ID plus the read name and stored under the reference of the mate.
//!
//! Query-name order: one map in memory, as Picard's
//! `MemoryBasedReadEndsForMarkDuplicatesMap` (mates are adjacent, so it
//! stays small); a key put twice replaces the first entry.
//!
//! Coordinate order: htsjdk's `CoordinateSortedPairInfoMap`. Memory holds
//! the entries whose mate is on the reference being read. An entry whose
//! mate is on a later reference is appended to that reference's file in
//! `TMP_DIR` (through a buffer of at most `MAX_RECORDS_IN_RAM` entries,
//! capped at 1,048,576), and the file is read back when the first paired
//! record of that reference arrives. Rules copied from htsjdk 5.0.0:
//! * A key put twice for the reference in memory is an error (`Putting
//!   value into PairInfoMap that already existed`); a key found twice when
//!   a file is read back is an error too (`Value was put into PairInfoMap
//!   more than once`), raised when that reference is reached.
//! * An entry for a reference already passed is never looked up again:
//!   htsjdk writes it to a file it never reads, bammetrics drops it. The
//!   entries left in memory when the read moves to a later reference are
//!   dropped the same way. Coordinate order makes references only grow; an
//!   input that breaks it fails when `OUTPUT` is written.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::PathBuf;

use super::read_ends::ReadEnds;
use crate::error::{Error, Result};
use crate::spill::{Spill, TempDir, get_u64, put_u64};

/// Most entries buffered for the reference files before they are written.
const MAX_BUFFERED: usize = 1 << 20;

pub enum Mates<'d> {
    Queryname(HashMap<(i32, Vec<u8>), ReadEnds>),
    Coordinate(ByReference<'d>),
}

impl<'d> Mates<'d> {
    pub fn new(queryname: bool, dir: &'d TempDir, records_in_ram: usize) -> Self {
        if queryname {
            return Self::Queryname(HashMap::new());
        }
        Self::Coordinate(ByReference {
            dir,
            current: None,
            in_ram: HashMap::new(),
            files: HashMap::new(),
            buffers: HashMap::new(),
            buffered: 0,
            max_buffered: records_in_ram.clamp(1, MAX_BUFFERED),
        })
    }

    /// The entry waiting under `key` on `reference` (the reference of the
    /// record now read), removed.
    pub fn remove(&mut self, reference: i32, key: &[u8]) -> Result<Option<ReadEnds>> {
        match self {
            Self::Queryname(map) => Ok(map.remove(&(reference, key.to_vec()))),
            Self::Coordinate(map) => {
                map.load(reference)?;
                Ok(map.in_ram.remove(key))
            }
        }
    }

    /// Let `ends` wait under `key` for a mate on `mate_reference`.
    pub fn put(&mut self, mate_reference: i32, key: Vec<u8>, ends: ReadEnds) -> Result<()> {
        match self {
            Self::Queryname(map) => {
                map.insert((mate_reference, key), ends);
                Ok(())
            }
            Self::Coordinate(map) => map.put(mate_reference, key, ends),
        }
    }
}

/// The coordinate-order map: one reference in memory, later ones on disk.
pub struct ByReference<'d> {
    dir: &'d TempDir,
    /// The reference whose entries are in memory.
    current: Option<i32>,
    in_ram: HashMap<Vec<u8>, ReadEnds>,
    /// The file of each later reference that has one.
    files: HashMap<i32, PathBuf>,
    /// Encoded entries of each later reference not yet in its file.
    buffers: HashMap<i32, (Vec<u8>, usize)>,
    buffered: usize,
    max_buffered: usize,
}

fn key_text(key: &[u8]) -> String {
    String::from_utf8_lossy(key).into_owned()
}

impl ByReference<'_> {
    fn put(&mut self, mate_reference: i32, key: Vec<u8>, ends: ReadEnds) -> Result<()> {
        if Some(mate_reference) == self.current {
            if self.in_ram.contains_key(&key) {
                return Err(Error::Metric(format!(
                    "Putting value into PairInfoMap that already existed. {mate_reference}: {}",
                    key_text(&key)
                )));
            }
            self.in_ram.insert(key, ends);
            return Ok(());
        }
        if self.current.is_some_and(|current| mate_reference < current) {
            return Ok(());
        }
        let (buffer, entries) = self.buffers.entry(mate_reference).or_default();
        put_u64(buffer, key.len() as u64);
        buffer.extend_from_slice(&key);
        ends.encode(buffer);
        *entries += 1;
        self.buffered += 1;
        if self.buffered >= self.max_buffered {
            self.flush()?;
        }
        Ok(())
    }

    /// Append every buffer to its reference's file.
    fn flush(&mut self) -> Result<()> {
        for (reference, (bytes, _)) in self.buffers.drain() {
            let path = match self.files.get(&reference) {
                Some(path) => path.clone(),
                None => {
                    let path = self.dir.file("mates")?;
                    self.files.insert(reference, path.clone());
                    path
                }
            };
            let error = |source| Error::Output {
                path: path.clone(),
                source,
            };
            File::options()
                .create(true)
                .append(true)
                .open(&path)
                .and_then(|mut file| file.write_all(&bytes))
                .map_err(error)?;
        }
        self.buffered = 0;
        Ok(())
    }

    /// `ensureSequenceLoaded`: make `reference` the one in memory, with the
    /// entries its file and buffer hold.
    fn load(&mut self, reference: i32) -> Result<()> {
        if self.current == Some(reference) {
            return Ok(());
        }
        self.in_ram.clear();
        self.current = Some(reference);
        let path = self.files.remove(&reference);
        let (buffer, entries) = self.buffers.remove(&reference).unwrap_or_default();
        self.buffered -= entries;
        self.drop_passed(reference);
        let loaded = self.read_back(reference, path.as_ref(), &buffer);
        if let Some(path) = &path {
            self.dir.remove(path);
        }
        loaded
    }

    /// Drop the waiting entries of references before `reference`: no
    /// record of theirs comes any more.
    fn drop_passed(&mut self, reference: i32) {
        let passed: Vec<i32> = self
            .files
            .keys()
            .chain(self.buffers.keys())
            .copied()
            .filter(|r| *r < reference)
            .collect();
        for r in passed {
            if let Some(path) = self.files.remove(&r) {
                self.dir.remove(&path);
            }
            if let Some((_, entries)) = self.buffers.remove(&r) {
                self.buffered -= entries;
            }
        }
    }

    /// Put the entries of `path`, then of `buffer`, in memory.
    fn read_back(&mut self, reference: i32, path: Option<&PathBuf>, buffer: &[u8]) -> Result<()> {
        let error = |source| match path {
            Some(path) => Error::Input {
                path: path.clone(),
                source,
            },
            None => Error::Record(source),
        };
        let file: Box<dyn Read> = match path {
            Some(path) => Box::new(File::open(path).map_err(error)?),
            None => Box::new(std::io::empty()),
        };
        let mut input = BufReader::with_capacity(1 << 18, file.chain(buffer));
        while let Some(length) = get_u64(&mut input).map_err(error)? {
            let mut key = vec![0; usize::try_from(length).unwrap_or(usize::MAX)];
            input.read_exact(&mut key).map_err(error)?;
            let ends = ReadEnds::decode(&mut input)
                .map_err(error)?
                .ok_or_else(|| error(std::io::ErrorKind::UnexpectedEof.into()))?;
            if self.in_ram.contains_key(&key) {
                return Err(Error::Metric(format!(
                    "Value was put into PairInfoMap more than once.  {reference}: {}",
                    key_text(&key)
                )));
            }
            self.in_ram.insert(key, ends);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ends(index: i64) -> ReadEnds {
        ReadEnds {
            read1_index: index,
            ..ReadEnds::default()
        }
    }

    #[test]
    fn mates_on_later_references_wait_on_disk() {
        let dir = TempDir::new(None);
        let mut mates = Mates::new(false, &dir, 2);
        assert_eq!(mates.remove(0, b"a").expect("remove"), None);
        for (i, name) in [b"a", b"b", b"c"].iter().enumerate() {
            mates.put(2, name.to_vec(), ends(i as i64)).expect("put");
        }
        mates.put(0, b"d".to_vec(), ends(9)).expect("put");
        assert!(dir.path().is_some(), "two entries filled the buffer");
        assert_eq!(
            mates.remove(1, b"d").expect("remove"),
            None,
            "d was dropped"
        );
        let b = mates
            .remove(2, b"b")
            .expect("remove")
            .expect("b waits on 2");
        assert_eq!(b.read1_index, 1);
        assert_eq!(
            mates
                .remove(2, b"c")
                .expect("remove")
                .map(|e| e.read1_index),
            Some(2)
        );
        let left = std::fs::read_dir(dir.path().expect("dir"))
            .expect("dir")
            .count();
        assert_eq!(left, 0, "the file of reference 2 was removed once read");
    }

    #[test]
    fn a_key_twice_fails_where_htsjdk_fails() {
        let dir = TempDir::new(None);
        let mut mates = Mates::new(false, &dir, 100);
        mates.remove(0, b"x").expect("load 0");
        mates.put(0, b"k".to_vec(), ends(0)).expect("put");
        let twice = mates.put(0, b"k".to_vec(), ends(1)).expect_err("in memory");
        assert!(
            twice
                .to_string()
                .starts_with("Putting value into PairInfoMap")
        );
        mates.put(1, b"k".to_vec(), ends(2)).expect("put");
        mates
            .put(1, b"k".to_vec(), ends(3))
            .expect("put, checked later");
        let loaded = mates.remove(1, b"k").expect_err("read back");
        assert_eq!(
            loaded.to_string(),
            "Value was put into PairInfoMap more than once.  1: k"
        );
    }

    #[test]
    fn query_name_order_replaces() {
        let dir = TempDir::new(None);
        let mut mates = Mates::new(true, &dir, 1);
        mates.put(0, b"k".to_vec(), ends(0)).expect("put");
        mates.put(0, b"k".to_vec(), ends(1)).expect("put");
        assert_eq!(
            mates
                .remove(0, b"k")
                .expect("remove")
                .map(|e| e.read1_index),
            Some(1)
        );
        assert!(dir.path().is_none());
    }
}
