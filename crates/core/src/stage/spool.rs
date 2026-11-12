//! The spool of a barrier that replays (`MarkDuplicates`) when its input
//! is not the source itself: the records it takes in its first read are
//! also written to a file in the chain's temporary directory, and its
//! second read replays that file (docs/07 section 2.3, replay or spool).
//!
//! The file is the record stream of a BAM without its header: each record
//! as `block_size` and its bytes, in BGZF blocks at level 3
//! ([`crate::sam_sort::SPILL_LEVEL`]), compressed on the chain's threads.
//! It is made by [`TempDir::create_file`] (mode 0600, a new name) and
//! removed when the replay is done, or with the directory.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::bgzf;

use super::{BATCH_RECORDS, Records};
use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::spill::TempDir;

/// A spool being written.
pub(super) struct SpoolWriter {
    path: PathBuf,
    out: bgzf::io::MultithreadedWriter<BufWriter<File>>,
}

impl SpoolWriter {
    pub(super) fn create(temp: &TempDir, threads: usize) -> Result<Self> {
        let (path, file) = temp.create_file("spool")?;
        let level = bgzf::io::writer::CompressionLevel::new(crate::sam_sort::SPILL_LEVEL)
            .unwrap_or(bgzf::io::writer::CompressionLevel::FAST);
        let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
        let out = bgzf::io::multithreaded_writer::Builder::default()
            .set_compression_level(level)
            .set_worker_count(workers)
            .build_from_writer(BufWriter::with_capacity(1 << 20, file));
        Ok(Self { path, out })
    }

    fn error(&self, source: std::io::Error) -> Error {
        Error::Output {
            path: self.path.clone(),
            source,
        }
    }

    pub(super) fn write(&mut self, batch: &[RawRecord]) -> Result<()> {
        for record in batch {
            record
                .write_to(&mut self.out)
                .map_err(|source| self.error(source))?;
        }
        Ok(())
    }

    /// Flush the last block; the file to replay.
    pub(super) fn finish(mut self) -> Result<PathBuf> {
        let mut file = self.out.finish().map_err(|source| Error::Output {
            path: self.path.clone(),
            source,
        })?;
        file.flush().map_err(|source| self.error(source))?;
        Ok(self.path)
    }
}

/// A spool read back in batches.
pub(super) struct SpoolRecords {
    path: PathBuf,
    reader: bgzf::io::MultithreadedReader<File>,
}

impl SpoolRecords {
    pub(super) fn open(path: &Path, threads: usize) -> Result<Self> {
        let file = File::open(path).map_err(|source| Error::Input {
            path: path.to_path_buf(),
            source,
        })?;
        let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
        Ok(Self {
            path: path.to_path_buf(),
            reader: bgzf::io::MultithreadedReader::with_worker_count(workers, file),
        })
    }
}

impl Records for SpoolRecords {
    fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool> {
        let mut filled = 0;
        while filled < BATCH_RECORDS {
            if filled == out.len() {
                out.push(RawRecord::default());
            }
            let more = out[filled]
                .read_from(&mut self.reader)
                .map_err(|source| Error::Input {
                    path: self.path.clone(),
                    source,
                })?;
            if !more {
                break;
            }
            filled += 1;
        }
        out.truncate(filled);
        Ok(filled > 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_spool_gives_back_its_records_and_is_private() {
        let temp = TempDir::new(None);
        let records: Vec<RawRecord> = (0..3000)
            .map(|i| crate::bam_raw::build(&format!("r{i}"), 0, 0, i, &[], &[], b""))
            .collect();
        let mut spool = SpoolWriter::create(&temp, 2).expect("create");
        for batch in records.chunks(700) {
            spool.write(batch).expect("write");
        }
        let path = spool.finish().expect("finish");
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mode = std::fs::metadata(&path).expect("meta").permissions().mode();
            assert_eq!(mode & 0o777, 0o600);
        }
        let mut back = SpoolRecords::open(&path, 2).expect("open");
        let mut all = Vec::new();
        let mut batch = Vec::new();
        while back.next_batch(&mut batch).expect("read") {
            assert!(batch.len() <= BATCH_RECORDS);
            all.append(&mut batch);
        }
        assert_eq!(all, records);
        temp.remove(&path);
        assert!(!path.exists());
    }
}
