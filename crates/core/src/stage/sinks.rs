//! The sinks of a chain: the tap adapter over a [`Collector`] and the
//! `write` step.

use std::io::{self, Read};
use std::path::{Path, PathBuf};

use noodles::{bam, sam};

use super::{Frame, Products, Sink};
use crate::bam_raw::RawRecord;
use crate::bam_writer::{self, BamWriter, Companions};
use crate::error::{Error, Result};
use crate::mark_duplicates;
use crate::metric::ToolOutput;
use crate::pass::Collector;

/// A collector that takes its first `left` records and ignores the rest:
/// a tap's `STOP_AFTER` on a read shared with other taps.
pub(super) struct Limited {
    pub(super) inner: Box<dyn Collector>,
    pub(super) left: u64,
}

impl Collector for Limited {
    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if self.left == 0 {
            return Ok(());
        }
        self.left -= 1;
        self.inner.accept(record)
    }

    fn finish(self: Box<Self>) -> Result<ToolOutput> {
        self.inner.finish()
    }
}

/// `block_size` and then the bytes of one record, the form
/// `bam::io::Reader::read_record` reads.
struct Framed<'a> {
    size: [u8; 4],
    body: &'a [u8],
    at: usize,
}

impl Read for Framed<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let (from, offset) = if self.at < 4 {
            (self.size.as_slice(), self.at)
        } else {
            (self.body, self.at - 4)
        };
        let n = buf.len().min(from.len() - offset);
        buf[..n].copy_from_slice(&from[offset..offset + n]);
        self.at += n;
        Ok(n)
    }
}

/// Fill `record` from the bytes of `raw`.
fn fill(record: &mut bam::Record, raw: &RawRecord) -> io::Result<()> {
    let body = raw.as_bytes();
    let size = u32::try_from(body.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        .to_le_bytes();
    let mut reader = bam::io::Reader::from(Framed { size, body, at: 0 });
    reader.read_record(record).map(|_| ())
}

/// The tap adapter: an existing [`Collector`] as a [`Sink`]. Each raw
/// record is copied into one reused `bam::Record` and accepted, so the
/// collector sees what it sees when it reads a BAM holding these records.
pub struct CollectorSink {
    collector: Box<dyn Collector>,
    record: bam::Record,
    /// Records still to take under `STOP_AFTER`; `None` takes all.
    left: Option<u64>,
}

impl CollectorSink {
    /// `stop_after` 0 takes every record.
    pub fn new(collector: Box<dyn Collector>, stop_after: u64) -> Self {
        Self {
            collector,
            record: bam::Record::default(),
            left: (stop_after > 0).then_some(stop_after),
        }
    }
}

// A collector is built on the chain's thread and moved to the tap's thread
// once, before its first record, as `pass` moves it to its worker.
impl Sink for CollectorSink {
    fn accept_batch(&mut self, batch: &[RawRecord]) -> Result<()> {
        for raw in batch {
            if self.left == Some(0) {
                break;
            }
            fill(&mut self.record, raw)?;
            self.collector.accept(&self.record)?;
            if let Some(left) = self.left.as_mut() {
                *left -= 1;
            }
        }
        Ok(())
    }

    fn done(&self) -> bool {
        self.left == Some(0)
    }

    fn finish(self: Box<Self>) -> Result<Products> {
        self.collector.finish().map(Products::metrics)
    }
}

/// The file the `write` step writes.
#[derive(Clone, Debug)]
pub struct WriteOptions {
    /// The BAM. Written under a temporary name and renamed when complete.
    pub path: PathBuf,
    /// BGZF level, 0 to 9. Picard's default is 5.
    pub compression_level: u8,
    /// A `.bai` beside a coordinate-sorted BAM.
    pub create_index: bool,
    /// `<path>.md5`.
    pub create_md5: bool,
}

impl WriteOptions {
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            compression_level: 5,
            create_index: false,
            create_md5: false,
        }
    }
}

/// The `write` step: a BAM of the records at its position.
struct WriteSink {
    writer: BamWriter,
    header: sam::Header,
    companions: Companions,
}

pub(super) fn write_sink(
    options: &WriteOptions,
    threads: usize,
    frame: &Frame,
) -> Result<Box<dyn Sink>> {
    let path = &options.path;
    let extension = |wanted: &str| {
        path.extension()
            .is_some_and(|e| e.eq_ignore_ascii_case(wanted))
    };
    if extension("sam") || extension("cram") {
        return Err(Error::Usage(format!(
            "write {}: bammetrics writes BAM only (SAM and CRAM output are not supported)",
            path.display()
        )));
    }
    if let Some(read) = frame
        .reads
        .iter()
        .find(|read| mark_duplicates::same_file(read, path))
    {
        return Err(Error::Usage(format!(
            "write {}: it is the input {}; writing it would destroy it",
            path.display(),
            read.display()
        )));
    }
    Ok(Box::new(WriteSink {
        writer: BamWriter::create(path, options.compression_level, threads, &frame.header)?,
        header: frame.header.clone(),
        companions: Companions {
            index: options.create_index,
            md5: options.create_md5,
        },
    }))
}

impl Sink for WriteSink {
    fn accept_batch(&mut self, batch: &[RawRecord]) -> Result<()> {
        batch
            .iter()
            .try_for_each(|record| self.writer.write(record))
    }

    fn finish(self: Box<Self>) -> Result<Products> {
        let files = bam_writer::finish(self.writer, &self.header, self.companions)?;
        Ok(Products {
            files,
            ..Products::default()
        })
    }
}
