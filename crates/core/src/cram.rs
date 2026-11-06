//! CRAM input (docs/04 section 5.9: BAM, SAM and CRAM with reference).
//!
//! A file whose first four bytes are `CRAM` is read as CRAM whatever its
//! name; every other file goes to the BAM reader. The records of a CRAM
//! reach the tools as `bam::Record`s: each container is decoded with the
//! `REFERENCE_SEQUENCE` FASTA on a worker thread and its records are
//! encoded as BAM records into one buffer. The buffers come back in file
//! order, so a tool sees the records of a CRAM in the order the BAM of the
//! same data gives them, and the `Accumulator` and `Collector` traits stay
//! as they are.
//!
//! Reference bases are read as htsjdk's `ReferenceSource` hands them to the
//! CRAM decoder: whole contigs, upper case. A container is decoded only
//! with the contig it names, and a contig is dropped once no container in
//! flight uses it, so at most two contigs are held at a time (a container
//! over several contigs holds every contig it touched). Errors come where
//! htsjdk raises them: a CRAM read without a reference fails at the first
//! slice that needs reference bases, and a slice whose reference MD5 does
//! not match the bases fails there.

use std::fs::File;
use std::io::{self, Chain, Cursor, Read, Seek, SeekFrom};
use std::num::NonZero;
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::sync::{Arc, Mutex};
use std::thread;

use noodles::core::Region;
use noodles::cram::io::reader::Container;
use noodles::sam::alignment::io::Write as _;
use noodles::{bam, cram, fasta, sam};

use crate::error::{Error, Result};
use crate::pass;

/// The CRAM file definition starts with these bytes.
pub const MAGIC: [u8; 4] = *b"CRAM";

/// A reader whose first bytes were looked at and are put back in front.
pub(crate) struct Sniffed<R> {
    pub(crate) cram: bool,
    pub(crate) reader: Chain<Cursor<Vec<u8>>, R>,
}

/// Read the first four bytes of `inner` (fewer for a shorter file) and tell
/// a CRAM by its magic. Nothing is lost: the bytes read come first again.
pub(crate) fn sniff<R: Read>(mut inner: R) -> io::Result<Sniffed<R>> {
    let mut head = Vec::with_capacity(MAGIC.len());
    (&mut inner)
        .take(MAGIC.len() as u64)
        .read_to_end(&mut head)?;
    Ok(Sniffed {
        cram: head == MAGIC,
        reader: Cursor::new(head).chain(inner),
    })
}

/// True when `path` is a regular file that starts with the CRAM magic. A
/// pipe or a device is never opened for this, since the bytes looked at
/// would be gone for the read that follows.
pub fn is_cram(path: &Path) -> bool {
    let Ok(file) = File::open(path) else {
        return false;
    };
    file.metadata().is_ok_and(|m| m.is_file()) && sniff(file).is_ok_and(|s| s.cram)
}

/// A FASTA with its `.fai`. A plain file gives a contig with one
/// `read_exact` and its line ends dropped in place, so a contig costs one
/// buffer of its size; a bgzipped one (`.gz`, `.bgz`) goes through the
/// noodles indexed reader.
enum Fasta {
    Plain {
        file: File,
        index: fasta::fai::Index,
    },
    Compressed(fasta::io::IndexedReader<fasta::io::BufReader<File>>),
}

impl Fasta {
    fn open(path: &Path) -> io::Result<Self> {
        if matches!(
            path.extension().and_then(|e| e.to_str()),
            Some("gz" | "bgz")
        ) {
            return fasta::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .map(Self::Compressed);
        }
        let mut index = path.as_os_str().to_owned();
        index.push(".fai");
        Ok(Self::Plain {
            index: fasta::fai::fs::read(index)?,
            file: File::open(path)?,
        })
    }

    /// The bases of contig `name` as the file has them, `None` when the
    /// index has no such contig.
    fn bases(&mut self, name: &[u8]) -> io::Result<Option<Vec<u8>>> {
        let (file, index) = match self {
            Self::Plain { file, index } => (file, index),
            Self::Compressed(reader) => {
                if !reader.index().as_ref().iter().any(|r| r.name() == name) {
                    return Ok(None);
                }
                let record = reader.query(&Region::new(name, ..))?;
                return Ok(Some(record.sequence().as_ref().to_vec()));
            }
        };
        let Some(record) = index.as_ref().iter().find(|r| r.name() == name) else {
            return Ok(None);
        };
        let length = record.length();
        let full_lines = length.saturating_sub(1) / record.line_base_count().get();
        let span = if length == 0 {
            0
        } else {
            full_lines * record.line_width().get() + length
                - full_lines * record.line_base_count().get()
        };
        let too_long = |e| io::Error::new(io::ErrorKind::InvalidData, e);
        let mut bases = vec![0; usize::try_from(span).map_err(too_long)?];
        file.seek(SeekFrom::Start(record.position()))?;
        file.read_exact(&mut bases)?;
        bases.retain(|b| !matches!(b, b'\n' | b'\r'));
        if bases.len() as u64 != length {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "the .fai line lengths do not match the FASTA",
            ));
        }
        Ok(Some(bases))
    }
}

/// The reference FASTA the decoder reads contigs from, shared by every
/// [`fasta::Repository`] of one read. `None` without `REFERENCE_SEQUENCE`.
#[derive(Clone)]
struct Bases {
    reader: Option<Arc<Mutex<Fasta>>>,
    path: PathBuf,
}

impl Bases {
    fn open(reference: Option<&Path>) -> Result<Self> {
        let Some(path) = reference else {
            return Ok(Self {
                reader: None,
                path: PathBuf::new(),
            });
        };
        let reader = Fasta::open(path).map_err(|e| Error::Reference {
            path: path.to_path_buf(),
            message: format!("cannot open with its .fai index: {e}"),
        })?;
        Ok(Self {
            reader: Some(Arc::new(Mutex::new(reader))),
            path: path.to_path_buf(),
        })
    }

    /// A repository with an empty cache over the same FASTA.
    fn repository(&self) -> fasta::Repository {
        fasta::Repository::new(self.clone())
    }

    /// Contig `name` upper case, as htsjdk's `ReferenceSource` gives it.
    fn contig(&self, name: &[u8]) -> io::Result<fasta::Record> {
        let shown = String::from_utf8_lossy(name);
        let Some(reader) = &self.reader else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "the records on {shown} need reference bases: a CRAM input needs \
                     REFERENCE_SEQUENCE (R)"
                ),
            ));
        };
        let mut reader = reader
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        let context = |e: io::Error| {
            io::Error::new(
                e.kind(),
                format!("{}: contig {shown}: {e}", self.path.display()),
            )
        };
        let Some(mut bases) = reader.bases(name).map_err(context)? else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: no contig {shown} in the reference",
                    self.path.display()
                ),
            ));
        };
        bases.make_ascii_uppercase();
        Ok(fasta::Record::new(
            fasta::record::Definition::new(name, None),
            fasta::record::Sequence::from(bases),
        ))
    }
}

impl fasta::repository::Adapter for Bases {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        Some(self.contig(name))
    }
}

/// One container's records as a BAM record stream (each record with its
/// `block_size`), or the error that stopped it.
type Decoded = io::Result<Vec<u8>>;

/// The decoded records of every container in file order: what a
/// `bam::io::Reader` reads after the header.
pub(crate) struct Stream {
    order: Receiver<Receiver<Decoded>>,
    current: Cursor<Vec<u8>>,
}

impl Read for Stream {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        loop {
            let n = self.current.read(buf)?;
            if n > 0 || buf.is_empty() {
                return Ok(n);
            }
            let Ok(next) = self.order.recv() else {
                return Ok(0);
            };
            let bytes = next
                .recv()
                .map_err(|_| io::Error::other("a CRAM decoder thread stopped"))??;
            self.current = Cursor::new(bytes);
        }
    }
}

/// A container to decode, the repository holding its contig, and where its
/// records go.
struct Job {
    container: Container,
    repository: fasta::Repository,
    done: SyncSender<Decoded>,
}

/// Read the containers of `reader` (positioned after the header) and hand
/// `consume` their records as a BAM record stream. Containers are decoded
/// on `workers` threads, at most two per worker queued ahead of `consume`.
/// When `consume` returns early (`STOP_AFTER`), reading stops.
pub(crate) fn decode<R, T>(
    mut reader: cram::io::Reader<Peek<R>>,
    header: &sam::Header,
    reference: Option<&Path>,
    workers: NonZero<usize>,
    consume: impl FnOnce(&mut bam::io::Reader<Stream>) -> Result<T>,
) -> Result<T>
where
    R: Read + Send,
{
    let bases = Bases::open(reference)?;
    let ahead = 2 * workers.get();
    let (order_tx, order) = mpsc::sync_channel::<Receiver<Decoded>>(ahead);
    let (job_tx, jobs) = mpsc::sync_channel::<Job>(ahead);
    let jobs = Mutex::new(jobs);
    thread::scope(|scope| {
        scope.spawn(move || dispatch(&mut reader, header, &bases, &order_tx, &job_tx));
        for _ in 0..workers.get() {
            // Returning on a failed spawn drops `order`, which ends the
            // dispatcher and with it the workers already started.
            thread::Builder::new()
                .stack_size(DECODER_STACK)
                .spawn_scoped(scope, || decode_jobs(&jobs, header))
                .map_err(Error::Record)?;
        }
        let mut records = bam::io::Reader::from(Stream {
            order,
            current: Cursor::new(Vec::new()),
        });
        let result = consume(&mut records);
        // Dropping the stream ends the dispatcher at its next send, and with
        // it the job queue, so every thread of the scope returns.
        drop(records);
        result
    })
}

/// A worker: decode queued containers until the queue closes. A panic in
/// the decoder (a reference shorter than a slice) becomes that
/// container's error.
fn decode_jobs(jobs: &Mutex<Receiver<Job>>, header: &sam::Header) {
    loop {
        let job = match jobs.lock() {
            Ok(queue) => queue.recv(),
            Err(_) => return,
        };
        let Ok(job) = job else { return };
        let decoded = panic::catch_unwind(AssertUnwindSafe(|| {
            decode_container(&job.container, header, &job.repository)
        }))
        .unwrap_or_else(|panic| {
            Err(io::Error::other(format!(
                "decoding a CRAM container failed: {}",
                pass::panic_text(panic.as_ref())
            )))
        });
        // The receiver is gone only after `consume` stopped early.
        let _ = job.done.send(decoded);
    }
}

/// The dispatcher: each container read is queued for a worker, its result
/// slot queued for the stream in the same order. A read error ends the
/// stream with that error.
fn dispatch<R: Read>(
    reader: &mut cram::io::Reader<Peek<R>>,
    header: &sam::Header,
    bases: &Bases,
    order: &SyncSender<Receiver<Decoded>>,
    jobs: &SyncSender<Job>,
) {
    let mut contig: Option<(i32, fasta::Repository)> = None;
    let mut several: Option<fasta::Repository> = None;
    loop {
        let mut container = Container::default();
        reader.get_mut().head.clear();
        let read = reader.read_container(&mut container);
        let (done, result) = mpsc::sync_channel(1);
        match read {
            Ok(0) => return,
            Ok(_) => {}
            Err(e) => {
                let _ = done.send(Err(e));
                let _ = order.send(result);
                return;
            }
        }
        if order.send(result).is_err() {
            return;
        }
        let repository = match reader.get_ref().reference_sequence_id() {
            Some(MANY) => several.get_or_insert_with(|| bases.repository()).clone(),
            Some(id) if id >= 0 => match &contig {
                Some((current, repository)) if *current == id => repository.clone(),
                _ => {
                    let repository = bases.repository();
                    preload(&repository, header, id);
                    contig = Some((id, repository.clone()));
                    repository
                }
            },
            // Unmapped records need no bases; the FASTA stays at hand anyway.
            _ => bases.repository(),
        };
        let job = Job {
            container,
            repository,
            done,
        };
        if jobs.send(job).is_err() {
            return;
        }
    }
}

/// Stack of a decoder thread. The CRAM 3.1 codecs recurse deeper than the
/// 2 MiB a spawned thread gets by default.
const DECODER_STACK: usize = 16 << 20;

/// The container header's reference sequence id for several contigs.
const MANY: i32 = -2;

/// The container length (4 bytes) and an ITF8 of at most 5 bytes.
const HEAD: usize = 9;

/// The input, keeping the first bytes read since `head` was cleared: the
/// start of a container header, whose reference sequence id noodles reads
/// but does not expose.
pub(crate) struct Peek<R> {
    inner: R,
    head: Vec<u8>,
}

impl<R> Peek<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            head: Vec::with_capacity(HEAD),
        }
    }

    /// The ITF8 after the `i32` length: the contig id, -1 for unmapped
    /// records, -2 for several contigs (CRAM 3.1 sections 2.3 and 7.1).
    fn reference_sequence_id(&self) -> Option<i32> {
        let bytes = self.head.get(4..)?;
        let b = |i: usize| bytes.get(i).copied().map(u32::from);
        let lead = bytes.first()?.leading_ones();
        let first = b(0)?;
        let value = match lead {
            0 => first,
            1 => (first & 0x7f) << 8 | b(1)?,
            2 => (first & 0x3f) << 16 | b(1)? << 8 | b(2)?,
            3 => (first & 0x1f) << 24 | b(1)? << 16 | b(2)? << 8 | b(3)?,
            _ => (first & 0x0f) << 28 | b(1)? << 20 | b(2)? << 12 | b(3)? << 4 | (b(4)? & 0x0f),
        };
        Some(value.cast_signed())
    }
}

impl<R: Read> Read for Peek<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        let wanted = HEAD.saturating_sub(self.head.len()).min(n);
        self.head.extend_from_slice(&buf[..wanted]);
        Ok(n)
    }
}

/// Load contig `id` once, before the workers of its containers ask for it
/// at the same time. A failure is left for the decoder to report when it
/// needs the bases.
fn preload(repository: &fasta::Repository, header: &sam::Header, id: i32) {
    let found = usize::try_from(id)
        .ok()
        .and_then(|i| header.reference_sequences().get_index(i));
    if let Some((name, _)) = found {
        let _ = repository.get(name);
    }
}

/// Every record of `container` encoded as a BAM record with its
/// `block_size`, in slice order.
fn decode_container(
    container: &Container,
    header: &sam::Header,
    repository: &fasta::Repository,
) -> io::Result<Vec<u8>> {
    let compression_header = container.compression_header()?;
    let mut writer = bam::io::Writer::from(Vec::new());
    for slice in container.slices() {
        let slice = slice?;
        let (core, external) = slice.decode_blocks()?;
        let records = slice
            .records(
                repository.clone(),
                header,
                &compression_header,
                &core,
                &external,
            )
            .map_err(checksum)?;
        for record in &records {
            writer.write_alignment_record(header, record)?;
        }
    }
    Ok(writer.into_inner())
}

/// The slice MD5 check in the words of htsjdk's `CRAMException`.
fn checksum(e: io::Error) -> io::Error {
    if e.to_string().contains("checksum mismatch") {
        io::Error::new(
            e.kind(),
            "the MD5 of the reference bases does not match the CRAM slice: the reference is \
             not the one the CRAM was written with",
        )
    } else {
        e
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sniffing_puts_the_bytes_back() {
        for data in [&b"CRAM\x03\x00rest"[..], b"\x1f\x8b\x08\x04", b"CR", b""] {
            let mut sniffed = sniff(data).expect("read");
            assert_eq!(sniffed.cram, data.starts_with(&MAGIC));
            let mut all = Vec::new();
            sniffed.reader.read_to_end(&mut all).expect("read");
            assert_eq!(all, data);
        }
    }

    #[test]
    fn the_container_contig_id_is_the_itf8_after_the_length() {
        let id = |itf8: &[u8]| {
            let mut peek = Peek::new(io::empty());
            peek.head = [&[9, 0, 0, 0][..], itf8, &[0; 5]].concat();
            peek.head.truncate(HEAD);
            peek.reference_sequence_id()
        };
        assert_eq!(id(&[0]), Some(0));
        assert_eq!(id(&[0x7f]), Some(127));
        assert_eq!(id(&[0x80, 200]), Some(200));
        assert_eq!(id(&[0xc1, 0x02, 0x03]), Some(0x01_0203));
        assert_eq!(id(&[0xff, 0xff, 0xff, 0xff, 0x0f]), Some(-1));
        assert_eq!(id(&[0xff, 0xff, 0xff, 0xff, 0x0e]), Some(MANY));
        assert_eq!(Peek::new(io::empty()).reference_sequence_id(), None);
    }

    #[test]
    fn a_missing_reference_fails_when_bases_are_needed() {
        let bases = Bases::open(None).expect("no reference opens");
        let err = bases.contig(b"chr1").expect_err("no reference");
        assert!(err.to_string().contains("REFERENCE_SEQUENCE (R)"), "{err}");
    }

    #[test]
    fn contigs_come_upper_case_and_unknown_names_fail() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/ref_100k.fa");
        let bases = Bases::open(Some(&path)).expect("opens");
        let chr1 = bases.contig(b"chr1").expect("chr1");
        assert_eq!(chr1.sequence().len(), 100_000);
        assert!(!chr1.sequence().as_ref().iter().any(u8::is_ascii_lowercase));
        let err = bases.contig(b"chrZ").expect_err("no chrZ");
        assert!(err.to_string().contains("no contig chrZ"), "{err}");
    }
}
