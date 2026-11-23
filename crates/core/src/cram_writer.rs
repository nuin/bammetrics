//! The CRAM file a tool or a chain writes, beside [`crate::bam_writer`] and
//! [`crate::sam_text::SamWriter`]: a header and [`RawRecord`]s into CRAM 3.0
//! containers, under a temporary name ([`AtomicFile`]) until [`finish`]
//! renames it into place with its MD5.
//!
//! The truth is the decoded file, never its bytes (docs/concordance.md,
//! tier B): htsjdk, samtools and bammetrics decode the CRAM to the header
//! and the records a BAM of the same stream holds.
//!
//! What is written, and why:
//! * **CRAM 3.0.** htsjdk 5.0.0 (Picard 3.5.0) writes 3.1. Every reader of
//!   3.1 reads 3.0, and htsjdk before 4.2 and htslib before 1.12 read 3.0
//!   alone. Blocks are gzip at `COMPRESSION_LEVEL`, the quality scores rANS
//!   4x8 of order 1.
//! * **The header as the BAM holds it**, byte for byte
//!   ([`crate::bam_writer::header_bytes`]). htsjdk adds no `M5` and no `UR`
//!   to an `@SQ` line and needs neither to read; a reader finds the bases
//!   through its own reference (`R=`, `samtools view -T`).
//! * **Containers of one slice, cut as htsjdk's `ContainerFactory` cuts
//!   them**: [`READS_PER_SLICE`] records on one reference sequence; at a
//!   change of reference the slice ends when it holds
//!   [`MIN_SINGLE_REFERENCE_SLICE`] records, else it goes on as a
//!   multi-reference slice of at most that many. Unplaced records fill
//!   slices of their own. A slice on one reference stores the MD5 of the
//!   upper-case reference bases it spans, cut at the contig's end.
//! * **Every record detached**: its own mate reference, mate start, mate
//!   flags and template length, as the BAM holds them. htsjdk, htslib and
//!   noodles each rebuild the fields of an attached mate by a rule of their
//!   own, so only a detached record decodes to the same fields in all three.
//! * **Every tag where the BAM has it, typed as the BAM types it**, `RG`,
//!   `MD` and `NM` among them. The read group data series is -1: a decoder
//!   that takes `RG` from it prints the tag last.
//! * **Bases as features against the reference**, upper case. A mismatch
//!   with a base outside `ACGTN` on either side is stored as the base, and
//!   so is every read base past the end of its contig, whatever it is.
//!   htsjdk leaves a read `N` there to the decoder, which reads the missing
//!   reference base as `N`. From such a CRAM samtools 1.25 warns `Ref pos
//!   outside of ref sequence boundary` and rebuilds another `MD` than from
//!   its own CRAM of the same record, and Picard 3.5.0 rebuilds an empty
//!   `MD` (seen on a 20-base read with 14 bases past the end, 4 of them
//!   `N`). With the base stored, samtools reads 28 shapes of such reads as
//!   it reads its own CRAM of them, without the warning
//!   (`tests/concordance/test_cram_writer.py`).
//! * **The qualities as an array**, `0xFF` bytes for a record without
//!   them, which samtools, Picard and this crate's reader each print `*`.
//!   A record without bases stores its CIGAR alone, an insertion or a clip
//!   as `N`s.
//!
//! Memory: the records of one container per compression thread, and one
//! contig of the reference. A slice on one reference loads that contig
//! whole, as htsjdk does; the records of a multi-reference slice read the
//! bases under each record through the `.fai`, so an unsorted stream never
//! holds a contig it does not stay on.
//!
//! What a CRAM cannot hold, so a decoder gives another record than the BAM:
//! * `=` and `X` operators, read back as `M`; neighbouring operators of one
//!   kind, read back as one; an operator of length 0.
//! * The CIGAR of a record flagged unmapped, and its mapping quality.
//! * CIGAR `*` on a mapped record with bases, read back as one `M`.
//! * A record flagged mapped without a reference sequence or a start: it
//!   has no reference to differ from, and the writer fails naming it.
//! * The `bin` field and the width of an integer tag, which a decoder
//!   chooses; the text of the record is the same.
//!
//! What a CRAM holds and a reader does not give back, from any writer's
//! file (`tests/concordance/test_cram_writer.py` holds each to htslib's own
//! CRAM of the same BAM):
//! * The mate reference and mate start of a record without the paired
//!   flag: htslib 1.25 reads the reference as `*`, htsjdk 5.0.0 both as
//!   `*` and 0. This crate's reader returns them.
//! * The unsigned subtype of a `B` array: htsjdk 5.0.0 and this crate's
//!   reader, which follows it, read `B:C`, `B:S` and `B:I` as `B:c`, `B:s`
//!   and `B:i`. htslib returns the subtype stored.
//! * A deletion that runs past the end of its contig: `samtools view`
//!   fails on the slice (`CRAM CIGAR extends beyond slice reference
//!   extents`) unless it is told `decode_md=0`.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use noodles::core::Position;
use noodles::cram::io::writer::raw::{self, Record, SliceContext};
use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::bam_writer::{self, Companions, Presorted};
use crate::cleanup::AtomicFile;
use crate::error::{Error, Result};
use crate::md5::Md5;
use crate::reference::Reference;

mod convert;

use convert::Window;

/// htsjdk `CRAMEncodingStrategy.DEFAULT_READS_PER_SLICE`; a container holds
/// one slice.
pub const READS_PER_SLICE: usize = 10_000;
/// htsjdk `DEFAULT_MINIMUM_SINGLE_REFERENCE_SLICE_THRESHOLD`: a slice of
/// fewer records does not end at a change of reference, and a
/// multi-reference slice ends at this many.
pub const MIN_SINGLE_REFERENCE_SLICE: usize = 1_000;

/// How a CRAM is written.
#[derive(Clone, Copy, Debug)]
pub struct CramOptions {
    /// The gzip level of the blocks, 0 to 9.
    pub compression_level: u8,
    /// The threads that encode containers.
    pub threads: usize,
}

impl Default for CramOptions {
    fn default() -> Self {
        Self {
            compression_level: 5,
            threads: 1,
        }
    }
}

/// What the records of the container being filled lie on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Slice {
    Single(usize),
    Unplaced,
    Many,
}

/// A container cut and converted, waiting for its encoder.
struct Cut {
    /// The records in the containers before it.
    counter: u64,
    context: Context,
    /// The alignment starts never go down, so each is stored as its
    /// difference from the one before.
    deltas: bool,
    records: Vec<Record>,
}

/// [`SliceContext`] before the MD5 of a single reference is known.
#[derive(Clone, Copy)]
enum Context {
    Single { id: usize, start: usize, end: usize },
    Unplaced,
    Many,
}

/// The CRAM output. Records must come in the header's order, as the BAM
/// writer asks ([`Presorted`]).
pub struct CramWriter {
    path: PathBuf,
    file: AtomicFile,
    out: BufWriter<File>,
    order: Presorted,
    options: raw::Options,
    threads: usize,
    reference: Option<Reference>,
    /// The name and the length of each reference sequence of the header.
    contigs: Vec<(Vec<u8>, usize)>,
    /// The records of the container being filled: the first `filled` are
    /// its records, the rest are buffers kept for the next ones.
    pending: Vec<RawRecord>,
    filled: usize,
    slice: Option<Slice>,
    /// Containers cut and not yet encoded, at most `threads`.
    cut: Vec<Cut>,
    /// The contig in memory, upper case: the one of the last slice on a
    /// single reference.
    loaded: Option<(usize, Vec<u8>)>,
    /// The records in the containers cut so far.
    counter: u64,
}

impl CramWriter {
    /// Start `path` under its temporary name and write `header`.
    /// `reference` is the FASTA the records are stored against, with its
    /// `.fai`; its contigs must be the header's. Without one the file holds
    /// unplaced records alone, as htsjdk's: the first record with a
    /// reference sequence fails.
    pub fn create(
        path: &Path,
        reference: Option<&Path>,
        header: &sam::Header,
        options: &CramOptions,
    ) -> Result<Self> {
        let reference = reference
            .map(|fasta| Reference::open(fasta, header))
            .transpose()?;
        let path = path.to_path_buf();
        let error = |source| Error::Output {
            path: path.clone(),
            source,
        };
        let file = AtomicFile::new(&path);
        let mut out = BufWriter::with_capacity(1 << 20, file.create().map_err(error)?);
        let encoding = raw::Options::gzip(u32::from(options.compression_level));
        raw::write_file_definition(&mut out, &encoding).map_err(error)?;
        raw::write_header_text(&mut out, &bam_writer::header_bytes(header)).map_err(error)?;
        Ok(Self {
            path,
            file,
            out,
            order: Presorted::of(header),
            options: encoding,
            threads: options.threads.max(1),
            reference,
            contigs: header
                .reference_sequences()
                .iter()
                .map(|(name, map)| (name.to_vec(), usize::from(map.length())))
                .collect(),
            pending: Vec::new(),
            filled: 0,
            slice: None,
            cut: Vec::new(),
            loaded: None,
            counter: 0,
        })
    }

    fn output(&self, source: io::Error) -> Error {
        Error::Output {
            path: self.path.clone(),
            source,
        }
    }

    /// What `record` asks of the slice it goes into: its reference
    /// sequence, none, or a multi-reference slice for a record with a
    /// reference sequence and no start on it (start 0, or past its end),
    /// which a slice on one reference cannot span.
    fn class(&self, record: &RawRecord) -> Result<Slice> {
        let Ok(id) = usize::try_from(record.reference_id()) else {
            return Ok(Slice::Unplaced);
        };
        let Some((name, length)) = self.contigs.get(id) else {
            return Err(Error::Record(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "record {}: reference index {id} is not in the header",
                    String::from_utf8_lossy(record.name())
                ),
            )));
        };
        if self.reference.is_none() {
            // htsjdk's message, from `CRAMReferenceRegion`.
            return Err(Error::Metric(format!(
                "A reference must be supplied that includes the reference sequence for {}.",
                String::from_utf8_lossy(name)
            )));
        }
        let start = usize::try_from(record.alignment_start()).unwrap_or(0);
        Ok(if (1..=*length).contains(&start) {
            Slice::Single(id)
        } else {
            Slice::Many
        })
    }

    /// Add `record`, which must follow the one before it in the header's
    /// order.
    pub fn write(&mut self, record: &RawRecord) -> Result<()> {
        self.order.check(record, &self.path)?;
        let class = self.class(record)?;
        // htsjdk `ContainerFactory.shouldEmitSlice`.
        let full = match (self.slice, class) {
            (None, _) => false,
            (Some(Slice::Single(_) | Slice::Unplaced), next) if self.slice == Some(next) => {
                self.filled >= READS_PER_SLICE
            }
            _ => self.filled >= MIN_SINGLE_REFERENCE_SLICE,
        };
        if full {
            self.cut_container()?;
        }
        self.slice = Some(match self.slice {
            None => class,
            Some(current) if current == class => current,
            Some(_) => Slice::Many,
        });
        match self.pending.get_mut(self.filled) {
            Some(buffer) => buffer.clone_from(record),
            None => self.pending.push(record.clone()),
        }
        self.filled += 1;
        Ok(())
    }

    /// The contig at `id` whole and upper case, in place of the one held.
    /// The containers that wait for their encoder are encoded first: they
    /// take the MD5 of their slice from the contig held.
    fn load(&mut self, id: usize) -> Result<()> {
        if self.loaded.as_ref().is_some_and(|(held, _)| *held == id) {
            return Ok(());
        }
        self.encode()?;
        self.loaded = None;
        let reference = self.reference.as_mut().expect("a placed record has one");
        let mut bases = reference.load(id)?;
        bases.make_ascii_uppercase();
        self.loaded = Some((id, bases));
        Ok(())
    }

    /// Convert the records of the container being filled and queue it.
    fn cut_container(&mut self) -> Result<()> {
        let Some(slice) = self.slice.take() else {
            return Ok(());
        };
        if let Slice::Single(id) = slice {
            self.load(id)?;
        }
        let mut records = Vec::with_capacity(self.filled);
        let mut deltas = matches!(slice, Slice::Single(_));
        let (mut first, mut last, mut previous) = (usize::MAX, 0, 0);
        let mut fetched;
        for record in &self.pending[..self.filled] {
            let start = usize::try_from(record.alignment_start()).unwrap_or(0);
            let window = match (usize::try_from(record.reference_id()), &self.loaded) {
                (Ok(id), Some((held, bases))) if id == *held => Window { bases, first: 1 },
                (Ok(id), _) if record.flags() & 0x4 == 0 && start > 0 => {
                    // A multi-reference slice: the bases under this record.
                    let length = self.contigs.get(id).map_or(0, |(_, length)| *length);
                    let span = crate::bam_raw::reference_length(record.cigar());
                    let end = (start + usize::try_from(span).unwrap_or(0)).min(length);
                    let reference = self.reference.as_mut().expect("a placed record has one");
                    fetched = reference.subsequence(id, start as i64, end as i64)?;
                    fetched.make_ascii_uppercase();
                    Window {
                        bases: &fetched,
                        first: start,
                    }
                }
                _ => Window {
                    bases: &[],
                    first: 1,
                },
            };
            let (converted, end) = convert::convert(record, window).map_err(Error::Record)?;
            deltas &= start >= previous;
            previous = start;
            first = first.min(start);
            last = last.max(end);
            records.push(converted);
        }
        let context = match slice {
            Slice::Single(id) => Context::Single {
                id,
                start: first,
                end: last.min(self.contigs[id].1),
            },
            Slice::Unplaced => Context::Unplaced,
            Slice::Many => Context::Many,
        };
        self.cut.push(Cut {
            counter: self.counter,
            context,
            deltas,
            records,
        });
        self.counter += self.filled as u64;
        self.filled = 0;
        if self.cut.len() >= self.threads {
            self.encode()?;
        }
        Ok(())
    }

    /// Encode the containers that wait, one thread each, and write them in
    /// their order.
    fn encode(&mut self) -> Result<()> {
        let cuts = std::mem::take(&mut self.cut);
        if cuts.is_empty() {
            return Ok(());
        }
        let options = &self.options;
        let contig = self.loaded.as_ref().map(|(_, bases)| bases.as_slice());
        let encoded: Vec<io::Result<Vec<u8>>> = if cuts.len() == 1 {
            vec![encode_container(options, contig, &cuts[0])]
        } else {
            std::thread::scope(|scope| {
                let workers: Vec<_> = cuts
                    .iter()
                    .map(|cut| scope.spawn(move || encode_container(options, contig, cut)))
                    .collect();
                workers
                    .into_iter()
                    .map(|worker| {
                        worker
                            .join()
                            .unwrap_or_else(|_| Err(io::Error::other("the CRAM encoder failed")))
                    })
                    .collect()
            })
        };
        for container in encoded {
            let bytes = container.map_err(|e| self.output(e))?;
            self.out.write_all(&bytes).map_err(|e| self.output(e))?;
        }
        Ok(())
    }

    /// Write the last containers and the end-of-file container. The file
    /// keeps its temporary name until [`finish`] or [`AtomicFile::commit`].
    pub fn close(mut self) -> Result<AtomicFile> {
        self.cut_container()?;
        self.encode()?;
        raw::write_eof(&mut self.out)
            .and_then(|()| self.out.flush())
            .map_err(|e| self.output(e))?;
        Ok(self.file)
    }
}

/// The bytes of the container of `cut`. `contig` is the contig in memory,
/// upper case: the one a slice on a single reference lies on.
fn encode_container(
    options: &raw::Options,
    contig: Option<&[u8]>,
    cut: &Cut,
) -> io::Result<Vec<u8>> {
    let invalid = |what: &str| io::Error::new(io::ErrorKind::InvalidData, what.to_owned());
    let context = match cut.context {
        Context::Unplaced => SliceContext::Unplaced,
        Context::Many => SliceContext::Many,
        Context::Single { id, start, end } => {
            let bases = contig
                .and_then(|bases| bases.get(start - 1..end))
                .ok_or_else(|| invalid("CRAM: a slice outside its reference sequence"))?;
            let mut md5 = Md5::new();
            md5.update(bases);
            SliceContext::Single {
                reference_sequence_id: id,
                alignment_start: Position::new(start)
                    .ok_or_else(|| invalid("CRAM: a slice without a start"))?,
                alignment_end: Position::new(end)
                    .ok_or_else(|| invalid("CRAM: a slice without an end"))?,
                reference_md5: md5.finish(),
            }
        }
    };
    raw::encode_container(options, cut.counter, context, cut.deltas, &cut.records)
}

/// Close `writer`, write `<CRAM>.md5` when `md5`, then rename the files
/// into place. Returns the final paths, the CRAM first.
pub fn finish(writer: CramWriter, md5: bool) -> Result<Vec<PathBuf>> {
    let threads = writer.threads;
    let output = writer.close()?;
    // No index: the header is not consulted.
    let companions = Companions { index: false, md5 };
    bam_writer::finish_file(
        output,
        &sam::Header::default(),
        companions,
        Path::new(""),
        threads,
    )
}

/// The usage error for `CREATE_INDEX` with a CRAM `output`, named as the
/// caller names it (`OUTPUT=out.cram`, `out.cram`). Picard 3.5.0 writes
/// `<CRAM>.bai`, a BAI over the slices that htsjdk alone reads; bammetrics
/// writes no index of a CRAM.
pub fn index_refusal(output: &str) -> Error {
    Error::Usage(format!(
        "{output}: no index is written beside a CRAM (Picard writes a .cram.bai that only \
         htsjdk reads): write it without one and run `samtools index` for a .crai"
    ))
}

#[cfg(test)]
mod tests;
