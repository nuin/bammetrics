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
//! * **CRAM 3.0.** htsjdk 5.0.0 (Picard 3.5.0) writes 3.1: its files start
//!   `CRAM 03 01`. 3.0 is the version a reader of 3.1 also reads; samtools
//!   1.25, Picard 3.5.0 and this crate read the files. Blocks are gzip at
//!   `COMPRESSION_LEVEL`, the quality scores rANS 4x8 of order 1.
//! * **The header as htsjdk's CRAM writer stores it**: the BAM header text
//!   ([`crate::bam_writer::header_bytes`]), one character per byte as
//!   htsjdk holds a header, encoded as UTF-8 ([`crate::java::latin1_utf8`]).
//!   `café` read from SAM text is held as `caf\xE9` and stored as UTF-8
//!   `caf\xC3\xA9`; read from a BAM it is held as `caf\xC3\xA9` and stored
//!   as `caf\xC3\x83\xC2\xA9`, as in Picard. [`crate::cram::read_header`]
//!   reads it back as UTF-8. htsjdk adds no `M5` and no `UR`
//!   to an `@SQ` line and needs neither to read; a reader finds the bases
//!   through its own reference (`R=`, `samtools view -T`).
//! * **Each reference sequence found in the FASTA by its name**, as
//!   htsjdk's `ReferenceSource` finds it. The FASTA may hold more contigs
//!   than the header, hold them in another order, and lack one no record
//!   lies on. A record on a contig the FASTA lacks fails with htsjdk's
//!   sentence. htsjdk also takes `chr1` of the FASTA for `1` of the header;
//!   this writer does not: samtools 1.25 reads no record of the CRAM Picard
//!   wrote that way (`ref '1' not present`). Where the `@SQ` length and the
//!   FASTA's differ, the bases both hold are stored against the reference
//!   and the bases past the shorter of the two as they are.
//! * **Containers of one slice, cut as htsjdk's `ContainerFactory` cuts
//!   them**: [`READS_PER_SLICE`] records on one reference sequence; at a
//!   change of reference the slice ends when it holds
//!   [`MIN_SINGLE_REFERENCE_SLICE`] records, else it goes on as a
//!   multi-reference slice of at most that many. Unplaced records fill
//!   slices of their own. A slice also ends at [`BYTES_PER_SLICE`] of BAM
//!   records, which htsjdk's does not: 10,000 reads of 100,000 bases are
//!   1.5 GB in one container. A slice on one reference stores the MD5 of
//!   the upper-case reference bases it spans, cut at the contig's end.
//! * **Every record detached**: its own mate reference, mate start, mate
//!   flags and template length, as the BAM holds them. htsjdk, htslib and
//!   noodles each rebuild the fields of an attached mate by a rule of their
//!   own, so only a detached record decodes to the same fields in all three.
//!   A mate reference index the header lacks fails the writer: no reader
//!   opens a CRAM that holds one.
//! * **Every tag where the BAM has it, typed as the BAM types it**, `RG`,
//!   `MD` and `NM` among them. The read group data series is -1: a decoder
//!   that takes `RG` from it prints the tag last. A `Z` or `H` text is
//!   stored as the bytes of the BAM, whatever they are: a UTF-8 letter, a
//!   control byte. htsjdk's CRAM holds an `H` as `B:c` and each byte of a
//!   `Z` over 127 as `?`.
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
//! Memory is bounded by bytes, whatever the read length and the thread
//! count:
//! * The container being filled: at most [`BYTES_PER_SLICE`] (8 MiB) of BAM
//!   records and one record. A mapped record without bases (SEQ `*`) also
//!   counts two bytes per read base of its CIGAR, the `N` and the missing
//!   quality it is stored with, and one whose CIGAR covers more than
//!   16 Mi read bases fails the writer
//!   (`convert::MAX_MISSING_SEQUENCE`): a 194-byte BAM of 4 such records
//!   with CIGAR `268435455M` cost 2,058 MiB before. 40 records of
//!   `16000000S1M` peak at 76 MiB with 1 thread under macOS with its
//!   large-allocation cache off (`MallocLargeCache=0`), and at 590 MiB with
//!   it on: libmalloc keeps the freed blocks of each container.
//! * The containers that wait for their encoders: one per thread, and at
//!   most `QUEUED_BYTES` (64 MiB) of BAM records together. 3,000 reads of
//!   50,000 bases (225 MB of BAM records) peak at 185 MiB with 8 threads,
//!   and at 777 MiB when a container held all of them
//!   (`benches/memgate.toml`, `clean-sam cram long reads`).
//! * The reference bases one slice spans, read through the `.fai`: 10 MiB
//!   of peak memory for 5 reads within 250 bases of a 250 Mbp contig. A
//!   slice that spans more than an eighth of its contig loads the contig
//!   whole, and the writer keeps it until such a slice loads another: 248
//!   MiB for 2 reads at the two ends of that contig, and for a stream that
//!   is not sorted within one contig. A stream whose wide slices alternate
//!   between two contigs loads a contig at each change. The records of a
//!   multi-reference slice read the bases under each record, so a stream
//!   sorted by name holds no contig.
//!
//! No step searches a list by name per record or per contig: 160,000
//! contigs with one read each are written in 0.5 s (19 s when each contig
//! was looked up by its name in the `.fai`).
//!
//! What a CRAM cannot hold, so a decoder gives another record than the BAM:
//! * `=` and `X` operators, read back as `M`; neighbouring operators of one
//!   kind, read back as one. An operator of length 0 is left out, so
//!   `5M0D5M` reads back `10M`. htsjdk leaves out `0I` and `0H` alone; its
//!   `10M0S` is a slice samtools 1.25 fails on.
//! * The CIGAR of a record flagged unmapped, and its mapping quality.
//!   htsjdk's writer fails on such a record (`9 > 0`); htslib writes it as
//!   this writer does.
//! * CIGAR `*` on a mapped record with bases, read back as one `M`.
//! * A record flagged mapped without a reference sequence or a start: it
//!   has no reference to differ from, and the writer fails naming it.
//! * The `bin` field and the width of an integer tag, which a decoder
//!   chooses; the text of the record is the same.
//!
//! A read that starts past the end of its contig is written, every base
//! stored: htsjdk's writer fails on it (`Requested start 100000 is beyond
//! the sequence length chr1`), htslib writes it.
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
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::map::reference_sequence::tag::{
    ALTERNATIVE_NAMES, ASSEMBLY_ID,
};

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
/// A slice also ends once its records hold this many bytes as BAM records
/// (each with its `block_size`), so a container of long reads is no larger
/// than one of short reads. htsjdk has no such limit; htslib ends a slice
/// at 5,000,000 bases, which is 7.5 MB of BAM bases and qualities.
pub const BYTES_PER_SLICE: usize = 8 << 20;
/// The containers that wait for their encoders hold at most this many
/// bytes of BAM records together, whatever the number of threads.
const QUEUED_BYTES: usize = 64 << 20;
/// A slice on one reference reads the bases it spans through the `.fai`.
/// One that spans more than one part in this many of its contig loads the
/// contig whole, and the writer keeps it for the slices after it.
const WHOLE_CONTIG_PART: usize = 8;

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

/// What [`SliceContext`] is made from.
#[derive(Clone, Copy)]
enum Context {
    Single {
        id: usize,
        start: usize,
        end: usize,
        /// Of the upper-case reference bases `start` to `end`.
        md5: [u8; 16],
    },
    Unplaced,
    Many,
}

/// One reference sequence of the header.
struct Contig {
    name: Vec<u8>,
    /// The bases a record can lie on: the `@SQ` length, or the FASTA's when
    /// that is shorter.
    length: usize,
    /// Its line in the `.fai`, found by name. `None` without a FASTA and
    /// for a name the FASTA lacks: a record on it fails.
    fasta: Option<usize>,
    /// For a name the FASTA lacks, the error of a record on it.
    absent: Option<String>,
}

/// htsjdk's message for a record on the reference sequence `name` at
/// `index` when the FASTA has no contig of that name (`CRAMReferenceRegion`,
/// with `SAMSequenceRecord.toString`).
fn absent_from_fasta(index: usize, name: &[u8], map: &Map<ReferenceSequence>) -> String {
    let field = |tag| map.other_fields().get(&tag).map(|value| value.to_string());
    format!(
        "A reference must be supplied (reference sequence SAMSequenceRecord(name={},length={},\
         dict_index={index},assembly={},alternate_names=[{}]) not found).",
        String::from_utf8_lossy(name),
        usize::from(map.length()),
        field(ASSEMBLY_ID).unwrap_or_else(|| "null".into()),
        field(ALTERNATIVE_NAMES).map_or_else(String::new, |names| names.replace(',', ", ")),
    )
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
    contigs: Vec<Contig>,
    /// The records of the container being filled, as a BAM holds them:
    /// each after its `block_size`. At most [`BYTES_PER_SLICE`] and one
    /// record.
    pending: Vec<u8>,
    /// The records in `pending`.
    filled: usize,
    /// The bytes of `pending` and the bytes its records hold once converted
    /// beyond them (`N`s and missing qualities of a mapped record without
    /// bases): what [`BYTES_PER_SLICE`] and `QUEUED_BYTES` count.
    weight: usize,
    slice: Option<Slice>,
    /// The first and the last reference base of the records in `pending`,
    /// kept while they lie on one reference.
    span: (usize, usize),
    /// Containers cut and not yet encoded, at most `threads` and
    /// [`QUEUED_BYTES`].
    cut: Vec<Cut>,
    /// The `weight` of the containers in `cut`.
    queued: usize,
    /// A contig in memory, upper case and cut at its `@SQ` length: the one
    /// of the last slice that spanned more than a part of its contig.
    loaded: Option<(usize, Vec<u8>)>,
    /// The record the bytes of `pending` are read back into.
    scratch: RawRecord,
    /// The records in the containers cut so far.
    counter: u64,
}

impl CramWriter {
    /// Start `path` under its temporary name and write `header`.
    /// `reference` is the FASTA the records are stored against, with its
    /// `.fai`. Each reference sequence of the header is found in it by
    /// name, as htsjdk's `ReferenceSource` finds it: the FASTA may hold
    /// more contigs, in another order, and may lack one no record lies on.
    /// Without a FASTA the file holds unplaced records alone, as htsjdk's:
    /// the first record on a reference sequence fails, and so does the
    /// first on one the FASTA lacks.
    pub fn create(
        path: &Path,
        reference: Option<&Path>,
        header: &sam::Header,
        options: &CramOptions,
    ) -> Result<Self> {
        let reference = reference.map(Reference::open_unchecked).transpose()?;
        let places = reference.as_ref().map(|fasta| fasta.places(header));
        let contigs = header
            .reference_sequences()
            .iter()
            .enumerate()
            .map(|(id, (name, map))| {
                let declared = usize::from(map.length());
                let place = places.as_ref().and_then(|places| places[id]);
                Contig {
                    name: name.to_vec(),
                    length: place.map_or(declared, |(_, length)| {
                        declared.min(usize::try_from(length).unwrap_or(usize::MAX))
                    }),
                    fasta: place.map(|(line, _)| line),
                    absent: (places.is_some() && place.is_none())
                        .then(|| absent_from_fasta(id, name, map)),
                }
            })
            .collect();
        let path = path.to_path_buf();
        let error = |source| Error::Output {
            path: path.clone(),
            source,
        };
        let file = AtomicFile::new(&path);
        let mut out = BufWriter::with_capacity(1 << 20, file.create().map_err(error)?);
        let encoding = raw::Options::gzip(u32::from(options.compression_level));
        raw::write_file_definition(&mut out, &encoding).map_err(error)?;
        let text = bam_writer::header_bytes(header);
        raw::write_header_text(&mut out, &crate::java::latin1_utf8(&text)).map_err(error)?;
        Ok(Self {
            path,
            file,
            out,
            order: Presorted::of(header),
            options: encoding,
            threads: options.threads.max(1),
            reference,
            contigs,
            pending: Vec::new(),
            filled: 0,
            weight: 0,
            slice: None,
            span: (0, 0),
            cut: Vec::new(),
            queued: 0,
            loaded: None,
            scratch: RawRecord::default(),
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
        let outside = |what: &str, id: usize| {
            Error::Record(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "record {}: {what} {id} is not in the header",
                    String::from_utf8_lossy(record.name())
                ),
            ))
        };
        // A reader of the file fails on a mate reference the header lacks
        // (samtools: `Requested unknown mate reference ID`).
        if let Ok(mate) = usize::try_from(record.mate_reference_id())
            && mate >= self.contigs.len()
        {
            return Err(outside("mate reference index", mate));
        }
        let Ok(id) = usize::try_from(record.reference_id()) else {
            return Ok(Slice::Unplaced);
        };
        let Some(contig) = self.contigs.get(id) else {
            return Err(outside("reference index", id));
        };
        if contig.fasta.is_none() {
            return Err(self.no_reference(id));
        }
        let start = usize::try_from(record.alignment_start()).unwrap_or(0);
        Ok(if (1..=contig.length).contains(&start) {
            Slice::Single(id)
        } else {
            Slice::Many
        })
    }

    /// htsjdk's message for a record on a reference sequence it has no
    /// bases for, from `CRAMReferenceRegion`: one sentence without a FASTA,
    /// another for a FASTA that lacks the name.
    fn no_reference(&self, id: usize) -> Error {
        let contig = self.contigs.get(id);
        let absent = contig.and_then(|contig| contig.absent.clone());
        Error::Metric(absent.unwrap_or_else(|| {
            format!(
                "A reference must be supplied that includes the reference sequence for {}.",
                contig.map_or_else(
                    || id.to_string(),
                    |contig| String::from_utf8_lossy(&contig.name).into_owned()
                )
            )
        }))
    }

    /// Add `record`, which must follow the one before it in the header's
    /// order.
    pub fn write(&mut self, record: &RawRecord) -> Result<()> {
        self.order.check(record, &self.path)?;
        let class = self.class(record)?;
        let missing = convert::missing_sequence_bytes(record).map_err(Error::Record)?;
        // htsjdk `ContainerFactory.shouldEmitSlice`, and the byte limit.
        let full = match (self.slice, class) {
            (None, _) => false,
            _ if self.weight >= BYTES_PER_SLICE => true,
            (Some(Slice::Single(_) | Slice::Unplaced), next) if self.slice == Some(next) => {
                self.filled >= READS_PER_SLICE
            }
            _ => self.filled >= MIN_SINGLE_REFERENCE_SLICE,
        };
        if full {
            self.cut_container()?;
        }
        if matches!(class, Slice::Single(_)) {
            let (first, last) = convert::extent(record).map_err(Error::Record)?;
            self.span = match self.slice {
                None => (first, last),
                Some(_) => (self.span.0.min(first), self.span.1.max(last)),
            };
        }
        self.slice = Some(match self.slice {
            None => class,
            Some(current) if current == class => current,
            Some(_) => Slice::Many,
        });
        let before = self.pending.len();
        record.write_to(&mut self.pending).map_err(Error::Record)?;
        self.weight += self.pending.len() - before + missing;
        self.filled += 1;
        Ok(())
    }

    /// The FASTA and the line of the contig at `id` in its `.fai`.
    /// [`CramWriter::class`] lets no record on a contig without one in.
    fn fasta(&mut self, id: usize) -> Result<(&mut Reference, usize)> {
        let line = self.contigs.get(id).and_then(|contig| contig.fasta);
        let Some(line) = line.filter(|_| self.reference.is_some()) else {
            return Err(self.no_reference(id));
        };
        match self.reference.as_mut() {
            Some(reference) => Ok((reference, line)),
            None => Err(Error::Metric("the CRAM writer has no reference".into())),
        }
    }

    /// The contig at `id` whole, upper case and cut at its `@SQ` length, in
    /// place of the one held.
    fn load(&mut self, id: usize) -> Result<()> {
        if self.loaded.as_ref().is_some_and(|(held, _)| *held == id) {
            return Ok(());
        }
        self.loaded = None;
        let length = self.contigs[id].length;
        let (reference, line) = self.fasta(id)?;
        let mut bases = reference.load(line)?;
        bases.truncate(length);
        bases.make_ascii_uppercase();
        self.loaded = Some((id, bases));
        Ok(())
    }

    /// The context of a slice on the contig at `id`, and the upper-case
    /// bases it spans: read through the `.fai`, or `None` when the contig
    /// is in memory ([`WHOLE_CONTIG_PART`]).
    fn single(&mut self, id: usize) -> Result<(Context, Option<Vec<u8>>)> {
        let length = self.contigs[id].length;
        let (start, end) = (self.span.0, self.span.1.min(length));
        let held = self.loaded.as_ref().is_some_and(|(held, _)| *held == id);
        let mut md5 = Md5::new();
        let fetched = if held || (end - start + 1).saturating_mul(WHOLE_CONTIG_PART) > length {
            self.load(id)?;
            let bases = self
                .loaded
                .as_ref()
                .and_then(|(_, bases)| bases.get(start - 1..end));
            let Some(bases) = bases else {
                return Err(self.short(id, end));
            };
            md5.update(bases);
            None
        } else {
            let (reference, line) = self.fasta(id)?;
            let mut bases = reference.subsequence(line, start as i64, end as i64)?;
            if bases.len() != end - start + 1 {
                return Err(self.short(id, end));
            }
            bases.make_ascii_uppercase();
            md5.update(&bases);
            Some(bases)
        };
        let context = Context::Single {
            id,
            start,
            end,
            md5: md5.finish(),
        };
        Ok((context, fetched))
    }

    /// The error of a FASTA whose lines of the contig at `id` end before
    /// the base `end`, which its `.fai` line holds.
    fn short(&self, id: usize, end: usize) -> Error {
        let message = format!(
            "contig {} has fewer than the {end} bases its .fai line gives it",
            String::from_utf8_lossy(&self.contigs[id].name)
        );
        match &self.reference {
            Some(reference) => reference.error(message),
            None => Error::Metric(message),
        }
    }

    /// Convert the records of the container being filled and queue it.
    fn cut_container(&mut self) -> Result<()> {
        let Some(slice) = self.slice.take() else {
            return Ok(());
        };
        let (context, spanned) = match slice {
            Slice::Single(id) => self.single(id)?,
            Slice::Unplaced => (Context::Unplaced, None),
            Slice::Many => (Context::Many, None),
        };
        let pending = std::mem::take(&mut self.pending);
        let mut record = std::mem::take(&mut self.scratch);
        let mut reader = pending.as_slice();
        let mut records = Vec::with_capacity(self.filled);
        let mut deltas = matches!(slice, Slice::Single(_));
        let mut previous = 0;
        let mut fetched;
        while record.read_from(&mut reader).map_err(Error::Record)? {
            let start = usize::try_from(record.alignment_start()).unwrap_or(0);
            let id = usize::try_from(record.reference_id()).ok();
            let window = match (context, &spanned, id, &self.loaded) {
                // A slice on one reference: the bases it spans.
                (Context::Single { start: first, .. }, Some(bases), _, _) => {
                    Window { bases, first }
                }
                (_, _, Some(id), Some((held, bases))) if id == *held => Window { bases, first: 1 },
                (_, _, Some(id), _) if record.flags() & 0x4 == 0 && start > 0 => {
                    // A multi-reference slice: the bases under this record.
                    let (_, last) = convert::extent(&record).map_err(Error::Record)?;
                    let end = last.min(self.contigs.get(id).map_or(0, |contig| contig.length));
                    let (reference, line) = self.fasta(id)?;
                    fetched = reference.subsequence(line, start as i64, end as i64)?;
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
            let (converted, _) = convert::convert(&record, window).map_err(Error::Record)?;
            deltas &= start >= previous;
            previous = start;
            records.push(converted);
        }
        self.cut.push(Cut {
            counter: self.counter,
            context,
            deltas,
            records,
        });
        self.counter += self.filled as u64;
        self.filled = 0;
        self.queued += std::mem::take(&mut self.weight);
        self.pending = pending;
        self.pending.clear();
        self.scratch = record;
        if self.cut.len() >= self.threads || self.queued >= QUEUED_BYTES {
            self.encode()?;
        }
        Ok(())
    }

    /// Encode the containers that wait, one thread each, and write them in
    /// their order.
    fn encode(&mut self) -> Result<()> {
        let cuts = std::mem::take(&mut self.cut);
        self.queued = 0;
        if cuts.is_empty() {
            return Ok(());
        }
        let options = &self.options;
        let encoded: Vec<io::Result<Vec<u8>>> = if cuts.len() == 1 {
            vec![encode_container(options, &cuts[0])]
        } else {
            std::thread::scope(|scope| {
                let workers: Vec<_> = cuts
                    .iter()
                    .map(|cut| scope.spawn(move || encode_container(options, cut)))
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

/// The bytes of the container of `cut`.
fn encode_container(options: &raw::Options, cut: &Cut) -> io::Result<Vec<u8>> {
    let invalid = |what: &str| io::Error::new(io::ErrorKind::InvalidData, what.to_owned());
    let context = match cut.context {
        Context::Unplaced => SliceContext::Unplaced,
        Context::Many => SliceContext::Many,
        Context::Single {
            id,
            start,
            end,
            md5,
        } => SliceContext::Single {
            reference_sequence_id: id,
            alignment_start: Position::new(start)
                .ok_or_else(|| invalid("CRAM: a slice without a start"))?,
            alignment_end: Position::new(end)
                .ok_or_else(|| invalid("CRAM: a slice without an end"))?,
            reference_md5: md5,
        },
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
