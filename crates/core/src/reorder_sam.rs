//! `ReorderSam`, mirroring Picard 3.5.0 (`picard.sam.ReorderSam`): the
//! records of a SAM, BAM or CRAM under the sequence dictionary of
//! `SEQUENCE_DICTIONARY`, matched to the input's by contig name.
//!
//! The dictionary comes from a `.dict`, a FASTA with its `.dict` beside it
//! (Picard fails on a FASTA without one), a SAM, BAM or CRAM header, an
//! interval list or a VCF's `##contig` lines, told by the name's extension
//! as `SAMSequenceDictionaryExtractor` tells them ([`DictionaryKind`]); any
//! other name is refused. A file read whole holds at most
//! [`bed_to_interval_list::MAX_DICTIONARY_BYTES`]. A dictionary with no
//! `@SQ` line is taken, as Picard takes it; a length below 1 or over
//! 2^31 - 1 is refused. Every input contig needs a contig of the
//! same name (`ALLOW_INCOMPLETE_DICT_CONCORDANCE` lets one go) and of the
//! same length (`ALLOW_CONTIG_LENGTH_DISCORDANCE` lets it differ, with
//! Picard's warning on standard error).
//!
//! Each record gets the new index of its reference and its mate's. A
//! record whose contig the new dictionary lacks becomes unplaced: start 0,
//! unmapped flag, CIGAR `*`, MAPQ 0; its other flags and `TLEN` stay. A
//! record whose mate's contig it lacks gets mate start 0, the mate
//! unmapped flag and no `MC` tag. Both changes encode the record again as
//! htsjdk does (tags in [`crate::bam_raw::Tags`] order, qualities whose
//! first byte is `0xFF` written as all `0xFF`). The output header is the
//! input's with the new dictionary.
//!
//! Picard writes through a writer made with `presorted` false, so a
//! `coordinate` input is sorted again under the new dictionary,
//! `queryname` and `duplicate` too, `unsorted` is written as read, and
//! `unknown` fails on a second record ([`record_edit::WriterOrder`]).
//! Picard reads an indexed input contig by contig in the order of the
//! input's dictionary, then the unplaced records: for a coordinate-sorted
//! input that is file order, which bammetrics reads in every case.
//!
//! In a chain the step sorts the stream when its order names a sort
//! ([`record_edit::sorting_stage`]) and streams an unsorted one.

use std::path::{Path, PathBuf};

use noodles::sam;

use crate::add_or_replace_read_groups::htsjdk_qualities;
use crate::bam_raw::RawRecord;
use crate::bed_to_interval_list::{self, Sequence};
use crate::error::{Error, Result};
use crate::header_merge::SortOrder;
use crate::metric::ToolOutput;
use crate::record_edit::{self, Outcome, RecordEdit, Reorder, Setup, WriterOrder};
use crate::registry::{self, Invocation, Opt, Tool};
use crate::spill::TempDir;
use crate::stage::{Frame, Stage};

pub const TOOL: Tool = Tool {
    name: "ReorderSam",
    summary: "no rows; OUTPUT gets every record of the input under the sequence dictionary of \
              SEQUENCE_DICTIONARY, references matched by name, a record on a contig the \
              dictionary lacks made unplaced, sorted again in the input's order: a BAM, SAM \
              text for a name ending in .sam, or a CRAM (with REFERENCE_SEQUENCE) for one ending \
              in .cram.",
    options: &[
        Opt::path(
            &["SEQUENCE_DICTIONARY", "SD"],
            "The dictionary of OUTPUT: a .dict, a FASTA with its .dict beside it, a SAM, BAM \
             or CRAM, an interval list or a VCF.",
        )
        .required(),
        Opt::bool(
            &["ALLOW_INCOMPLETE_DICT_CONCORDANCE", "S"],
            "Let input contigs the dictionary lacks go: their records become unplaced \
             (default false).",
        ),
        Opt::bool(
            &["ALLOW_CONTIG_LENGTH_DISCORDANCE", "U"],
            "Map a contig to the dictionary's contig of the same name and another length \
             (default false).",
        ),
        record_edit::OUTPUT_OPTION,
        record_edit::COMPRESSION_OPTION,
        record_edit::INDEX_OPTION,
        record_edit::MD5_OPTION,
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records held in memory before a sorted run spills to TMP_DIR, 1 or more \
             (default 500000).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the spilled runs of the sort (default: the system temporary \
             directory).",
        )
        .repeated(),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// The options the edit takes.
#[derive(Clone, Debug)]
struct Options {
    dictionary: PathBuf,
    allow_incomplete: bool,
    allow_length: bool,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let dictionary = inv
            .get("SEQUENCE_DICTIONARY")
            .map(PathBuf::from)
            .ok_or_else(|| Error::Usage("SEQUENCE_DICTIONARY is required".to_owned()))?;
        Ok(Self {
            dictionary,
            allow_incomplete: inv.flag("ALLOW_INCOMPLETE_DICT_CONCORDANCE", false)?,
            allow_length: inv.flag("ALLOW_CONTIG_LENGTH_DISCORDANCE", false)?,
        })
    }
}

/// htsjdk's `FileExtensions.FASTA`.
const FASTA_EXTENSIONS: [&str; 10] = [
    ".fasta",
    ".fasta.gz",
    ".fas",
    ".fas.gz",
    ".fa",
    ".fa.gz",
    ".fna",
    ".fna.gz",
    ".txt",
    ".txt.gz",
];

/// The file types of `SAMSequenceDictionaryExtractor.TYPE`, in the order
/// `forFile` tries them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DictionaryKind {
    /// A FASTA: the `.dict` beside it, its extension replaced.
    Fasta(usize),
    Dictionary,
    Cram,
    Sam,
    Vcf,
    IntervalList,
}

impl DictionaryKind {
    /// `TYPE.forFile`: the first type with an extension the name ends with,
    /// case as written.
    fn of(name: &str) -> Option<Self> {
        let ends = |extensions: &[&str]| extensions.iter().any(|e| name.ends_with(e));
        if let Some(e) = FASTA_EXTENSIONS.iter().find(|e| name.ends_with(*e)) {
            return Some(Self::Fasta(e.len()));
        }
        [
            (&[".dict"][..], Self::Dictionary),
            (&[".cram"], Self::Cram),
            (&[".sam", ".bam"], Self::Sam),
            (&[".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"], Self::Vcf),
            (&[".interval_list", ".interval_list.gz"], Self::IntervalList),
        ]
        .into_iter()
        .find_map(|(extensions, kind)| ends(extensions).then_some(kind))
    }
}

fn not_found(path: &Path, message: String) -> Error {
    Error::Input {
        path: path.to_path_buf(),
        source: std::io::Error::new(std::io::ErrorKind::NotFound, message),
    }
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`.
fn read_sequences(path: &Path) -> Result<Vec<Sequence>> {
    let text = path.to_string_lossy();
    let Some(kind) = DictionaryKind::of(&text) else {
        return Err(Error::Input {
            path: path.to_path_buf(),
            source: std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "Cannot figure out type of file {} from extension. Current implementation \
                     understands the following types: FASTA, DICTIONARY, CRAM, SAM, VCF, \
                     INTERVAL_LIST",
                    path.display()
                ),
            ),
        });
    };
    match kind {
        DictionaryKind::Fasta(extension) => {
            if !path.is_file() {
                return Err(not_found(
                    path,
                    format!("Reference file {} not found", path.display()),
                ));
            }
            let dict = PathBuf::from(format!("{}.dict", &text[..text.len() - extension]));
            if !dict.is_file() {
                return Err(not_found(
                    path,
                    format!(
                        "Could not find dictionary next to reference file {}",
                        path.display()
                    ),
                ));
            }
            bed_to_interval_list::read_sequences(&dict)
        }
        DictionaryKind::Vcf => crate::vcf::read_dictionary(path)?
            .into_iter()
            .map(|contig| {
                let length = contig.length.ok_or_else(|| Error::Input {
                    path: path.to_path_buf(),
                    source: std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("##contig={} has no length", contig.name),
                    ),
                })?;
                Ok(Sequence {
                    name: contig.name,
                    length,
                    other_fields: Vec::new(),
                })
            })
            .collect(),
        DictionaryKind::Dictionary
        | DictionaryKind::Cram
        | DictionaryKind::Sam
        | DictionaryKind::IntervalList => bed_to_interval_list::read_sequences(path),
    }
}

/// The dictionary of `path` as header `@SQ` records, each with the fields
/// the file gives in its order.
fn read_dictionary(path: &Path) -> Result<sam::header::ReferenceSequences> {
    let sequences = read_sequences(path)?;
    let mut text = String::new();
    for sequence in &sequences {
        if sequence.length < 1 || sequence.length > i64::from(i32::MAX) {
            return Err(Error::Input {
                path: path.to_path_buf(),
                source: std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "SEQUENCE_DICTIONARY @SQ SN:{} LN:{}: a BAM header holds a length of 1 \
                         to {}",
                        sequence.name,
                        sequence.length,
                        i32::MAX
                    ),
                ),
            });
        }
        text.push_str("@SQ\tSN:");
        text.push_str(&sequence.name);
        text.push_str(&format!("\tLN:{}", sequence.length));
        for field in &sequence.other_fields {
            text.push('\t');
            text.push_str(field);
        }
        text.push('\n');
    }
    let header: sam::Header = text.parse().map_err(|e| Error::Input {
        path: path.to_path_buf(),
        source: std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("its sequence dictionary does not parse: {e}"),
        ),
    })?;
    Ok(header.reference_sequences().clone())
}

/// `buildSequenceDictionaryMap`: the new index of each input reference, -1
/// for one the new dictionary lacks when that is allowed.
fn sequence_map(
    options: &Options,
    reads: &sam::header::ReferenceSequences,
    reference: &sam::header::ReferenceSequences,
) -> Result<Vec<i32>> {
    let mut map = vec![None; reads.len()];
    for (new, (name, sequence)) in reference.iter().enumerate() {
        let Some((old, _, read)) = reads.get_full(name) else {
            continue;
        };
        if read.length() != sequence.length() {
            let message = format!(
                "Discordant contig lengths: read {name} LN={}, ref {name} LN={}",
                read.length(),
                sequence.length()
            );
            if !options.allow_length {
                return Err(Error::Metric(message));
            }
            eprintln!("WARNING ReorderSam {message}");
        }
        map[old] = Some(
            i32::try_from(new)
                .map_err(|_| Error::Metric(format!("{name}: more than 2^31 sequences")))?,
        );
    }
    map.into_iter()
        .zip(reads.keys())
        .map(|(new, name)| match new {
            Some(new) => Ok(new),
            None if options.allow_incomplete => Ok(-1),
            None => Err(Error::Metric(format!(
                "New reference sequence does not contain a matching contig for {name}"
            ))),
        })
        .collect()
}

const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;

/// The edit: the new dictionary in the header, each record's references
/// moved to it.
struct Reordering {
    dictionary: sam::header::ReferenceSequences,
    map: Vec<i32>,
}

impl Reordering {
    fn new_index(&self, id: i32) -> i32 {
        usize::try_from(id)
            .ok()
            .and_then(|at| self.map.get(at).copied())
            .unwrap_or(-1)
    }
}

impl RecordEdit for Reordering {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        let mut header = input.clone();
        *header.reference_sequences_mut() = self.dictionary.clone();
        Ok(header)
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        let (reference, mate) = (record.reference_id(), record.mate_reference_id());
        let (new_reference, new_mate) = (self.new_index(reference), self.new_index(mate));
        record.set_reference_ids(new_reference, new_mate);
        let lost = reference != -1 && new_reference == -1;
        let mate_lost = mate != -1 && new_mate == -1;
        if lost {
            record.set_alignment_start(0);
            record.set_flags(record.flags() | UNMAPPED);
            record.set_cigar(&[]).map_err(Error::Record)?;
            record.set_mapping_quality(0);
        }
        if mate_lost {
            record.set_mate_alignment_start(0);
            record.set_flags(record.flags() | MATE_UNMAPPED);
        }
        if lost || mate_lost {
            let mut tags = record.tags().map_err(Error::Record)?;
            if mate_lost {
                tags.remove(*b"MC");
            }
            record.set_tags(&tags).map_err(Error::Record)?;
            htsjdk_qualities(record);
        }
        Ok(
            if (new_reference, new_mate) == (reference, mate) && !lost && !mate_lost {
                Outcome::Same
            } else {
                Outcome::Changed
            },
        )
    }

    fn rewrites(&self) -> bool {
        true
    }
}

/// The edit's builder: the dictionary read once, matched against the
/// header of the records read.
fn build(
    options: Options,
    dictionary: sam::header::ReferenceSequences,
) -> impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>> {
    move |setup| {
        let map = sequence_map(&options, setup.header.reference_sequences(), &dictionary)?;
        Ok(Box::new(Reordering { dictionary, map }))
    }
}

/// The writer Picard makes: `presorted` false in the input's order.
const WRITER: WriterOrder = WriterOrder::Header {
    order: None,
    presorted_if_same: false,
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let options = Options::from_invocation(inv)?;
    let mut opts = record_edit::Options::from_invocation(inv)?;
    opts.reads = vec![options.dictionary.clone()];
    record_edit::check_output(&opts)?;
    let dictionary = read_dictionary(&options.dictionary)?;
    let (limits, tmp_dir) = record_edit::sort_limits(inv)?;
    opts.reorder = Reorder {
        regroup: None,
        writer: WRITER,
        limits,
        tmp_dir,
    };
    record_edit::run(&opts, build(options, dictionary))?;
    Ok(ToolOutput::NONE)
}

/// The usage errors of the step, asked before the first read of a chain:
/// an option of the written file, a `MAX_RECORDS_IN_RAM` below 1.
pub fn check_stage(inv: &Invocation) -> Result<()> {
    record_edit::check_stage(inv)?;
    Options::from_invocation(inv)?;
    record_edit::stage_limits(inv, 1).map(|_| ())
}

/// `ReorderSam` as a step of a chain: the edit on each batch, then the
/// sort the tool's writer makes in the stream's order; an unsorted stream
/// streams on.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let options = Options::from_invocation(inv)?;
    let dictionary = read_dictionary(&options.dictionary)?;
    let writer = WRITER.against(SortOrder::of(&frame.header));
    let stage = record_edit::stage(inv, frame, build(options, dictionary))?;
    if writer == WriterOrder::Presorted {
        return Ok(stage);
    }
    let Stage::Transform(transform) = stage else {
        return Err(Error::Metric(
            "ReorderSam: the edit of a step is a transform".to_owned(),
        ));
    };
    let reorder = Reorder {
        writer,
        limits: record_edit::stage_limits(inv, 1)?,
        ..Reorder::default()
    };
    record_edit::sorting_stage(inv, frame, temp, &reorder, transform)
}

#[cfg(test)]
mod tests;
