//! `FilterSamReads`, mirroring Picard 3.5.0 (`picard.sam.FilterSamReads`
//! over htsjdk 5.0.0's `FilteringSamIterator` and filters).
//!
//! `FILTER` picks the records of `INPUT` that `OUTPUT` gets:
//!
//! * `includeReadList`, `excludeReadList`: the records whose name is (is
//!   not) in `READ_LIST_FILE` (`ReadNameFilter`). Each line that is not
//!   blank gives the text before its first whitespace character (`\s+`), so
//!   a line that starts with a space or a tab names the empty name, which
//!   no record has; no comment syntax. The file is decoded as UTF-8 and a
//!   name matches a record name whose bytes are its characters (a name with
//!   a character above U+00FF matches none). A name ending in `.gz`,
//!   `.gzip` or `.bfq` is inflated first; an empty file fails; at most
//!   [`MAX_READ_LIST_BYTES`] bytes and lines of at most
//!   [`MAX_READ_LIST_LINE_BYTES`] are read. The names go into a hash set.
//! * `includeTagValues`, `excludeTagValues`: the records whose `TAG` holds
//!   (does not hold) one of `TAG_VALUE` (`TagFilter`). htsjdk compares the
//!   tag's Java object with the strings of `TAG_VALUE`, so only a string
//!   (`Z`) value can match: an `i`, `f`, `A`, `H` or `B` tag matches no
//!   value. A `TAG` of other than two characters fails at the first record.
//! * `includePairedIntervals`: primary paired records whose alignment, or
//!   whose mate's (start and `MC`), shares a base with an interval of
//!   `INTERVAL_LIST` (`IntervalKeepPairFilter`), so both ends of a pair are
//!   kept or dropped together. Unpaired, secondary and supplementary
//!   records are dropped; a mapped mate without `MC` fails. The header's
//!   order, after `SORT_ORDER`, must be `coordinate`.
//! * `includeAligned`, `excludeAligned`: aligned (unmapped) records, with
//!   the pairs of a queryname-sorted input judged together: a first-of-pair
//!   record and the record after it (which must be its mate, second of the
//!   pair) are kept when both are aligned (`includeAligned`) or when either
//!   is unmapped (`excludeAligned`); a second-of-pair record that does not
//!   follow its first is kept. The records must come in queryname order.
//! * `includeJavascript` is refused (exit 2): bammetrics runs no
//!   JavaScript, and `JAVASCRIPT_FILE` is refused with it (docs/07 Q6).
//!
//! The output header is the input's with `SO` set to `SORT_ORDER` (default
//! the input's). The writer is presorted when the order is the input's
//! (htsjdk's order check) and sorts otherwise (`MAX_RECORDS_IN_RAM`,
//! `TMP_DIR`). With `WRITE_READS_FILES` (default false) the SAM line of
//! every record read goes to `<INPUT basename>.reads` and of every record
//! written to `<OUTPUT basename>.reads`, both beside `OUTPUT`; a name equal
//! to the other's holds the output's lines. No `@PG` line is added.
//!
//! In a chain the step leaves `filtered` and streams, or holds and sorts
//! its records when `SORT_ORDER` differs from the stream's order.
//! `WRITE_READS_FILES` is a usage error there.

use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use noodles::sam;

use crate::bam_raw::{HtsjdkTags, RawRecord, TagIter, TagValue};
use crate::cleanup::AtomicFile;
use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::header_merge::{HeaderTag, SortOrder, set_header_tag};
use crate::interval_ops::IntervalList;
use crate::interval_ops::ops::OverlapTree;
use crate::mark_duplicates::same_file;
use crate::metric::ToolOutput;
use crate::pass;
use crate::record_edit::{
    self, Mode, Outcome, Pipeline, RecordEdit, Reorder, Restore, Setup, WriterOrder,
};
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_text::{parse_cigar, write_record};
use crate::spill::TempDir;
use crate::stage::{Frame, Stage};

/// The most bytes of `READ_LIST_FILE` read, after inflation.
pub const MAX_READ_LIST_BYTES: u64 = 1 << 30;
/// The longest line of `READ_LIST_FILE`.
pub const MAX_READ_LIST_LINE_BYTES: usize = 16 << 20;

pub const TOOL: Tool = Tool {
    name: "FilterSamReads",
    summary: "no rows; OUTPUT gets the records of INPUT that FILTER keeps (aligned or unmapped \
              reads, a read name list, tag values, pairs over intervals), in SORT_ORDER: a BAM, \
              or SAM text for a name ending in .sam.",
    options: &[
        Opt::text(
            &["FILTER"],
            "includeAligned, excludeAligned, includeReadList, excludeReadList, \
             includeTagValues, excludeTagValues or includePairedIntervals (includeJavascript \
             is refused).",
        )
        .required(),
        Opt::path(
            &["READ_LIST_FILE", "RLF"],
            "Read names, one per line, for FILTER=includeReadList or excludeReadList.",
        ),
        Opt::path(
            &["INTERVAL_LIST", "IL"],
            "Picard interval list for FILTER=includePairedIntervals.",
        ),
        Opt::text(
            &["TAG", "T"],
            "The two-character tag of FILTER=includeTagValues or excludeTagValues.",
        ),
        Opt::text(
            &["TAG_VALUE", "TV"],
            "A value of TAG to keep or drop; repeat for several. Only a string (Z) tag matches.",
        )
        .repeated(),
        Opt::text(
            &["SORT_ORDER", "SO"],
            "The order of OUTPUT: unsorted, queryname, coordinate, duplicate or unknown \
             (default: the input's).",
        ),
        Opt::bool(
            &["WRITE_READS_FILES"],
            "Write the SAM line of each record read to <INPUT basename>.reads and of each \
             record written to <OUTPUT basename>.reads, beside OUTPUT (default false).",
        ),
        Opt::rejected(&["JAVASCRIPT_FILE", "JS"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        record_edit::OUTPUT_OPTION,
        record_edit::COMPRESSION_OPTION,
        record_edit::INDEX_OPTION,
        record_edit::MD5_OPTION,
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records held in memory before a sorted run spills to TMP_DIR when SORT_ORDER sorts, \
             1 or more (default 500000).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the spilled runs (default: the system temporary directory).",
        )
        .repeated(),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST: u16 = 0x40;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

/// `FilterSamReads.Filter`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Filter {
    IncludeAligned,
    ExcludeAligned,
    IncludeReadList,
    ExcludeReadList,
    IncludePairedIntervals,
    IncludeTagValues,
    ExcludeTagValues,
}

const FILTER_NAMES: &str = "includeAligned, excludeAligned, includeReadList, excludeReadList, \
                            includeJavascript, includePairedIntervals, includeTagValues, \
                            excludeTagValues";

impl Filter {
    fn parse(text: &str) -> Result<Self> {
        Ok(match text {
            "includeAligned" => Self::IncludeAligned,
            "excludeAligned" => Self::ExcludeAligned,
            "includeReadList" => Self::IncludeReadList,
            "excludeReadList" => Self::ExcludeReadList,
            "includePairedIntervals" => Self::IncludePairedIntervals,
            "includeTagValues" => Self::IncludeTagValues,
            "excludeTagValues" => Self::ExcludeTagValues,
            "includeJavascript" => {
                return Err(Error::Usage(
                    "FILTER=includeJavascript: bammetrics runs no JavaScript, and \
                     JAVASCRIPT_FILE is refused (docs/07_full_parity_program.md Q6)"
                        .to_owned(),
                ));
            }
            other => {
                return Err(Error::Usage(format!(
                    "FILTER={other}: expected one of {FILTER_NAMES}"
                )));
            }
        })
    }

    fn name(self) -> &'static str {
        match self {
            Self::IncludeAligned => "includeAligned",
            Self::ExcludeAligned => "excludeAligned",
            Self::IncludeReadList => "includeReadList",
            Self::ExcludeReadList => "excludeReadList",
            Self::IncludePairedIntervals => "includePairedIntervals",
            Self::IncludeTagValues => "includeTagValues",
            Self::ExcludeTagValues => "excludeTagValues",
        }
    }

    /// `FilteringSamIterator(iterator, filter, true)`: the pairs of a
    /// queryname-sorted input judged together.
    fn by_pair(self) -> bool {
        matches!(self, Self::IncludeAligned | Self::ExcludeAligned)
    }
}

/// The options both faces share, checked as Picard's
/// `customCommandLineValidation` checks them.
#[derive(Clone, Debug)]
struct FilterOptions {
    filter: Filter,
    read_list: Option<PathBuf>,
    interval_list: Option<PathBuf>,
    tag: Option<String>,
    tag_values: Vec<String>,
    order: Option<SortOrder>,
    write_reads: bool,
}

/// `checkInputs`: an option given with a `FILTER` that does not read it,
/// or missing for one that does.
fn check_input(filter: Filter, takes: &[Filter], given: bool, name: &str) -> Result<()> {
    let wanted = takes.contains(&filter);
    if wanted && !given {
        return Err(Error::Usage(format!(
            "{name} must be specified when using FILTER={}, but it was null.",
            filter.name()
        )));
    }
    if !wanted && given {
        let from: Vec<&str> = takes.iter().map(|f| f.name()).collect();
        return Err(Error::Usage(format!(
            "{name} may only be specified when using FILTER from {}, FILTER value: {}",
            from.join(","),
            filter.name()
        )));
    }
    Ok(())
}

impl FilterOptions {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let filter = Filter::parse(
            inv.get("FILTER")
                .ok_or_else(|| Error::Usage("Option 'FILTER' is required.".to_owned()))?,
        )?;
        let read_list = inv.get("READ_LIST_FILE").map(PathBuf::from);
        let interval_list = inv.get("INTERVAL_LIST").map(PathBuf::from);
        let tag = inv.get("TAG").map(str::to_owned);
        use Filter::*;
        check_input(
            filter,
            &[IncludeReadList, ExcludeReadList],
            read_list.is_some(),
            "READ_LIST_FILE",
        )?;
        check_input(
            filter,
            &[IncludePairedIntervals],
            interval_list.is_some(),
            "INTERVAL_LIST",
        )?;
        check_input(
            filter,
            &[IncludeTagValues, ExcludeTagValues],
            tag.is_some(),
            "TAG",
        )?;
        let order = inv
            .get("SORT_ORDER")
            .filter(|value| *value != "null")
            .map(|value| {
                SortOrder::parse(value).ok_or_else(|| {
                    Error::Usage(format!("'{value}' is not a valid value for SortOrder."))
                })
            })
            .transpose()?;
        Ok(Self {
            filter,
            read_list,
            interval_list,
            tag,
            tag_values: inv.all("TAG_VALUE").to_vec(),
            order,
            write_reads: inv.flag("WRITE_READS_FILES", false)?,
        })
    }

    /// The files the filter reads beside `INPUT`.
    fn reads(&self) -> Vec<PathBuf> {
        self.read_list
            .iter()
            .chain(&self.interval_list)
            .cloned()
            .collect()
    }

    /// The order of the output header: `SORT_ORDER`, else the input's.
    fn output_order(&self, input: &sam::Header) -> SortOrder {
        self.order.unwrap_or_else(|| SortOrder::of(input))
    }

    /// The rule, its files read; `header` is the input's.
    fn rule(&self, header: &sam::Header) -> Result<Rule> {
        Ok(match self.filter {
            Filter::IncludeAligned => Rule::Aligned(true),
            Filter::ExcludeAligned => Rule::Aligned(false),
            Filter::IncludeReadList | Filter::ExcludeReadList => Rule::Names {
                names: read_names(self.read_list.as_deref().expect("checked"))?,
                include: self.filter == Filter::IncludeReadList,
            },
            Filter::IncludeTagValues | Filter::ExcludeTagValues => Rule::Tag {
                tag: java_tag(self.tag.as_deref().expect("checked")),
                values: self
                    .tag_values
                    .iter()
                    .filter_map(|value| latin1(value))
                    .collect(),
                include: self.filter == Filter::IncludeTagValues,
                tags: HtsjdkTags::default(),
            },
            Filter::IncludePairedIntervals => Rule::Intervals(Intervals::new(
                self.interval_list.as_deref().expect("checked"),
                header,
            )?),
        })
    }
}

/// The bytes a Java string has as a record field read from a BAM, where
/// htsjdk makes one character per byte: `None` when a character is above
/// U+00FF, which no record field can hold.
fn latin1(text: &str) -> Option<Vec<u8>> {
    text.chars()
        .map(|c| u8::try_from(u32::from(c)).ok())
        .collect()
}

/// `TAG` as `SAMTag.makeBinaryTag` takes it: two UTF-16 units, each a
/// byte of a record's tag (a unit above `0xFF` names no tag of any record).
/// Any other length fails when a record is filtered.
#[derive(Clone, Debug)]
enum JavaTag {
    Bytes([u8; 2]),
    NoRecordHasIt,
    WrongLength(String),
}

fn java_tag(text: &str) -> JavaTag {
    let units: Vec<u16> = text.encode_utf16().collect();
    match units.as_slice() {
        [a, b] => match (u8::try_from(*a), u8::try_from(*b)) {
            (Ok(a), Ok(b)) => JavaTag::Bytes([a, b]),
            _ => JavaTag::NoRecordHasIt,
        },
        _ => JavaTag::WrongLength(text.to_owned()),
    }
}

/// `ReadNameFilter`'s set of `READ_LIST_FILE`, as [`latin1`] bytes.
fn read_names(path: &Path) -> Result<HashSet<Vec<u8>>> {
    let in_error = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = File::open(path).map_err(in_error)?;
    if file.metadata().map_err(in_error)?.len() == 0 {
        return Err(in_error(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} has length 0", path.display()),
        )));
    }
    let gzip = path
        .file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| [".gz", ".gzip", ".bfq"].iter().any(|e| name.ends_with(e)));
    let inner: Box<dyn Read> = if gzip {
        Box::new(crate::gzip::GzipReader::new(BufReader::new(file)))
    } else {
        Box::new(file)
    };
    let mut reader = BufReader::new(inner.take(MAX_READ_LIST_BYTES + 1));
    let mut names = HashSet::new();
    let mut line = Vec::new();
    let mut total: u64 = 0;
    loop {
        line.clear();
        let read = (&mut reader)
            .take(MAX_READ_LIST_LINE_BYTES as u64 + 1)
            .read_until(b'\n', &mut line)
            .map_err(in_error)?;
        if read == 0 {
            break;
        }
        total += read as u64;
        if total > MAX_READ_LIST_BYTES {
            return Err(in_error(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("the read list holds more than {MAX_READ_LIST_BYTES} bytes, the most read"),
            )));
        }
        if line.len() > MAX_READ_LIST_LINE_BYTES {
            return Err(in_error(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("a line of more than {MAX_READ_LIST_LINE_BYTES} bytes"),
            )));
        }
        let text = line.strip_suffix(b"\n").unwrap_or(&line);
        // `readLine` also ends a line at `\r`; the empty lines this makes
        // are blank, which `ReadNameFilter` skips.
        for piece in text.split(|&b| b == b'\r') {
            if let Some(name) = first_name(piece) {
                names.insert(name);
            }
        }
    }
    Ok(names)
}

/// `line.trim().isEmpty()` then `line.split("\\s+")[0]` of one line, as
/// [`latin1`] bytes; `None` for a blank line or a name no record has.
fn first_name(line: &[u8]) -> Option<Vec<u8>> {
    if line.iter().all(|&b| b <= b' ') {
        return None;
    }
    let end = line
        .iter()
        .position(|b| matches!(b, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r'))
        .unwrap_or(line.len());
    if end == 0 {
        return None;
    }
    latin1(&crate::java::decode_utf8(&line[..end]))
}

/// `INTERVAL_LIST` in `IntervalKeepPairFilter`'s `OverlapDetector`, its
/// contigs found by name from the input's references.
struct Intervals {
    tree: OverlapTree,
    /// Per reference of the input header, the contig of the list.
    contigs: Vec<Option<usize>>,
}

impl Intervals {
    fn new(path: &Path, header: &sam::Header) -> Result<Self> {
        let list = IntervalList::read(path)?;
        let index = list.header.name_index();
        let contigs = header
            .reference_sequences()
            .keys()
            .map(|name| {
                index
                    .get(&String::from_utf8_lossy(name).into_owned())
                    .copied()
            })
            .collect();
        Ok(Self {
            tree: OverlapTree::new(&list.intervals, list.header.sequences.len()),
            contigs,
        })
    }

    /// `hasOverlaps(reference, start, end)`.
    fn overlaps(&self, reference: i32, start: i32, end: i32) -> bool {
        let Some(contig) = usize::try_from(reference)
            .ok()
            .and_then(|at| self.contigs.get(at).copied().flatten())
        else {
            return false;
        };
        let mut found = false;
        self.tree
            .for_each_overlap(contig, start, end, |_| found = true);
        found
    }
}

/// One record's test.
enum Rule {
    /// `AlignedFilter(include)`.
    Aligned(bool),
    Names {
        names: HashSet<Vec<u8>>,
        include: bool,
    },
    Tag {
        tag: JavaTag,
        values: Vec<Vec<u8>>,
        include: bool,
        tags: HtsjdkTags,
    },
    Intervals(Intervals),
}

fn failure(record: &RawRecord, what: &str) -> Error {
    Error::Metric(format!(
        "{}: record {}: {what}",
        TOOL.name,
        String::from_utf8_lossy(record.name())
    ))
}

impl Rule {
    /// `filterOut(record)`: true drops it.
    fn out(&mut self, record: &RawRecord) -> Result<bool> {
        let flags = record.flags();
        Ok(match self {
            Self::Aligned(include) => (flags & UNMAPPED != 0) == *include,
            Self::Names { names, include } => names.contains(record.name()) != *include,
            Self::Tag {
                tag,
                values,
                include,
                tags,
            } => {
                let found = match tag {
                    JavaTag::WrongLength(text) => {
                        return Err(Error::Metric(format!(
                            "String tag does not have length() == 2: {text}"
                        )));
                    }
                    JavaTag::NoRecordHasIt => false,
                    JavaTag::Bytes(name) => string_tag(record, *name, tags)?
                        .is_some_and(|value| values.contains(&value)),
                };
                found != *include
            }
            Self::Intervals(intervals) => intervals.out(record)?,
        })
    }

    /// `filterOut(first, second)` of `AlignedFilter`.
    fn pair_out(&self, first: &RawRecord, second: &RawRecord) -> bool {
        let unmapped = |r: &RawRecord| r.flags() & UNMAPPED != 0;
        match self {
            Self::Aligned(true) => unmapped(first) || unmapped(second),
            Self::Aligned(false) => !(unmapped(first) || unmapped(second)),
            _ => false,
        }
    }
}

/// `getAttribute(tag)` when it is a `String`: the value of the field
/// htsjdk's tag list keeps ([`HtsjdkTags`]).
fn string_tag(record: &RawRecord, name: [u8; 2], tags: &mut HtsjdkTags) -> Result<Option<Vec<u8>>> {
    tags.try_load(TagIter::new(record.tag_bytes()).map(|field| field.map(|(n, _)| n)))
        .map_err(Error::Record)?;
    let Some(at) = tags.position(name) else {
        return Ok(None);
    };
    let field = TagIter::new(record.tag_bytes())
        .nth(at)
        .transpose()
        .map_err(Error::Record)?;
    Ok(match field {
        Some((_, TagValue::String(text))) => Some(text),
        _ => None,
    })
}

impl Intervals {
    /// `IntervalKeepPairFilter.filterOut`.
    fn out(&self, record: &RawRecord) -> Result<bool> {
        let flags = record.flags();
        if flags & (SECONDARY | SUPPLEMENTARY) != 0 || flags & PAIRED == 0 {
            return Ok(true);
        }
        if flags & UNMAPPED == 0 {
            let start = record.alignment_start();
            let end = start
                .wrapping_add(record.reference_length())
                .wrapping_sub(1);
            if self.overlaps(record.reference_id(), start, end) {
                return Ok(false);
            }
        }
        if flags & MATE_UNMAPPED != 0 {
            return Ok(true);
        }
        let start = record.mate_alignment_start();
        let end = start
            .wrapping_add(mate_reference_length(record)?)
            .wrapping_sub(1);
        Ok(!self.overlaps(record.mate_reference_id(), start, end))
    }
}

/// The reference length of the mate's CIGAR in `MC`
/// (`SAMUtils.getMateAlignmentEnd`).
fn mate_reference_length(record: &RawRecord) -> Result<i32> {
    match record.tag(*b"MC").map_err(Error::Record)? {
        None => Err(failure(record, "Mate CIGAR (Tag MC) not found")),
        Some(TagValue::String(text)) => {
            let cigar = parse_cigar(&text).map_err(|e| failure(record, &e))?;
            Ok(crate::bam_raw::reference_length(cigar.into_iter()))
        }
        Some(_) => Err(failure(record, "Value for tag MC is not a String")),
    }
}

/// `SamPairUtil.assertMate(first, second)`.
fn assert_mate(first: &RawRecord, second: &RawRecord) -> Result<()> {
    let line = |r: &RawRecord| String::from_utf8_lossy(r.name()).into_owned();
    let error = |message: String| Err(Error::Metric(format!("{}: {message}", TOOL.name)));
    if first.flags() & PAIRED == 0 {
        return error(format!(
            "First record is not a paired read: {}",
            line(first)
        ));
    }
    if second.flags() & PAIRED == 0 {
        return error(format!(
            "Second record is not a paired read: {}",
            line(second)
        ));
    }
    if first.flags() & FIRST == 0 {
        return error(format!(
            "First record is not marked as first of pair: {}",
            line(first)
        ));
    }
    if second.flags() & SECOND == 0 {
        return error(format!(
            "Second record is not marked as second of pair: {}",
            line(second)
        ));
    }
    if first.name() != second.name() {
        return error(format!(
            "Read names do not match: {} : {}",
            line(first),
            line(second)
        ));
    }
    Ok(())
}

/// `FilteringSamIterator` over the records as they come: each record is
/// [`Filtering::take`]n, and the records it keeps wait in `ready`.
struct Filtering {
    rule: Rule,
    by_pair: bool,
    /// A first-of-pair record waiting for the record after it.
    held: Option<RawRecord>,
    ready: VecDeque<RawRecord>,
    /// The last name read, for the queryname order a pair filter asserts.
    previous: Option<Vec<u8>>,
}

impl Filtering {
    fn new(rule: Rule, by_pair: bool) -> Self {
        Self {
            rule,
            by_pair,
            held: None,
            ready: VecDeque::new(),
            previous: None,
        }
    }

    /// The `AssertingIterator` that `assertSorted(queryname)` puts on the
    /// reader: names never go down (`fileOrderCompare`).
    fn check_order(&mut self, record: &RawRecord) -> Result<()> {
        if let Some(previous) = &self.previous {
            if previous.as_slice() > record.name() {
                return Err(Error::Metric(format!(
                    "{}: Record {} should come after {} when sorting with queryname ordering.",
                    TOOL.name,
                    String::from_utf8_lossy(previous),
                    String::from_utf8_lossy(record.name())
                )));
            }
            self.previous.as_mut().expect("some").clear();
        }
        self.previous
            .get_or_insert_with(Vec::new)
            .extend_from_slice(record.name());
        Ok(())
    }

    fn take(&mut self, record: RawRecord) -> Result<()> {
        if !self.by_pair {
            if !self.rule.out(&record)? {
                self.ready.push_back(record);
            }
            return Ok(());
        }
        self.check_order(&record)?;
        if let Some(first) = self.held.take() {
            assert_mate(&first, &record)?;
            if self.rule.pair_out(&first, &record) {
                return Ok(());
            }
            self.ready.push_back(first);
        }
        let flags = record.flags();
        if flags & PAIRED != 0 && flags & FIRST != 0 {
            self.held = Some(record);
            return Ok(());
        }
        self.alone(record)
    }

    /// A record judged without the record after it.
    fn alone(&mut self, record: RawRecord) -> Result<()> {
        let flags = record.flags();
        let kept = (self.by_pair && flags & PAIRED != 0 && flags & SECOND != 0)
            || !self.rule.out(&record)?;
        if kept {
            self.ready.push_back(record);
        }
        Ok(())
    }

    /// After the last record: a first-of-pair record with no record after
    /// it is judged alone.
    fn end(&mut self) -> Result<()> {
        match self.held.take() {
            Some(record) => self.alone(record),
            None => Ok(()),
        }
    }
}

/// The edit of the chain step: [`Filtering`] one record behind, so each
/// call hands on at most one record and `finish` the rest.
struct FilterEdit {
    filtering: Filtering,
    order: Option<SortOrder>,
}

/// The input header with `SO` set to `order` when one is given.
fn with_order(input: &sam::Header, order: Option<SortOrder>) -> sam::Header {
    let mut header = input.clone();
    if let Some(order) = order {
        set_header_tag(&mut header, HeaderTag::SortOrder, order.name());
    }
    header
}

impl RecordEdit for FilterEdit {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        Ok(with_order(input, self.order))
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        self.filtering.take(std::mem::take(record))?;
        match self.filtering.ready.pop_front() {
            Some(next) => {
                *record = next;
                Ok(Outcome::Same)
            }
            None => Ok(Outcome::Dropped),
        }
    }

    fn finish(&mut self, out: &mut Vec<RawRecord>) -> Result<()> {
        self.filtering.end()?;
        out.extend(self.filtering.ready.drain(..));
        Ok(())
    }

    fn rewrites(&self) -> bool {
        true
    }
}

/// `IOUtil.basename(file)` + `.reads` beside `output`.
fn reads_path(output: &Path, file: &Path) -> PathBuf {
    use std::os::unix::ffi::OsStrExt;
    let name = file.file_name().map_or(&[][..], |name| name.as_bytes());
    let base = match name.iter().rposition(|&b| b == b'.') {
        Some(dot) if dot > 0 => &name[..dot],
        _ => name,
    };
    let mut reads = base.to_vec();
    reads.extend_from_slice(b".reads");
    let reads = PathBuf::from(std::ffi::OsStr::from_bytes(&reads));
    match output.parent() {
        Some(dir) => dir.join(reads),
        None => reads,
    }
}

/// A `.reads` file being written: `SAMRecord.toString`, the SAM line, of
/// each record.
struct ReadsText {
    file: AtomicFile,
    out: BufWriter<File>,
    references: Vec<Vec<u8>>,
    line: Vec<u8>,
}

impl ReadsText {
    fn create(path: &Path, header: &sam::Header) -> Result<Self> {
        let file = AtomicFile::new(path);
        let out = BufWriter::new(file.create().map_err(|source| Error::Output {
            path: path.to_path_buf(),
            source,
        })?);
        Ok(Self {
            file,
            out,
            references: header
                .reference_sequences()
                .keys()
                .map(|name| name.to_vec())
                .collect(),
            line: Vec::new(),
        })
    }

    fn error(&self, source: io::Error) -> Error {
        Error::Output {
            path: self.file.path().to_path_buf(),
            source,
        }
    }

    fn write(&mut self, record: &RawRecord) -> Result<()> {
        self.line.clear();
        write_record(&mut self.line, record, &self.references).map_err(Error::Record)?;
        self.out
            .write_all(&self.line)
            .map_err(|source| self.error(source))
    }

    fn commit(self) -> Result<()> {
        let path = self.file.path().to_path_buf();
        let error = |source| Error::Output {
            path: path.clone(),
            source,
        };
        let file = self.out.into_inner().map_err(|e| error(e.into_error()))?;
        file.sync_all().map_err(error)?;
        drop(file);
        self.file.commit().map_err(error)
    }
}

/// A stand-alone run.
struct Run {
    filter: FilterOptions,
    edit: record_edit::Options,
}

impl Run {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let filter = FilterOptions::from_invocation(inv)?;
        let mut edit = record_edit::Options::from_invocation(inv)?;
        edit.reads = filter.reads();
        let (limits, tmp_dir) = record_edit::sort_limits(inv)?;
        edit.reorder = Reorder {
            writer: WriterOrder::Header {
                order: filter.order,
                presorted_if_same: true,
            },
            limits,
            tmp_dir,
            ..Reorder::default()
        };
        Ok(Self { filter, edit })
    }

    /// The `.reads` files: the input's (left out when it is the output's
    /// name) and the output's.
    fn reads_files(&self) -> (Option<PathBuf>, Option<PathBuf>) {
        if !self.filter.write_reads {
            return (None, None);
        }
        let input = reads_path(&self.edit.output, &self.edit.input);
        let output = reads_path(&self.edit.output, &self.edit.output);
        ((input != output).then_some(input), Some(output))
    }

    /// A `.reads` file may be no file the run reads, nor `OUTPUT`.
    fn check_reads_files(&self) -> Result<()> {
        let (input, output) = self.reads_files();
        let read = std::iter::once(&self.edit.input)
            .chain(&self.edit.reference)
            .chain(&self.edit.reads)
            .chain(std::iter::once(&self.edit.output));
        for named in read {
            for written in input.iter().chain(&output) {
                if written == named || same_file(named, written) {
                    return Err(Error::Usage(format!(
                        "WRITE_READS_FILES: {} is {}; writing it would destroy it",
                        written.display(),
                        named.display()
                    )));
                }
            }
        }
        Ok(())
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    run(&Run::from_invocation(inv)?)?;
    Ok(ToolOutput::NONE)
}

/// The tool on its own: one read of `INPUT`, the kept records written to
/// `OUTPUT` (and the `.reads` files) under temporary names, renamed when
/// every file is complete.
fn run(run: &Run) -> Result<()> {
    record_edit::check_output(&run.edit)?;
    run.check_reads_files()?;
    let edit = &run.edit;
    let temp = TempDir::new(edit.reorder.tmp_dir.as_deref(), "filter-sam-reads");
    pass::with_record_stream_kind(
        &edit.input,
        edit.reference.as_deref(),
        edit.threads,
        edit.stringency,
        |header, reader, kind| {
            record_edit::check_sort_order(header, edit.stringency, &edit.input)?;
            let order = run.filter.output_order(header);
            if run.filter.filter == Filter::IncludePairedIntervals && order != SortOrder::Coordinate
            {
                return Err(Error::Metric(
                    "Input must be coordinate sorted to use includePairedIntervals".to_owned(),
                ));
            }
            let mut filtering =
                Filtering::new(run.filter.rule(header)?, run.filter.filter.by_pair());
            let out_header = record_edit::with_version(with_order(header, run.filter.order), true);
            let cram = kind == Kind::Cram;
            let mut restore = Restore::new(cram, Mode::Rewrite, edit.reference.as_deref(), header)?;
            let (input_reads, output_reads) = run.reads_files();
            let mut input_reads = input_reads
                .map(|path| ReadsText::create(&path, header))
                .transpose()?;
            let mut output_reads = output_reads
                .map(|path| ReadsText::create(&path, &out_header))
                .transpose()?;
            let mut pipeline = Pipeline::new(edit, SortOrder::of(header), &out_header, &temp, None);
            if let Some(reads) = output_reads.as_mut() {
                pipeline.observe(Box::new(move |record| reads.write(record)));
            }
            let references = header.reference_sequences().len();
            let mut record = RawRecord::default();
            while record
                .read_from(reader)
                .map_err(|e| pass::read_error(&edit.input, cram, e))?
            {
                record_edit::check_reference_ids(&record, references)?;
                restore.apply(&mut record)?;
                if let Some(reads) = input_reads.as_mut() {
                    reads.write(&record)?;
                }
                filtering.take(std::mem::take(&mut record))?;
                while let Some(mut kept) = filtering.ready.pop_front() {
                    pipeline.push(&mut kept)?;
                }
            }
            filtering.end()?;
            while let Some(mut kept) = filtering.ready.pop_front() {
                pipeline.push(&mut kept)?;
            }
            let out = pipeline.finish()?;
            out.finish(&out_header, edit)?;
            for reads in input_reads.into_iter().chain(output_reads) {
                reads.commit()?;
            }
            Ok(())
        },
    )
}

/// The usage errors of the step, asked before the first read of a chain:
/// an option of the written file, `WRITE_READS_FILES`,
/// `SORT_ORDER=unknown`, and the `FILTER` checks of the tool.
pub fn check_stage(inv: &Invocation) -> Result<()> {
    record_edit::check_stage(inv)?;
    let options = FilterOptions::from_invocation(inv)?;
    if options.write_reads {
        return Err(Error::Usage(
            "WRITE_READS_FILES: the .reads files sit beside OUTPUT, and a chain step hands its \
             records on (wave W14 of docs/07_full_parity_program.md)"
                .to_owned(),
        ));
    }
    if options.order == Some(SortOrder::Unknown) {
        return Err(Error::Usage(
            "SORT_ORDER=unknown in a chain: Picard sorts unknown with no comparator and fails \
             on any input whose order is not unknown"
                .to_owned(),
        ));
    }
    if matches!(
        options.order,
        Some(SortOrder::Queryname | SortOrder::Coordinate | SortOrder::Duplicate)
    ) {
        record_edit::stage_limits(inv, 1)?;
    }
    Ok(())
}

/// `FilterSamReads` as a step of a chain: a transform on the stream at its
/// position, or a barrier that sorts when `SORT_ORDER` differs from the
/// stream's order.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let options = FilterOptions::from_invocation(inv)?;
    let order = options.output_order(&frame.header);
    if options.filter == Filter::IncludePairedIntervals && order != SortOrder::Coordinate {
        return Err(Error::Metric(
            "Input must be coordinate sorted to use includePairedIntervals".to_owned(),
        ));
    }
    let rule = options.rule(&frame.header)?;
    let by_pair = options.filter.by_pair();
    let filter_order = options.order;
    let build = move |_: &Setup| -> Result<Box<dyn RecordEdit>> {
        Ok(Box::new(FilterEdit {
            filtering: Filtering::new(rule, by_pair),
            order: filter_order,
        }))
    };
    let writer = WriterOrder::Header {
        order: options.order,
        presorted_if_same: true,
    }
    .against(SortOrder::of(&frame.header));
    let Stage::Transform(transform) = record_edit::stage(inv, frame, build)? else {
        return Err(Error::Metric(format!(
            "{}: the edit of a step is a transform",
            TOOL.name
        )));
    };
    if writer == WriterOrder::Presorted {
        return record_edit::presorted_stage(frame, transform, TOOL.name);
    }
    let reorder = Reorder {
        writer,
        limits: record_edit::stage_limits(inv, 1)?,
        ..Reorder::default()
    };
    record_edit::sorting_stage(inv, frame, temp, &reorder, transform)
}

#[cfg(test)]
mod tests;
