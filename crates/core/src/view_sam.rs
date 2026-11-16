//! `ViewSam`, mirroring Picard 3.5.0 (`picard.sam.ViewSam`).
//!
//! The header and every record of a SAM, BAM or CRAM as SAM text on
//! standard output, as htsjdk's `SAMTextHeaderCodec` and `SAMTextWriter`
//! write them ([`crate::sam_text`]): the header re-encoded with its own
//! `VN`, each record's tags in htsjdk's order with every integer as `i`.
//! A mapped CRAM record gets the `MD` and `NM` tags htsjdk rebuilds while
//! it decodes one.
//!
//! `ALIGNMENT_STATUS` (`Aligned`, `Unaligned`, `All`) keeps records by
//! their unmapped flag, `PF_STATUS` (`PF`, `NonPF`, `All`) by their
//! QC-fail flag. `HEADER_ONLY` prints the header and reads no record.
//! `RECORDS_ONLY` prints no header. Both together are an error.
//!
//! `INTERVAL_LIST` keeps the records that overlap an interval of the list
//! after `IntervalList.uniqued()`. Picard has two paths and they differ:
//!
//! - With an index beside the input (a `.bai` or `.csi` of a BAM, a
//!   `.crai` of a CRAM, found as `SamFiles.findIndex` finds them) it
//!   queries the index: a record matches when its alignment overlaps an
//!   interval, and an unmapped record placed with its mate matches by its
//!   start. bammetrics opens no index: it reads the file whole and applies
//!   the same test (`BAMQueryMultipleIntervalsIteratorFilter`), which gives
//!   the same records for a coordinate-sorted file.
//! - Without an index (SAM text, a pipe, a BAM without one) it scans with
//!   `IntervalFilter`: reading stops at the first record without a
//!   reference or past the last interval, and an unmapped record never
//!   matches, because its alignment end is 0.
//!
//! The intervals must be in the order of the input's `@SQ` lines, as both
//! Picard paths require.
//!
//! [`view`] streams. The registry `run` ([`TOOL`]) serves Python and
//! `--format`: it holds the text and one row per record ([`COLUMNS`]), at
//! most [`MAX_TEXT_BYTES`] of text.

use std::io::{self, Write};
use std::path::{Path, PathBuf};

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::cleanup::AtomicFile;
use crate::cram;
use crate::error::{Error, Result};
use crate::intervals;
use crate::mark_duplicates::{ContigCache, restore_md_and_nm};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::pass;
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_text;

pub const TOOL: Tool = Tool {
    name: "ViewSam",
    summary: "the header and each record of a SAM, BAM or CRAM as SAM text, as htsjdk writes \
              it; one row per record, the 11 SAM fields and the tags.",
    options: &[
        Opt::text(
            &["ALIGNMENT_STATUS"],
            "Aligned, Unaligned or All: the records printed, by their unmapped flag (default All).",
        ),
        Opt::text(
            &["PF_STATUS"],
            "PF, NonPF or All: the records printed, by their QC-fail flag (default All).",
        ),
        Opt::bool(
            &["HEADER_ONLY"],
            "Print the header and no record (default false).",
        ),
        Opt::bool(
            &["RECORDS_ONLY"],
            "Print the records and no header (default false).",
        ),
        Opt::path(
            &["INTERVAL_LIST"],
            "An interval list: only records that overlap one of its intervals are printed.",
        ),
        Opt::rejected(&["OUTPUT", "O"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

/// The columns of a row: the mandatory SAM fields, then every tag field
/// joined by tabs.
pub const COLUMNS: [&str; 12] = [
    "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL",
    "TAGS",
];

/// The SAM text the registry `run` holds at most, with its rows.
pub const MAX_TEXT_BYTES: usize = 256 << 20;

const UNMAPPED: u16 = 0x4;
const QC_FAIL: u16 = 0x200;

/// A three-way record filter on one flag.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Keep {
    All,
    /// Records with the flag set.
    Set,
    /// Records with the flag clear.
    Clear,
}

impl Keep {
    fn takes(self, flags: u16, flag: u16) -> bool {
        match self {
            Keep::All => true,
            Keep::Set => flags & flag != 0,
            Keep::Clear => flags & flag == 0,
        }
    }
}

/// The options of one run.
#[derive(Clone, Debug)]
pub struct Options {
    /// By the unmapped flag: `Aligned` keeps it clear.
    alignment: Keep,
    /// By the QC-fail flag: `PF` keeps it clear.
    pf: Keep,
    pub header_only: bool,
    pub records_only: bool,
    pub interval_list: Option<PathBuf>,
}

impl Options {
    /// Picard's enum values are case-sensitive.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let choice = |name: &str, clear: &str, set: &str| match inv.get(name) {
            None | Some("All") => Ok(Keep::All),
            Some(value) if value == clear => Ok(Keep::Clear),
            Some(value) if value == set => Ok(Keep::Set),
            Some(value) => Err(Error::Usage(format!(
                "{name}={value}: expected {clear}, {set} or All"
            ))),
        };
        let options = Self {
            alignment: choice("ALIGNMENT_STATUS", "Aligned", "Unaligned")?,
            pf: choice("PF_STATUS", "PF", "NonPF")?,
            header_only: inv.flag("HEADER_ONLY", false)?,
            records_only: inv.flag("RECORDS_ONLY", false)?,
            interval_list: inv
                .get("INTERVAL_LIST")
                .filter(|path| *path != "null")
                .map(PathBuf::from),
        };
        if options.header_only && options.records_only {
            return Err(Error::Usage(
                "Cannot specify both HEADER_ONLY=true and RECORDS_ONLY=true.".to_owned(),
            ));
        }
        Ok(options)
    }
}

/// True when htsjdk finds an index for `path` (`SamFiles.findIndex`): a
/// BAM or CRAM in a regular file with `x.bai` or `x.csi` (for `x.bam`),
/// `x.crai` or `x.cram.crai` (for `x.cram`), `NAME.bai` or `NAME.csi`
/// beside it, or the same names beside the file a symbolic link names.
fn has_index(path: &Path) -> bool {
    let beside = |file: &Path| {
        let Some(name) = file.file_name().and_then(|name| name.to_str()) else {
            return false;
        };
        let mut names = Vec::new();
        if let Some(stem) = name.strip_suffix(".bam") {
            names.extend([format!("{stem}.bai"), format!("{stem}.csi")]);
        } else if let Some(stem) = name.strip_suffix(".cram") {
            names.extend([format!("{stem}.crai"), format!("{name}.crai")]);
        }
        names.extend([format!("{name}.bai"), format!("{name}.csi")]);
        names
            .into_iter()
            .any(|index| file.with_file_name(index).is_file())
    };
    let regular = std::fs::metadata(path).is_ok_and(|m| m.is_file()) && !path.starts_with("/dev");
    regular
        && !cram::is_sam(path)
        && (beside(path) || path.canonicalize().is_ok_and(|real| beside(&real)))
}

/// What the interval filter says of a record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Step {
    Print,
    Skip,
    /// No later record can match.
    Stop,
}

/// One interval on the input's references; `reference` is -1 for a contig
/// the input's header lacks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Region {
    reference: i32,
    start: i32,
    end: i32,
}

/// The `INTERVAL_LIST` filter over records in file order.
struct IntervalFilter {
    regions: Vec<Region>,
    at: usize,
    /// True for Picard's index query, false for its scan.
    query: bool,
}

impl IntervalFilter {
    /// The uniqued intervals of `path` on the references of `header`.
    fn load(path: &Path, header: &sam::Header, query: bool) -> Result<Self> {
        let (names, list) = intervals::read_interval_list(path)?;
        let unique =
            intervals::uniqued(&list.iter().map(|named| named.interval).collect::<Vec<_>>());
        let references = header.reference_sequences();
        let regions: Vec<Region> = unique
            .iter()
            .map(|interval| Region {
                reference: references
                    .get_index_of(names[interval.contig].as_bytes())
                    .map_or(-1, |index| index as i32),
                start: interval.start as i32,
                end: interval.end as i32,
            })
            .collect();
        // `IntervalUtil.assertOrderedNonOverlapping` on the scan,
        // `QueryInterval.assertIntervalsOptimized` on the query.
        if let Some(at) = regions
            .windows(2)
            .position(|pair| pair[0].reference > pair[1].reference)
        {
            let show = |interval: &intervals::Interval| {
                format!(
                    "{}:{}-{}",
                    names[interval.contig], interval.start, interval.end
                )
            };
            return Err(Error::Intervals {
                path: path.to_path_buf(),
                message: format!(
                    "Intervals not in order: {}; {}: the interval list orders its contigs \
                     differently from the input's header",
                    show(&unique[at]),
                    show(&unique[at + 1])
                ),
            });
        }
        Ok(Self {
            regions,
            at: 0,
            query,
        })
    }

    fn step(&mut self, record: &RawRecord) -> Step {
        if self.query {
            self.query_step(record)
        } else {
            self.scan_step(record)
        }
    }

    /// `BAMQueryMultipleIntervalsIteratorFilter.compareToFilter`, records
    /// overlapping.
    fn query_step(&mut self, record: &RawRecord) -> Step {
        let reference = record.reference_id();
        let start = record.alignment_start();
        let unmapped = record.flags() & UNMAPPED != 0;
        let end = if unmapped && start != 0 {
            start
        } else {
            alignment_end(record)
        };
        while let Some(region) = self.regions.get(self.at) {
            // An end of 0 or less means the end of the reference.
            let region_end = if region.end <= 0 {
                i32::MAX
            } else {
                region.end
            };
            if region.reference < reference || (region.reference == reference && region_end < start)
            {
                self.at += 1;
            } else if region.reference > reference || end < region.start {
                return Step::Skip;
            } else {
                return Step::Print;
            }
        }
        Step::Stop
    }

    /// `SamRecordIntervalIteratorFactory.StopAfterFilteringIterator` over
    /// `IntervalFilter`.
    fn scan_step(&mut self, record: &RawRecord) -> Step {
        let reference = record.reference_id();
        let start = record.alignment_start();
        let (stop_reference, stop_position) = self
            .regions
            .last()
            .map_or((-1, -1), |last| (last.reference, last.end));
        if reference == -1
            || reference > stop_reference
            || (reference == stop_reference && start > stop_position)
        {
            return Step::Stop;
        }
        while self.regions.get(self.at).is_some_and(|region| {
            region.reference < reference || (region.reference == reference && region.end < start)
        }) {
            self.at += 1;
        }
        match self.regions.get(self.at) {
            Some(region)
                if region.reference == reference && region.start <= alignment_end(record) =>
            {
                Step::Print
            }
            _ => Step::Skip,
        }
    }
}

/// `SAMRecord.getAlignmentEnd`: 0 for an unmapped record.
fn alignment_end(record: &RawRecord) -> i32 {
    if record.flags() & UNMAPPED != 0 {
        return 0;
    }
    record
        .alignment_start()
        .wrapping_add(record.reference_length())
        .wrapping_sub(1)
}

/// Where the text of a view goes.
pub trait Sink {
    /// The header lines, each ended by `\n`.
    fn header(&mut self, text: &[u8]) -> Result<()>;
    /// One record line with its `\n`.
    fn record(&mut self, line: &[u8]) -> Result<()>;
}

/// A [`Sink`] that writes the text as it comes.
pub struct Stream<W> {
    out: W,
    /// The name an error gives the output.
    path: PathBuf,
}

impl<W: Write> Stream<W> {
    pub fn new(out: W, path: &Path) -> Self {
        Self {
            out,
            path: path.to_path_buf(),
        }
    }

    fn write(&mut self, bytes: &[u8]) -> Result<()> {
        self.out.write_all(bytes).map_err(|source| Error::Output {
            path: self.path.clone(),
            source,
        })
    }

    pub fn finish(mut self) -> Result<W> {
        self.out.flush().map_err(|source| Error::Output {
            path: self.path.clone(),
            source,
        })?;
        Ok(self.out)
    }
}

impl<W: Write> Sink for Stream<W> {
    fn header(&mut self, text: &[u8]) -> Result<()> {
        self.write(text)
    }

    fn record(&mut self, line: &[u8]) -> Result<()> {
        self.write(line)
    }
}

/// `ViewSam.writeSamText`: the header and the records `inv` selects, into
/// `sink`.
pub fn view(inv: &Invocation, sink: &mut dyn Sink) -> Result<()> {
    let options = Options::from_invocation(inv)?;
    let stringency = inv.stringency()?;
    if options.header_only {
        let header = pass::read_header_at(&inv.input, stringency)?;
        return sink.header(&sam_text::header_text(&header, true));
    }
    let indexed = options.interval_list.is_some() && has_index(&inv.input);
    pass::with_record_stream_at(
        &inv.input,
        inv.reference(),
        inv.threads,
        stringency,
        |header, reader, cram| {
            let mut filter = options
                .interval_list
                .as_deref()
                .map(|path| IntervalFilter::load(path, header, indexed))
                .transpose()?;
            if !options.records_only {
                sink.header(&sam_text::header_text(header, true))?;
            }
            let references: Vec<Vec<u8>> = header
                .reference_sequences()
                .keys()
                .map(|name| name.to_vec())
                .collect();
            // htsjdk rebuilds `MD` and `NM` while it decodes a CRAM record.
            let mut contigs = match (cram, inv.reference()) {
                (true, Some(fasta)) => Some(ContigCache::new(Reference::open(fasta, header)?)),
                _ => None,
            };
            let mut record = RawRecord::default();
            let mut line = Vec::new();
            while record
                .read_from(reader)
                .map_err(|e| pass::read_error(&inv.input, cram, e))?
            {
                match filter.as_mut().map_or(Step::Print, |f| f.step(&record)) {
                    Step::Stop => break,
                    Step::Skip => continue,
                    Step::Print => {}
                }
                let flags = record.flags();
                if !options.alignment.takes(flags, UNMAPPED) || !options.pf.takes(flags, QC_FAIL) {
                    continue;
                }
                if let Some(contigs) = contigs.as_mut() {
                    let mut tags = record.tags().map_err(Error::Record)?;
                    restore_md_and_nm(&record, &mut tags, contigs)?;
                    record.set_tags(&tags).map_err(Error::Record)?;
                }
                line.clear();
                sam_text::write_record(&mut line, &record, &references).map_err(Error::Record)?;
                sink.record(&line)?;
            }
            Ok(())
        },
    )
}

/// A [`Sink`] that keeps the text and one row per record, up to
/// [`MAX_TEXT_BYTES`].
#[derive(Default)]
struct Held {
    text: Vec<u8>,
    rows: Vec<Vec<Value>>,
}

impl Held {
    fn push(&mut self, bytes: &[u8]) -> Result<()> {
        if self.text.len() + bytes.len() > MAX_TEXT_BYTES {
            return Err(Error::Metric(format!(
                "the SAM text is over {MAX_TEXT_BYTES} bytes, the most held as rows; the \
                 view-sam verb without --format and the Picard face print any size"
            )));
        }
        self.text.extend_from_slice(bytes);
        Ok(())
    }
}

impl Sink for Held {
    fn header(&mut self, text: &[u8]) -> Result<()> {
        self.push(text)
    }

    fn record(&mut self, line: &[u8]) -> Result<()> {
        self.push(line)?;
        let line = line.strip_suffix(b"\n").unwrap_or(line);
        let text = |field: &[u8]| String::from_utf8_lossy(field).into_owned();
        let mut fields = line.splitn(COLUMNS.len(), |&b| b == b'\t');
        let mut row = Vec::with_capacity(COLUMNS.len());
        for column in COLUMNS {
            let field = fields.next().unwrap_or_default();
            let number = matches!(column, "FLAG" | "POS" | "MAPQ" | "PNEXT" | "TLEN")
                .then(|| text(field).parse::<i64>().ok())
                .flatten();
            row.push(number.map_or_else(|| Value::Str(text(field)), Value::Int));
        }
        self.rows.push(row);
        Ok(())
    }
}

/// The registry's `run`: the text and its rows in memory.
fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let mut held = Held::default();
    view(inv, &mut held)?;
    let text = String::from_utf8_lossy(&held.text).into_owned();
    Ok(MetricsOutput::text(&COLUMNS, held.rows, text).into())
}

/// Print the view of `inv` on standard output as it is read. The status
/// is Picard's exit code: 1 when standard output is closed before the
/// last record (`PrintStream.checkError`), else 0.
pub fn print(inv: &Invocation) -> Result<u8> {
    let stdout = io::stdout().lock();
    let mut stream = Stream::new(
        io::BufWriter::with_capacity(1 << 16, stdout),
        Path::new("-"),
    );
    let closed = |e: &Error| matches!(e, Error::Output { source, .. } if source.kind() == io::ErrorKind::BrokenPipe);
    match view(inv, &mut stream).and_then(|()| stream.finish().map(drop)) {
        Ok(()) => Ok(0),
        Err(e) if closed(&e) => Ok(1),
        Err(e) => Err(e),
    }
}

/// Write the view of `inv` to the file at `path`, under a temporary name
/// until it is complete ([`AtomicFile`]).
pub fn write_file(inv: &Invocation, path: &Path) -> Result<()> {
    let output = |source| Error::Output {
        path: path.to_path_buf(),
        source,
    };
    let file = AtomicFile::new(path);
    let out = io::BufWriter::with_capacity(1 << 20, file.create().map_err(output)?);
    let mut stream = Stream::new(out, path);
    view(inv, &mut stream)?;
    drop(stream.finish()?);
    file.commit().map_err(output)
}

/// The Picard face: parse, then [`print`].
pub fn run_cli(args: &[String]) -> Result<u8> {
    let inv = TOOL.parse(args)?;
    Options::from_invocation(&inv)?;
    for note in &inv.notes {
        eprintln!("bammetrics: {note}");
    }
    print(&inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_raw::{Fields, Tags};

    fn record(flags: u16, reference_id: i32, start: i32, cigar: &[(u32, u8)]) -> RawRecord {
        RawRecord::from_fields(&Fields {
            name: b"r",
            flags,
            reference_id,
            alignment_start: start,
            mapping_quality: 0,
            cigar,
            mate_reference_id: -1,
            mate_alignment_start: 0,
            template_length: 0,
            bases: b"",
            qualities: None,
            tags: &Tags::default(),
        })
        .expect("a record")
    }

    fn filter(query: bool) -> IntervalFilter {
        let region = |reference, start, end| Region {
            reference,
            start,
            end,
        };
        IntervalFilter {
            regions: vec![region(0, 100, 200), region(1, 50, 60)],
            at: 0,
            query,
        }
    }

    #[test]
    fn the_index_query_and_the_scan_keep_overlapping_records() {
        for query in [true, false] {
            let mut filter = filter(query);
            let steps: Vec<Step> = [
                record(0, 0, 80, &[(10, 0)]),
                record(0, 0, 91, &[(10, 0)]),
                record(0, 0, 200, &[(10, 0)]),
                record(0, 0, 201, &[(10, 0)]),
                record(0, 1, 40, &[(5, 0), (20, 3), (5, 0)]),
                record(0, 1, 60, &[(1, 0)]),
            ]
            .iter()
            .map(|record| filter.step(record))
            .collect();
            assert_eq!(
                steps,
                [
                    Step::Skip,
                    Step::Print,
                    Step::Print,
                    Step::Skip,
                    Step::Print,
                    Step::Print
                ],
                "query {query}"
            );
        }
    }

    #[test]
    fn only_the_index_query_keeps_an_unmapped_record_placed_with_its_mate() {
        let placed = record(0x4, 0, 150, &[]);
        assert_eq!(filter(true).step(&placed), Step::Print);
        assert_eq!(filter(false).step(&placed), Step::Skip);
    }

    #[test]
    fn the_scan_stops_after_the_last_interval_and_at_a_record_without_a_reference() {
        let past = record(0, 1, 61, &[(10, 0)]);
        let unplaced = record(0x4, -1, 0, &[]);
        assert_eq!(filter(false).step(&past), Step::Stop);
        assert_eq!(filter(false).step(&unplaced), Step::Stop);
        // The query stops once every interval is behind the record.
        assert_eq!(filter(true).step(&past), Step::Stop);
        assert_eq!(filter(true).step(&unplaced), Step::Skip);
        let mut none = IntervalFilter {
            regions: Vec::new(),
            at: 0,
            query: false,
        };
        assert_eq!(none.step(&record(0, 0, 1, &[(1, 0)])), Step::Stop);
    }

    #[test]
    fn options_take_picards_values_in_picards_case() {
        let parse = |args: &[&str]| {
            let args: Vec<String> = std::iter::once("I=x.bam")
                .chain(args.iter().copied())
                .map(str::to_owned)
                .collect();
            TOOL.parse(&args)
                .and_then(|inv| Options::from_invocation(&inv))
        };
        let options = parse(&["ALIGNMENT_STATUS=Unaligned", "PF_STATUS=PF"]).expect("options");
        assert_eq!((options.alignment, options.pf), (Keep::Set, Keep::Clear));
        assert!(Keep::Set.takes(0x4, UNMAPPED) && !Keep::Set.takes(0, UNMAPPED));
        assert!(Keep::Clear.takes(0, QC_FAIL) && !Keep::Clear.takes(0x200, QC_FAIL));
        assert!(matches!(parse(&["PF_STATUS=pf"]), Err(Error::Usage(_))));
        assert!(matches!(
            parse(&["HEADER_ONLY=true", "RECORDS_ONLY=true"]),
            Err(Error::Usage(text)) if text.starts_with("Cannot specify both")
        ));
        assert!(matches!(parse(&["O=x.sam"]), Err(Error::Usage(_))));
        assert!(
            parse(&["INTERVAL_LIST=null"])
                .expect("null")
                .interval_list
                .is_none()
        );
    }

    #[test]
    fn held_text_has_a_limit_and_rows_have_the_twelve_columns() {
        let mut held = Held::default();
        held.header(b"@HD\tVN:1.6\n").expect("header");
        held.record(b"r1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tXA:i:1\tXB:Z:a b\n")
            .expect("record");
        assert_eq!(held.rows.len(), 1);
        assert_eq!(held.rows[0].len(), COLUMNS.len());
        assert!(matches!(held.rows[0][1], Value::Int(4)));
        assert!(matches!(&held.rows[0][11], Value::Str(tags) if tags == "XA:i:1\tXB:Z:a b"));
        held.text = vec![0; MAX_TEXT_BYTES];
        assert!(held.header(b"x").is_err());
    }
}
