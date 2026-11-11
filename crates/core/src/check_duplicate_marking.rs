//! `CheckDuplicateMarking`, mirroring Picard 3.5.0
//! (`picard.sam.markduplicates.CheckDuplicateMarking`).
//!
//! Every record of a queryname must carry the same duplicate flag (0x400).
//! The first checked record of a name sets the flag the rest are held to;
//! each later record of the same name with the other flag is one bad
//! record. Whether the marking is right is not checked.
//!
//! The header's `SO` decides the read. `SO:queryname` is read in file
//! order, and only adjacent records are one name: a name that comes back
//! later starts a new group, as Picard's single `currentReadName` does.
//! Any other sort order (coordinate, unsorted, unknown, no `@HD`) is sorted
//! first, as Picard sorts it into temp space with
//! `SAMRecordQueryNameComparator`: name as `String.compareTo` (byte order),
//! then paired before unpaired, first of pair before second, forward
//! before reverse, primary before secondary, secondary after supplementary,
//! records without `HI` before those with one and `HI` ascending, and file
//! order between records that compare equal (the sort is stable, as
//! `SortingCollection` is). bammetrics keeps the sort keys in memory and
//! writes no temporary file, so `TMP_DIR` and `MAX_RECORDS_IN_RAM` change
//! nothing. The sort reads `HI` as `SAMRecord.getIntegerAttribute` does: a
//! float is truncated, a string, character or array is an error (`Value for
//! tag HI is not Number`), a value past `Integer` is an error too. A
//! `SO:queryname` file never reads `HI`.
//!
//! htsjdk decodes a read name byte by byte (Latin-1), so a name holding
//! bytes above 127 prints each byte as one character; [`latin1`] does the
//! same for the log, `OUTPUT` and the rows.
//!
//! Picard's `ProgressLogger` lines (`Checked.` every 1,000,000 checked
//! records, `Read into sorter` every 1,000,000 records sorted) carry
//! timings and are not printed.
//!
//! `MODE` filters the records checked: `ALL` every one; `PRIMARY_ONLY`
//! drops secondary and supplementary records; `PRIMARY_MAPPED_ONLY` drops
//! unmapped records too; `PRIMARY_PROPER_PAIR_ONLY` keeps proper pairs and
//! fails at an unpaired record, as `SAMRecord.getProperPairFlag` throws
//! `Inappropriate call if not paired read`. A dropped record does not
//! change the current name.
//!
//! Picard logs `Reads with queryname NAME have different duplicate flags
//! (at CONTIG:START)` at WARNING for the first 100 bad records (`CONTIG` is
//! `null` for an unmapped record, placed or not; `START` is the alignment
//! start, 0 without one), then `Further warnings will be suppressed.`, and
//! ends with `Found N records that do not agree on their duplicate flag.`
//! at ERROR (exit 1) or `All records' duplicate markings agree.` at INFO
//! (exit 0). `OUTPUT`, when given, gets the queryname of each bad record on
//! its own line, so a name repeats when several of its records disagree,
//! and is an empty file when all agree. An `OUTPUT` that is the `INPUT`
//! file is refused before the read. [`run_cli`] prints the log lines
//! without the timestamp field and returns the exit code.

use std::cmp::Ordering;

use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::{bam, sam};

use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::pass::Accumulator;
use crate::registry::{Invocation, Opt, Tool, write_files};
use crate::sam_order::QuerynameKey;

pub const TOOL: Tool = Tool {
    name: "CheckDuplicateMarking",
    summary: "one row per record whose duplicate flag differs from the first record of its \
              queryname (read_name, contig, start); no row when every queryname agrees. OUTPUT \
              gets the queryname of each such record, one per line. The command line exits 1 \
              when any is found.",
    options: &[
        Opt::output(&["OUTPUT", "O"])
            .with_help("The bad querynames, one line per record that disagrees."),
        Opt::text(
            &["MODE"],
            "ALL, PRIMARY_ONLY, PRIMARY_MAPPED_ONLY or PRIMARY_PROPER_PAIR_ONLY: which records \
             of a queryname are checked (default ALL).",
        ),
        Opt::text(
            &["VERBOSITY"],
            "ERROR drops the INFO and WARNING log lines, WARNING the INFO lines (default INFO).",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::ignored(
            &["TMP_DIR"],
            "a file without SO:queryname is sorted in memory and no temporary file is written",
        ),
        Opt::ignored(
            &["MAX_RECORDS_IN_RAM"],
            "a file without SO:queryname is sorted in memory whatever the value",
        ),
    ],
    run: run_tool,
    collector: None,
};

/// The columns of the rows [`TOOL`] returns, one per bad record.
pub const COLUMNS: &[&str] = &["read_name", "contig", "start"];

/// Picard's `NUM_WARNINGS`: bad records logged one by one.
const NUM_WARNINGS: usize = 100;

/// A read name as htsjdk's `BAMRecord` decodes it: each byte is one
/// character, so `0xC3 0xA9` is `Ã©` and `0xFF` is `ÿ`.
pub fn latin1(bytes: &[u8]) -> String {
    bytes.iter().map(|&b| char::from(b)).collect()
}

/// Picard's `Mode`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    All,
    PrimaryOnly,
    PrimaryMappedOnly,
    PrimaryProperPairOnly,
}

impl Mode {
    /// The enum constant, or Picard's parser error.
    fn parse(text: &str) -> Result<Self> {
        match text {
            "ALL" => Ok(Self::All),
            "PRIMARY_ONLY" => Ok(Self::PrimaryOnly),
            "PRIMARY_MAPPED_ONLY" => Ok(Self::PrimaryMappedOnly),
            "PRIMARY_PROPER_PAIR_ONLY" => Ok(Self::PrimaryProperPairOnly),
            other => Err(Error::Usage(format!(
                "'{other}' is not a valid value for Mode."
            ))),
        }
    }
}

/// htsjdk log levels, lowest first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Level {
    Error,
    Warning,
    Info,
}

impl Level {
    fn name(self) -> &'static str {
        match self {
            Self::Error => "ERROR",
            Self::Warning => "WARNING",
            Self::Info => "INFO",
        }
    }
}

#[derive(Clone, Debug)]
struct Options {
    mode: Mode,
    /// The highest level printed.
    verbosity: Level,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let mode = inv.get("MODE").map_or(Ok(Mode::All), Mode::parse)?;
        let verbosity = match inv.get("VERBOSITY") {
            None | Some("INFO" | "DEBUG") => Level::Info,
            Some("WARNING") => Level::Warning,
            Some("ERROR") => Level::Error,
            Some(other) => {
                return Err(Error::Usage(format!(
                    "VERBOSITY={other}: one of ERROR, WARNING, INFO, DEBUG"
                )));
            }
        };
        Ok(Self { mode, verbosity })
    }
}

/// One record whose duplicate flag differs from the first checked record
/// of its queryname.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BadRecord {
    pub read_name: String,
    /// `SAMRecord.getContig`: `None` for an unmapped record.
    pub contig: Option<String>,
    /// `SAMRecord.getStart`: the alignment start, 0 without one.
    pub start: i64,
}

impl BadRecord {
    /// Picard's warning for this record.
    fn warning(&self) -> String {
        format!(
            "Reads with queryname {} have different duplicate flags (at {}:{})",
            self.read_name,
            self.contig.as_deref().unwrap_or("null"),
            self.start
        )
    }

    fn row(&self) -> Vec<Value> {
        vec![
            Value::Str(self.read_name.clone()),
            self.contig.clone().into(),
            Value::Int(self.start),
        ]
    }
}

/// One run's result.
#[derive(Clone, Debug)]
pub struct Check {
    pub bad: Vec<BadRecord>,
    /// True when the header was not `SO:queryname` and the records were
    /// sorted.
    pub sorted: bool,
    /// Picard's log lines in order: the level and the message, already
    /// filtered by `VERBOSITY`.
    pub log: Vec<(&'static str, String)>,
    /// The rows of [`COLUMNS`] with the `OUTPUT` text (one queryname per
    /// bad record) as the main file.
    pub output: ToolOutput,
}

impl Check {
    /// Picard's exit code: 1 when a record disagrees.
    pub fn status(&self) -> u8 {
        u8::from(!self.bad.is_empty())
    }

    /// The closing log line: `Found N records ...` or `All records' ...`.
    pub fn summary(&self) -> String {
        if self.bad.is_empty() {
            "All records' duplicate markings agree.".to_owned()
        } else {
            format!(
                "Found {} records that do not agree on their duplicate flag.",
                self.bad.len()
            )
        }
    }

    /// The summary when `VERBOSITY` shows it: an agreeing summary is
    /// `INFO`, a disagreeing one `ERROR`. `None` under `VERBOSITY=ERROR`
    /// or `WARNING` when every name agrees.
    pub fn summary_shown(&self) -> Option<String> {
        let summary = self.summary();
        self.log
            .last()
            .filter(|(_, message)| *message == summary)
            .map(|_| summary)
    }

    /// Each log line as Picard prints it, without the timestamp field:
    /// `WARNING\tCheckDuplicateMarking\t<message>`.
    pub fn print_log(&self) {
        for (level, message) in &self.log {
            eprintln!("{level}\tCheckDuplicateMarking\t{message}");
        }
    }

    fn build(bad: Vec<BadRecord>, sorted: bool, verbosity: Level) -> Self {
        let mut lines: Vec<(Level, String)> = Vec::new();
        if sorted {
            lines.push((
                Level::Info,
                "Input file isn't queryname sorted. Sorting into temp space.".to_owned(),
            ));
        }
        for record in bad.iter().take(NUM_WARNINGS) {
            lines.push((Level::Warning, record.warning()));
        }
        if bad.len() >= NUM_WARNINGS {
            lines.push((
                Level::Warning,
                "Further warnings will be suppressed.".to_owned(),
            ));
        }
        let text: String = bad.iter().map(|b| format!("{}\n", b.read_name)).collect();
        let rows = bad.iter().map(BadRecord::row).collect();
        let mut check = Self {
            bad,
            sorted,
            log: Vec::new(),
            output: ToolOutput::file(MetricsOutput::text(COLUMNS, rows, text)),
        };
        let level = if check.bad.is_empty() {
            Level::Info
        } else {
            Level::Error
        };
        lines.push((level, check.summary()));
        check.log = lines
            .into_iter()
            .filter(|(level, _)| *level <= verbosity)
            .map(|(level, message)| (level.name(), message))
            .collect();
        check
    }
}

/// The fields `SAMRecordQueryNameComparator` and the check read from a
/// record, kept when the file has to be sorted.
#[derive(Clone, Debug)]
struct Key {
    name: Vec<u8>,
    flags: sam::alignment::record::Flags,
    contig: Option<usize>,
    start: i64,
    hit_index: Option<i64>,
}

impl Key {
    /// The fields of `record`; `HI` is read only for a file that is
    /// sorted, as only `SAMRecordQueryNameComparator` asks for it.
    fn read(record: &bam::Record, sorted: bool) -> Result<Self> {
        let hit_index = match record.data().get(&Tag::HIT_INDEX) {
            Some(Ok(value)) if sorted => Some(hit_index(&value)?),
            Some(Err(e)) => return Err(e.into()),
            _ => None,
        };
        Ok(Self {
            // htsjdk keeps `*` as the name of a record without one.
            name: record.name().map_or_else(|| b"*".to_vec(), |n| n.to_vec()),
            flags: record.flags(),
            contig: record.reference_sequence_id().transpose()?,
            start: record
                .alignment_start()
                .transpose()?
                .map_or(0, |p| usize::from(p) as i64),
            hit_index,
        })
    }

    /// The fields `SAMRecordQueryNameComparator` reads.
    fn order_key(&self) -> QuerynameKey<'_> {
        QuerynameKey::new(&self.name, self.flags.bits(), self.hit_index)
    }

    /// `SAMRecordQueryNameComparator.compare` ([`QuerynameKey::compare`]).
    fn compare(&self, other: &Self) -> Ordering {
        self.order_key().compare(&other.order_key())
    }
}

/// `SAMRecord.getIntegerAttribute("HI")`: an integer as it is, a float
/// through `Number.longValue` (truncated), anything else a
/// `RuntimeException` naming the Java class; a value outside `Integer` is
/// an error too.
fn hit_index(value: &sam::alignment::record::data::field::Value<'_>) -> Result<i64> {
    use sam::alignment::record::data::field::Value;
    use sam::alignment::record::data::field::value::Array;
    let number = match value {
        Value::Float(f) => Some(*f as i64),
        other => other.as_int(),
    };
    let class = match value {
        Value::Character(_) => "java.lang.Character",
        Value::String(_) => "java.lang.String",
        Value::Hex(_) | Value::Array(Array::Int8(_) | Array::UInt8(_)) => "[B",
        Value::Array(Array::Int16(_) | Array::UInt16(_)) => "[S",
        Value::Array(Array::Int32(_) | Array::UInt32(_)) => "[I",
        Value::Array(Array::Float(_)) => "[F",
        _ => "",
    };
    match number {
        Some(n) if i64::from(i32::MIN) <= n && n <= i64::from(i32::MAX) => Ok(n),
        Some(_) => Err(Error::Metric(
            "Value for tag HI is not in Integer range".to_owned(),
        )),
        None => Err(Error::Metric(format!(
            "Value for tag HI is not Number: class {class}"
        ))),
    }
}

/// `checkDuplicateMarkingsInIterator` over the records as they come
/// (`SO:queryname`), or over the sorted keys.
struct Checker {
    mode: Mode,
    contigs: Vec<String>,
    sorted: bool,
    /// The records waiting to be sorted, when the file is not queryname
    /// sorted.
    pending: Vec<Key>,
    current_name: Option<Vec<u8>>,
    current_duplicate: bool,
    bad: Vec<BadRecord>,
}

impl Checker {
    fn new(header: &sam::Header, mode: Mode) -> Self {
        let order = header
            .header()
            .and_then(|h| h.other_fields().get(&SORT_ORDER))
            .map(ToString::to_string);
        let contigs = header
            .reference_sequences()
            .keys()
            .map(|name| name.to_string())
            .collect();
        Self {
            mode,
            contigs,
            sorted: order.as_deref() != Some("queryname"),
            pending: Vec::new(),
            current_name: None,
            current_duplicate: false,
            bad: Vec::new(),
        }
    }

    /// The MODE filter, then `checkAndTallyRecordDuplicateMarking`.
    fn check(&mut self, key: &Key) -> Result<()> {
        let flags = key.flags;
        if self.mode != Mode::All && (flags.is_secondary() || flags.is_supplementary()) {
            return Ok(());
        }
        if self.mode == Mode::PrimaryMappedOnly && flags.is_unmapped() {
            return Ok(());
        }
        if self.mode == Mode::PrimaryProperPairOnly {
            if !flags.is_segmented() {
                return Err(Error::Metric(
                    "Inappropriate call if not paired read".to_owned(),
                ));
            }
            if !flags.is_properly_segmented() {
                return Ok(());
            }
        }
        if self.current_name.as_deref() != Some(key.name.as_slice()) {
            self.current_name = Some(key.name.clone());
            self.current_duplicate = flags.is_duplicate();
            return Ok(());
        }
        if flags.is_duplicate() != self.current_duplicate {
            self.bad.push(BadRecord {
                read_name: latin1(&key.name),
                contig: if flags.is_unmapped() {
                    None
                } else {
                    Some(key.contig.map_or_else(
                        || "*".to_owned(),
                        |id| {
                            self.contigs
                                .get(id)
                                .cloned()
                                .unwrap_or_else(|| id.to_string())
                        },
                    ))
                },
                start: key.start,
            });
        }
        Ok(())
    }
}

impl Accumulator for Checker {
    type Output = (Vec<BadRecord>, bool);

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let key = Key::read(record, self.sorted)?;
        if self.sorted {
            self.pending.push(key);
            Ok(())
        } else {
            self.check(&key)
        }
    }

    fn finish(mut self) -> Result<Self::Output> {
        let mut pending = std::mem::take(&mut self.pending);
        pending.sort_by(Key::compare);
        for key in &pending {
            self.check(key)?;
        }
        Ok((self.bad, self.sorted))
    }
}

/// Parse the options and check `INPUT`. `OUTPUT` must be writable before
/// the read starts, as `IOUtil.assertFileIsWritable` requires, and must
/// not be the `INPUT` file (Picard truncates it under its own reader).
pub fn check(inv: &Invocation) -> Result<Check> {
    let options = Options::from_invocation(inv)?;
    if let Some(output) = &inv.output {
        crate::compare_metrics::check_writable(&output.to_string_lossy())?;
        if same_file(output, &inv.input) {
            return Err(Error::Usage(format!(
                "OUTPUT {} is the INPUT file; writing it would destroy it",
                output.display()
            )));
        }
    }
    let (bad, sorted) =
        inv.read_whole(&inv.input, |header| Ok(Checker::new(header, options.mode)))?;
    Ok(Check::build(bad, sorted, options.verbosity))
}

/// True when both paths name one existing file, links and `..` resolved.
fn same_file(a: &std::path::Path, b: &std::path::Path) -> bool {
    match (a.canonicalize(), b.canonicalize()) {
        (Ok(a), Ok(b)) => a == b,
        _ => false,
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    Ok(check(inv)?.output)
}

/// The Picard face: parse, check, print the log lines on standard error,
/// write `OUTPUT` when given (Picard writes nothing without it). The
/// status is Picard's exit code: 0 when every queryname agrees, 1
/// otherwise.
pub fn run_cli(args: &[String]) -> Result<u8> {
    let inv = TOOL.parse(args)?;
    for note in &inv.notes {
        eprintln!("bammetrics: {note}");
    }
    let check = check(&inv)?;
    check.print_log();
    if inv.output.is_some() {
        write_files(&inv, &check.output)?;
    }
    Ok(check.status())
}

#[cfg(test)]
mod tests {
    use super::*;
    use sam::alignment::record::Flags;

    fn key(name: &str, flags: u16, hit_index: Option<i64>) -> Key {
        Key {
            name: name.as_bytes().to_vec(),
            flags: Flags::from_bits_retain(flags),
            contig: Some(0),
            start: 1,
            hit_index,
        }
    }

    /// The order the probe of Picard 3.5.0 showed, pair by pair.
    #[test]
    fn comparator_orders_as_htsjdk_does() {
        let less = |a: Key, b: Key| {
            assert_eq!(a.compare(&b), Ordering::Less);
            assert_eq!(b.compare(&a), Ordering::Greater);
        };
        less(key("N3", 0, None), key("g1", 0, None));
        less(key("n1", 0, None), key("n10", 0, None));
        less(key("n10", 0, None), key("n2", 0, None));
        less(key("g", 99, None), key("g", 147, None));
        less(key("g", 73, None), key("g", 0, None));
        less(key("g", 0, None), key("g", 16, None));
        less(key("g", 0, None), key("g", 256, None));
        less(key("g", 0, None), key("g", 2048, None));
        less(key("g", 2048, None), key("g", 256, None));
        less(key("g", 0, Some(1)), key("g", 0, Some(2)));
        less(key("g", 0, None), key("g", 0, Some(1)));
        assert_eq!(
            key("g", 0, None).compare(&key("g", 1024, None)),
            Ordering::Equal
        );
    }

    #[test]
    fn log_lines_follow_the_verbosity() {
        let bad = vec![BadRecord {
            read_name: "r2".to_owned(),
            contig: None,
            start: 30,
        }];
        let all = Check::build(bad.clone(), true, Level::Info);
        assert_eq!(
            all.log,
            vec![
                (
                    "INFO",
                    "Input file isn't queryname sorted. Sorting into temp space.".to_owned()
                ),
                (
                    "WARNING",
                    "Reads with queryname r2 have different duplicate flags (at null:30)"
                        .to_owned()
                ),
                (
                    "ERROR",
                    "Found 1 records that do not agree on their duplicate flag.".to_owned()
                ),
            ]
        );
        assert_eq!(all.status(), 1);
        let errors = Check::build(bad, true, Level::Error);
        assert_eq!(errors.log.len(), 1);
        let agree = Check::build(Vec::new(), false, Level::Warning);
        assert!(agree.log.is_empty());
        assert_eq!(agree.status(), 0);
        assert_eq!(
            agree.output.main.as_ref().and_then(|m| m.text.as_deref()),
            Some("")
        );
    }

    #[test]
    fn the_hundredth_bad_record_ends_the_warnings() {
        let bad: Vec<BadRecord> = (0..NUM_WARNINGS + 5)
            .map(|i| BadRecord {
                read_name: format!("r{i}"),
                contig: Some("chr1".to_owned()),
                start: i as i64,
            })
            .collect();
        let check = Check::build(bad, false, Level::Info);
        let warnings: Vec<&str> = check
            .log
            .iter()
            .filter(|(level, _)| *level == "WARNING")
            .map(|(_, m)| m.as_str())
            .collect();
        assert_eq!(warnings.len(), NUM_WARNINGS + 1);
        assert_eq!(
            warnings[NUM_WARNINGS],
            "Further warnings will be suppressed."
        );
        assert!(warnings[NUM_WARNINGS - 1].contains("r99 "));
        let text = check.output.main.as_ref().and_then(|m| m.text.clone());
        assert_eq!(text.map(|t| t.lines().count()), Some(NUM_WARNINGS + 5));
    }

    #[test]
    fn names_decode_byte_by_byte_as_htsjdk_does() {
        assert_eq!(latin1(b"z1"), "z1");
        assert_eq!(latin1("é1".as_bytes()), "\u{c3}\u{a9}1");
        assert_eq!(latin1(b"\xff"), "\u{ff}");
    }

    #[test]
    fn hit_index_reads_as_get_integer_attribute_does() {
        use sam::alignment::record::data::field::Value;
        assert_eq!(hit_index(&Value::Int32(2)).ok(), Some(2));
        assert_eq!(hit_index(&Value::UInt8(7)).ok(), Some(7));
        assert_eq!(hit_index(&Value::Float(1.5)).ok(), Some(1));
        assert_eq!(
            hit_index(&Value::String(b"abc".into()))
                .unwrap_err()
                .to_string(),
            "Value for tag HI is not Number: class java.lang.String"
        );
        assert_eq!(
            hit_index(&Value::Character(b'x')).unwrap_err().to_string(),
            "Value for tag HI is not Number: class java.lang.Character"
        );
        assert_eq!(
            hit_index(&Value::UInt32(u32::MAX)).unwrap_err().to_string(),
            "Value for tag HI is not in Integer range"
        );
    }

    #[test]
    fn a_queryname_file_never_reads_hi() {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::RecordBuf;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;

        let mut buf = RecordBuf::default();
        buf.data_mut()
            .insert(Tag::HIT_INDEX, BufValue::from("abc".to_owned()));
        let header = sam::Header::default();
        let mut writer = bam::io::Writer::new(Vec::new());
        writer.write_header(&header).expect("header");
        writer
            .write_alignment_record(&header, &buf)
            .expect("encode");
        let bytes = writer.into_inner().finish().expect("finish");
        let mut reader = bam::io::Reader::new(bytes.as_slice());
        reader.read_header().expect("header");
        let mut record = bam::Record::default();
        reader.read_record(&mut record).expect("decode");
        assert!(Key::read(&record, false).is_ok());
        let err = Key::read(&record, true).unwrap_err().to_string();
        assert_eq!(
            err,
            "Value for tag HI is not Number: class java.lang.String"
        );
    }

    #[test]
    fn an_output_that_is_the_input_is_refused_before_the_read() {
        let bam = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/duplicate_marking.bam"
        );
        let args: Vec<String> = [format!("I={bam}"), format!("O={bam}")].to_vec();
        let inv = TOOL.parse(&args).expect("parses");
        let err = check(&inv).unwrap_err();
        assert!(matches!(err, Error::Usage(_)), "{err}");
        assert!(
            err.to_string()
                .ends_with("is the INPUT file; writing it would destroy it")
        );
        assert_eq!(std::fs::metadata(bam).expect("fixture").len(), 363);
    }

    #[test]
    fn the_verb_summary_follows_the_verbosity() {
        let agree = Check::build(Vec::new(), false, Level::Warning);
        assert_eq!(agree.summary_shown(), None);
        let agree = Check::build(Vec::new(), false, Level::Info);
        assert_eq!(
            agree.summary_shown().as_deref(),
            Some("All records' duplicate markings agree.")
        );
        let bad = vec![BadRecord {
            read_name: "r2".to_owned(),
            contig: None,
            start: 30,
        }];
        let disagree = Check::build(bad, false, Level::Error);
        assert_eq!(
            disagree.summary_shown().as_deref(),
            Some("Found 1 records that do not agree on their duplicate flag.")
        );
    }

    #[test]
    fn mode_parses_picard_constants_only() {
        assert_eq!(Mode::parse("PRIMARY_ONLY").ok(), Some(Mode::PrimaryOnly));
        let error = Mode::parse("primary_only").unwrap_err().to_string();
        assert_eq!(error, "'primary_only' is not a valid value for Mode.");
    }
}
