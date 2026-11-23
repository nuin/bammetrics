//! `CompareSAMs`, mirroring Picard 3.5.0 (`picard.sam.CompareSAMs` over
//! `picard.sam.util.SamComparison`).
//!
//! Two SAM, BAM or CRAM files, the left and the right, are compared: first
//! their headers, then their primary records. Picard prints each header
//! difference on standard output as three lines (`X differs.`, `File 1:
//! A`, `File 2: B`), writes one `SamComparisonMetric` row to `OUTPUT` when
//! given, prints `SAM files match.` (exit 0) or `SAM files differ.` (exit
//! 1). [`run_cli`] is that face. Picard takes the two files as positional
//! arguments; the registry row holds them as `INPUT` given twice, which the
//! verb and Python use, and [`picard_arguments`] turns Picard's form into.
//!
//! # Header
//!
//! Compared in this order (`compare_sams/header.rs`): the `@HD` version
//! (`1.6` for a header without `@HD`), `CR`, `SO`; the sequence
//! dictionaries (the count, then per index the name, the length and,
//! unless `LENIENT_HEADER`, `SP`, `AS`, `M5` and `UR`); the read groups (the
//! count, then per index `ID`, `SM`, `LB`, `DS`, `PU`, `PI`, `CN`, `DT`,
//! `PL`); and, unless `LENIENT_HEADER`, the program records (the count,
//! then per index `ID`, `VN`, `CL`). `GO`, `@CO` lines and every other
//! field are not compared.
//!
//! # Records
//!
//! Secondary and supplementary records are skipped. Two records are one
//! read when their names are equal and both are unpaired, or both first, or
//! both second of their pair (`PrimaryAlignmentKey`). The walk depends on
//! the sort order, which both headers must share (`Sort Order differs.`
//! and `Cannot compare alignments if sort orders differ.` otherwise, with
//! every count 0):
//!
//! - `coordinate`: refused when the dictionaries differ. The left records
//!   of one place (reference index and start; a record without a reference
//!   has index -1) are paired by key with the right records of that place.
//!   A record without a partner at its place waits in a map, and the
//!   waiting records of both files are paired by key at the end. A key
//!   given twice at one place keeps its last record.
//! - `queryname`: a merge by key, the name in `String.compareTo` order.
//!   When the right file ends before the left, Picard 3.5.0 throws; so
//!   does this.
//! - `unsorted`, `duplicate` and no `SO`: every left record is held by
//!   key, the last record of a key given twice, then each right record
//!   looks its partner up.
//! - `unknown`, or an `SO` value htsjdk does not know: Picard throws
//!   `Unrecognized sort order (unknown) found.`
//!
//! A pair counts as `UNMAPPED_BOTH`, `UNMAPPED_LEFT` or `UNMAPPED_RIGHT` by
//! the unmapped flags, else as `MAPPINGS_MATCH` when the reference names
//! and the starts are equal, else as `MAPPINGS_DIFFER`. The strand is never
//! compared: Picard's `alignmentsMatch` compares the left record's strand
//! with itself (docs/04 Q-018, kept; KNOWN_DIFFERENCES.md). With
//! `LENIENT_LOW_MQ_ALIGNMENT` a pair whose mapping qualities are both at
//! most `LOW_MQ_THRESHOLD` matches wherever it maps (Picard's help says
//! "below"; the code tests `<=`). With `LENIENT_UNKNOWN_MQ_ALIGNMENT` a
//! pair whose qualities are both 255 matches. A record without a partner
//! counts in `MISSING_RIGHT` (a left record) or `MISSING_LEFT`.
//!
//! A pair whose duplicate flags differ counts in
//! `DUPLICATE_MARKINGS_DIFFER`, or under `LENIENT_DUP` goes through
//! `compare_sams/duplicates.rs`. `COMPARE_MQ` adds a histogram of pairs by
//! `LEFT,RIGHT` mapping quality, keyed by that text.
//!
//! `ARE_EQUAL` is true when the headers are equal, no record is missing,
//! differs or is unmapped on one side, and `DUPLICATE_MARKINGS_DIFFER` is
//! 0.
//!
//! # Memory
//!
//! The queryname walk holds two records. The other two hold waiting
//! records, as Picard's maps do, and the worst case is every primary
//! record of both files (coordinate: two files that share no read at any
//! place) or of the left file (unsorted). An entry is 48 bytes and the
//! read name, where Picard holds the whole `SAMRecord`; under
//! `LENIENT_DUP` it also holds the record cut to what the duplicate order
//! reads (`records::slim`). Measured on 8.4 million unsorted records: 233
//! bytes a record at the peak, 418 under `LENIENT_DUP`
//! (docs/benchmarks.md). Two files in the same order hold the records of
//! one place. The `LENIENT_DUP` sorts hold `MAX_RECORDS_IN_RAM` keys each
//! and spill to `TMP_DIR`.
//!
//! A BAM or a CRAM is not validated (`VALIDATION_STRINGENCY` applies to
//! SAM text), as for every tool. One check htsjdk makes at any stringency
//! is kept: a reference or mate reference index outside the dictionary
//! fails the read, on a secondary or supplementary record too.

mod alignments;
mod duplicates;
mod header;
mod records;
#[cfg(test)]
mod tests;

use std::io::Write;
use std::path::PathBuf;

use crate::compare_metrics::JavaFile;
use crate::error::{Error, Result};
use crate::metric::{
    HistogramSection, HistogramSeries, MetricRow, MetricsOutput, ToolOutput, Value,
};
use crate::metrics_file::KeyClass;
use crate::pass;
use crate::registry::{self, Invocation, Opt, Tool, Ty, write_files};
use crate::spill::TempDir;

use alignments::{Comparer, Input};
use duplicates::Lenient;
use header::Say;

/// Picard's `MAX_RECORDS_IN_RAM` default, the size of a sorted run.
const DEFAULT_RECORDS_IN_RAM: i64 = 500_000;

/// A `LENIENT_DUP` run also spills at this many bytes of keys.
const SORT_BYTES: usize = 256 << 20;

/// The histogram labels Picard writes under `COMPARE_MQ`.
pub const KEY_LABEL: &str = "BIN";
pub const SERIES_LABEL: &str = "VALUE";

pub const TOOL: Tool = Tool {
    name: "CompareSAMs",
    summary: "one SamComparisonMetric row for two SAM, BAM or CRAM files: the pairs of primary \
              records that map to the same place, differ, are unmapped or missing on a side, \
              the duplicate flags that differ, and ARE_EQUAL. The command line prints each \
              header difference and exits 1 when the files differ.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The two files, the left then the right (Picard's two positional arguments).",
        )
        .repeated()
        .required(),
        Opt::output(&["OUTPUT", "O"]).with_help("The SamComparisonMetric metrics file."),
        Opt::bool(
            &["LENIENT_HEADER"],
            "Do not compare SP, AS, M5 and UR of the @SQ lines nor the @PG lines (default \
             false).",
        ),
        Opt::bool(
            &["LENIENT_DUP"],
            "Duplicate flags may differ where the two files chose another representative of \
             one duplicate set (default false).",
        ),
        Opt::bool(
            &["LENIENT_LOW_MQ_ALIGNMENT"],
            "A pair whose mapping qualities are both at most LOW_MQ_THRESHOLD matches wherever \
             it maps (default false).",
        ),
        Opt::bool(
            &["LENIENT_UNKNOWN_MQ_ALIGNMENT"],
            "A pair whose mapping qualities are both 255 matches wherever it maps (default \
             false).",
        ),
        Opt::int(
            &["LOW_MQ_THRESHOLD"],
            "The mapping quality LENIENT_LOW_MQ_ALIGNMENT tests against (default 3).",
        ),
        Opt::bool(
            &["COMPARE_MQ"],
            "Add a histogram of pairs by LEFT,RIGHT mapping quality (default false).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the sorted runs of LENIENT_DUP (default: the system temporary \
             directory).",
        )
        .repeated(),
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records each LENIENT_DUP sort holds in memory before a run spills to TMP_DIR \
             (default 500000).",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
        registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// The options of one comparison.
#[derive(Clone, Debug)]
pub struct Options {
    left: PathBuf,
    right: PathBuf,
    lenient_header: bool,
    lenient_dup: bool,
    lenient_low_mq: bool,
    lenient_unknown_mq: bool,
    low_mq_threshold: i32,
    compare_mq: bool,
    records_in_ram: usize,
    tmp_dir: Option<PathBuf>,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let [left, right] = inv.all("INPUT") else {
            return Err(Error::Usage(format!(
                "CompareSAMs: exactly two input files are compared, got {}",
                inv.all("INPUT").len()
            )));
        };
        let lenient_dup = inv.flag("LENIENT_DUP", false)?;
        let records: i64 = inv.parse("MAX_RECORDS_IN_RAM", DEFAULT_RECORDS_IN_RAM)?;
        if lenient_dup && records < 1 {
            return Err(Error::Usage(format!(
                "MAX_RECORDS_IN_RAM={records}: LENIENT_DUP sorts need 1 or more"
            )));
        }
        Ok(Self {
            left: PathBuf::from(left),
            right: PathBuf::from(right),
            lenient_header: inv.flag("LENIENT_HEADER", false)?,
            lenient_dup,
            lenient_low_mq: inv.flag("LENIENT_LOW_MQ_ALIGNMENT", false)?,
            lenient_unknown_mq: inv.flag("LENIENT_UNKNOWN_MQ_ALIGNMENT", false)?,
            low_mq_threshold: inv.parse("LOW_MQ_THRESHOLD", 3)?,
            compare_mq: inv.flag("COMPARE_MQ", false)?,
            records_in_ram: usize::try_from(records.max(1)).unwrap_or(usize::MAX),
            tmp_dir: inv.all("TMP_DIR").first().map(PathBuf::from),
        })
    }
}

/// The counts of `SamComparisonMetric`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub mappings_match: i64,
    pub mappings_differ: i64,
    pub unmapped_both: i64,
    pub unmapped_left: i64,
    pub unmapped_right: i64,
    pub missing_left: i64,
    pub missing_right: i64,
    pub duplicate_markings_differ: i64,
}

impl Counts {
    /// `allVisitedAlignmentsEqual`.
    fn all_visited_alignments_equal(&self) -> bool {
        self.missing_left == 0
            && self.missing_right == 0
            && self.mappings_differ == 0
            && self.unmapped_left == 0
            && self.unmapped_right == 0
    }
}

/// `SamComparisonMetric`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    /// `File.getAbsolutePath` of each input.
    pub left_file: String,
    pub right_file: String,
    pub counts: Counts,
    pub are_equal: bool,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.sam.SamComparisonMetric";
    const COLUMNS: &'static [&'static str] = &[
        "LEFT_FILE",
        "RIGHT_FILE",
        "MAPPINGS_MATCH",
        "MAPPINGS_DIFFER",
        "UNMAPPED_BOTH",
        "UNMAPPED_LEFT",
        "UNMAPPED_RIGHT",
        "MISSING_LEFT",
        "MISSING_RIGHT",
        "DUPLICATE_MARKINGS_DIFFER",
        "ARE_EQUAL",
    ];

    fn values(&self) -> Vec<Value> {
        let counts = &self.counts;
        vec![
            Value::Str(self.left_file.clone()),
            Value::Str(self.right_file.clone()),
            counts.mappings_match.into(),
            counts.mappings_differ.into(),
            counts.unmapped_both.into(),
            counts.unmapped_left.into(),
            counts.unmapped_right.into(),
            counts.missing_left.into(),
            counts.missing_right.into(),
            counts.duplicate_markings_differ.into(),
            Value::Bool(self.are_equal),
        ]
    }
}

/// One base-11 digit per decimal digit of a mapping quality, three in all:
/// 0 for no digit, else the digit plus 1. A shorter number is a prefix, so
/// the codes order as the decimal texts order.
fn mq_code(quality: u8) -> i64 {
    let text = quality.to_string();
    (0..3).fold(0, |code, at| {
        let digit = text
            .as_bytes()
            .get(at)
            .map_or(0, |d| i64::from(d - b'0') + 1);
        code * 11 + digit
    })
}

/// The histogram key of the text `LEFT,RIGHT` ([`KeyClass::MqPair`]). Two
/// keys compare as their texts compare in Java: `,` sorts before every
/// digit, as the missing digit does here.
pub fn mq_pair_key(left: u8, right: u8) -> i64 {
    mq_code(left) * 1331 + mq_code(right)
}

/// The text of a [`mq_pair_key`]: `LEFT,RIGHT`.
pub fn mq_pair_text(key: i64) -> String {
    let number = |code: i64| -> String {
        [code / 121, code / 11 % 11, code % 11]
            .into_iter()
            .filter(|digit| (1..=10).contains(digit))
            .map(|digit| char::from(b'0' + (digit - 1) as u8))
            .collect()
    };
    format!(
        "{},{}",
        number(key.div_euclid(1331)),
        number(key.rem_euclid(1331))
    )
}

/// One comparison's result.
#[derive(Clone, Debug)]
pub struct Comparison {
    pub row: Row,
    /// The metrics file: the row and, under `COMPARE_MQ`, the histogram.
    pub output: ToolOutput,
    /// The most records the waiting maps held at once.
    pub peak_pending: usize,
}

impl Comparison {
    /// Picard's exit code: 0 when the files are equal, else 1.
    pub fn status(&self) -> u8 {
        u8::from(!self.row.are_equal)
    }

    /// The last line Picard prints.
    pub fn verdict(&self) -> &'static str {
        if self.row.are_equal {
            "SAM files match."
        } else {
            "SAM files differ."
        }
    }
}

/// Both inputs open on their headers.
struct Opened<'r> {
    left: Input<'r>,
    right: Input<'r>,
}

/// The comparison of two opened inputs: the counts, the histogram, the
/// verdict.
fn compare_opened(
    options: &Options,
    inv: &Invocation,
    opened: Opened<'_>,
    say: Say<'_>,
) -> Result<Comparison> {
    let Opened { left, right } = opened;
    let side = |input: &Input<'_>, header| header::Side {
        header,
        one_byte_chars: !input.cram && !crate::cram::is_sam(input.path),
    };
    let (left_side, right_side) = (side(&left, left.header), side(&right, right.header));
    let headers = header::compare(&left_side, &right_side, options.lenient_header, say)?;
    let orders = (left_side.sort_order(), right_side.sort_order());
    let temp = TempDir::new(options.tmp_dir.as_deref(), "compare-sams");
    let lenient = options.lenient_dup.then(|| {
        let bytes = inv.sort_memory.unwrap_or(SORT_BYTES);
        Lenient::new(
            &temp,
            options.records_in_ram,
            bytes,
            left.header,
            right.header,
        )
    });
    let mut comparer = Comparer::new(options, lenient, left.header, right.header);
    let alignments =
        comparer.alignments(orders, headers.dictionaries_differ, (left, right), say)?;
    let mut counts = comparer.counts;
    if let Some(lenient) = comparer.lenient.take() {
        counts.duplicate_markings_differ += lenient.count()?;
    }
    let row = Row {
        left_file: JavaFile::new(&options.left.to_string_lossy()).absolute,
        right_file: JavaFile::new(&options.right.to_string_lossy()).absolute,
        counts,
        are_equal: headers.equal && alignments && counts.duplicate_markings_differ == 0,
    };
    let histogram = options.compare_mq.then(|| HistogramSection {
        key_label: KEY_LABEL,
        key_class: KeyClass::MqPair,
        series: vec![HistogramSeries {
            label: SERIES_LABEL.to_owned(),
            histogram: comparer.mapping_qualities,
        }],
    });
    Ok(Comparison {
        output: MetricsOutput::from_rows(std::slice::from_ref(&row), histogram).into(),
        row,
        peak_pending: comparer.peak_pending,
    })
}

/// Compare the two `INPUT` files of `inv`. Each line Picard prints on
/// standard output before its verdict goes to `say` as it is found.
pub fn compare(inv: &Invocation, say: Say<'_>) -> Result<Comparison> {
    let options = Options::from_invocation(inv)?;
    if let Some(output) = &inv.output {
        registry::check_writable([output.as_path()], &inv.inputs)?;
    }
    let stringency = inv.stringency()?;
    let (reference, threads) = (inv.reference(), inv.threads);
    pass::with_record_stream_at(
        &options.left,
        reference,
        threads,
        stringency,
        |left_header, left_reader, left_cram| {
            pass::with_record_stream_at(
                &options.right,
                reference,
                threads,
                stringency,
                |right_header, right_reader, right_cram| {
                    let opened = Opened {
                        left: Input {
                            path: &options.left,
                            reader: left_reader,
                            cram: left_cram,
                            header: left_header,
                        },
                        right: Input {
                            path: &options.right,
                            reader: right_reader,
                            cram: right_cram,
                            header: right_header,
                        },
                    };
                    compare_opened(&options, inv, opened, say)
                },
            )
        },
    )
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    Ok(compare(inv, &mut |_| {})?.output)
}

/// Picard's command line as the registry reads it: each positional
/// argument (a word without `=` that is no option's value) becomes
/// `INPUT=WORD`. Picard has no `INPUT` option on this tool, so one given by
/// name is refused as Picard refuses it.
pub fn picard_arguments(args: &[String]) -> Result<Vec<String>> {
    let is_bool = |name: &str| {
        TOOL.table()
            .iter()
            .any(|opt| opt.names.contains(&name) && opt.ty == Ty::Bool)
    };
    let mut out = Vec::with_capacity(args.len());
    let mut rest = args.iter().peekable();
    while let Some(arg) = rest.next() {
        let body = arg.trim_start_matches('-');
        let name = body.split_once('=').map_or(body, |(name, _)| name);
        let option = arg.starts_with('-') || body.contains('=');
        if option && matches!(name, "INPUT" | "I") {
            return Err(Error::Usage(format!(
                "CompareSAMs: Unrecognized option: {name} (the two files are positional \
                 arguments)"
            )));
        }
        if !option {
            out.push(format!("INPUT={arg}"));
            continue;
        }
        out.push(arg.clone());
        let bare_bool = is_bool(body) && rest.peek().is_none_or(|next| next.starts_with('-'));
        if arg.starts_with('-') && !body.contains('=') && !bare_bool {
            out.extend(rest.next().cloned());
        }
    }
    Ok(out)
}

/// The Picard face: the header differences and the verdict on standard
/// output, `OUTPUT` when given (Picard writes no file without it), and
/// Picard's exit code. A closed standard output ends the printing, not the
/// run.
pub fn run_cli(args: &[String]) -> Result<u8> {
    let inv = TOOL.parse(&picard_arguments(args)?)?;
    for note in &inv.notes {
        eprintln!("bammetrics: {note}");
    }
    let mut print = |line: &str| {
        let _ = writeln!(std::io::stdout().lock(), "{line}");
    };
    let comparison = compare(&inv, &mut print)?;
    if inv.output.is_some() {
        write_files(&inv, &comparison.output)?;
    }
    print(comparison.verdict());
    Ok(comparison.status())
}
