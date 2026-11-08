//! `BedToIntervalList`, mirroring Picard 3.5.0 `BedToIntervalList.doWork`
//! with htsjdk's `BEDCodec`, `IntervalList.sorted`, `IntervalList.uniqued`
//! and `IntervalListWriter`.
//!
//! The BED is read line by line: blank lines and lines starting with `#`,
//! `track` or `browser` are skipped; fields are split on a tab or a run of
//! spaces; a line with one field is skipped; the start is the second field
//! plus one; the end is the third field, or the start without one; the
//! name is the fourth field with its double quotes removed (empty is no
//! name); the fifth field must parse as a float for the sixth to be read
//! as the strand, where only a leading `-` is negative. Each interval is
//! checked against the dictionary the way Picard checks it, and a
//! zero-length one is dropped unless `KEEP_LENGTH_ZERO_INTERVALS`.
//!
//! The dictionary comes from `SEQUENCE_DICTIONARY`: a `.dict`, a SAM, BAM
//! or CRAM header (a CRAM told by its magic), an interval list, or a FASTA with a `.dict` next to it (its
//! `.fai` when there is none, which then yields `SN` and `LN` only). The
//! file starts with `@HD VN:1.6 SO:coordinate` and the `@SQ` lines with
//! their fields in the order the dictionary gave them, then one line per
//! interval: contig, start, end, `+` or `-`, name or `.`.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use noodles::{bam, bgzf};

use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};

/// The columns of the rows returned to Python, one per interval written.
pub const COLUMNS: [&str; 5] = ["contig", "start", "end", "strand", "name"];

/// `SAMFileHeader.CURRENT_VERSION` in htsjdk 5.0.0.
const HEADER_VERSION: &str = "1.6";

pub const TOOL: Tool = Tool {
    name: "BedToIntervalList",
    summary: "one row per interval written (contig, start, end, strand, name); OUTPUT gets \
              the interval list.",
    options: &[
        Opt::path(&["INPUT", "I"], "The BED file to convert.").required(),
        crate::registry::REFERENCE_NOT_DECODED,
        Opt::path(
            &["SEQUENCE_DICTIONARY", "SD"],
            "The @SQ lines: a .dict, BAM, SAM, interval list, or FASTA with a .dict.",
        )
        .required(),
        Opt::bool(&["SORT"], "Sort the intervals (default true)."),
        Opt::bool(&["UNIQUE"], "Merge identical intervals (default false)."),
        Opt::bool(
            &["DROP_MISSING_CONTIGS"],
            "Drop intervals on contigs the dictionary lacks (default false).",
        ),
        Opt::bool(
            &["KEEP_LENGTH_ZERO_INTERVALS"],
            "Keep intervals of length zero (default false).",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub sequence_dictionary: PathBuf,
    pub sort: bool,
    pub unique: bool,
    pub drop_missing_contigs: bool,
    pub keep_length_zero_intervals: bool,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let sequence_dictionary = inv.get("SEQUENCE_DICTIONARY").ok_or_else(|| {
            Error::Usage("BedToIntervalList: SEQUENCE_DICTIONARY (SD) is required".to_owned())
        })?;
        Ok(Self {
            sequence_dictionary: sequence_dictionary.into(),
            sort: inv.flag("SORT", true)?,
            unique: inv.flag("UNIQUE", false)?,
            drop_missing_contigs: inv.flag("DROP_MISSING_CONTIGS", false)?,
            keep_length_zero_intervals: inv.flag("KEEP_LENGTH_ZERO_INTERVALS", false)?,
        })
    }
}

/// One `@SQ` line: name, length, and the other fields as written.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub length: i64,
    pub other_fields: Vec<String>,
}

impl Sequence {
    /// `SAMTextHeaderCodec.getSQLine`: `SN`, `LN`, then the other fields in
    /// the order they were read.
    fn line(&self) -> String {
        let mut fields = vec![
            "@SQ".to_owned(),
            format!("SN:{}", self.name),
            format!("LN:{}", self.length),
        ];
        fields.extend(self.other_fields.iter().cloned());
        fields.join("\t")
    }
}

/// `SAMSequenceDictionaryExtractor.extractDictionary` for the file types
/// a metrics harness meets.
pub fn read_dictionary(path: &Path) -> Result<Vec<Sequence>> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let name = path.to_string_lossy().to_ascii_lowercase();
    let text = if crate::cram::is_cram(path) {
        let mut reader =
            noodles::cram::io::Reader::new(BufReader::new(File::open(path).map_err(input)?));
        let mut header = reader.header_reader();
        header.read_magic_number().map_err(input)?;
        header.read_format_version().map_err(input)?;
        header.read_file_id().map_err(input)?;
        let mut container = header.container_reader().map_err(input)?;
        let mut text = String::new();
        container
            .raw_sam_header_reader()
            .map_err(input)?
            .read_to_string(&mut text)
            .map_err(input)?;
        text
    } else if name.ends_with(".bam") {
        let mut reader = bam::io::Reader::new(File::open(path).map_err(input)?);
        let mut header = reader.header_reader();
        header.read_magic_number().map_err(input)?;
        let mut raw = header.raw_sam_header_reader().map_err(input)?;
        let mut text = String::new();
        raw.read_to_string(&mut text).map_err(input)?;
        text
    } else if [".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz"]
        .iter()
        .any(|e| name.ends_with(e))
    {
        return fasta_dictionary(path);
    } else {
        std::fs::read_to_string(path).map_err(input)?
    };
    let sequences = parse_sq_lines(&text).map_err(|message| Error::Intervals {
        path: path.to_path_buf(),
        message,
    })?;
    if sequences.is_empty() {
        return Err(Error::Intervals {
            path: path.to_path_buf(),
            message: "no @SQ line: cannot extract a sequence dictionary".to_owned(),
        });
    }
    Ok(sequences)
}

/// `ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence`:
/// the `.dict` with the FASTA's extension replaced, else the `.fai`.
fn fasta_dictionary(path: &Path) -> Result<Vec<Sequence>> {
    let name = path.to_string_lossy();
    let stem = [".fa.gz", ".fasta.gz", ".fa", ".fasta", ".fna"]
        .iter()
        .find_map(|e| name.strip_suffix(e))
        .unwrap_or(&name)
        .to_owned();
    let dict = PathBuf::from(format!("{stem}.dict"));
    if dict.exists() {
        return read_dictionary(&dict);
    }
    let fai = PathBuf::from(format!("{name}.fai"));
    let text = std::fs::read_to_string(&fai).map_err(|source| Error::Input {
        path: fai.clone(),
        source,
    })?;
    let mut sequences = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let mut fields = line.split('\t');
        let (Some(name), Some(length)) = (fields.next(), fields.next()) else {
            return Err(Error::Intervals {
                path: fai,
                message: format!("bad .fai line: {line}"),
            });
        };
        sequences.push(Sequence {
            name: name.to_owned(),
            length: length.parse().map_err(|_| Error::Intervals {
                path: fai.clone(),
                message: format!("bad .fai length: {line}"),
            })?,
            other_fields: Vec::new(),
        });
    }
    Ok(sequences)
}

/// The `@SQ` lines of a SAM-style header text.
fn parse_sq_lines(text: &str) -> std::result::Result<Vec<Sequence>, String> {
    let mut sequences = Vec::new();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("@SQ\t") else {
            continue;
        };
        let (mut name, mut length, mut other) = (None, None, Vec::new());
        for field in rest.split('\t') {
            if let Some(v) = field.strip_prefix("SN:") {
                name = Some(v.to_owned());
            } else if let Some(v) = field.strip_prefix("LN:") {
                length = Some(
                    v.parse::<i64>()
                        .map_err(|_| format!("bad LN in @SQ line: {line}"))?,
                );
            } else {
                other.push(field.to_owned());
            }
        }
        match (name, length) {
            (Some(name), Some(length)) => sequences.push(Sequence {
                name,
                length,
                other_fields: other,
            }),
            _ => return Err(format!("@SQ line without SN and LN: {line}")),
        }
    }
    Ok(sequences)
}

/// An interval as `BedToIntervalList` builds it: 1-based closed, with the
/// dictionary index of its contig for sorting.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Interval {
    pub contig: usize,
    pub start: i64,
    pub end: i64,
    pub negative_strand: bool,
    pub name: Option<String>,
}

/// `BEDCodec.SPLIT_PATTERN`: a tab, or a run of spaces, ends a field.
fn bed_fields(line: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let bytes = line.as_bytes();
    let (mut start, mut i) = (0, 0);
    while i < bytes.len() {
        match bytes[i] {
            b'\t' => {
                fields.push(&line[start..i]);
                i += 1;
                start = i;
            }
            b' ' => {
                fields.push(&line[start..i]);
                while i < bytes.len() && bytes[i] == b' ' {
                    i += 1;
                }
                start = i;
            }
            _ => i += 1,
        }
    }
    fields.push(&line[start..]);
    fields
}

/// A decoded BED line: contig, 1-based start, end, name (`None` when
/// empty), negative strand.
type Feature = (String, i64, i64, Option<String>, bool);

/// `BEDCodec.decode`: `None` for a line that yields no feature.
fn decode_bed(line: &str) -> Result<Option<Feature>> {
    if line.trim().is_empty()
        || line.starts_with('#')
        || line.starts_with("track")
        || line.starts_with("browser")
    {
        return Ok(None);
    }
    let fields = bed_fields(line);
    if fields.len() < 2 {
        return Ok(None);
    }
    let number = |field: &str| -> Result<i64> {
        field.parse::<i32>().map(i64::from).map_err(|_| {
            Error::Metric(format!("For input string: \"{field}\" in BED line: {line}"))
        })
    };
    let start = number(fields[1])? + 1;
    let end = match fields.get(2) {
        Some(f) => number(f)?,
        None => start,
    };
    let name = fields
        .get(3)
        .map(|n| n.replace('"', ""))
        .filter(|n| !n.is_empty());
    let mut negative = false;
    if let Some(score) = fields.get(4) {
        // A score that is not a float ends the line's parsing there.
        if score.trim().parse::<f32>().is_ok()
            && let Some(strand) = fields.get(5)
        {
            negative = strand.trim().starts_with('-');
        }
    }
    Ok(Some((fields[0].to_owned(), start, end, name, negative)))
}

/// `BedToIntervalList.doWork` up to the sort: every feature checked against
/// the dictionary. A feature Picard rejects is an [`Error::Metric`], an
/// error in the data (exit 1 on every face), never a usage error.
pub fn read_bed(text: &str, dictionary: &[Sequence], opts: &Options) -> Result<Vec<Interval>> {
    let mut intervals = Vec::new();
    for line in text.lines() {
        let Some((contig, start, end, name, negative_strand)) = decode_bed(line)? else {
            continue;
        };
        let Some(index) = dictionary.iter().position(|s| s.name == contig) else {
            if opts.drop_missing_contigs {
                continue;
            }
            return Err(Error::Metric(format!(
                "Sequence '{contig}' was not found in the sequence dictionary"
            )));
        };
        let length = dictionary[index].length;
        if start < 1 {
            return Err(Error::Metric(format!(
                "Start on sequence '{contig}' was less than one: {start}"
            )));
        } else if length < start {
            return Err(Error::Metric(format!(
                "Start on sequence '{contig}' was past the end: {length} < {start}"
            )));
        } else if (end == 0 && start != 1) || end < 0 {
            return Err(Error::Metric(format!(
                "End on sequence '{contig}' was less than one: {end}"
            )));
        } else if length < end {
            return Err(Error::Metric(format!(
                "End on sequence '{contig}' was past the end: {length} < {end}"
            )));
        } else if end < start - 1 {
            return Err(Error::Metric(format!(
                "On sequence '{contig}', end < start-1: {end} <= {start}"
            )));
        }
        if start == end + 1 && !opts.keep_length_zero_intervals {
            continue;
        }
        intervals.push(Interval {
            contig: index,
            start,
            end,
            negative_strand,
            name,
        });
    }
    Ok(intervals)
}

/// `IntervalCoordinateComparator`: contig, start, end, `+` before `-`, then
/// a missing name before any name, names by code unit.
pub fn sort(intervals: &mut [Interval]) {
    intervals.sort_by(|a, b| {
        (a.contig, a.start, a.end, a.negative_strand, &a.name).cmp(&(
            b.contig,
            b.start,
            b.end,
            b.negative_strand,
            &b.name,
        ))
    });
}

/// `IntervalList.uniqued()`: sort, then merge runs of intervals that
/// overlap or abut (`IntervalMergerIterator` with `combineAbuttingIntervals`,
/// without `enforceSameStrands`). A merged interval keeps the strand of its
/// first member and joins the distinct names of its members with `|`, in
/// order, skipping missing ones (`IntervalList.merge`).
pub fn uniqued(intervals: &[Interval]) -> Vec<Interval> {
    let mut sorted = intervals.to_vec();
    sort(&mut sorted);
    let mut merged: Vec<(Interval, Vec<String>)> = Vec::new();
    for next in sorted {
        match merged.last_mut() {
            Some((current, names))
                if current.contig == next.contig && next.start <= current.end + 1 =>
            {
                current.end = current.end.max(next.end);
                if let Some(name) = next.name
                    && !names.contains(&name)
                {
                    names.push(name);
                }
            }
            _ => {
                let names = next.name.iter().cloned().collect();
                merged.push((next, names));
            }
        }
    }
    merged
        .into_iter()
        .map(|(mut interval, names)| {
            interval.name = (!names.is_empty()).then(|| names.join("|"));
            interval
        })
        .collect()
}

/// `IntervalListWriter`: the header, then one line per interval.
pub fn write_interval_list(dictionary: &[Sequence], intervals: &[Interval]) -> String {
    let mut text = format!("@HD\tVN:{HEADER_VERSION}\tSO:coordinate\n");
    for sequence in dictionary {
        text.push_str(&sequence.line());
        text.push('\n');
    }
    for interval in intervals {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            dictionary[interval.contig].name,
            interval.start,
            interval.end,
            if interval.negative_strand { '-' } else { '+' },
            interval.name.as_deref().unwrap_or(".")
        ));
    }
    text
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let dictionary = read_dictionary(&opts.sequence_dictionary)?;
    let bed = read_bed_file(&inv.input)?;
    let mut intervals = read_bed(&bed, &dictionary, &opts)?;
    if opts.sort {
        sort(&mut intervals);
    }
    if opts.unique {
        intervals = uniqued(&intervals);
    }
    let text = write_interval_list(&dictionary, &intervals);
    let rows = intervals
        .iter()
        .map(|i| {
            vec![
                Value::Str(dictionary[i.contig].name.clone()),
                Value::Int(i.start),
                Value::Int(i.end),
                Value::Str(if i.negative_strand { "-" } else { "+" }.to_owned()),
                i.name.clone().into(),
            ]
        })
        .collect();
    Ok(MetricsOutput::text(&COLUMNS, rows, text).into())
}

/// The BED text, inflated when the name ends in `.gz`.
fn read_bed_file(path: &Path) -> Result<String> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = File::open(path).map_err(input)?;
    let mut text = String::new();
    if path.to_string_lossy().ends_with(".gz") {
        bgzf::io::Reader::new(file)
            .read_to_string(&mut text)
            .map_err(input)?;
    } else {
        std::io::BufReader::new(file)
            .read_to_string(&mut text)
            .map_err(input)?;
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dictionary() -> Vec<Sequence> {
        parse_sq_lines(
            "@HD\tVN:1.5\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\tM5:abc\tUR:file:/x.fa\n@SQ\tSN:chr2\tLN:500\n",
        )
        .expect("parses")
    }

    fn options() -> Options {
        Options {
            sequence_dictionary: PathBuf::new(),
            sort: true,
            unique: false,
            drop_missing_contigs: false,
            keep_length_zero_intervals: false,
        }
    }

    #[test]
    fn sq_lines_keep_their_other_fields_in_order() {
        let dictionary = dictionary();
        assert_eq!(dictionary.len(), 2);
        assert_eq!(
            dictionary[0].line(),
            "@SQ\tSN:chr1\tLN:1000\tM5:abc\tUR:file:/x.fa"
        );
        assert_eq!(dictionary[1].line(), "@SQ\tSN:chr2\tLN:500");
        assert!(parse_sq_lines("@SQ\tSN:chr1\n").is_err());
    }

    #[test]
    fn bed_fields_split_on_tabs_and_space_runs() {
        assert_eq!(bed_fields("a\tb  c\t\td"), ["a", "b", "c", "", "d"]);
        assert_eq!(bed_fields("a b"), ["a", "b"]);
    }

    #[test]
    fn bed_lines_decode_like_the_codec() {
        let (contig, start, end, name, negative) = decode_bed("chr1\t9\t20\t\"n1\"\t0\t-")
            .expect("ok")
            .expect("feature");
        assert_eq!(
            (contig.as_str(), start, end, name.as_deref(), negative),
            ("chr1", 10, 20, Some("n1"), true)
        );
        // A score that is not a float stops the parse before the strand.
        let (.., negative) = decode_bed("chr1\t9\t20\tx\tscore\t-")
            .expect("ok")
            .expect("feature");
        assert!(!negative);
        let (_, start, end, name, _) = decode_bed("chr1 9").expect("ok").expect("feature");
        assert_eq!((start, end, name), (10, 10, None));
        assert!(decode_bed("# comment").expect("ok").is_none());
        assert!(decode_bed("track name=x").expect("ok").is_none());
        assert!(decode_bed("chr1").expect("ok").is_none());
        assert!(decode_bed("chr1\tx\t5").is_err());
    }

    #[test]
    fn intervals_are_checked_against_the_dictionary() {
        let dictionary = dictionary();
        let bed = "chr2\t0\t10\tb\nchr1\t5\t5\tzero\nchr1\t0\t100\ta\n";
        let intervals = read_bed(bed, &dictionary, &options()).expect("ok");
        assert_eq!(intervals.len(), 2);
        let keep = Options {
            keep_length_zero_intervals: true,
            ..options()
        };
        assert_eq!(read_bed(bed, &dictionary, &keep).expect("ok").len(), 3);
        assert!(read_bed("chrX\t0\t10\n", &dictionary, &options()).is_err());
        let drop = Options {
            drop_missing_contigs: true,
            ..options()
        };
        assert!(
            read_bed("chrX\t0\t10\n", &dictionary, &drop)
                .expect("ok")
                .is_empty()
        );
        assert!(read_bed("chr1\t0\t1001\n", &dictionary, &options()).is_err());
        assert!(read_bed("chr1\t1000\t1000\n", &dictionary, &options()).is_err());
        assert!(read_bed("chr1\t10\t5\n", &dictionary, &options()).is_err());
    }

    #[test]
    fn sorted_and_uniqued_output() {
        let dictionary = dictionary();
        let bed = "chr2\t0\t10\tb\nchr1\t50\t60\tc\t0\t-\nchr1\t0\t100\ta\nchr1\t100\t120\t\n";
        let mut intervals = read_bed(bed, &dictionary, &options()).expect("ok");
        sort(&mut intervals);
        let text = write_interval_list(&dictionary, &intervals);
        assert_eq!(
            text,
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\tM5:abc\tUR:file:/x.fa\n\
             @SQ\tSN:chr2\tLN:500\nchr1\t1\t100\t+\ta\nchr1\t51\t60\t-\tc\nchr1\t101\t120\t+\t.\n\
             chr2\t1\t10\t+\tb\n"
        );
        let unique = uniqued(&intervals);
        assert_eq!(unique.len(), 2);
        assert_eq!(unique[0].name.as_deref(), Some("a|c"));
        assert_eq!(
            (unique[0].start, unique[0].end, unique[0].negative_strand),
            (1, 120, false)
        );
    }
}
