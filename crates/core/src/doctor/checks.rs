//! Each check of `doctor` judged from the [`Facts`] read once. The
//! thresholds are starting values (docs/doctor.md lists them).

use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use super::input::{BgzfEnd, Line};
use super::records::Stats;
use super::reference;
use super::{Check, Eof, Facts, Options, Sex, SheetRow, Status};
use crate::cram::Kind;
use crate::genome_build::{self, BUILD_NAMES, Build, Detection};

mod fasta;

/// The checks, in the order a report lists them.
pub const NAMES: [&str; 19] = [
    "file",
    "eof",
    "header",
    "references",
    "build",
    "read_groups",
    "sort_order",
    "records",
    "record_read_groups",
    "paired",
    "mapped",
    "duplicates",
    "read_length",
    "qualities",
    "md_nm",
    "reference",
    "index",
    "sample_sheet",
    "expected_sex",
];

/// Under this fraction of mapped primary records, `mapped` warns.
const MAPPED_WARN: f64 = 0.5;

/// Phred+64 suspicion: every quality at least this (`@` read as Phred+33)
/// and some quality above [`PHRED64_MAX`].
const PHRED64_MIN: u8 = 31;
const PHRED64_MAX: u8 = 45;

/// Items named in one reason at most.
const SHOWN: usize = 3;

/// `@PG` IDs and names of tools that mark or remove duplicates, lower case.
const DUPLICATE_MARKERS: [&str; 6] = [
    "markduplicates",
    "markdup",
    "samblaster",
    "sambamba",
    "bammarkduplicates",
    "dedup",
];

/// `n` with a comma every three digits.
pub(crate) fn thousands(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn percent(part: u64, whole: u64) -> String {
    format!("{:.1}%", 100.0 * part as f64 / whole.max(1) as f64)
}

/// At most [`SHOWN`] items joined, with the count of the rest.
fn shown<I: IntoIterator<Item = String>>(items: I) -> String {
    let items: Vec<String> = items.into_iter().collect();
    let mut text = items
        .iter()
        .take(SHOWN)
        .cloned()
        .collect::<Vec<_>>()
        .join(", ");
    if items.len() > SHOWN {
        let _ = write!(text, " and {} more", items.len() - SHOWN);
    }
    text
}

/// The reference names records point at: the `@SQ` lines, else a BAM's
/// binary list.
/// The first `@SQ` contig whose name or length the `.fai` of `fasta`
/// does not have, as a reason not to decode a CRAM with it; `None` when
/// every contig is there or the `.fai` cannot be read (the `reference`
/// check reports that).
pub(crate) fn fasta_disagrees(fasta: &Path, lines: &[Line]) -> Option<String> {
    let entries = reference::read_fai(&reference::fai_path(fasta)).ok()?;
    let by_name: HashMap<&str, u64> = entries
        .iter()
        .map(|e| (e.name.as_str(), e.length))
        .collect();
    lines.iter().filter(|l| l.kind == "SQ").find_map(|l| {
        let name = l.get("SN")?;
        let length: u64 = l.get("LN")?.parse().ok()?;
        (by_name.get(name) != Some(&length)).then(|| {
            format!("the records were not decoded: {name} differs between @SQ and the -R FASTA")
        })
    })
}

pub(crate) fn reference_names(lines: &[Line], binary: &[(String, u64)]) -> Vec<String> {
    let text: Vec<String> = lines
        .iter()
        .filter(|l| l.kind == "SQ")
        .map(|l| l.get("SN").unwrap_or_default().to_owned())
        .collect();
    if text.is_empty() {
        binary.iter().map(|(name, _)| name.clone()).collect()
    } else {
        text
    }
}

/// The dictionary the checks use: `@SQ` names and lengths (a length that
/// does not parse as 0), else a BAM's binary list.
fn dictionary(facts: &Facts) -> Vec<(String, u64)> {
    let text: Vec<(String, u64)> = facts
        .lines
        .iter()
        .filter(|l| l.kind == "SQ")
        .map(|l| {
            let length = l.get("LN").and_then(|v| v.parse().ok()).unwrap_or(0);
            (l.get("SN").unwrap_or_default().to_owned(), length)
        })
        .collect();
    if !text.is_empty() {
        return text;
    }
    match &facts.header {
        Some(Ok(header)) => header.binary.clone(),
        _ => Vec::new(),
    }
}

/// Every check of `path`, in [`NAMES`] order.
pub(crate) fn judge(path: &Path, facts: &Facts, options: &Options) -> Vec<Check> {
    let mut out = Vec::new();
    let Some(kind) = facts.kind else {
        let value = match facts.short {
            Some(0) => "empty file",
            Some(_) => "under 4 bytes",
            None => "not BAM, CRAM or SAM",
        };
        out.push(
            Check::new("file", Status::Fail, value)
                .because(facts.kind_error.clone().unwrap_or_default()),
        );
        skip_rest(&mut out, "the file is not BAM, CRAM or SAM");
        return out;
    };
    out.push(file(kind, facts));
    out.push(eof(facts.eof));
    let header = match &facts.header {
        Some(Ok(header)) => header,
        Some(Err(e)) => {
            out.push(Check::new("header", Status::Fail, "unreadable").because(e.clone()));
            skip_rest(&mut out, "the header could not be read");
            return out;
        }
        None => {
            skip_rest(&mut out, "the header could not be read");
            return out;
        }
    };
    out.push(header_lines(&facts.lines));
    let dict = dictionary(facts);
    let stats = facts.stats.as_ref();
    out.push(references(&facts.lines, &header.binary, &dict, stats));
    let detection = genome_build::detect(dict.iter().map(|(n, l)| (n.as_str(), *l)));
    out.push(build(&dict, &detection, options.expect_build));
    out.push(read_groups(&facts.lines));
    let so = facts
        .lines
        .iter()
        .find(|l| l.kind == "HD")
        .and_then(|l| l.get("SO"))
        .unwrap_or("")
        .to_owned();
    out.push(sort_order(&so, stats));
    out.push(records(facts, options.records, options.reference.is_some()));
    let none = if options.records == 0 {
        "--records 0 reads no record"
    } else {
        "no record was read"
    };
    match stats.filter(|s| s.read > 0) {
        Some(s) => {
            out.push(record_read_groups(&facts.lines, s));
            out.push(paired(s));
            out.push(mapped(s, dict.is_empty()));
            out.push(duplicates(&facts.lines, s));
            out.push(read_length(s));
            out.push(qualities(s));
            out.push(md_nm(kind, s));
        }
        None => {
            for name in [
                "record_read_groups",
                "paired",
                "mapped",
                "duplicates",
                "read_length",
                "qualities",
                "md_nm",
            ] {
                out.push(Check::new(name, Status::Skip, "").because(none));
            }
        }
    }
    out.push(reference_check(kind, facts, &dict, options));
    out.push(index(path, kind, &so));
    out.push(sample_sheet(
        path,
        &facts.lines,
        &detection,
        options.sheet.as_ref(),
    ));
    out.push(expected_sex(&dict, options));
    out
}

/// Skip every check not yet in `out`.
fn skip_rest(out: &mut Vec<Check>, reason: &str) {
    for name in NAMES {
        if !out.iter().any(|c| c.name == name) {
            out.push(Check::new(name, Status::Skip, "").because(reason));
        }
    }
}

fn file(kind: Kind, facts: &Facts) -> Check {
    let mut value = super::kind_name(kind, facts);
    if let Some(version) = facts
        .lines
        .iter()
        .find(|l| l.kind == "HD")
        .and_then(|l| l.get("VN"))
    {
        let _ = write!(value, ", @HD VN:{version}");
    }
    Check::new("file", Status::Pass, value)
}

fn eof(eof: Option<Eof>) -> Check {
    let check = |status, value: &str, reason: &str| {
        Check::new("eof", status, value).because(reason.to_owned())
    };
    match eof {
        Some(Eof::Bgzf(BgzfEnd::Marker)) => check(
            Status::Pass,
            "end-of-file block present",
            "the file ends with the 28-byte BGZF end-of-file block",
        ),
        Some(Eof::Bgzf(BgzfEnd::WholeBlock)) => check(
            Status::Fail,
            "no end-of-file block",
            "the last BGZF block is whole but the 28-byte end-of-file block is missing: the \
             writer stopped before closing the file, or the file was cut at a block boundary",
        ),
        Some(Eof::Bgzf(BgzfEnd::Cut)) => check(
            Status::Fail,
            "cut short",
            "the last BGZF block is incomplete: the file was cut short (an interrupted copy or \
             a full disk)",
        ),
        Some(Eof::Bgzf(BgzfEnd::Trailing(n))) => {
            Check::new("eof", Status::Fail, "bytes after the end-of-file block").because(format!(
                "{} after the 28-byte end-of-file block start no BGZF block: readers that read \
             past the block fail on them (something was appended to the file)",
                plural(n as usize, "byte")
            ))
        }
        Some(Eof::Cram(true)) => check(
            Status::Pass,
            "EOF container present",
            "the file ends with the CRAM 3 end-of-file container",
        ),
        Some(Eof::Cram(false)) => check(
            Status::Fail,
            "no EOF container",
            "the file does not end with the CRAM 3 end-of-file container: it was cut short or \
             its writer stopped",
        ),
        Some(Eof::Text(true)) => check(Status::Pass, "ends with a line end", ""),
        Some(Eof::Text(false)) => check(
            Status::Warn,
            "no line end at the end",
            "the last line has no line end: the file may be cut short",
        ),
        Some(Eof::Gzip) | None => check(
            Status::Skip,
            "",
            "plain gzip has no end-of-file marker to look for",
        ),
    }
}

fn header_lines(lines: &[Line]) -> Check {
    let count = |kind: &str| lines.iter().filter(|l| l.kind == kind).count();
    let value = format!(
        "{} lines: {} @HD, {} @SQ, {} @RG, {} @PG, {} @CO",
        lines.len(),
        count("HD"),
        count("SQ"),
        count("RG"),
        count("PG"),
        count("CO")
    );
    Check::new("header", Status::Pass, value)
}

fn references(
    lines: &[Line],
    binary: &[(String, u64)],
    dict: &[(String, u64)],
    stats: Option<&Stats>,
) -> Check {
    let sq: Vec<&Line> = lines.iter().filter(|l| l.kind == "SQ").collect();
    let total: u64 = dict.iter().map(|(_, l)| *l).sum();
    let value = format!(
        "{}, {} bases",
        plural(dict.len(), "contig"),
        thousands(total)
    );
    let fail =
        |reason: String| Check::new("references", Status::Fail, value.clone()).because(reason);
    let mut problems = Vec::new();
    let mut seen = HashSet::new();
    for (i, line) in sq.iter().enumerate() {
        let number = i + 1;
        match line.get("SN") {
            None => problems.push(format!("@SQ line {number} has no SN")),
            Some(name) if !seen.insert(name) => problems.push(format!("{name} appears twice")),
            Some(_) => {}
        }
        if line.get("LN").and_then(|v| v.parse::<u64>().ok()).is_none() {
            problems.push(format!("@SQ line {number} has no LN number"));
        }
    }
    if !problems.is_empty() {
        return fail(shown(problems));
    }
    if !sq.is_empty() && !binary.is_empty() && dict != binary {
        let at = dict
            .iter()
            .zip(binary)
            .position(|(a, b)| a != b)
            .unwrap_or(dict.len().min(binary.len()));
        let what = |list: &[(String, u64)]| {
            list.get(at)
                .map_or("nothing".to_owned(), |(n, l)| format!("{n} of {l}"))
        };
        return fail(format!(
            "the BAM's binary reference list differs from the @SQ lines at contig {}: {} in \
             the list, {} in the text",
            at + 1,
            what(binary),
            what(dict)
        ));
    }
    if dict.is_empty() {
        if stats.is_some_and(|s| s.placed > 0) {
            return fail("records are placed on contigs, but the header names none".to_owned());
        }
        return Check::new("references", Status::Pass, value)
            .because("no @SQ line: an unaligned file");
    }
    if let Some(s) = stats.filter(|s| s.beyond_dictionary > 0) {
        return fail(format!(
            "{} on a contig number the header lacks (it names {})",
            plural(s.beyond_dictionary as usize, "record"),
            plural(dict.len(), "contig")
        ));
    }
    let empty: Vec<String> = dict
        .iter()
        .filter(|(_, length)| *length == 0)
        .map(|(name, _)| name.clone())
        .collect();
    if !empty.is_empty() {
        return Check::new("references", Status::Warn, value).because(format!(
            "a zero-length contig: {} (samtools refuses such a header)",
            shown(empty)
        ));
    }
    let reason = if sq.is_empty() {
        "no @SQ line; the BAM's binary reference list names them"
    } else {
        ""
    };
    Check::new("references", Status::Pass, value).because(reason)
}

/// What `detection` says against `expected`: `None` when they agree.
fn build_mismatch(detection: &Detection, expected: Build) -> Option<String> {
    match detection {
        Detection::Known { assembly, chr, .. } => {
            if *assembly != expected.assembly {
                Some(format!(
                    "expected {expected}, found {}",
                    assembly.name(*chr)
                ))
            } else if expected.chr.is_some_and(|want| want != *chr) {
                let names = if *chr { "with" } else { "without" };
                Some(format!(
                    "expected {expected}, found {} (contig names {names} chr)",
                    assembly.name(*chr)
                ))
            } else {
                None
            }
        }
        Detection::Mixed(_) => Some(format!("expected {expected}, found a mix of builds")),
        Detection::Unknown => Some(format!("expected {expected}, found no known build")),
    }
}

fn build(dict: &[(String, u64)], detection: &Detection, expected: Option<Build>) -> Check {
    if dict.is_empty() {
        return Check::new("build", Status::Skip, "").because("no contig in the header");
    }
    let c = genome_build::contigs(dict.iter().map(|(n, l)| (n.as_str(), *l)));
    let composition = format!(
        "{} primary, {} mitochondrial, {} alt, {} decoy, {} HLA, {} other contigs",
        c.primary, c.mito, c.alt, c.decoy, c.hla, c.other
    );
    let (status, value, mut reason) = match detection {
        Detection::Known { assembly, chr, .. } => {
            (Status::Pass, assembly.name(*chr).to_owned(), composition)
        }
        Detection::Mixed(said) => (
            Status::Fail,
            "mixed".to_owned(),
            format!(
                "primary chromosomes of different builds: {}",
                said.join("; ")
            ),
        ),
        Detection::Unknown => (
            Status::Pass,
            "unknown".to_owned(),
            format!(
                "no primary chromosome has a GRCh37, GRCh38 or T2T-CHM13 length; {composition}"
            ),
        ),
    };
    match expected.map(|e| (e, build_mismatch(detection, e))) {
        Some((_, Some(mismatch))) => Check::new("build", Status::Fail, value).because(mismatch),
        Some((e, None)) => {
            reason = format!("{e} as expected; {reason}");
            Check::new("build", status, value).because(reason)
        }
        None => Check::new("build", status, value).because(reason),
    }
}

fn read_groups(lines: &[Line]) -> Check {
    let rg: Vec<&Line> = lines.iter().filter(|l| l.kind == "RG").collect();
    if rg.is_empty() {
        return Check::new("read_groups", Status::Warn, "0 @RG lines").because(
            "no @RG line: tools that report by sample, library or read group, and GATK, need one",
        );
    }
    let mut samples: Vec<&str> = rg.iter().filter_map(|l| l.get("SM")).collect();
    samples.sort_unstable();
    samples.dedup();
    let value = format!(
        "{}, {}: {}",
        plural(rg.len(), "@RG line"),
        plural(samples.len(), "sample"),
        shown(samples.iter().map(|s| (*s).to_owned()))
    );
    let mut counts: HashMap<&str, usize> = HashMap::new();
    for id in rg.iter().filter_map(|l| l.get("ID")) {
        *counts.entry(id).or_default() += 1;
    }
    let mut duplicated: Vec<(&str, usize)> = counts.into_iter().filter(|(_, n)| *n > 1).collect();
    duplicated.sort_unstable();
    let missing = |tag: &str| rg.iter().filter(|l| l.get(tag).is_none()).count();
    let mut fails = Vec::new();
    if !duplicated.is_empty() {
        fails.push(format!(
            "duplicate IDs: {}",
            shown(duplicated.iter().map(|(id, n)| format!("{id} ({n} lines)")))
        ));
    }
    for tag in ["ID", "SM"] {
        let n = missing(tag);
        if n > 0 {
            fails.push(format!("{n} lines without {tag}"));
        }
    }
    if !fails.is_empty() {
        return Check::new("read_groups", Status::Fail, value).because(fails.join("; "));
    }
    let warns: Vec<String> = ["LB", "PL"]
        .into_iter()
        .filter_map(|tag| {
            let n = missing(tag);
            (n > 0).then(|| format!("{n} lines without {tag}"))
        })
        .collect();
    if !warns.is_empty() {
        return Check::new("read_groups", Status::Warn, value).because(format!(
            "{}: MarkDuplicates groups by LB, base recalibration needs PL",
            warns.join("; ")
        ));
    }
    Check::new("read_groups", Status::Pass, value).because("every line has ID, SM, LB and PL")
}

fn sort_order(so: &str, stats: Option<&Stats>) -> Check {
    let shown_so = if so.is_empty() { "absent" } else { so };
    let value = format!("header {shown_so}");
    let Some(s) = stats.filter(|s| s.read > 0) else {
        return Check::new("sort_order", Status::Skip, value).because("no record was read");
    };
    if s.read < 2 {
        return Check::new("sort_order", Status::Pass, value)
            .because("one record: no order to check");
    }
    let n = thousands(s.read);
    match so {
        "coordinate" => match &s.coordinate_break {
            Some(b) => Check::new("sort_order", Status::Fail, value)
                .because(format!("the header says coordinate, but {}", b.what)),
            None => Check::new("sort_order", Status::Pass, value)
                .because(format!("the first {n} records are in coordinate order")),
        },
        "queryname" => match (&s.name_break, &s.natural_break) {
            (Some(b), Some(_)) => Check::new("sort_order", Status::Fail, value).because(format!(
                "the header says queryname, but {} (in byte order and in natural order)",
                b.what
            )),
            (None, _) => Check::new("sort_order", Status::Pass, value).because(format!(
                "the first {n} records are in name order (byte order, as Picard sorts)"
            )),
            (Some(_), None) => Check::new("sort_order", Status::Pass, value).because(format!(
                "the first {n} records are in name order (natural order, as samtools sorts; \
                 Picard expects byte order)"
            )),
        },
        _ if s.coordinate_break.is_none() && s.placed > 0 => {
            Check::new("sort_order", Status::Warn, value).because(format!(
                "the first {n} records are in coordinate order but the header does not say so: \
                 tools that need sorted input refuse the file (SortSam or samtools sort sets \
                 SO:coordinate)"
            ))
        }
        _ => Check::new("sort_order", Status::Pass, value)
            .because(format!("the first {n} records are not in coordinate order")),
    }
}

fn records(facts: &Facts, limit: u64, reference_given: bool) -> Check {
    if limit == 0 {
        return Check::new("records", Status::Skip, "").because("--records 0 reads no record");
    }
    if let Some(why) = &facts.records_skipped {
        return Check::new("records", Status::Skip, "").because(why.clone());
    }
    if let Some(e) = &facts.records_error {
        return Check::new("records", Status::Fail, "0 records read").because(e.clone());
    }
    let Some(s) = &facts.stats else {
        return Check::new("records", Status::Skip, "").because("no record was read");
    };
    let value = format!("{} records read", thousands(s.read));
    if let Some(e) = &s.error {
        if s.read == 0 && !reference_given && needs_reference(e) {
            return Check::new("records", Status::Skip, value)
                .because("the CRAM records need reference bases: give -R to read them");
        }
        return Check::new("records", Status::Fail, value).because(e.clone());
    }
    if s.read == 0 {
        return Check::new("records", Status::Warn, value).because("the file holds no record");
    }
    let reason = if s.read < limit {
        "every record of the file".to_owned()
    } else {
        format!("the first {}", thousands(limit))
    };
    Check::new("records", Status::Pass, value).because(reason)
}

fn record_read_groups(lines: &[Line], s: &Stats) -> Check {
    let header_has_rg = lines.iter().any(|l| l.kind == "RG");
    let n = thousands(s.read);
    if s.unknown_rg_records > 0 {
        let ids = shown(s.unknown_rg.iter().map(|(id, k)| format!("{id} ({k})")));
        return Check::new(
            "record_read_groups",
            Status::Fail,
            format!(
                "{} records with an RG the header lacks",
                s.unknown_rg_records
            ),
        )
        .because(format!("RG not in the header: {ids}"));
    }
    if header_has_rg && s.without_rg > 0 {
        return Check::new(
            "record_read_groups",
            Status::Warn,
            format!("{} records without RG", s.without_rg),
        )
        .because(format!("of the first {n}; the header has @RG lines"));
    }
    if !header_has_rg {
        return Check::new("record_read_groups", Status::Pass, "no RG tag")
            .because("no @RG line and no RG tag");
    }
    Check::new("record_read_groups", Status::Pass, format!("{n} records"))
        .because("every RG is in the header")
}

fn paired(s: &Stats) -> Check {
    if s.primary == 0 {
        return Check::new("paired", Status::Skip, "").because("no primary record");
    }
    let value = percent(s.paired, s.primary);
    match s.paired {
        0 => Check::new("paired", Status::Pass, value).because("single-end"),
        p if p == s.primary => Check::new("paired", Status::Pass, value).because("paired-end"),
        _ => Check::new("paired", Status::Warn, value).because(format!(
            "paired and single-end records mixed among {} primary records",
            thousands(s.primary)
        )),
    }
}

fn mapped(s: &Stats, unaligned: bool) -> Check {
    if s.primary == 0 {
        return Check::new("mapped", Status::Skip, "").because("no primary record");
    }
    let value = percent(s.mapped, s.primary);
    if unaligned && s.mapped == 0 {
        return Check::new("mapped", Status::Pass, value).because("an unaligned file");
    }
    let of = thousands(s.primary);
    if (s.mapped as f64) < MAPPED_WARN * s.primary as f64 {
        return Check::new("mapped", Status::Warn, value).because(format!(
            "under 50% of {of} primary records mapped (a starting value)"
        ));
    }
    Check::new("mapped", Status::Pass, value).because(format!("of {of} primary records"))
}

fn duplicates(lines: &[Line], s: &Stats) -> Check {
    if s.mapped == 0 {
        return Check::new("duplicates", Status::Skip, "")
            .because("no mapped primary record: duplicates are marked on aligned records");
    }
    let marker = lines
        .iter()
        .filter(|l| l.kind == "PG")
        .filter_map(|l| {
            let named = |tag| l.get(tag).map(str::to_ascii_lowercase);
            [named("ID"), named("PN")]
                .into_iter()
                .flatten()
                .find(|name| DUPLICATE_MARKERS.iter().any(|m| name.contains(m)))
                .map(|_| l.get("ID").or(l.get("PN")).unwrap_or_default().to_owned())
        })
        .next();
    let value = format!("{} flagged", percent(s.duplicate, s.primary));
    let of = thousands(s.primary);
    match (s.duplicate, marker) {
        (0, None) => Check::new("duplicates", Status::Warn, value).because(format!(
            "no duplicate flag among {of} primary records and no @PG of a duplicate marker: \
             duplicates were probably not marked"
        )),
        (0, Some(pg)) => Check::new("duplicates", Status::Pass, value).because(format!(
            "no flag among {of} primary records; @PG {pg} ran (duplicates may have been removed)"
        )),
        (_, Some(pg)) => {
            Check::new("duplicates", Status::Pass, value).because(format!("@PG {pg} ran"))
        }
        (_, None) => Check::new("duplicates", Status::Pass, value)
            .because("duplicates are marked (no @PG names the tool)"),
    }
}

fn read_length(s: &Stats) -> Check {
    if s.with_bases == 0 {
        return Check::new("read_length", Status::Warn, "no bases")
            .because("no primary record has bases (SEQ is *)");
    }
    let mean = s.bases as f64 / s.with_bases as f64;
    Check::new("read_length", Status::Pass, format!("mean {mean:.1}")).because(format!(
        "over {} primary records with bases",
        thousands(s.with_bases)
    ))
}

fn qualities(s: &Stats) -> Check {
    if s.with_qualities == 0 {
        return Check::new("qualities", Status::Warn, "no qualities")
            .because("no primary record has qualities (QUAL is *)");
    }
    let value = format!("Q{} to Q{}", s.quality_min, s.quality_max);
    if s.quality_min >= PHRED64_MIN && s.quality_max > PHRED64_MAX {
        return Check::new("qualities", Status::Fail, value).because(format!(
            "Phred+64 suspected: no quality under Q{PHRED64_MIN} and some over Q{PHRED64_MAX}; \
             qualities read from Phred+64 FASTQ as Phred+33 are 31 too high"
        ));
    }
    Check::new("qualities", Status::Pass, value).because("Phred+33")
}

fn md_nm(kind: Kind, s: &Stats) -> Check {
    if kind == Kind::Cram {
        return Check::new("md_nm", Status::Skip, "")
            .because("a CRAM stores no MD or NM; readers rebuild them from the reference");
    }
    if s.mapped_primary == 0 {
        return Check::new("md_nm", Status::Skip, "").because("no mapped primary record");
    }
    let value = format!(
        "MD {}, NM {}",
        percent(s.with_md, s.mapped_primary),
        percent(s.with_nm, s.mapped_primary)
    );
    if s.with_md == s.mapped_primary && s.with_nm == s.mapped_primary {
        return Check::new("md_nm", Status::Pass, value).because(format!(
            "on every one of {} mapped primary records",
            thousands(s.mapped_primary)
        ));
    }
    Check::new("md_nm", Status::Warn, value).because(
        "some mapped records lack MD or NM: tools that need them (variant callers, NM-based \
         filters) recompute them or skip the records (samtools calmd adds them)",
    )
}

/// True for the error of a CRAM decoded without the reference it needs.
fn needs_reference(error: &str) -> bool {
    error.contains("REFERENCE_SEQUENCE")
}

/// `n` and `noun`, with an `s` unless `n` is 1.
fn plural(n: usize, noun: &str) -> String {
    if n == 1 {
        format!("1 {noun}")
    } else {
        format!("{} {noun}s", thousands(n as u64))
    }
}

/// The error text of the record read, when there is one.
fn record_error(facts: &Facts) -> Option<&str> {
    facts
        .records_error
        .as_deref()
        .or_else(|| facts.stats.as_ref().and_then(|s| s.error.as_deref()))
}

fn reference_check(kind: Kind, facts: &Facts, dict: &[(String, u64)], options: &Options) -> Check {
    let Some(fasta) = &options.reference else {
        if kind != Kind::Cram {
            return Check::new("reference", Status::Skip, "").because("no -R given");
        }
        let url = facts
            .lines
            .iter()
            .find(|l| l.kind == "SQ")
            .and_then(|l| l.get("UR"));
        let hint = url.map_or_else(String::new, |u| format!(" (the first @SQ names UR:{u})"));
        if record_error(facts).is_some_and(needs_reference) {
            return Check::new("reference", Status::Warn, "needed").because(format!(
                "the CRAM slices need reference bases: give -R, and the run needs \
                 REFERENCE_SEQUENCE{hint}"
            ));
        }
        return match &facts.stats {
            Some(s) if s.read > 0 && s.error.is_none() => {
                Check::new("reference", Status::Pass, "not needed").because(format!(
                    "the first {} records decode without a reference",
                    thousands(s.read)
                ))
            }
            _ => Check::new("reference", Status::Skip, "")
                .because("no -R given and no record decoded"),
        };
    };
    fasta::against_fasta(fasta, facts, dict, options)
}

fn file_name(path: &Path) -> String {
    path.file_name().map_or_else(
        || path.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

/// The index beside a CRAM: `x.cram.crai`, then `x.crai`.
fn crai(path: &Path) -> Option<PathBuf> {
    let name = path.file_name()?.to_str()?;
    let mut candidates = vec![format!("{name}.crai")];
    if let Some(stem) = name.strip_suffix(".cram") {
        candidates.push(format!("{stem}.crai"));
    }
    candidates
        .into_iter()
        .map(|c| path.with_file_name(c))
        .find(|p| p.is_file())
}

fn index(path: &Path, kind: Kind, so: &str) -> Check {
    let found = match kind {
        Kind::Bam => crate::bam_index_stats::find_index(path),
        Kind::Cram => crai(path),
        Kind::Sam | Kind::GzipSam => {
            return Check::new("index", Status::Skip, "").because("SAM text is read whole");
        }
    };
    if so != "coordinate" {
        return Check::new("index", Status::Skip, "")
            .because("an index needs a coordinate-sorted file");
    }
    let Some(index) = found else {
        let make = if kind == Kind::Cram {
            "samtools index"
        } else {
            "bammetrics build-bam-index"
        };
        return Check::new("index", Status::Warn, "none").because(format!(
            "no index beside the file: region queries need one ({make} writes it)"
        ));
    };
    let modified = |p: &Path| std::fs::metadata(p).and_then(|m| m.modified()).ok();
    let name = file_name(&index);
    match (modified(&index), modified(path)) {
        (Some(i), Some(f)) if i < f => {
            let seconds = f.duration_since(i).map_or(0, |d| d.as_secs());
            Check::new("index", Status::Warn, name).because(format!(
                "the index is {} s older than the file: it may describe an earlier version \
                 (rewrite it), or a copy or download wrote the two in that order",
                thousands(seconds)
            ))
        }
        _ => Check::new("index", Status::Pass, name).because("not older than the file"),
    }
}

/// The `bam` column resolved against the sheet's directory.
fn sheet_path(sheet: &Path, named: &str) -> PathBuf {
    let named = Path::new(named);
    if named.is_absolute() {
        return named.to_path_buf();
    }
    sheet
        .parent()
        .map_or_else(|| named.to_path_buf(), |dir| dir.join(named))
}

/// True when the sheet's `bam` and `path` name one file: by device and
/// inode when both exist, else by resolved path.
fn same_file(a: &Path, b: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(x), Ok(y)) = (std::fs::metadata(a), std::fs::metadata(b)) {
            return x.dev() == y.dev() && x.ino() == y.ino();
        }
    }
    a.canonicalize()
        .ok()
        .zip(b.canonicalize().ok())
        .is_some_and(|(x, y)| x == y)
}

fn present(value: Option<&String>) -> Option<&str> {
    value
        .map(|v| v.trim())
        .filter(|v| !v.is_empty() && !v.eq_ignore_ascii_case("na"))
}

fn sample_sheet(
    path: &Path,
    lines: &[Line],
    detection: &Detection,
    row: Option<&SheetRow>,
) -> Check {
    let Some(row) = row else {
        return Check::new("sample_sheet", Status::Skip, "").because("no --sample-sheet given");
    };
    let sheet = file_name(&row.sheet);
    let value = format!("{} in {sheet}", row.sample);
    if row.rows != 1 {
        let what = if row.rows == 0 {
            "no row names the sample".to_owned()
        } else {
            format!("{} rows name the sample", row.rows)
        };
        return Check::new("sample_sheet", Status::Fail, value).because(what);
    }
    let mut compared = vec!["SM"];
    let mut problems = Vec::new();
    let samples: Vec<&str> = lines
        .iter()
        .filter(|l| l.kind == "RG")
        .filter_map(|l| l.get("SM"))
        .collect();
    if samples.is_empty() {
        problems.push("no @RG line names a sample".to_owned());
    } else if !samples.contains(&row.sample.as_str()) {
        let mut distinct = samples.clone();
        distinct.sort_unstable();
        distinct.dedup();
        problems.push(format!(
            "SM is {}, the sheet says {}",
            shown(distinct.iter().map(|s| (*s).to_owned())),
            row.sample
        ));
    }
    if let Some(name) = present(row.build.as_ref()) {
        compared.push("build");
        match Build::parse(name) {
            None => problems.push(format!("build {name} is none of {BUILD_NAMES}")),
            Some(b) => problems.extend(build_mismatch(detection, b)),
        }
    }
    if let Some(named) = present(row.bam.as_ref()) {
        compared.push("bam");
        if !same_file(&sheet_path(&row.sheet, named), path) {
            problems.push(format!("the sheet names {named}, not this file"));
        }
    }
    if let Some(library) = present(row.library.as_ref()) {
        compared.push("library");
        let libraries: Vec<&str> = lines
            .iter()
            .filter(|l| l.kind == "RG")
            .filter_map(|l| l.get("LB"))
            .collect();
        if !libraries.contains(&library) {
            problems.push(format!(
                "no @RG has LB {library} (LB: {})",
                if libraries.is_empty() {
                    "none".to_owned()
                } else {
                    shown(libraries.iter().map(|s| (*s).to_owned()))
                }
            ));
        }
    }
    if let Some(sex) = present(row.sex.as_ref())
        && !sex.eq_ignore_ascii_case("unknown")
        && Sex::parse(sex).is_none()
    {
        problems.push(format!("sex {sex} is neither male nor female"));
    }
    if problems.is_empty() {
        return Check::new("sample_sheet", Status::Pass, value)
            .because(format!("{} agree", compared.join(", ")));
    }
    Check::new("sample_sheet", Status::Fail, value).because(problems.join("; "))
}

fn expected_sex(dict: &[(String, u64)], options: &Options) -> Check {
    let sheet = options
        .sheet
        .as_ref()
        .and_then(|row| present(row.sex.as_ref()))
        .and_then(Sex::parse);
    let Some(sex) = options.expect_sex.or(sheet) else {
        return Check::new("expected_sex", Status::Skip, "")
            .because("no --expect-sex and no sex in the sample sheet");
    };
    let value = format!("{} expected", sex.name());
    let has_y = dict
        .iter()
        .any(|(name, _)| matches!(name.as_str(), "Y" | "chrY"));
    let limit = "doctor reads too few records to infer sex; bammetrics sex-check does";
    if sex == Sex::Male && !dict.is_empty() && !has_y {
        return Check::new("expected_sex", Status::Warn, value)
            .because(format!("the header has no Y contig; {limit}"));
    }
    Check::new("expected_sex", Status::Pass, value).because(format!("not inferred: {limit}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thousands_groups_digits() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1000), "1,000");
        assert_eq!(thousands(3_099_922_541), "3,099,922,541");
    }

    #[test]
    fn shown_names_the_rest_by_count() {
        let items = (1..=5).map(|i| i.to_string());
        assert_eq!(shown(items), "1, 2, 3 and 2 more");
    }
}
