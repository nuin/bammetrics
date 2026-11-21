//! `CollectArraysVariantCallingMetrics`: an arrays VCF as `GtcToVcf` writes
//! it into three files under the prefix `OUTPUT`, as Picard 3.5.0 writes
//! them:
//!
//! - `.arrays_variant_calling_detail_metrics`: one
//!   `ArraysVariantCallingDetailMetrics` row per sample, in Java `HashMap`
//!   order of the sample names (`Collectors.groupingBy` in
//!   `ArraysCallingMetricAccumulator.Result.merge`);
//! - `.arrays_variant_calling_summary_metrics`: one
//!   `ArraysVariantCallingSummaryMetrics` row, the sum of the detail rows;
//! - `.arrays_control_code_summary_metrics`: one
//!   `ArraysControlCodesSummaryMetrics` row per control of
//!   `ArraysControlInfo.CONTROL_INFO`, read from its header line
//!   (`##DNP(High)=CONTROL|CATEGORY|RED|GREEN`).
//!
//! `INPUT` needs the index htsjdk opens (`VCFFileReader(INPUT, true)`): a
//! `.tbi` beside a name ending in `.gz`, `.gzip`, `.bgz` or `.bgzf`, else a
//! Tribble `.idx`. The index is not read: a record counts once when it
//! overlaps 1 to the length of its `##contig` line (a line without a
//! length is length 0), which is what the `VariantProcessor` index queries
//! return. `NUM_PROCESSORS` changes nothing.
//!
//! The detail rows follow a `HashMap` filled with the sample names in
//! reverse header order. A table in which more than 8 names share a bin
//! turns that bin into a tree, whose order is not modelled: such a VCF is
//! refused. At most [`MAX_SAMPLES`] samples are read. A tab, line feed or
//! carriage return in a string cell is refused before any file is
//! written, as htsjdk's `MetricsFile` refuses it.
//!
//! Per sample and record (`updateDetailMetric`): every record adds to
//! `NUM_ASSAYS`. A record that is not filtered, or whose filters include
//! `DUPE`, adds to `NUM_NON_FILTERED_ASSAYS`; its genotype adds to
//! `NUM_CALLS` when called (also to `NUM_AUTOCALL_CALLS` unless its `GTA`
//! is `./.`) and to `NUM_NO_CALLS` otherwise, whatever the genotype; a SNP
//! record adds to `NUM_SNPS` (and to `NUM_IN_DB_SNP` at a dbSNP SNP
//! position), an indel record to `NUM_INDELS`; the sample adds to
//! `NUM_SINGLETONS` when it is the record's singleton sample, and a het or
//! hom-var genotype to the hidden counts of `HET_PCT` and
//! `HET_HOMVAR_RATIO`. Any other filtered record adds to
//! `NUM_FILTERED_ASSAYS`, and to `NUM_ZEROED_OUT_ASSAYS` with the filter
//! `ZEROED_OUT_ASSAY`. Rates divide in `float`; `HET_PCT` and
//! `HET_HOMVAR_RATIO` in `double`.
//!
//! The header gives every detail field but the counts. Missing
//! `sampleAlias`, `arrayType`, `autocallGender`, `autocallVersion`,
//! `autocallDate`, `imagingDate`, `extendedIlluminaManifestVersion`,
//! `clusterFile`, `p95Green`, `p95Red`, `scannerName` or a control line
//! fails the run; `pipelineVersion`, `analysisVersionNumber`,
//! `expectedGender`, `fingerprintGender`, `gtcCallRate`, `zcallVersion`
//! and `zcallThresholds` may be absent. A header key given twice keeps its
//! first line. The dates are read by [`super::header_date`].
//!
//! The VCF is read under [`MAX_VCF_LINE_BYTES`] per line and
//! [`MAX_VCF_LINES`] lines.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles::vcf::header::record::value::Collection;

use super::content_error;
use super::header_date::{self, Pattern};
use crate::error::{Error, Result};
use crate::java::{self, hash_map_order_checked, string_hash_code};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};
use crate::variant_calling_metrics::{output_path, singleton_sample};
use crate::vcf::{
    self, DbSnpSets, Dictionary, Genotype, GenotypeType, Variant, VariantType, VcfReader,
};

pub const DETAIL_EXTENSION: &str = "arrays_variant_calling_detail_metrics";
pub const SUMMARY_EXTENSION: &str = "arrays_variant_calling_summary_metrics";
pub const CONTROL_EXTENSION: &str = "arrays_control_code_summary_metrics";

/// The longest line of the input VCF, in bytes.
pub const MAX_VCF_LINE_BYTES: usize = 8 << 20;
/// The most lines of the input VCF.
pub const MAX_VCF_LINES: u64 = 100_000_000;
/// The most samples of the input VCF.
pub const MAX_SAMPLES: usize = 50_000;

/// `ArraysControlInfo.CONTROL_INFO`: the header key of each control, in
/// the order of the control-code file.
pub const CONTROLS: [&str; 23] = [
    "DNP(High)",
    "DNP(Bgnd)",
    "Biotin(High)",
    "Biotin(Bgnd)",
    "Extension(A)",
    "Extension(T)",
    "Extension(C)",
    "Extension(G)",
    "TargetRemoval",
    "Hyb(High)",
    "Hyb(Medium)",
    "Hyb(Low)",
    "String(PM)",
    "String(MM)",
    "NSB(Bgnd)Red",
    "NSB(Bgnd)Purple",
    "NSB(Bgnd)Blue",
    "NSB(Bgnd)Green",
    "NP(A)",
    "NP(T)",
    "NP(C)",
    "NP(G)",
    "Restore",
];

const SUMMARY_COLUMNS: [&str; 15] = [
    "NUM_ASSAYS",
    "NUM_NON_FILTERED_ASSAYS",
    "NUM_FILTERED_ASSAYS",
    "NUM_ZEROED_OUT_ASSAYS",
    "NUM_SNPS",
    "NUM_INDELS",
    "NUM_CALLS",
    "NUM_AUTOCALL_CALLS",
    "NUM_NO_CALLS",
    "NUM_IN_DB_SNP",
    "NOVEL_SNPS",
    "PCT_DBSNP",
    "CALL_RATE",
    "AUTOCALL_CALL_RATE",
    "NUM_SINGLETONS",
];

const DETAIL_COLUMNS: [&str; 38] = [
    "CHIP_WELL_BARCODE",
    "SAMPLE_ALIAS",
    "ANALYSIS_VERSION",
    "CHIP_TYPE",
    "AUTOCALL_PF",
    "AUTOCALL_DATE",
    "IMAGING_DATE",
    "IS_ZCALLED",
    "GTC_CALL_RATE",
    "AUTOCALL_GENDER",
    "FP_GENDER",
    "REPORTED_GENDER",
    "GENDER_CONCORDANCE_PF",
    "HET_PCT",
    "CLUSTER_FILE_NAME",
    "P95_GREEN",
    "P95_RED",
    "AUTOCALL_VERSION",
    "ZCALL_VERSION",
    "EXTENDED_MANIFEST_VERSION",
    "HET_HOMVAR_RATIO",
    "SCANNER_NAME",
    "PIPELINE_VERSION",
    "NUM_ASSAYS",
    "NUM_NON_FILTERED_ASSAYS",
    "NUM_FILTERED_ASSAYS",
    "NUM_ZEROED_OUT_ASSAYS",
    "NUM_SNPS",
    "NUM_INDELS",
    "NUM_CALLS",
    "NUM_AUTOCALL_CALLS",
    "NUM_NO_CALLS",
    "NUM_IN_DB_SNP",
    "NOVEL_SNPS",
    "PCT_DBSNP",
    "CALL_RATE",
    "AUTOCALL_CALL_RATE",
    "NUM_SINGLETONS",
];

/// The `@MergeByAdding` counts of one sample, or of all of them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub assays: i64,
    pub non_filtered: i64,
    pub filtered: i64,
    pub zeroed_out: i64,
    pub snps: i64,
    pub indels: i64,
    pub calls: i64,
    pub autocall_calls: i64,
    pub no_calls: i64,
    pub in_db_snp: i64,
    pub singletons: i64,
    /// `numHets`, not printed.
    pub hets: i64,
    /// `numHomVar`, not printed.
    pub hom_vars: i64,
}

impl Counts {
    fn add(&mut self, other: &Counts) {
        self.assays += other.assays;
        self.non_filtered += other.non_filtered;
        self.filtered += other.filtered;
        self.zeroed_out += other.zeroed_out;
        self.snps += other.snps;
        self.indels += other.indels;
        self.calls += other.calls;
        self.autocall_calls += other.autocall_calls;
        self.no_calls += other.no_calls;
        self.in_db_snp += other.in_db_snp;
        self.singletons += other.singletons;
        self.hets += other.hets;
        self.hom_vars += other.hom_vars;
    }

    /// `AUTOCALL_CALL_RATE`, a `float`.
    fn autocall_call_rate(&self) -> f32 {
        self.autocall_calls as f32 / self.non_filtered as f32
    }

    /// The summary columns after `calculateDerivedFields`.
    fn values(&self) -> Vec<Value> {
        let rate = |num: i64, den: i64| f64::from(num as f32 / den as f32);
        vec![
            self.assays.into(),
            self.non_filtered.into(),
            self.filtered.into(),
            self.zeroed_out.into(),
            self.snps.into(),
            self.indels.into(),
            self.calls.into(),
            self.autocall_calls.into(),
            self.no_calls.into(),
            self.in_db_snp.into(),
            (self.snps - self.in_db_snp).into(),
            rate(self.in_db_snp, self.snps).into(),
            rate(self.calls, self.non_filtered).into(),
            f64::from(self.autocall_call_rate()).into(),
            self.singletons.into(),
        ]
    }
}

/// One `ArraysVariantCallingSummaryMetrics` row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SummaryRow(pub Counts);

impl MetricRow for SummaryRow {
    const CLASS: &'static str =
        "picard.arrays.CollectArraysVariantCallingMetrics$ArraysVariantCallingSummaryMetrics";
    const COLUMNS: &'static [&'static str] = &SUMMARY_COLUMNS;

    fn values(&self) -> Vec<Value> {
        self.0.values()
    }
}

/// The detail fields the header gives, the same for every sample.
#[derive(Clone, Debug, PartialEq)]
pub struct HeaderFields {
    pub sample_alias: String,
    pub analysis_version: Option<i64>,
    pub chip_type: String,
    pub autocall_date: String,
    pub imaging_date: String,
    pub zcall_thresholds: Option<String>,
    pub gtc_call_rate: Option<f64>,
    pub autocall_gender: Sex,
    pub fp_gender: Sex,
    pub reported_gender: Sex,
    pub cluster_file: String,
    pub p95_green: i64,
    pub p95_red: i64,
    pub autocall_version: String,
    pub zcall_version: Option<String>,
    pub extended_manifest_version: String,
    pub scanner_name: String,
    pub pipeline_version: Option<String>,
}

/// One `ArraysVariantCallingDetailMetrics` row.
#[derive(Clone, Debug, PartialEq)]
pub struct DetailRow {
    pub chip_well_barcode: String,
    /// One copy for every sample.
    pub header: Arc<HeaderFields>,
    pub counts: Counts,
    /// `CALL_RATE_PF_THRESHOLD`.
    pub threshold: f64,
}

impl MetricRow for DetailRow {
    const CLASS: &'static str =
        "picard.arrays.CollectArraysVariantCallingMetrics$ArraysVariantCallingDetailMetrics";
    const COLUMNS: &'static [&'static str] = &DETAIL_COLUMNS;

    fn values(&self) -> Vec<Value> {
        let h = &self.header;
        let c = &self.counts;
        let autocall_pf = f64::from(c.autocall_call_rate()) > self.threshold;
        let mut values = vec![
            self.chip_well_barcode.as_str().into(),
            h.sample_alias.as_str().into(),
            h.analysis_version.into(),
            h.chip_type.as_str().into(),
            autocall_pf.into(),
            h.autocall_date.as_str().into(),
            h.imaging_date.as_str().into(),
            h.zcall_thresholds
                .as_deref()
                .is_some_and(|t| !t.is_empty())
                .into(),
            h.gtc_call_rate.into(),
            h.autocall_gender.symbol().into(),
            h.fp_gender.symbol().into(),
            h.reported_gender.symbol().into(),
            gender_concordance(h.reported_gender, h.fp_gender, h.autocall_gender).into(),
            (c.hets as f64 / c.calls as f64).into(),
            h.cluster_file.as_str().into(),
            h.p95_green.into(),
            h.p95_red.into(),
            h.autocall_version.as_str().into(),
            h.zcall_version.as_deref().into(),
            h.extended_manifest_version.as_str().into(),
            (c.hets as f64 / c.hom_vars as f64).into(),
            h.scanner_name.as_str().into(),
            h.pipeline_version.as_deref().into(),
        ];
        values.extend(c.values());
        values
    }
}

/// `ArraysControlCodesSummaryMetrics`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ControlRow {
    pub control: String,
    pub category: String,
    pub red: i64,
    pub green: i64,
}

impl MetricRow for ControlRow {
    const CLASS: &'static str =
        "picard.arrays.CollectArraysVariantCallingMetrics$ArraysControlCodesSummaryMetrics";
    const COLUMNS: &'static [&'static str] = &["CONTROL", "CATEGORY", "RED", "GREEN"];

    fn values(&self) -> Vec<Value> {
        vec![
            self.control.as_str().into(),
            self.category.as_str().into(),
            self.red.into(),
            self.green.into(),
        ]
    }
}

/// `picard.pedigree.Sex`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    Unknown,
    NotReported,
}

impl Sex {
    const ALL: [(Sex, &'static str, &'static str); 4] = [
        (Sex::Male, "Male", "M"),
        (Sex::Female, "Female", "F"),
        (Sex::Unknown, "Unknown", "U"),
        (Sex::NotReported, "NotReported", "N"),
    ];

    /// `Sex.fromString`: the one value whose symbol or name equals `text`
    /// ignoring case.
    pub fn from_text(text: &str) -> std::result::Result<Sex, String> {
        let found: Vec<Sex> = Self::ALL
            .iter()
            .filter(|(_, name, symbol)| {
                java_equals_ignore_case(text, name) || java_equals_ignore_case(text, symbol)
            })
            .map(|(sex, _, _)| *sex)
            .collect();
        match found.as_slice() {
            [one] => Ok(*one),
            _ => Err(format!("Unrecognized Sex string: {text}")),
        }
    }

    pub fn symbol(self) -> &'static str {
        Self::ALL
            .iter()
            .find(|(sex, _, _)| *sex == self)
            .map_or("U", |(_, _, symbol)| symbol)
    }
}

/// `String.equalsIgnoreCase`: the same UTF-16 length and, unit by unit,
/// equal, equal upper-cased or equal lower-cased (one-to-one case maps).
fn java_equals_ignore_case(a: &str, b: &str) -> bool {
    let single = |mut chars: std::char::ToUppercase| match (chars.next(), chars.next()) {
        (Some(c), None) => Some(c),
        _ => None,
    };
    let lower = |c: char| {
        let mut chars = c.to_lowercase();
        match (chars.next(), chars.next()) {
            (Some(l), None) => l,
            _ => c,
        }
    };
    a.encode_utf16().count() == b.encode_utf16().count()
        && a.chars().zip(b.chars()).all(|(x, y)| {
            let upper = |c: char| single(c.to_uppercase()).unwrap_or(c);
            x == y || upper(x) == upper(y) || lower(x) == lower(y)
        })
}

/// `getSexConcordance(reported, fingerprint, autocall)`: false when all
/// three are unknown or not reported; true when two or three say female
/// and none male, or two or three male and none female.
pub fn gender_concordance(reported: Sex, fingerprint: Sex, autocall: Sex) -> bool {
    let all = [reported, fingerprint, autocall];
    let count = |sex: Sex| all.iter().filter(|s| **s == sex).count();
    if count(Sex::Unknown) + count(Sex::NotReported) == 3 {
        return false;
    }
    let (female, male) = (count(Sex::Female), count(Sex::Male));
    (female > 1 && male == 0) || (male > 1 && female == 0)
}

/// `CollectArraysVariantCallingMetrics` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CollectArraysVariantCallingMetrics",
    summary: "the input is an arrays VCF as GtcToVcf writes it, with its index. OUTPUT is a \
              prefix: PREFIX.arrays_variant_calling_detail_metrics (one row per sample), \
              PREFIX.arrays_variant_calling_summary_metrics and \
              PREFIX.arrays_control_code_summary_metrics; Python returns one row per file: \
              file and rows.",
    options: &[
        Opt::output(&["OUTPUT", "O"])
            .with_help("The output prefix of the detail, summary and control-code files.")
            .required(),
        crate::registry::REFERENCE_NOT_DECODED,
        crate::registry::STRINGENCY_NO_RECORDS,
        Opt::path(
            &["INPUT", "I"],
            "The arrays VCF to read, plain with a .idx or bgzip with a .tbi.",
        )
        .required(),
        Opt::path(&["DBSNP"], "A dbSNP VCF.").required(),
        Opt::float(
            &["CALL_RATE_PF_THRESHOLD"],
            "AUTOCALL_PF is Y when AUTOCALL_CALL_RATE is above this; over 0, at most 1 \
             (default 0.98).",
        ),
        Opt::int(
            &["NUM_PROCESSORS"],
            "Accepted; the output does not depend on it (default 0).",
        ),
        Opt::path(
            &["SEQUENCE_DICTIONARY", "SD"],
            "A .dict, VCF, BAM, interval list or FASTA for the dbSNP contigs (default the \
             INPUT ##contig lines).",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

/// The options of a run.
pub struct Options {
    pub prefix: PathBuf,
    pub dbsnp: PathBuf,
    pub sequence_dictionary: Option<PathBuf>,
    pub call_rate_pf_threshold: f64,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let usage = |message: String| Error::Usage(message);
        let prefix = inv
            .output
            .clone()
            .ok_or_else(|| usage("Option 'OUTPUT' is required.".to_owned()))?;
        let dbsnp = inv
            .get("DBSNP")
            .map(PathBuf::from)
            .ok_or_else(|| usage("Option 'DBSNP' is required.".to_owned()))?;
        if let Some(value) = inv.get("NUM_PROCESSORS") {
            java::parse_integer(value, 32)
                .map_err(|e| usage(format!("NUM_PROCESSORS={value}: {e}")))?;
        }
        let threshold = match inv.get("CALL_RATE_PF_THRESHOLD") {
            Some(value) => java::parse_double(value)
                .map_err(|e| usage(format!("CALL_RATE_PF_THRESHOLD={value}: {e}")))?,
            None => 0.98,
        };
        // `customCommandLineValidation`: NaN passes both comparisons.
        if threshold <= 0.0 || threshold > 1.0 {
            return Err(usage(
                "The parameter CALL_RATE_PF_THRESHOLD must be > 0 and <= 1.0".to_owned(),
            ));
        }
        Ok(Self {
            prefix,
            dbsnp,
            sequence_dictionary: inv.get("SEQUENCE_DICTIONARY").map(PathBuf::from),
            call_rate_pf_threshold: threshold,
        })
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let report = collect(&inv.input, &opts)?;
    let prefix = &opts.prefix;
    let files = [
        (
            DETAIL_EXTENSION,
            MetricsOutput::from_rows(&report.details, None),
        ),
        (
            SUMMARY_EXTENSION,
            MetricsOutput::from_rows(std::slice::from_ref(&report.summary), None),
        ),
        (
            CONTROL_EXTENSION,
            MetricsOutput::from_rows(&report.controls, None),
        ),
    ];
    for (extension, file) in &files {
        check_cells(&inv.input, extension, file)?;
    }
    Ok(ToolOutput {
        main: None,
        extra: files
            .into_iter()
            .map(|(extension, file)| (output_path(prefix, extension), file))
            .collect(),
    })
}

/// htsjdk's `MetricsFile` refuses a string cell holding a tab or a line
/// end (`Supplied String contains illegal character`): the whole run
/// fails before any file is written.
fn check_cells(input: &Path, extension: &str, file: &MetricsOutput) -> Result<()> {
    for (i, row) in file.rows.iter().enumerate() {
        for (column, value) in file.columns.iter().zip(row) {
            let Value::Str(text) = value else { continue };
            let found = text.chars().find_map(|c| match c {
                '\t' => Some("tab"),
                '\n' => Some("line feed"),
                '\r' => Some("carriage return"),
                _ => None,
            });
            if let Some(character) = found {
                return Err(content_error(
                    input,
                    format!(
                        "Picard text cannot hold the {character} in {extension} row {i}, \
                         column {column}: {text:?} (Picard: Supplied String contains illegal \
                         character)"
                    ),
                ));
            }
        }
    }
    Ok(())
}

/// The rows of the three files.
pub struct Report {
    pub details: Vec<DetailRow>,
    pub summary: SummaryRow,
    pub controls: Vec<ControlRow>,
}

/// The index `VCFFileReader(file, true)` finds: `<file>.tbi` for a
/// block-compressed name, else `<file>.idx`.
pub fn has_index(path: &Path) -> bool {
    let name = path.to_string_lossy();
    let lower = name.to_lowercase();
    let compressed = [".gz", ".gzip", ".bgz", ".bgzf"]
        .iter()
        .any(|ext| lower.ends_with(ext));
    (compressed && Path::new(&format!("{name}.tbi")).exists())
        || Path::new(&format!("{name}.idx")).exists()
}

/// `doWork` up to the rows.
pub fn collect(input: &Path, opts: &Options) -> Result<Report> {
    std::fs::metadata(input).map_err(|source| Error::Input {
        path: input.to_path_buf(),
        source,
    })?;
    if !has_index(input) {
        return Err(content_error(
            input,
            "no index: Picard opens INPUT with its index (.tbi beside a bgzip VCF, .idx beside \
             a plain one) and fails with An index is required, but none found"
                .to_owned(),
        ));
    }
    let mut reader = VcfReader::open_capped(input, MAX_VCF_LINE_BYTES, MAX_VCF_LINES)?;
    let header_dictionary: Dictionary = reader.dictionary().ok_or_else(|| {
        content_error(
            input,
            "no ##contig line (Picard's VariantProcessor fails with a NullPointerException)"
                .to_owned(),
        )
    })?;
    let dictionary = match &opts.sequence_dictionary {
        Some(path) => vcf::read_dictionary(path)?,
        None => header_dictionary.clone(),
    };
    // `DbSnpBitSetUtil` compares the dictionaries only when the dbSNP file
    // has `##contig` lines.
    let dbsnp = DbSnpSets::load_lenient(&opts.dbsnp, &dictionary)?;
    let header = reader.header();
    let controls = CONTROLS
        .iter()
        .map(|name| control(header, name))
        .collect::<std::result::Result<Vec<_>, String>>()
        .map_err(|message| content_error(input, message))?;
    let fields = Arc::new(header_fields(header).map_err(|message| content_error(input, message))?);
    let samples = reader.samples();
    if samples.len() > MAX_SAMPLES {
        return Err(content_error(
            input,
            format!(
                "{} samples, over {MAX_SAMPLES}, the most this tool reads",
                samples.len()
            ),
        ));
    }
    let order = detail_order(&samples).map_err(|message| content_error(input, message))?;
    let mut counts = vec![Counts::default(); samples.len()];
    // A contig line without a length is htsjdk's length 0: no query
    // interval overlaps it.
    let lengths: HashMap<String, i64> = header_dictionary
        .into_iter()
        .map(|c| (c.name, c.length.unwrap_or(0)))
        .collect();
    while let Some(vc) = reader.next_variant()? {
        match lengths.get(&vc.contig) {
            Some(&length) if vc.end >= 1 && vc.start <= length => {}
            _ => continue,
        }
        let record = reader.record();
        let filter_field = record.filters();
        let filters: &str = filter_field.as_ref();
        let autocall_no_calls = gta_no_calls(record, samples.len());
        let singleton = singleton_sample(&vc);
        let kind = vc.kind();
        for (i, (genotype, sample)) in vc.genotypes.iter().zip(counts.iter_mut()).enumerate() {
            let site = Site {
                vc: &vc,
                kind,
                filters,
                singleton: singleton == Some(i),
                autocall_no_call: autocall_no_calls[i],
                dbsnp: &dbsnp,
            };
            update(sample, genotype, &site);
        }
    }
    let mut total = Counts::default();
    for sample in &counts {
        total.add(sample);
    }
    let details = order
        .into_iter()
        .map(|i| DetailRow {
            chip_well_barcode: samples[i].clone(),
            header: Arc::clone(&fields),
            counts: counts[i].clone(),
            threshold: opts.call_rate_pf_threshold,
        })
        .collect();
    Ok(Report {
        details,
        summary: SummaryRow(total),
        controls,
    })
}

/// The sample indices in the order Picard writes the detail rows: a
/// `HashMap` keyed by sample name, filled in reverse header order.
fn detail_order(samples: &[String]) -> std::result::Result<Vec<usize>, String> {
    let last = samples.len().saturating_sub(1);
    let reversed: Vec<i32> = samples.iter().rev().map(|s| string_hash_code(s)).collect();
    let order = hash_map_order_checked(&reversed).ok_or_else(|| {
        "more than 8 sample names fall in one HashMap bin, which Java turns into a tree whose \
         order is not modelled"
            .to_owned()
    })?;
    Ok(order.into_iter().map(|i| last - i).collect())
}

/// What one record says for every sample.
struct Site<'a> {
    vc: &'a Variant,
    /// `getType`, once per record.
    kind: VariantType,
    filters: &'a str,
    singleton: bool,
    autocall_no_call: bool,
    dbsnp: &'a DbSnpSets,
}

/// `updateDetailMetric`.
fn update(counts: &mut Counts, genotype: &Genotype, site: &Site) {
    let vc = site.vc;
    let has_filter = |name: &str| site.filters.split(';').any(|f| f == name);
    counts.assays += 1;
    if vc.filtered && !has_filter("DUPE") {
        counts.filtered += 1;
        if has_filter("ZEROED_OUT_ASSAY") {
            counts.zeroed_out += 1;
        }
        return;
    }
    counts.non_filtered += 1;
    let kind = genotype.kind();
    if matches!(kind, GenotypeType::NoCall | GenotypeType::Unavailable) {
        counts.no_calls += 1;
    } else {
        counts.calls += 1;
        if !site.autocall_no_call {
            counts.autocall_calls += 1;
        }
    }
    if site.kind == VariantType::Snp {
        counts.snps += 1;
        if site.dbsnp.snps.is_site(&vc.contig, vc.start) {
            counts.in_db_snp += 1;
        }
    } else if site.kind == VariantType::Indel {
        counts.indels += 1;
    }
    if site.singleton {
        counts.singletons += 1;
    }
    match kind {
        GenotypeType::Het => counts.hets += 1,
        GenotypeType::HomVar => counts.hom_vars += 1,
        _ => {}
    }
}

/// Per sample, whether its `GTA` value is `./.` (the genotype string
/// stands in for a `GTA` that is absent or `.`, and is never `./.` for a
/// called genotype).
fn gta_no_calls(record: &noodles::vcf::Record, samples: usize) -> Vec<bool> {
    let raw = record.samples();
    let keys = raw.keys();
    let Some(index) = keys.iter().position(|k| k == "GTA") else {
        return vec![false; samples];
    };
    let mut out: Vec<bool> = raw
        .iter()
        .map(|sample| {
            let text: &str = sample.as_ref();
            text.split(':').nth(index) == Some("./.")
        })
        .collect();
    out.resize(samples, false);
    out
}

/// `getOtherHeaderLine(key).getValue()`: the first unstructured line of the
/// key.
fn header_value<'h>(
    header: &'h noodles::vcf::Header,
    key: &str,
) -> std::result::Result<Option<&'h str>, String> {
    match header.get(key) {
        None => Ok(None),
        Some(Collection::Unstructured(values)) => Ok(values.first().map(String::as_str)),
        Some(Collection::Structured(_)) => Err(format!(
            "the ##{key} header line holds a <...> map, where a plain value is read"
        )),
    }
}

fn required<'h>(
    header: &'h noodles::vcf::Header,
    key: &str,
) -> std::result::Result<&'h str, String> {
    header_value(header, key)?
        .ok_or_else(|| format!("Input VCF file is missing header line of type '{key}'"))
}

fn integer(key: &str, text: &str) -> std::result::Result<i64, String> {
    java::parse_integer(text, 32).map_err(|e| format!("##{key}: {e}"))
}

/// The `ArraysControlInfo` of one control's header line.
fn control(header: &noodles::vcf::Header, name: &str) -> std::result::Result<ControlRow, String> {
    let value = required(header, name)?;
    let parts = java_split(value, '|');
    let [control, category, red, green, ..] = parts.as_slice() else {
        return Err(format!(
            "##{name}={value}: fewer than 4 fields split at | (Picard: \
             ArrayIndexOutOfBoundsException)"
        ));
    };
    Ok(ControlRow {
        control: (*control).to_owned(),
        category: (*category).to_owned(),
        red: integer(name, red)?,
        green: integer(name, green)?,
    })
}

/// `String.split` at one character: no empty field at the end.
fn java_split(text: &str, separator: char) -> Vec<&str> {
    let mut parts: Vec<&str> = text.split(separator).collect();
    while parts.len() > 1 && parts.last() == Some(&"") {
        parts.pop();
    }
    if parts == [""] && !text.is_empty() {
        parts.clear();
    }
    parts
}

/// `ArraysCallingMetricAccumulator.setup`, in its order.
fn header_fields(header: &noodles::vcf::Header) -> std::result::Result<HeaderFields, String> {
    let optional = |key: &str| header_value(header, key).map(|v| v.map(str::to_owned));
    let sex = |key: &str, text: &str| Sex::from_text(text).map_err(|e| format!("##{key}: {e}"));
    let optional_sex = |key: &str| match header_value(header, key)? {
        Some(text) => sex(key, text),
        None => Ok(Sex::NotReported),
    };
    let sample_alias = required(header, "sampleAlias")?.to_owned();
    let pipeline_version = optional("pipelineVersion")?;
    let analysis_version = header_value(header, "analysisVersionNumber")?
        .map(|v| integer("analysisVersionNumber", v))
        .transpose()?;
    let chip_type = required(header, "arrayType")?.to_owned();
    let reported_gender = optional_sex("expectedGender")?;
    let fp_gender = optional_sex("fingerprintGender")?;
    let gtc_call_rate = header_value(header, "gtcCallRate")?
        .map(|v| java::parse_double(v).map_err(|e| format!("##gtcCallRate: {e}")))
        .transpose()?;
    let autocall_gender = sex("autocallGender", required(header, "autocallGender")?)?;
    let autocall_version = required(header, "autocallVersion")?.to_owned();
    let date = |key: &str, pattern| {
        header_date::iso_date(required(header, key)?, pattern).map_err(|e| format!("##{key}: {e}"))
    };
    let autocall_date = date("autocallDate", Pattern::Autocall)?;
    let imaging_date = date("imagingDate", Pattern::Imaging)?;
    Ok(HeaderFields {
        sample_alias,
        analysis_version,
        chip_type,
        autocall_date,
        imaging_date,
        extended_manifest_version: required(header, "extendedIlluminaManifestVersion")?.to_owned(),
        zcall_version: optional("zcallVersion")?,
        zcall_thresholds: optional("zcallThresholds")?,
        cluster_file: required(header, "clusterFile")?.to_owned(),
        p95_green: integer("p95Green", required(header, "p95Green")?)?,
        p95_red: integer("p95Red", required(header, "p95Red")?)?,
        scanner_name: required(header, "scannerName")?.to_owned(),
        gtc_call_rate,
        autocall_gender,
        fp_gender,
        reported_gender,
        autocall_version,
        pipeline_version,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sex_reads_names_and_symbols_without_case() {
        assert_eq!(Sex::from_text("female"), Ok(Sex::Female));
        assert_eq!(Sex::from_text("m"), Ok(Sex::Male));
        assert_eq!(Sex::from_text("NOTREPORTED"), Ok(Sex::NotReported));
        assert_eq!(Sex::from_text("u"), Ok(Sex::Unknown));
        assert!(Sex::from_text("X").is_err());
        assert!(Sex::from_text("").is_err());
        // KELVIN SIGN lower-cases to k, as in Java.
        assert_eq!(Sex::from_text("Un\u{212A}nown"), Ok(Sex::Unknown));
    }

    #[test]
    fn concordance_needs_two_agreeing_and_none_against() {
        use Sex::*;
        assert!(gender_concordance(Female, Female, Female));
        assert!(gender_concordance(Female, NotReported, Female));
        assert!(!gender_concordance(Female, Male, Female));
        assert!(!gender_concordance(Female, NotReported, Unknown));
        assert!(!gender_concordance(NotReported, Unknown, NotReported));
        assert!(gender_concordance(Male, Male, Unknown));
    }

    #[test]
    fn detail_rows_follow_a_map_filled_in_reverse() {
        let names = |list: &[&str]| list.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>();
        // Both names fall in bin 14 of 16: the later header name comes first.
        assert_eq!(detail_order(&names(&["S1", "."])).unwrap(), [1, 0]);
        assert_eq!(detail_order(&names(&[".", "S1"])).unwrap(), [1, 0]);
        assert_eq!(detail_order(&[]).unwrap(), Vec::<usize>::new());
        // 16 names with one hash code (Aa and BB) make a tree bin.
        let same: Vec<String> = (0..16u32)
            .map(|i| {
                (0..4)
                    .map(|bit| if i >> bit & 1 == 1 { "Aa" } else { "BB" })
                    .collect()
            })
            .collect();
        assert!(detail_order(&same).is_err());
    }

    #[test]
    fn split_drops_trailing_empty_fields() {
        assert_eq!(java_split("a|b|1|2||", '|'), ["a", "b", "1", "2"]);
        assert_eq!(java_split("|b|1|2", '|'), ["", "b", "1", "2"]);
        assert_eq!(java_split("abc", '|'), ["abc"]);
        assert!(java_split("||", '|').is_empty());
        assert_eq!(java_split("", '|'), [""]);
    }
}
