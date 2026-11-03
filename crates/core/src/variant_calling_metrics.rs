//! `CollectVariantCallingMetrics` and `AccumulateVariantCallingMetrics`
//! (plan rows P7), over the VCF layer of `vcf.rs`.
//!
//! The collector is `CallingMetricAccumulator`: a record counts when it is
//! a variant (its allele list holds an alternate) and some genotype is not
//! hom-ref. A filtered record adds to `FILTERED_SNPS` or `FILTERED_INDELS`
//! by its type and nothing else. A passing record first adds to
//! `NUM_SINGLETONS` when exactly one sample carries exactly one variant
//! chromosome, then goes to one bucket: biallelic SNP (dbSNP membership by
//! contig and start, transition or transversion, the het samples' `AD`
//! into the reference bias), multiallelic SNP, simple indel (insertion or
//! deletion, dbSNP by the indel bit set) or complex indel. MNP, symbolic
//! and mixed records fall in no bucket. The summary row sees every record
//! without a genotype; a sample's detail row sees the records where its
//! genotype is not hom-ref, skips its no-calls, and also counts `GQ == 0`,
//! het and hom-var genotypes at unfiltered sites. `PCT_DBSNP` and
//! `PCT_DBSNP_INDELS` divide in `float`, the other ratios in `double`; a
//! Ti/Tv or ins/del ratio stays 0 on a zero denominator.
//!
//! Detail rows come out in Java `HashMap` order of the sample names
//! (`CollectionUtil.DefaultingMap`, then `Collectors.groupingBy` in
//! `Result.merge`), so `THREAD_COUNT` cannot change the file. `GVCF_INPUT`
//! replaces each record by its single sample's sub-context before counting.
//!
//! `AccumulateVariantCallingMetrics` reads shard files back, inverts the
//! printed ratios into the hidden counts (`invertFromRatio`, lossy where a
//! ratio printed 0 or `?`), adds the counts and prints the merged rows,
//! samples in `HashMap` order of first appearance.

use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use crate::java::hash_map_order;
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};
use crate::vcf::{self, DbSnpSets, Dictionary, Genotype, Regions, Variant, VariantType, VcfReader};

pub const DETAIL_EXTENSION: &str = "variant_calling_detail_metrics";
pub const SUMMARY_EXTENSION: &str = "variant_calling_summary_metrics";

const SUMMARY_COLUMNS: [&str; 20] = [
    "TOTAL_SNPS",
    "NUM_IN_DB_SNP",
    "NOVEL_SNPS",
    "FILTERED_SNPS",
    "PCT_DBSNP",
    "DBSNP_TITV",
    "NOVEL_TITV",
    "TOTAL_INDELS",
    "NOVEL_INDELS",
    "FILTERED_INDELS",
    "PCT_DBSNP_INDELS",
    "NUM_IN_DB_SNP_INDELS",
    "DBSNP_INS_DEL_RATIO",
    "NOVEL_INS_DEL_RATIO",
    "TOTAL_MULTIALLELIC_SNPS",
    "NUM_IN_DB_SNP_MULTIALLELIC",
    "TOTAL_COMPLEX_INDELS",
    "NUM_IN_DB_SNP_COMPLEX_INDELS",
    "SNP_REFERENCE_BIAS",
    "NUM_SINGLETONS",
];

const DETAIL_COLUMNS: [&str; 25] = [
    "SAMPLE_ALIAS",
    "HET_HOMVAR_RATIO",
    "PCT_GQ0_VARIANTS",
    "TOTAL_GQ0_VARIANTS",
    "TOTAL_HET_DEPTH",
    "TOTAL_SNPS",
    "NUM_IN_DB_SNP",
    "NOVEL_SNPS",
    "FILTERED_SNPS",
    "PCT_DBSNP",
    "DBSNP_TITV",
    "NOVEL_TITV",
    "TOTAL_INDELS",
    "NOVEL_INDELS",
    "FILTERED_INDELS",
    "PCT_DBSNP_INDELS",
    "NUM_IN_DB_SNP_INDELS",
    "DBSNP_INS_DEL_RATIO",
    "NOVEL_INS_DEL_RATIO",
    "TOTAL_MULTIALLELIC_SNPS",
    "NUM_IN_DB_SNP_MULTIALLELIC",
    "TOTAL_COMPLEX_INDELS",
    "NUM_IN_DB_SNP_COMPLEX_INDELS",
    "SNP_REFERENCE_BIAS",
    "NUM_SINGLETONS",
];

/// The `@MergeByAdding` fields of `VariantCallingSummaryMetrics`, printed
/// and hidden.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Counts {
    pub total_snps: i64,
    pub num_in_db_snp: i64,
    pub novel_snps: i64,
    pub filtered_snps: i64,
    pub total_indels: i64,
    pub novel_indels: i64,
    pub filtered_indels: i64,
    pub num_in_db_snp_indels: i64,
    pub total_multiallelic_snps: f64,
    pub num_in_db_snp_multiallelic: f64,
    pub total_complex_indels: f64,
    pub num_in_db_snp_complex_indels: f64,
    pub num_singletons: i64,
    pub ref_allele_obs: i64,
    pub alt_allele_obs: i64,
    pub novel_deletions: i64,
    pub novel_insertions: i64,
    pub novel_transitions: i64,
    pub novel_transversions: i64,
    pub db_snp_deletions: i64,
    pub db_snp_insertions: i64,
    pub db_snp_transitions: i64,
    pub db_snp_transversions: i64,
}

impl Counts {
    /// `MergeableMetricBase.merge`: every field added.
    pub fn add(&mut self, other: &Counts) {
        self.total_snps += other.total_snps;
        self.num_in_db_snp += other.num_in_db_snp;
        self.novel_snps += other.novel_snps;
        self.filtered_snps += other.filtered_snps;
        self.total_indels += other.total_indels;
        self.novel_indels += other.novel_indels;
        self.filtered_indels += other.filtered_indels;
        self.num_in_db_snp_indels += other.num_in_db_snp_indels;
        self.total_multiallelic_snps += other.total_multiallelic_snps;
        self.num_in_db_snp_multiallelic += other.num_in_db_snp_multiallelic;
        self.total_complex_indels += other.total_complex_indels;
        self.num_in_db_snp_complex_indels += other.num_in_db_snp_complex_indels;
        self.num_singletons += other.num_singletons;
        self.ref_allele_obs += other.ref_allele_obs;
        self.alt_allele_obs += other.alt_allele_obs;
        self.novel_deletions += other.novel_deletions;
        self.novel_insertions += other.novel_insertions;
        self.novel_transitions += other.novel_transitions;
        self.novel_transversions += other.novel_transversions;
        self.db_snp_deletions += other.db_snp_deletions;
        self.db_snp_insertions += other.db_snp_insertions;
        self.db_snp_transitions += other.db_snp_transitions;
        self.db_snp_transversions += other.db_snp_transversions;
    }

    /// `calculateDerivedFields`: the printed summary values.
    fn summary_values(&self) -> Vec<Value> {
        let ratio = |num: i64, den: i64| {
            if den > 0 {
                num as f64 / den as f64
            } else {
                0.0
            }
        };
        let pct_dbsnp = f64::from(self.num_in_db_snp as f32 / self.total_snps as f32);
        let pct_dbsnp_indels =
            f64::from(self.num_in_db_snp_indels as f32 / self.total_indels as f32);
        let bias = self.ref_allele_obs as f64 / (self.ref_allele_obs + self.alt_allele_obs) as f64;
        vec![
            self.total_snps.into(),
            self.num_in_db_snp.into(),
            (self.total_snps - self.num_in_db_snp).into(),
            self.filtered_snps.into(),
            pct_dbsnp.into(),
            ratio(self.db_snp_transitions, self.db_snp_transversions).into(),
            ratio(self.novel_transitions, self.novel_transversions).into(),
            self.total_indels.into(),
            (self.total_indels - self.num_in_db_snp_indels).into(),
            self.filtered_indels.into(),
            pct_dbsnp_indels.into(),
            self.num_in_db_snp_indels.into(),
            ratio(self.db_snp_insertions, self.db_snp_deletions).into(),
            ratio(self.novel_insertions, self.novel_deletions).into(),
            self.total_multiallelic_snps.into(),
            self.num_in_db_snp_multiallelic.into(),
            self.total_complex_indels.into(),
            self.num_in_db_snp_complex_indels.into(),
            bias.into(),
            self.num_singletons.into(),
        ]
    }
}

/// One `VariantCallingSummaryMetrics` row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SummaryRow(pub Counts);

impl MetricRow for SummaryRow {
    const CLASS: &'static str =
        "picard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics";
    const COLUMNS: &'static [&'static str] = &SUMMARY_COLUMNS;

    fn values(&self) -> Vec<Value> {
        self.0.summary_values()
    }
}

/// One `VariantCallingDetailMetrics` row: a sample's counts.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DetailRow {
    pub sample: String,
    pub counts: Counts,
    pub total_gq0_variants: i64,
    pub num_hets: i64,
    pub num_hom_var: i64,
}

impl DetailRow {
    fn add(&mut self, other: &DetailRow) {
        self.counts.add(&other.counts);
        self.total_gq0_variants += other.total_gq0_variants;
        self.num_hets += other.num_hets;
        self.num_hom_var += other.num_hom_var;
    }
}

impl MetricRow for DetailRow {
    const CLASS: &'static str =
        "picard.vcf.CollectVariantCallingMetrics$VariantCallingDetailMetrics";
    const COLUMNS: &'static [&'static str] = &DETAIL_COLUMNS;

    fn values(&self) -> Vec<Value> {
        let het_homvar = self.num_hets as f64 / self.num_hom_var as f64;
        let pct_gq0 = self.total_gq0_variants as f64 / (self.num_hets + self.num_hom_var) as f64;
        let mut values = vec![
            self.sample.clone().into(),
            het_homvar.into(),
            pct_gq0.into(),
            self.total_gq0_variants.into(),
            (self.counts.ref_allele_obs + self.counts.alt_allele_obs).into(),
        ];
        values.extend(self.counts.summary_values());
        values
    }
}

/// `CallingMetricAccumulator`: the summary and one detail per sample.
pub struct CallingMetrics {
    dbsnp: DbSnpSets,
    summary: Counts,
    details: Vec<DetailRow>,
}

impl CallingMetrics {
    pub fn new(dbsnp: DbSnpSets, samples: &[String]) -> Self {
        Self {
            dbsnp,
            summary: Counts::default(),
            details: samples
                .iter()
                .map(|s| DetailRow {
                    sample: s.clone(),
                    ..DetailRow::default()
                })
                .collect(),
        }
    }

    /// `accumulate`.
    pub fn accept(&mut self, vc: &Variant) {
        if vc.kind() == VariantType::NoVariation || vc.genotypes.iter().all(Genotype::is_hom_ref) {
            return;
        }
        let singleton = singleton_sample(vc);
        let mut bias = (0, 0);
        update(
            &mut self.summary,
            None,
            vc,
            singleton.is_some(),
            &self.dbsnp,
            &mut bias,
        );
        for (i, genotype) in vc.genotypes.iter().enumerate() {
            if genotype.is_hom_ref() {
                continue;
            }
            let detail = &mut self.details[i];
            update(
                &mut detail.counts,
                Some(genotype),
                vc,
                singleton == Some(i),
                &self.dbsnp,
                &mut bias,
            );
            // The GQ 0 count does not skip a no-call genotype.
            if !vc.filtered {
                if genotype.gq == 0 {
                    detail.total_gq0_variants += 1;
                }
                if genotype.is_het() {
                    detail.num_hets += 1;
                } else if genotype.is_hom_var() {
                    detail.num_hom_var += 1;
                }
            }
        }
        self.summary.ref_allele_obs += bias.0;
        self.summary.alt_allele_obs += bias.1;
    }

    /// The summary row and the detail rows in Picard's file order.
    pub fn finish(self) -> (SummaryRow, Vec<DetailRow>) {
        (SummaryRow(self.summary), hash_map_ordered(self.details))
    }
}

/// `getSingletonSample`: the one sample with exactly one variant
/// chromosome, when no other sample has any.
fn singleton_sample(vc: &Variant) -> Option<usize> {
    // The first two het or hom-var genotypes: a het counts one variant
    // chromosome, a hom-var two; the sum must be exactly 1.
    let mut sum = 0;
    let mut sample = None;
    let mut seen = 0;
    for (i, g) in vc.genotypes.iter().enumerate() {
        let chromosomes = if g.is_het() {
            1
        } else if g.is_hom_var() {
            2
        } else {
            continue;
        };
        sample = Some(i);
        sum += chromosomes;
        seen += 1;
        if seen == 2 {
            break;
        }
    }
    if sum == 1 { sample } else { None }
}

/// `updateSummaryMetric` for the summary (`genotype` None) or a detail
/// (`genotype` the sample's). `bias` collects the het `AD` observations the
/// summary adds once per record.
fn update(
    metric: &mut Counts,
    genotype: Option<&Genotype>,
    vc: &Variant,
    has_singleton: bool,
    dbsnp: &DbSnpSets,
    bias: &mut (i64, i64),
) {
    if genotype.is_some_and(Genotype::is_no_call) {
        return;
    }
    let kind = vc.kind();
    if vc.filtered {
        if kind == VariantType::Snp {
            metric.filtered_snps += 1;
        } else if kind == VariantType::Indel {
            metric.filtered_indels += 1;
        }
        return;
    }
    if has_singleton {
        metric.num_singletons += 1;
    }
    if vc.is_biallelic() && kind == VariantType::Snp {
        let in_dbsnp = dbsnp.snps.is_site(&vc.contig, vc.start);
        let transition = vc.is_transition();
        metric.total_snps += 1;
        match (in_dbsnp, transition) {
            (true, true) => {
                metric.num_in_db_snp += 1;
                metric.db_snp_transitions += 1;
            }
            (true, false) => {
                metric.num_in_db_snp += 1;
                metric.db_snp_transversions += 1;
            }
            (false, true) => metric.novel_transitions += 1,
            (false, false) => metric.novel_transversions += 1,
        }
        if let Some(g) = genotype
            && g.is_het()
            && let Some(ad) = &g.ad
        {
            // `alleleDepths[0]` and `[1]`: an AD with one value is Picard's
            // ArrayIndexOutOfBoundsException.
            let (r, a) = (ad[0], ad[1]);
            metric.ref_allele_obs += r;
            metric.alt_allele_obs += a;
            bias.0 += r;
            bias.1 += a;
        }
    } else if kind == VariantType::Snp && vc.alternates().len() > 1 {
        metric.total_multiallelic_snps += 1.0;
        if dbsnp.snps.is_site(&vc.contig, vc.start) {
            metric.num_in_db_snp_multiallelic += 1.0;
        }
    } else if kind == VariantType::Indel && !vc.is_complex_indel() {
        let in_dbsnp = dbsnp.indels.is_site(&vc.contig, vc.start);
        let insertion = vc.is_simple_insertion();
        metric.total_indels += 1;
        match (in_dbsnp, insertion) {
            (true, true) => {
                metric.num_in_db_snp_indels += 1;
                metric.db_snp_insertions += 1;
            }
            (true, false) => {
                metric.num_in_db_snp_indels += 1;
                metric.db_snp_deletions += 1;
            }
            (false, true) => metric.novel_insertions += 1,
            (false, false) => metric.novel_deletions += 1,
        }
    } else if vc.is_complex_indel() {
        metric.total_complex_indels += 1.0;
        if dbsnp.indels.is_site(&vc.contig, vc.start) {
            metric.num_in_db_snp_complex_indels += 1.0;
        }
    }
}

/// Detail rows in the iteration order of a Java `HashMap` keyed by sample
/// name, filled in the given order.
fn hash_map_ordered(details: Vec<DetailRow>) -> Vec<DetailRow> {
    let names: Vec<&str> = details.iter().map(|d| d.sample.as_str()).collect();
    let order = hash_map_order(&names);
    let mut slots: Vec<Option<DetailRow>> = details.into_iter().map(Some).collect();
    order
        .into_iter()
        .map(|i| slots[i].take().expect("each index once"))
        .collect()
}

/// `CollectVariantCallingMetrics` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CollectVariantCallingMetrics",
    options: &[
        Opt::value(&["DBSNP"]),
        Opt::value(&["TARGET_INTERVALS", "TI"]),
        Opt::value(&["SEQUENCE_DICTIONARY", "SD"]),
        Opt::value(&["GVCF_INPUT"]),
        Opt::value(&["THREAD_COUNT"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

/// `AccumulateVariantCallingMetrics` in the tool registry.
pub const ACCUMULATE_TOOL: Tool = Tool {
    name: "AccumulateVariantCallingMetrics",
    options: &[Opt::rejected(&["ASSUME_SORTED", "AS"])],
    run: run_accumulate,
    collector: None,
};

pub struct Options {
    pub prefix: PathBuf,
    pub dbsnp: PathBuf,
    pub target_intervals: Option<PathBuf>,
    pub sequence_dictionary: Option<PathBuf>,
    pub gvcf_input: bool,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let required = |name: &str, value: Option<&str>| {
            value
                .map(str::to_owned)
                .ok_or_else(|| Error::Usage(format!("Option '{name}' is required.")))
        };
        let prefix = required("OUTPUT", inv.output.as_deref().and_then(Path::to_str))?;
        let threads: i64 = inv.parse("THREAD_COUNT", 1)?;
        if threads < 1 {
            return Err(Error::Usage(
                "Multithreading value must exceed 0.".to_owned(),
            ));
        }
        Ok(Self {
            prefix: PathBuf::from(prefix),
            dbsnp: PathBuf::from(required("DBSNP", inv.get("DBSNP"))?),
            target_intervals: inv.get("TARGET_INTERVALS").map(PathBuf::from),
            sequence_dictionary: inv.get("SEQUENCE_DICTIONARY").map(PathBuf::from),
            gvcf_input: inv.flag("GVCF_INPUT", false)?,
        })
    }
}

/// `<prefix>.<extension>`, the `OUTPUT.getAbsolutePath() + "."` rule.
pub fn output_path(prefix: &Path, extension: &str) -> PathBuf {
    PathBuf::from(format!("{}.{extension}", prefix.display()))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let (summary, details) = collect(&inv.input, &opts)?;
    Ok(files(&opts.prefix, &summary, &details))
}

/// Both files of the tool.
fn files(prefix: &Path, summary: &SummaryRow, details: &[DetailRow]) -> ToolOutput {
    ToolOutput {
        main: None,
        extra: vec![
            (
                output_path(prefix, DETAIL_EXTENSION),
                MetricsOutput::from_rows(details, None),
            ),
            (
                output_path(prefix, SUMMARY_EXTENSION),
                MetricsOutput::from_rows(std::slice::from_ref(summary), None),
            ),
        ],
    }
}

/// `doWork` up to the rows: the dictionary, the target intervals, the dbSNP
/// bit sets, then every record of `input` on a dictionary contig (and
/// overlapping a target interval, when given).
pub fn collect(input: &Path, opts: &Options) -> Result<(SummaryRow, Vec<DetailRow>)> {
    let mut reader = VcfReader::open(input)?;
    let samples = reader.samples();
    let header_dictionary: Dictionary = reader.dictionary().ok_or_else(|| Error::Intervals {
        path: input.to_path_buf(),
        message: "no ##contig line (Picard's VariantProcessor fails with a NullPointerException)"
            .to_owned(),
    })?;
    let dictionary = match &opts.sequence_dictionary {
        Some(path) => vcf::read_dictionary(path)?,
        None => header_dictionary.clone(),
    };
    let regions = opts
        .target_intervals
        .as_deref()
        .map(|p| Regions::read(p).map(|r| r.uniqued()))
        .transpose()?;
    let dbsnp = DbSnpSets::load(&opts.dbsnp, Some(&dictionary), regions.as_ref())?;
    if opts.gvcf_input && samples.len() != 1 {
        return Err(Error::Usage(format!(
            "Expected to have exactly 1 sample in a GVCF, found {}",
            samples.len()
        )));
    }
    let mut metrics = CallingMetrics::new(dbsnp, &samples);
    while let Some(vc) = reader.next_variant()? {
        if !header_dictionary.iter().any(|c| c.name == vc.contig) {
            continue;
        }
        if let Some(regions) = &regions
            && !regions.overlaps(&vc.contig, vc.start, vc.end)
        {
            continue;
        }
        if opts.gvcf_input {
            metrics.accept(&vc.sub_context(0));
        } else {
            metrics.accept(&vc);
        }
    }
    Ok(metrics.finish())
}

/// `Math.round(double)`.
fn java_round(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

/// `invertFromRatio`: Y from X/Y and X+Y; 0 when the ratio is NaN.
fn invert_from_ratio(sum: i64, ratio: f64) -> i64 {
    if ratio.is_nan() {
        0
    } else {
        java_round(sum as f64 / (ratio + 1.0))
    }
}

/// A row of a shard file by column name.
struct Fields<'a> {
    columns: &'a [String],
    row: &'a [String],
}

impl Fields<'_> {
    fn text(&self, name: &str) -> Result<&str> {
        self.columns
            .iter()
            .position(|c| c == name)
            .map(|i| self.row[i].as_str())
            .ok_or_else(|| Error::Metric(format!("column {name} is missing")))
    }

    fn long(&self, name: &str) -> Result<i64> {
        vcf::parse_metric_long(self.text(name)?)
    }

    fn double(&self, name: &str) -> Result<f64> {
        vcf::parse_metric_double(self.text(name)?)
    }

    /// `VariantCallingSummaryMetrics.calculateFromDerivedFields`.
    fn counts(&self, total_het_depth: i64) -> Result<Counts> {
        let mut c = Counts {
            total_snps: self.long("TOTAL_SNPS")?,
            num_in_db_snp: self.long("NUM_IN_DB_SNP")?,
            novel_snps: self.long("NOVEL_SNPS")?,
            filtered_snps: self.long("FILTERED_SNPS")?,
            total_indels: self.long("TOTAL_INDELS")?,
            novel_indels: self.long("NOVEL_INDELS")?,
            filtered_indels: self.long("FILTERED_INDELS")?,
            num_in_db_snp_indels: self.long("NUM_IN_DB_SNP_INDELS")?,
            total_multiallelic_snps: self.double("TOTAL_MULTIALLELIC_SNPS")?,
            num_in_db_snp_multiallelic: self.double("NUM_IN_DB_SNP_MULTIALLELIC")?,
            total_complex_indels: self.double("TOTAL_COMPLEX_INDELS")?,
            num_in_db_snp_complex_indels: self.double("NUM_IN_DB_SNP_COMPLEX_INDELS")?,
            num_singletons: self.long("NUM_SINGLETONS")?,
            ..Counts::default()
        };
        c.db_snp_transversions = invert_from_ratio(c.num_in_db_snp, self.double("DBSNP_TITV")?);
        c.db_snp_transitions = c.num_in_db_snp - c.db_snp_transversions;
        c.novel_transversions = invert_from_ratio(c.novel_snps, self.double("NOVEL_TITV")?);
        c.novel_transitions = c.novel_snps - c.novel_transversions;
        c.db_snp_deletions =
            invert_from_ratio(c.num_in_db_snp_indels, self.double("DBSNP_INS_DEL_RATIO")?);
        c.db_snp_insertions = c.num_in_db_snp_indels - c.db_snp_deletions;
        c.novel_deletions = invert_from_ratio(c.novel_indels, self.double("NOVEL_INS_DEL_RATIO")?);
        c.novel_insertions = c.novel_indels - c.novel_deletions;
        let bias = self.double("SNP_REFERENCE_BIAS")?;
        c.ref_allele_obs = if bias.is_nan() {
            0
        } else {
            java_round(total_het_depth as f64 * bias)
        };
        c.alt_allele_obs = total_het_depth - c.ref_allele_obs;
        Ok(c)
    }

    /// `VariantCallingDetailMetrics.calculateFromDerivedFields`.
    fn detail(&self) -> Result<DetailRow> {
        let total_snps = self.long("TOTAL_SNPS")?;
        let num_hom_var = invert_from_ratio(total_snps, self.double("HET_HOMVAR_RATIO")?);
        let total_het_depth = self.long("TOTAL_HET_DEPTH")?;
        Ok(DetailRow {
            sample: self.text("SAMPLE_ALIAS")?.to_owned(),
            counts: self.counts(total_het_depth)?,
            total_gq0_variants: self.long("TOTAL_GQ0_VARIANTS")?,
            num_hets: total_snps - num_hom_var,
            num_hom_var,
        })
    }
}

fn run_accumulate(inv: &Invocation) -> Result<ToolOutput> {
    let prefix = inv
        .output
        .clone()
        .ok_or_else(|| Error::Usage("Option 'OUTPUT' is required.".to_owned()))?;
    let inputs: Vec<PathBuf> = inv.all("INPUT").iter().map(PathBuf::from).collect();
    let (summary, details) = accumulate(&inputs)?;
    Ok(files(&prefix, &summary, &details))
}

/// `AccumulateVariantCallingMetrics.doWork`: each shard's detail rows
/// merged by sample (first appearance order into a `HashMap`), its summary
/// merged into one.
pub fn accumulate(prefixes: &[PathBuf]) -> Result<(SummaryRow, Vec<DetailRow>)> {
    let mut details: Vec<DetailRow> = Vec::new();
    let mut summary = Counts::default();
    for prefix in prefixes {
        let detail_path = output_path(prefix, DETAIL_EXTENSION);
        let (columns, rows) = vcf::read_metrics_rows(&detail_path)?;
        let mut total_het_depth = 0;
        for row in &rows {
            let fields = Fields {
                columns: &columns,
                row,
            };
            let detail = fields.detail()?;
            total_het_depth += detail.counts.ref_allele_obs + detail.counts.alt_allele_obs;
            match details.iter_mut().find(|d| d.sample == detail.sample) {
                Some(existing) => existing.add(&detail),
                None => {
                    let mut fresh = DetailRow {
                        sample: detail.sample.clone(),
                        ..DetailRow::default()
                    };
                    fresh.add(&detail);
                    details.push(fresh);
                }
            }
        }
        let summary_path = output_path(prefix, SUMMARY_EXTENSION);
        let (columns, rows) = vcf::read_metrics_rows(&summary_path)?;
        if rows.len() != 1 {
            return Err(Error::Metric(format!(
                "Expected 1 row in the summary metrics file but saw {}",
                rows.len()
            )));
        }
        let fields = Fields {
            columns: &columns,
            row: &rows[0],
        };
        summary.add(&fields.counts(total_het_depth)?);
    }
    Ok((SummaryRow(summary), hash_map_ordered(details)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::Allele;

    fn genotype(alleles: &[Option<usize>], gq: i64, ad: Option<Vec<i64>>) -> Genotype {
        Genotype {
            alleles: alleles.to_vec(),
            gq,
            dp: -1,
            ad,
            filtered: false,
        }
    }

    fn site(reference: &str, alts: &[&str], genotypes: Vec<Genotype>, filtered: bool) -> Variant {
        let mut alleles = vec![Allele {
            bases: reference.to_owned(),
            symbolic: false,
        }];
        alleles.extend(alts.iter().map(|a| Allele {
            bases: (*a).to_owned(),
            symbolic: a.starts_with('<'),
        }));
        Variant {
            contig: "c1".to_owned(),
            start: 100,
            end: 100 + reference.len() as i64 - 1,
            alleles,
            filtered,
            genotypes,
        }
    }

    #[test]
    fn singletons_and_buckets() {
        let samples = ["a".to_owned(), "b".to_owned()];
        let mut m = CallingMetrics::new(DbSnpSets::default(), &samples);
        let het = |ad| genotype(&[Some(0), Some(1)], 0, ad);
        let hom_ref = || genotype(&[Some(0), Some(0)], 30, None);
        // A het singleton SNP (transition) with AD in sample a.
        m.accept(&site(
            "A",
            &["G"],
            vec![het(Some(vec![3, 5])), hom_ref()],
            false,
        ));
        // Both samples variant: no singleton; transversion.
        m.accept(&site(
            "A",
            &["C"],
            vec![het(None), het(Some(vec![1, 1]))],
            false,
        ));
        // A filtered indel: only the filtered column, and no singleton.
        m.accept(&site("A", &["AT"], vec![het(None), hom_ref()], true));
        // A no-call sample is skipped; the other is a hom-var singleton? No:
        // hom-var counts 2 chromosomes, so no singleton.
        m.accept(&site(
            "AT",
            &["A"],
            vec![
                genotype(&[None, None], 0, None),
                genotype(&[Some(1), Some(1)], 0, None),
            ],
            false,
        ));
        // Every genotype hom-ref: excluded.
        m.accept(&site("A", &["G"], vec![hom_ref(), hom_ref()], false));
        // Multiallelic SNP, complex indel, MNP and mixed.
        m.accept(&site(
            "A",
            &["G", "T"],
            vec![genotype(&[Some(1), Some(2)], 0, None), hom_ref()],
            false,
        ));
        m.accept(&site("AT", &["GCC"], vec![het(None), hom_ref()], false));
        m.accept(&site("AT", &["GC"], vec![het(None), hom_ref()], false));
        m.accept(&site("A", &["G", "AT"], vec![het(None), hom_ref()], false));
        let (summary, details) = m.finish();
        let s = &summary.0;
        assert_eq!(
            (s.total_snps, s.novel_transitions, s.novel_transversions),
            (2, 1, 1)
        );
        assert_eq!(
            (s.filtered_indels, s.total_indels, s.novel_deletions),
            (1, 1, 1)
        );
        assert_eq!(
            (s.total_multiallelic_snps, s.total_complex_indels),
            (1.0, 1.0)
        );
        // Singletons: the first SNP, the multiallelic SNP, the complex indel,
        // the MNP and the mixed site.
        assert_eq!(s.num_singletons, 5);
        assert_eq!((s.ref_allele_obs, s.alt_allele_obs), (4, 6));
        let a = details.iter().find(|d| d.sample == "a").expect("a");
        let b = details.iter().find(|d| d.sample == "b").expect("b");
        assert_eq!((a.counts.num_singletons, b.counts.num_singletons), (5, 0));
        // a's GQ 0 count includes its no-call at the deletion.
        assert_eq!((a.num_hets, a.num_hom_var, a.total_gq0_variants), (6, 0, 7));
        assert_eq!((b.num_hets, b.num_hom_var, b.total_gq0_variants), (1, 1, 2));
        assert_eq!(a.counts.ref_allele_obs + a.counts.alt_allele_obs, 8);
        assert_eq!(a.values()[1].picard_text(), "?");
        assert_eq!(summary.values()[4].picard_text(), "0");
        assert_eq!(summary.values()[18].picard_text(), "0.4");
    }

    #[test]
    fn pct_dbsnp_divides_in_float() {
        let counts = Counts {
            total_snps: 11,
            num_in_db_snp: 6,
            ..Counts::default()
        };
        // 6/11 in float is 0.54545456, in double 0.5454545454.
        assert_eq!(counts.summary_values()[4].picard_text(), "0.545455");
        let counts = Counts {
            total_snps: 3,
            num_in_db_snp: 1,
            ..Counts::default()
        };
        assert_eq!(counts.summary_values()[4].picard_text(), "0.333333");
    }

    #[test]
    fn ratio_inversion_matches_picard() {
        assert_eq!(invert_from_ratio(10, f64::NAN), 0);
        assert_eq!(invert_from_ratio(10, 1.5), 4);
        assert_eq!(invert_from_ratio(7, 0.0), 7);
        assert_eq!(java_round(2.5), 3);
        assert_eq!(java_round(-2.5), -2);
    }

    #[test]
    fn detail_rows_come_in_hash_map_order() {
        let details: Vec<DetailRow> = ["L2", "L1", "Unknown", "L3"]
            .iter()
            .map(|s| DetailRow {
                sample: (*s).to_owned(),
                ..DetailRow::default()
            })
            .collect();
        let names: Vec<String> = hash_map_ordered(details)
            .into_iter()
            .map(|d| d.sample)
            .collect();
        assert_eq!(names, ["Unknown", "L1", "L2", "L3"]);
    }
}
