//! `CollectOxoGMetrics`: Picard 3.5.0's count of 8-oxoguanine damage, one
//! `CpcgMetrics` row per sequence context and library.
//!
//! Picard walks the BAM with a `SamLocusIterator` ([`crate::locus`],
//! [`Setup::OxoG`]): secondary and duplicate records are dropped, so are
//! records under `MINIMUM_MAPPING_QUALITY`, non-PF records when
//! `INCLUDE_NON_PF_READS=false`, and, unless both insert size bounds are
//! 0, unpaired records and records whose `|TLEN|` is outside the bounds.
//! Only covered loci come out. At each locus, in order:
//! * a position `DB_SNP` marks (any record type, start to end) is skipped;
//! * so is one within `CONTEXT_SIZE` bases of either contig end, and one
//!   whose reference base, upper-cased, is not C or G;
//! * the context is the `2 * CONTEXT_SIZE + 1` reference bases around it,
//!   upper-cased, reverse-complemented at a G; a context outside the set
//!   (an N in it, or one `CONTEXTS` leaves out) is skipped;
//! * each library's calculator for the context reads every base at the
//!   locus: quality from `OQ` when `USE_OQ` and the read has one, else
//!   `QUAL`; a base under `MINIMUM_QUALITY_SCORE` or of another library is
//!   passed over; a reference base (C or G) or alternate base (A at a C, T
//!   at a G) counts as oxidized when the base as sequenced (complemented on
//!   the reverse strand) is G or T in read 1 or C or A in read 2, and as a
//!   control otherwise. A calculator with any count adds one site;
//! * the site counts towards `STOP_AFTER` whether or not a base counted,
//!   and the walk stops once `STOP_AFTER` sites are reached (at the first
//!   site when it is 0 or less).
//!
//! Rows come out in the `HashMap` order of the contexts ([`crate::java`]),
//! each context's libraries in `HashSet` order; `SAMPLE_ALIAS` joins every
//! `@RG` sample in `HashSet` order. A read group without `SM` or `LB` is
//! `UnknownSample` or `UnknownLibrary`. A counted base on a record without a
//! header read group is an error, as Picard's `NullPointerException`.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use noodles::sam;

use crate::error::{Error, Result};
use crate::intervals::{self, uniqued};
use crate::java::{StringHashSet, hash_map_order_checked, string_hash_code};
use crate::levels::read_groups;
use crate::locus::{BaseEntry, Classic, Engine, Finished, Settings, Setup, Sink, check_sort_order};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool};
use crate::sequence::{complement, reverse_complement};
use crate::vcf::{Contig, DbSnpBitSet, Dictionary};

const UNKNOWN_LIBRARY: &str = "UnknownLibrary";
const UNKNOWN_SAMPLE: &str = "UnknownSample";
/// The largest `CONTEXT_SIZE` generated when `CONTEXTS` is empty: 65,536
/// contexts. From 5 on (1,048,576 contexts) some `HashMap` bins of
/// Picard's calculator map turn into trees and its row order follows them
/// ([`hash_order`]).
const MAX_GENERATED_CONTEXT_SIZE: i32 = 4;

/// `CollectOxoGMetrics$CpcgMetrics`, one per context and library.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub sample_alias: String,
    pub library: String,
    pub context: String,
    pub total_sites: i64,
    pub total_bases: i64,
    pub ref_nonoxo_bases: i64,
    pub ref_oxo_bases: i64,
    pub ref_total_bases: i64,
    pub alt_nonoxo_bases: i64,
    pub alt_oxo_bases: i64,
    pub oxidation_error_rate: f64,
    pub oxidation_q: f64,
    pub c_ref_ref_bases: i64,
    pub g_ref_ref_bases: i64,
    pub c_ref_alt_bases: i64,
    pub g_ref_alt_bases: i64,
    pub c_ref_oxo_error_rate: f64,
    pub c_ref_oxo_q: f64,
    pub g_ref_oxo_error_rate: f64,
    pub g_ref_oxo_q: f64,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.analysis.CollectOxoGMetrics$CpcgMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "SAMPLE_ALIAS",
        "LIBRARY",
        "CONTEXT",
        "TOTAL_SITES",
        "TOTAL_BASES",
        "REF_NONOXO_BASES",
        "REF_OXO_BASES",
        "REF_TOTAL_BASES",
        "ALT_NONOXO_BASES",
        "ALT_OXO_BASES",
        "OXIDATION_ERROR_RATE",
        "OXIDATION_Q",
        "C_REF_REF_BASES",
        "G_REF_REF_BASES",
        "C_REF_ALT_BASES",
        "G_REF_ALT_BASES",
        "C_REF_OXO_ERROR_RATE",
        "C_REF_OXO_Q",
        "G_REF_OXO_ERROR_RATE",
        "G_REF_OXO_Q",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.sample_alias.as_str().into(),
            self.library.as_str().into(),
            self.context.as_str().into(),
            self.total_sites.into(),
            self.total_bases.into(),
            self.ref_nonoxo_bases.into(),
            self.ref_oxo_bases.into(),
            self.ref_total_bases.into(),
            self.alt_nonoxo_bases.into(),
            self.alt_oxo_bases.into(),
            self.oxidation_error_rate.into(),
            self.oxidation_q.into(),
            self.c_ref_ref_bases.into(),
            self.g_ref_ref_bases.into(),
            self.c_ref_alt_bases.into(),
            self.g_ref_alt_bases.into(),
            self.c_ref_oxo_error_rate.into(),
            self.c_ref_oxo_q.into(),
            self.g_ref_oxo_error_rate.into(),
            self.g_ref_oxo_q.into(),
        ]
    }
}

/// `Math.max(double, double)`: NaN when either is NaN.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// `Counts`: one locus, one library.
#[derive(Clone, Copy, Debug, Default)]
struct Counts {
    control_a: i64,
    oxidated_a: i64,
    control_c: i64,
    oxidated_c: i64,
}

impl Counts {
    fn total(&self) -> i64 {
        self.control_c + self.oxidated_c + self.control_a + self.oxidated_a
    }
}

/// `Calculator`: the sums of one context and library, split by reference
/// base.
#[derive(Clone, Debug, Default)]
struct Calculator {
    sites: i64,
    ref_c: Counts,
    ref_g: Counts,
}

impl Calculator {
    /// `Calculator.finish`.
    fn finish(&self, sample_alias: &str, library: &str, context: &str) -> Row {
        let (c, g) = (&self.ref_c, &self.ref_g);
        let total_bases = c.total() + g.total();
        let ref_oxo_bases = c.oxidated_c + g.oxidated_c;
        let ref_nonoxo_bases = c.control_c + g.control_c;
        let alt_nonoxo_bases = c.control_a + g.control_a;
        let alt_oxo_bases = c.oxidated_a + g.oxidated_a;
        let oxidation_error_rate =
            (alt_oxo_bases - alt_nonoxo_bases).max(1) as f64 / total_bases as f64;
        let c_ref_ref_bases = c.control_c + c.oxidated_c;
        let g_ref_ref_bases = g.control_c + g.oxidated_c;
        let c_ref_alt_bases = c.control_a + c.oxidated_a;
        let g_ref_alt_bases = g.control_a + g.oxidated_a;
        let c_rate = c_ref_alt_bases as f64 / (c_ref_alt_bases + c_ref_ref_bases) as f64;
        let g_rate = g_ref_alt_bases as f64 / (g_ref_alt_bases + g_ref_ref_bases) as f64;
        let c_ref_oxo_error_rate = java_max(c_rate - g_rate, 1e-10);
        let g_ref_oxo_error_rate = java_max(g_rate - c_rate, 1e-10);
        Row {
            sample_alias: sample_alias.to_owned(),
            library: library.to_owned(),
            context: context.to_owned(),
            total_sites: self.sites,
            total_bases,
            ref_nonoxo_bases,
            ref_oxo_bases,
            ref_total_bases: ref_oxo_bases + ref_nonoxo_bases,
            alt_nonoxo_bases,
            alt_oxo_bases,
            oxidation_error_rate,
            oxidation_q: -10.0 * oxidation_error_rate.log10(),
            c_ref_ref_bases,
            g_ref_ref_bases,
            c_ref_alt_bases,
            g_ref_alt_bases,
            c_ref_oxo_error_rate,
            c_ref_oxo_q: -10.0 * c_ref_oxo_error_rate.log10(),
            g_ref_oxo_error_rate,
            g_ref_oxo_q: -10.0 * g_ref_oxo_error_rate.log10(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub reference: PathBuf,
    pub intervals: Option<PathBuf>,
    pub db_snp: Option<PathBuf>,
    pub minimum_quality_score: i32,
    pub minimum_mapping_quality: i32,
    pub minimum_insert_size: i32,
    pub maximum_insert_size: i32,
    pub include_non_pf_reads: bool,
    pub use_oq: bool,
    pub context_size: i32,
    /// `CONTEXTS` as given, repeats included.
    pub contexts: Vec<String>,
    pub stop_after: i32,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let reference = inv.get("REFERENCE_SEQUENCE").ok_or_else(|| {
            Error::Usage("CollectOxoGMetrics: REFERENCE_SEQUENCE is required".to_owned())
        })?;
        let opts = Self {
            reference: reference.into(),
            intervals: inv.get("INTERVALS").map(PathBuf::from),
            db_snp: inv.get("DB_SNP").map(PathBuf::from),
            minimum_quality_score: inv.parse("MINIMUM_QUALITY_SCORE", 20)?,
            minimum_mapping_quality: inv.parse("MINIMUM_MAPPING_QUALITY", 30)?,
            minimum_insert_size: inv.parse("MINIMUM_INSERT_SIZE", 60)?,
            maximum_insert_size: inv.parse("MAXIMUM_INSERT_SIZE", 600)?,
            include_non_pf_reads: inv.flag("INCLUDE_NON_PF_READS", true)?,
            use_oq: inv.flag("USE_OQ", true)?,
            context_size: inv.parse("CONTEXT_SIZE", 1)?,
            contexts: inv.list("CONTEXTS", Vec::new(), |s| Ok(s.to_owned()))?,
            stop_after: inv.parse("STOP_AFTER", i32::MAX)?,
        };
        opts.validate()?;
        Ok(opts)
    }

    /// `customCommandLineValidation`, then the two failures Picard meets
    /// later: a negative `CONTEXT_SIZE` (a negative array size) and one too
    /// large to generate.
    fn validate(&self) -> Result<()> {
        let size = 1 + 2 * i64::from(self.context_size);
        let mut messages = Vec::new();
        for context in &self.contexts {
            let chars: Vec<char> = context.chars().collect();
            if chars.len() as i64 != size {
                messages.push(format!(
                    "Context {context} is not {size} long as implied by CONTEXT_SIZE={}",
                    self.context_size
                ));
            } else if chars[chars.len() / 2] != 'C' {
                messages.push(format!(
                    "Middle base of context sequence {context} must be C"
                ));
            }
        }
        if self.minimum_insert_size < 0 {
            messages.push("MINIMUM_INSERT_SIZE cannot be negative".to_owned());
        }
        if self.maximum_insert_size < 0 {
            messages.push("MAXIMUM_INSERT_SIZE cannot be negative".to_owned());
        }
        if self.maximum_insert_size < self.minimum_insert_size {
            messages.push("MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE".to_owned());
        }
        if !messages.is_empty() {
            return Err(Error::Usage(messages.join("\n")));
        }
        if self.contexts.is_empty() && self.context_size < 0 {
            return Err(Error::Usage(format!(
                "CONTEXT_SIZE={}: cannot be negative (Picard fails with a NegativeArraySizeException)",
                self.context_size
            )));
        }
        if self.contexts.is_empty() && self.context_size > MAX_GENERATED_CONTEXT_SIZE {
            return Err(Error::Usage(format!(
                "CONTEXT_SIZE={}: at most {MAX_GENERATED_CONTEXT_SIZE} without CONTEXTS \
                 (Picard's row order then depends on HashMap tree bins, which bammetrics \
                 does not reproduce)",
                self.context_size
            )));
        }
        Ok(())
    }

    /// `InsertSizeFilter` exists unless both bounds are 0.
    fn insert_size(&self) -> Option<(i32, i32)> {
        (self.minimum_insert_size > 0 || self.maximum_insert_size > 0)
            .then_some((self.minimum_insert_size, self.maximum_insert_size))
    }

    /// The contexts in `HashMap` iteration order: `CONTEXTS` as a
    /// `HashSet`, or `makeContextStrings` (every k-mer with C in the
    /// middle, added in `generateAllKmers` order, A C G T at each position).
    /// The `ListMap` of calculators takes them in that order and has the
    /// same table size, so it iterates them the same way.
    fn context_order(&self) -> Result<Vec<String>> {
        let keys: Vec<String> = if self.contexts.is_empty() {
            let k = usize::try_from(self.context_size).unwrap_or(0);
            let flanks = 2 * k as u32;
            (0..4u64.pow(flanks))
                .map(|n| {
                    let mut kmer = String::with_capacity(2 * k + 1);
                    for i in (0..flanks).rev() {
                        kmer.push(['A', 'C', 'G', 'T'][((n >> (2 * i)) & 3) as usize]);
                    }
                    kmer.insert(k, 'C');
                    kmer
                })
                .collect()
        } else {
            let mut seen = HashSet::new();
            self.contexts
                .iter()
                .filter(|c| seen.insert(c.as_str()))
                .cloned()
                .collect()
        };
        let set = hash_order(&keys)?;
        hash_order(&set)
    }
}

/// `keys` (no repeats, in insertion order) in the order a `HashSet` or
/// `HashMap` of them iterates. A table where a bin turns into a tree is an
/// error: Picard's row order then follows the tree, which is not modelled.
fn hash_order(keys: &[String]) -> Result<Vec<String>> {
    let hashes: Vec<i32> = keys.iter().map(|k| string_hash_code(k)).collect();
    let order = hash_map_order_checked(&hashes).ok_or_else(|| {
        Error::Usage(format!(
            "{} contexts: Picard's row order then depends on HashMap tree bins, which \
             bammetrics does not reproduce; give fewer CONTEXTS or a smaller CONTEXT_SIZE",
            keys.len()
        ))
    })?;
    Ok(order.into_iter().map(|i| keys[i].clone()).collect())
}

/// `new HashSet<>()` filled with each `@RG` value, `unknown` for a missing
/// one, in header order; the keys in iteration order.
fn hash_set_order<'a>(values: impl Iterator<Item = Option<&'a str>>, unknown: &str) -> Vec<String> {
    let mut set = StringHashSet::new();
    for value in values {
        set.insert(value.unwrap_or(unknown));
    }
    set.iter().into_iter().map(str::to_owned).collect()
}

/// The locus loop of `doWork` with every calculator.
struct OxoGSink {
    minimum_quality_score: i32,
    use_oq: bool,
    context_size: usize,
    stop_after: i32,
    contig_names: Vec<String>,
    db_snp: Option<DbSnpBitSet>,
    /// Context bytes to its index in output order.
    contexts: HashMap<Vec<u8>, usize>,
    library_of_group: Vec<usize>,
    libraries: usize,
    /// `[context * libraries + library]`.
    calculators: Vec<Calculator>,
    /// Per-library counts of the current locus.
    scratch: Vec<Counts>,
    window: Vec<u8>,
    sites: i32,
    stopped: bool,
    /// Libraries in `HashSet` order, the row order within a context.
    library_names: Vec<String>,
    sample_alias: String,
}

impl OxoGSink {
    /// One row per calculator: contexts in the order given, libraries
    /// within each.
    fn rows(&self, contexts: &[String]) -> Vec<Row> {
        let mut rows = Vec::with_capacity(self.calculators.len());
        for (c, context) in contexts.iter().enumerate() {
            for (l, library) in self.library_names.iter().enumerate() {
                let calculator = &self.calculators[c * self.libraries + l];
                rows.push(calculator.finish(&self.sample_alias, library, context));
            }
        }
        rows
    }

    /// `computeAlleleFraction` for every library at once: the same tests
    /// in the same order, so a base that makes Picard throw makes this fail.
    fn count(&mut self, entries: &[BaseEntry], reference_base: u8) -> Result<()> {
        let alternate = if reference_base == b'C' { b'A' } else { b'T' };
        self.scratch.iter_mut().for_each(|c| *c = Counts::default());
        for entry in entries {
            let read = &entry.read;
            let original = if self.use_oq {
                read.original_quality(entry.offset)?
            } else {
                None
            };
            let quality = match original {
                Some(q) => q,
                None => read.quality(entry.offset)?,
            };
            if i32::from(quality) < self.minimum_quality_score {
                continue;
            }
            let library = read
                .read_group
                .and_then(|g| self.library_of_group.get(g).copied())
                .ok_or_else(|| {
                    read.malformed(
                        "no read group in the header (Picard fails with a NullPointerException)",
                    )
                })?;
            let base = read.base(entry.offset)?;
            let as_read = if read.flags.is_reverse_complemented() {
                complement(base)
            } else {
                base
            };
            let read_two = read.flags.is_segmented() && read.flags.is_last_segment();
            let counts = &mut self.scratch[library];
            if base == reference_base {
                match (as_read, read_two) {
                    (b'G', false) | (b'C', true) => counts.oxidated_c += 1,
                    (b'G', true) | (b'C', false) => counts.control_c += 1,
                    _ => {}
                }
            } else if base == alternate {
                match (as_read, read_two) {
                    (b'T', false) | (b'A', true) => counts.oxidated_a += 1,
                    (b'T', true) | (b'A', false) => counts.control_a += 1,
                    _ => {}
                }
            }
        }
        Ok(())
    }
}

impl Sink<BaseEntry> for OxoGSink {
    fn locus(
        &mut self,
        contig: usize,
        position: i64,
        entries: &[BaseEntry],
        reference: &[u8],
        _reference_n: bool,
        _counter: i64,
    ) -> Result<()> {
        if let Some(db_snp) = &self.db_snp
            && db_snp.is_site(&self.contig_names[contig], position)
        {
            return Ok(());
        }
        let k = self.context_size as i64;
        if position <= k || position > reference.len() as i64 - k {
            return Ok(());
        }
        let index = (position - 1) as usize;
        let base = reference[index].to_ascii_uppercase();
        if base != b'C' && base != b'G' {
            return Ok(());
        }
        self.window.clear();
        self.window.extend(
            reference[index - self.context_size..=index + self.context_size]
                .iter()
                .map(u8::to_ascii_uppercase),
        );
        if base == b'G' {
            self.window = reverse_complement(&self.window);
        }
        let Some(&context) = self.contexts.get(&self.window) else {
            return Ok(());
        };
        self.count(entries, base)?;
        for (library, counts) in self.scratch.iter().enumerate() {
            if counts.total() > 0 {
                let calculator = &mut self.calculators[context * self.libraries + library];
                calculator.sites += 1;
                let sums = if base == b'C' {
                    &mut calculator.ref_c
                } else {
                    &mut calculator.ref_g
                };
                sums.control_a += counts.control_a;
                sums.oxidated_a += counts.oxidated_a;
                sums.control_c += counts.control_c;
                sums.oxidated_c += counts.oxidated_c;
            }
        }
        self.sites = self.sites.wrapping_add(1);
        if self.sites >= self.stop_after {
            self.stopped = true;
        }
        Ok(())
    }

    fn done(&self) -> bool {
        self.stopped
    }
}

/// The BAM header's `@SQ` lines as a dictionary.
fn dictionary(header: &sam::Header) -> Dictionary {
    header
        .reference_sequences()
        .iter()
        .map(|(name, contig)| Contig {
            name: name.to_string(),
            length: Some(contig.length().get() as i64),
        })
        .collect()
}

/// `doWork`: the reference, the read groups, the calculators, dbSNP, the
/// intervals, the locus walk, then one row per calculator.
fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let contexts = opts.context_order()?;
    let Finished { sink, .. } = inv.read_whole(&inv.input, |header| {
        let reference = Reference::open(&opts.reference, header)?;
        check_sort_order(header)?;
        if header.read_groups().is_empty() {
            return Err(Error::Usage(
                "This analysis requires a read group entry in the alignment file header".to_owned(),
            ));
        }
        let groups = read_groups(header);
        let samples = hash_set_order(groups.iter().map(|g| g.sample.as_deref()), UNKNOWN_SAMPLE);
        let libraries =
            hash_set_order(groups.iter().map(|g| g.library.as_deref()), UNKNOWN_LIBRARY);
        let library_of_group = groups
            .iter()
            .map(|g| {
                let library = g.library.as_deref().unwrap_or(UNKNOWN_LIBRARY);
                libraries.iter().position(|l| l == library).unwrap_or(0)
            })
            .collect();
        let db_snp = match &opts.db_snp {
            Some(path) => Some(DbSnpBitSet::load_all(
                path,
                Some(&dictionary(header)),
                None,
            )?),
            None => None,
        };
        let unique = match &opts.intervals {
            Some(path) => Some(uniqued(&intervals::read_files(&[path], header)?)),
            None => None,
        };
        let sink = OxoGSink {
            minimum_quality_score: opts.minimum_quality_score,
            use_oq: opts.use_oq,
            context_size: usize::try_from(opts.context_size).unwrap_or(0),
            stop_after: opts.stop_after,
            contig_names: header
                .reference_sequences()
                .keys()
                .map(|n| n.to_string())
                .collect(),
            db_snp,
            contexts: contexts
                .iter()
                .enumerate()
                .map(|(i, c)| (c.as_bytes().to_vec(), i))
                .collect(),
            library_of_group,
            libraries: libraries.len(),
            calculators: vec![Calculator::default(); contexts.len() * libraries.len()],
            scratch: vec![Counts::default(); libraries.len()],
            window: Vec::new(),
            sites: 0,
            stopped: false,
            library_names: libraries,
            sample_alias: samples.join(","),
        };
        Ok(Engine::new(
            Classic,
            sink,
            reference,
            header,
            unique.as_deref(),
            Settings {
                minimum_mapping_quality: opts.minimum_mapping_quality,
                setup: Setup::OxoG {
                    insert_size: opts.insert_size(),
                    include_non_pf: opts.include_non_pf_reads,
                },
                locus_accumulation_cap: None,
                stop_after: 0,
            },
        ))
    })?;
    Ok(MetricsOutput::from_rows(&sink.rows(&contexts), None).into())
}

pub const TOOL: Tool = Tool {
    name: "CollectOxoGMetrics",
    summary: "one CpcgMetrics row per sequence context and library: reference and alternate \
              base counts at C and G sites split by read and strand, the oxidation error \
              rate and its Phred score.",
    options: &[
        Opt::path(
            &["REFERENCE_SEQUENCE", "R"],
            "Reference FASTA with a .fai index.",
        )
        .required(),
        Opt::path(
            &["INTERVALS"],
            "An interval list (or a BED): only these loci count.",
        ),
        Opt::path(
            &["DB_SNP"],
            "A VCF of known variants: every position a record spans is skipped.",
        ),
        Opt::int(
            &["MINIMUM_QUALITY_SCORE", "Q"],
            "Bases below this quality (OQ when USE_OQ) are not counted (default 20).",
        ),
        Opt::int(
            &["MINIMUM_MAPPING_QUALITY", "MQ"],
            "Reads below this MAPQ are skipped (default 30).",
        ),
        Opt::int(
            &["MINIMUM_INSERT_SIZE", "MIN_INS"],
            "Reads with |TLEN| below this are skipped; with both bounds 0 unpaired reads \
             count too (default 60).",
        ),
        Opt::int(
            &["MAXIMUM_INSERT_SIZE", "MAX_INS"],
            "Reads with |TLEN| above this are skipped (default 600).",
        ),
        Opt::bool(
            &["INCLUDE_NON_PF_READS", "NON_PF"],
            "Count reads that fail vendor quality checks (default true).",
        ),
        Opt::bool(
            &["USE_OQ"],
            "Test base quality on the OQ tag when a read has one (default true).",
        ),
        Opt::int(
            &["CONTEXT_SIZE"],
            "Reference bases on each side of the C or G (default 1).",
        ),
        Opt::text(
            &["CONTEXTS"],
            "A context to report, C in the middle; each adds one, null clears, none reports \
             every context.",
        )
        .repeated(),
        Opt::int(
            &["STOP_AFTER"],
            "Stop after this many C or G sites in a reported context (default 2147483647).",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

#[cfg(test)]
mod tests {
    use super::*;

    fn options() -> Options {
        Options {
            reference: PathBuf::from("ref.fa"),
            intervals: None,
            db_snp: None,
            minimum_quality_score: 20,
            minimum_mapping_quality: 30,
            minimum_insert_size: 60,
            maximum_insert_size: 600,
            include_non_pf_reads: true,
            use_oq: true,
            context_size: 1,
            contexts: Vec::new(),
            stop_after: i32::MAX,
        }
    }

    #[test]
    fn generated_contexts_have_c_in_the_middle() {
        let mut opts = options();
        opts.context_size = 0;
        assert_eq!(opts.context_order().unwrap(), ["C"]);
        opts.context_size = 1;
        let order = opts.context_order().unwrap();
        assert_eq!(order.len(), 16);
        assert!(order.iter().all(|c| c.len() == 3 && &c[1..2] == "C"));
        opts.context_size = 2;
        assert_eq!(opts.context_order().unwrap().len(), 256);
    }

    #[test]
    fn repeated_contexts_count_once() {
        let mut opts = options();
        opts.contexts = vec!["ACG".into(), "TCA".into(), "ACG".into()];
        let order = opts.context_order().unwrap();
        assert_eq!(order.len(), 2);
        assert!(order.contains(&"ACG".to_owned()) && order.contains(&"TCA".to_owned()));
    }

    #[test]
    fn validation_follows_picard() {
        let mut opts = options();
        opts.contexts = vec!["AC".into(), "AAG".into()];
        let message = opts.validate().unwrap_err().to_string();
        assert!(message.contains("Context AC is not 3 long as implied by CONTEXT_SIZE=1"));
        assert!(message.contains("Middle base of context sequence AAG must be C"));
        let mut opts = options();
        opts.minimum_insert_size = 700;
        assert!(opts.validate().is_err());
        let mut opts = options();
        opts.context_size = -1;
        assert!(opts.validate().is_err());
        let mut opts = options();
        opts.minimum_insert_size = 0;
        opts.maximum_insert_size = 0;
        assert_eq!(opts.insert_size(), None);
    }

    #[test]
    fn finish_matches_picard_arithmetic() {
        let calculator = Calculator {
            sites: 3,
            ref_c: Counts {
                control_a: 1,
                oxidated_a: 4,
                control_c: 50,
                oxidated_c: 45,
            },
            ref_g: Counts {
                control_a: 0,
                oxidated_a: 2,
                control_c: 40,
                oxidated_c: 58,
            },
        };
        let row = calculator.finish("S", "L", "ACA");
        assert_eq!(row.total_bases, 200);
        assert_eq!(row.ref_oxo_bases, 103);
        assert_eq!(row.ref_nonoxo_bases, 90);
        assert_eq!(row.alt_oxo_bases, 6);
        assert_eq!(row.alt_nonoxo_bases, 1);
        assert_eq!(row.oxidation_error_rate, 5.0 / 200.0);
        assert_eq!(row.c_ref_alt_bases, 5);
        assert!(row.g_ref_oxo_error_rate == 1e-10);

        let empty = Calculator::default().finish("S", "L", "ACA");
        assert!(empty.oxidation_error_rate.is_infinite());
        assert!(empty.c_ref_oxo_error_rate.is_nan() && empty.g_ref_oxo_q.is_nan());
    }
}
