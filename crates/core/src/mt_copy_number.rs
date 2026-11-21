//! `CollectMtCopyNumberMetrics` (docs/06 section 6): mitochondrial DNA
//! copy number per cell, `2 * MT_MEAN_DEPTH / AUTOSOMAL_MEAN_DEPTH`. Not a
//! Picard tool.
//!
//! One pass. A record counts when it is primary (not secondary, not
//! supplementary), mapped, not a duplicate, not QC-failed and at or above
//! `MINIMUM_MAPPING_QUALITY` (a missing MAPQ, 255, passes), the filters of
//! `sex_check.rs`. Its counted bases are the bases of its alignment blocks
//! (`blocks.rs`: `M`, `=` and `X`; soft clips, insertions and deletions add
//! nothing) at or above `MINIMUM_BASE_QUALITY` that fall in the territory:
//! the mitochondrial contig and the autosomes, or with `TARGET_INTERVALS`
//! the mitochondrial contig and the merged targets on the autosomes, less
//! the merged `EXCLUDE_INTERVALS` in both cases. A record without `QUAL`
//! counts no base unless `MINIMUM_BASE_QUALITY` is 0. One row per
//! accumulation unit (`levels.rs`), empty units included.

use std::path::PathBuf;

use noodles::{bam, sam};

use crate::blocks::alignment_blocks;
use crate::error::{Error, Result};
use crate::genome_build;
use crate::intervals::{self, Interval, OverlapIndex};
use crate::levels::{Labels, Level, Levels, parse_levels};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::pass::{self, Accumulator};
use crate::registry::{Invocation, Opt, Tool};

/// The metrics class string. Not a Picard class: readers key on the last
/// name.
pub const CLASS: &str = "bammetrics.qc.MtCopyNumberMetrics";

pub const DEFAULT_MINIMUM_MAPPING_QUALITY: u8 = 20;
pub const DEFAULT_MINIMUM_BASE_QUALITY: u8 = 20;

/// A BAM stores a missing `QUAL` (`*`) as 0xFF bytes; a record read from
/// SAM text may hold none.
const MISSING_QUALITY: u8 = 0xFF;

/// One row: a unit's counts against the shared territory.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub mt_contig: String,
    pub mt_length: i64,
    pub mt_bases: i64,
    pub mt_mean_depth: f64,
    pub autosomal_length: i64,
    pub autosomal_bases: i64,
    pub autosomal_mean_depth: f64,
    pub mt_copy_number: f64,
    pub reads_used: i64,
    pub labels: Labels,
}

impl MetricRow for Row {
    const CLASS: &'static str = CLASS;
    const COLUMNS: &'static [&'static str] = &[
        "MT_CONTIG",
        "MT_LENGTH",
        "MT_BASES",
        "MT_MEAN_DEPTH",
        "AUTOSOMAL_LENGTH",
        "AUTOSOMAL_BASES",
        "AUTOSOMAL_MEAN_DEPTH",
        "MT_COPY_NUMBER",
        "READS_USED",
        "SAMPLE",
        "LIBRARY",
        "READ_GROUP",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.mt_contig.clone().into(),
            self.mt_length.into(),
            self.mt_bases.into(),
            self.mt_mean_depth.into(),
            self.autosomal_length.into(),
            self.autosomal_bases.into(),
            self.autosomal_mean_depth.into(),
            self.mt_copy_number.into(),
            self.reads_used.into(),
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]
    }
}

/// The counts of one accumulation unit.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub mt_bases: i64,
    pub autosomal_bases: i64,
    pub reads_used: i64,
}

/// `bases / length`, NaN for an empty territory.
fn mean_depth(bases: i64, length: i64) -> f64 {
    if length == 0 {
        f64::NAN
    } else {
        bases as f64 / length as f64
    }
}

/// The row of one unit. `MT_COPY_NUMBER` is NaN (`?`) when the autosomal
/// depth is 0 or undefined.
pub fn row(mt_contig: &str, mt_length: i64, autosomal_length: i64, counts: Counts) -> Row {
    let mt_mean_depth = mean_depth(counts.mt_bases, mt_length);
    let autosomal_mean_depth = mean_depth(counts.autosomal_bases, autosomal_length);
    let mt_copy_number = if autosomal_mean_depth > 0.0 {
        2.0 * mt_mean_depth / autosomal_mean_depth
    } else {
        f64::NAN
    };
    Row {
        mt_contig: mt_contig.to_owned(),
        mt_length,
        mt_bases: counts.mt_bases,
        mt_mean_depth,
        autosomal_length,
        autosomal_bases: counts.autosomal_bases,
        autosomal_mean_depth,
        mt_copy_number,
        reads_used: counts.reads_used,
        labels: Labels::default(),
    }
}

/// The options of one run.
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub mt_contig: Option<String>,
    pub autosomes: Vec<String>,
    pub minimum_mapping_quality: u8,
    pub minimum_base_quality: u8,
    pub target_intervals: Vec<PathBuf>,
    pub exclude_intervals: Vec<PathBuf>,
    pub levels: Vec<Level>,
}

/// An option of 0 to 255.
fn quality(inv: &Invocation, name: &str, default: u8) -> Result<u8> {
    let value: i64 = inv.parse(name, i64::from(default))?;
    u8::try_from(value).map_err(|_| Error::Usage(format!("{name}={value}: 0 to 255")))
}

fn paths(inv: &Invocation, name: &str) -> Vec<PathBuf> {
    inv.all(name).iter().map(PathBuf::from).collect()
}

impl Options {
    pub(crate) fn from_invocation(inv: &Invocation) -> Result<Self> {
        Ok(Self {
            mt_contig: inv.get("MT_CONTIG").map(str::to_owned),
            autosomes: inv.all("AUTOSOMES").to_vec(),
            minimum_mapping_quality: quality(
                inv,
                "MINIMUM_MAPPING_QUALITY",
                DEFAULT_MINIMUM_MAPPING_QUALITY,
            )?,
            minimum_base_quality: quality(
                inv,
                "MINIMUM_BASE_QUALITY",
                DEFAULT_MINIMUM_BASE_QUALITY,
            )?,
            target_intervals: paths(inv, "TARGET_INTERVALS"),
            exclude_intervals: paths(inv, "EXCLUDE_INTERVALS"),
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
        })
    }
}

/// What a contig is to the tool.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Role {
    Mito,
    Autosome,
}

/// The header index of the mitochondrial contig and the role of each
/// contig. A missing contig is an error that names the option to set.
fn roles(opts: &Options, header: &sam::Header) -> Result<(usize, Vec<Option<Role>>)> {
    let names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();
    // The header's index map: one hash lookup per name, so a long AUTOSOMES
    // list costs one probe per value.
    let index_of = |name: &str| header.reference_sequences().get_index_of(name.as_bytes());
    let mt = match &opts.mt_contig {
        Some(name) => index_of(name).ok_or_else(|| {
            Error::Usage(format!("MT_CONTIG={name}: the header has no such contig"))
        })?,
        None => {
            let found: Vec<usize> = names
                .iter()
                .enumerate()
                .filter(|(_, n)| genome_build::is_mito(n))
                .map(|(i, _)| i)
                .collect();
            match found.as_slice() {
                [one] => *one,
                [] => {
                    return Err(Error::Usage(format!(
                        "the header has no mitochondrial contig ({}); set MT_CONTIG",
                        genome_build::MITO_NAMES.join(", ")
                    )));
                }
                several => {
                    let list: Vec<&str> = several.iter().map(|&i| names[i].as_str()).collect();
                    return Err(Error::Usage(format!(
                        "the header has several mitochondrial contigs ({}); set MT_CONTIG",
                        list.join(", ")
                    )));
                }
            }
        }
    };
    let mut roles: Vec<Option<Role>> = vec![None; names.len()];
    roles[mt] = Some(Role::Mito);
    if opts.autosomes.is_empty() {
        for (index, name) in names.iter().enumerate() {
            if index != mt && genome_build::is_autosome(name) {
                roles[index] = Some(Role::Autosome);
            }
        }
    } else {
        for name in &opts.autosomes {
            let index = index_of(name).ok_or_else(|| {
                Error::Usage(format!("AUTOSOMES={name}: the header has no such contig"))
            })?;
            if index == mt {
                return Err(Error::Usage(format!(
                    "AUTOSOMES={name}: that is the mitochondrial contig"
                )));
            }
            roles[index] = Some(Role::Autosome);
        }
    }
    if !roles.contains(&Some(Role::Autosome)) {
        return Err(Error::Usage(
            "the header has no autosome (chr1 to chr22 or 1 to 22); set AUTOSOMES".to_owned(),
        ));
    }
    Ok((mt, roles))
}

/// The bases of `include` outside `exclude`, both `uniqued`: one merge walk
/// (`interval_ops::subtract_disjoint`).
pub use crate::interval_ops::subtract_disjoint as subtract;

/// `CollectMtCopyNumberMetrics` as a pass [`Accumulator`].
pub struct MtCopyNumber {
    roles: Vec<Option<Role>>,
    territory: OverlapIndex,
    mt_contig: String,
    mt_length: i64,
    autosomal_length: i64,
    minimum_mapping_quality: u8,
    minimum_base_quality: u8,
    levels: Levels<Counts>,
}

impl MtCopyNumber {
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        let (mt, roles) = roles(opts, header)?;
        let lengths: Vec<i64> = header
            .reference_sequences()
            .values()
            .map(|contig| contig.length().get() as i64)
            .collect();
        let whole = |contig: usize| Interval {
            contig,
            start: 1,
            end: lengths[contig],
        };
        let mut include = vec![whole(mt)];
        if opts.target_intervals.is_empty() {
            include.extend(
                (0..roles.len())
                    .filter(|&i| roles[i] == Some(Role::Autosome))
                    .map(whole),
            );
        } else {
            include.extend(
                intervals::read_files(&opts.target_intervals, header)?
                    .into_iter()
                    .filter(|iv| roles[iv.contig] == Some(Role::Autosome)),
            );
        }
        let include = intervals::uniqued(&include);
        let exclude = intervals::uniqued(&intervals::read_files(&opts.exclude_intervals, header)?);
        let territory = subtract(&include, &exclude);
        let (mut mt_length, mut autosomal_length) = (0, 0);
        for interval in &territory {
            if interval.contig == mt {
                mt_length += interval.length();
            } else {
                autosomal_length += interval.length();
            }
        }
        let mt_contig = header
            .reference_sequences()
            .keys()
            .nth(mt)
            .map(ToString::to_string)
            .unwrap_or_default();
        Ok(Self {
            territory: OverlapIndex::new(territory, roles.len(), 0),
            roles,
            mt_contig,
            mt_length,
            autosomal_length,
            minimum_mapping_quality: opts.minimum_mapping_quality,
            minimum_base_quality: opts.minimum_base_quality,
            levels: Levels::new(header, &opts.levels, Counts::default()),
        })
    }

    /// The contig and role of a record the filters keep.
    fn counted(&self, record: &bam::Record) -> Result<Option<(usize, Role)>> {
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            return Ok(None);
        }
        if record
            .mapping_quality()
            .is_some_and(|q| u8::from(q) < self.minimum_mapping_quality)
        {
            return Ok(None);
        }
        let Some(contig) = record.reference_sequence_id().transpose()? else {
            return Ok(None);
        };
        Ok(self
            .roles
            .get(contig)
            .copied()
            .flatten()
            .map(|role| (contig, role)))
    }

    /// The record's aligned bases in the territory at or above the base
    /// quality limit.
    fn bases(&self, record: &bam::Record, contig: usize) -> Result<i64> {
        let qualities = record.quality_scores();
        let qualities = qualities.as_bytes();
        let missing = qualities.first().is_none_or(|&q| q == MISSING_QUALITY);
        if missing && self.minimum_base_quality > 0 {
            return Ok(0);
        }
        let all_pass = self.minimum_base_quality == 0;
        let intervals = self.territory.intervals();
        let mut bases = 0;
        for block in alignment_blocks(record)? {
            let end = block.reference_end();
            for interval in &intervals[self.territory.overlaps(contig, block.reference_start, end)]
            {
                let (from, to) = (
                    interval.start.max(block.reference_start),
                    interval.end.min(end),
                );
                if all_pass {
                    bases += to - from + 1;
                    continue;
                }
                let first = block.read_start - 1 + (from - block.reference_start) as usize;
                let count = (to - from + 1) as usize;
                let window = qualities
                    .get(first..first + count)
                    .ok_or(Error::MalformedRecord {
                        name: crate::record::record_name(record),
                        what: "fewer base qualities than the CIGAR aligns",
                    })?;
                bases += window
                    .iter()
                    .filter(|&&q| q >= self.minimum_base_quality)
                    .count() as i64;
            }
        }
        Ok(bases)
    }
}

impl Accumulator for MtCopyNumber {
    type Output = Vec<Row>;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let Some((contig, role)) = self.counted(record)? else {
            return Ok(());
        };
        let bases = self.bases(record, contig)?;
        if bases == 0 {
            return Ok(());
        }
        self.levels.route(record, |unit| {
            match role {
                Role::Mito => unit.mt_bases += bases,
                Role::Autosome => unit.autosomal_bases += bases,
            }
            unit.reads_used += 1;
        });
        Ok(())
    }

    fn finish(self) -> Result<Vec<Row>> {
        Ok(self
            .levels
            .into_buckets()
            .into_iter()
            .map(|bucket| Row {
                labels: bucket.labels,
                ..row(
                    &self.mt_contig,
                    self.mt_length,
                    self.autosomal_length,
                    bucket.state,
                )
            })
            .collect())
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectMtCopyNumberMetrics",
    summary: "one row of the 12 MtCopyNumberMetrics columns per accumulation unit: mean depth \
              over the mitochondrial contig and over the autosomes from aligned bases, and \
              MT_COPY_NUMBER, 2 x their ratio. Not a Picard tool (docs/06).",
    options: &[
        Opt::text(
            &["MT_CONTIG"],
            "The mitochondrial contig (default the one of chrM, MT, M or chrMT in the header).",
        ),
        Opt::text(
            &["AUTOSOMES"],
            "An autosome contig name; repeat for each (default chr1 to chr22 or 1 to 22).",
        )
        .repeated(),
        Opt::int(
            &["MINIMUM_MAPPING_QUALITY", "MQ"],
            "Reads below this MAPQ are not counted (default 20).",
        ),
        Opt::int(
            &["MINIMUM_BASE_QUALITY", "Q"],
            "Aligned bases below this quality are not counted (default 20).",
        ),
        Opt::path(
            &["TARGET_INTERVALS", "TI"],
            "Capture targets, an interval list or BED; repeat for several. The autosomal \
             depth is then over the merged targets only; the mitochondrial contig stays whole.",
        )
        .repeated(),
        Opt::path(
            &["EXCLUDE_INTERVALS"],
            "Regions left out of both depths (NUMTs, masked bases), an interval list or BED; \
             repeat for several.",
        )
        .repeated(),
        crate::registry::LEVEL_OPTION,
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn pass::Collector>> {
    let opts = Options::from_invocation(inv)?;
    Ok(pass::boxed(MtCopyNumber::new(&opts, header)?, |rows| {
        Ok(MetricsOutput::from_rows(&rows, None).into())
    }))
}

#[cfg(test)]
mod tests;
