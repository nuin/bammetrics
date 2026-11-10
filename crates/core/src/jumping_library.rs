//! `CollectJumpingLibraryMetrics`, mirroring Picard 3.5.0
//! (`CollectJumpingLibraryMetrics`, `JumpingLibraryMetrics`,
//! `Histogram.trimByTailLimit`, `DuplicationMetrics.estimateLibrarySize`).
//!
//! One row over every `INPUT`, read in turn. Rules copied from the source:
//! * Every record must be paired (`getFirstOfPairFlag` throws otherwise);
//!   only first-of-pair records count, secondary, supplementary, QC-failed
//!   and duplicate records among them.
//! * A first read that is unmapped with a mapped mate, or mapped with an
//!   unmapped mate, is a `FRAGMENT`. A first read unmapped with its mate
//!   unmapped and no reference index ends the read of that file (the
//!   unaligned records at the end of a coordinate-sorted BAM); one with a
//!   reference index is skipped.
//! * A pair is skipped when the read's mapping quality, or its `MQ` tag
//!   when present, is below `MINIMUM_MAPPING_QUALITY`.
//! * `|TLEN|` above the greater of `CHIMERA_KB_MIN` and the outward-facing
//!   mode is a supersized chimera; otherwise a same-strand pair is a tandem
//!   chimera; otherwise a pair on two contigs is a cross-contig chimera;
//!   otherwise `SamPairUtil.getPairOrientation` ([`crate::record`]) files
//!   the pair as RF (a jump) or FR (a non-jump), with its `|TLEN|` in that
//!   orientation's histogram and its duplicate flag counted.
//! * The outward-facing mode comes from a first pass that samples the first
//!   `50000 / INPUT.size()` RF pairs of each file (same tests as above but
//!   the chimera one, the pass over a file ending at its unaligned tail).
//!   Without a sample the mode is 0.
//! * Mean and SD of each orientation use the histogram after
//!   [`Histogram::trim_by_tail_limit`]: every bin up to the mode, then the
//!   consecutive bins holding at least mode count / `TAIL_LIMIT` pairs. An
//!   empty histogram gives a mean of `?` and an SD of `-0`.
//! * The library size is `estimateLibrarySize(pairs, pairs - duplicates)`
//!   when both pairs and duplicates are above 0, else 0; the percentages
//!   divide by the pairs of that orientation, or by jumps + non-jumps +
//!   chimeras, and are 0 without any.
//! * Every `INPUT` must be coordinate sorted by its header.

use std::path::Path;

use noodles::sam::alignment::record::data::field::Tag;
use noodles::{bam, sam};

use crate::duplication::estimate_library_size;
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::pass::{self, Accumulator};
use crate::record::{Orientation, orientation, record_name};
use crate::registry::{Invocation, Opt, Tool};

/// `SAMPLE_FOR_MODE`: outward-facing pairs sampled over all inputs for the
/// mode that sets the chimera size.
const SAMPLE_FOR_MODE: usize = 50_000;

#[derive(Clone, Debug)]
pub struct Options {
    pub minimum_mapping_quality: i64,
    pub tail_limit: i64,
    pub chimera_kb_min: i64,
}

impl Options {
    /// The options of a parsed `CollectJumpingLibraryMetrics` command line.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let defaults = Self::default();
        Ok(Self {
            minimum_mapping_quality: inv
                .parse("MINIMUM_MAPPING_QUALITY", defaults.minimum_mapping_quality)?,
            tail_limit: inv.parse("TAIL_LIMIT", defaults.tail_limit)?,
            chimera_kb_min: inv.parse("CHIMERA_KB_MIN", defaults.chimera_kb_min)?,
        })
    }
}

impl Default for Options {
    fn default() -> Self {
        Self {
            minimum_mapping_quality: 0,
            tail_limit: 10_000,
            chimera_kb_min: 100_000,
        }
    }
}

/// `picard.analysis.JumpingLibraryMetrics`.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub jump_pairs: i64,
    pub jump_duplicate_pairs: i64,
    pub jump_duplicate_pct: f64,
    pub jump_library_size: i64,
    pub jump_mean_insert_size: f64,
    pub jump_stdev_insert_size: f64,
    pub nonjump_pairs: i64,
    pub nonjump_duplicate_pairs: i64,
    pub nonjump_duplicate_pct: f64,
    pub nonjump_library_size: i64,
    pub nonjump_mean_insert_size: f64,
    pub nonjump_stdev_insert_size: f64,
    pub chimeric_pairs: i64,
    pub fragments: i64,
    pub pct_jumps: f64,
    pub pct_nonjumps: f64,
    pub pct_chimeras: f64,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.analysis.JumpingLibraryMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "JUMP_PAIRS",
        "JUMP_DUPLICATE_PAIRS",
        "JUMP_DUPLICATE_PCT",
        "JUMP_LIBRARY_SIZE",
        "JUMP_MEAN_INSERT_SIZE",
        "JUMP_STDEV_INSERT_SIZE",
        "NONJUMP_PAIRS",
        "NONJUMP_DUPLICATE_PAIRS",
        "NONJUMP_DUPLICATE_PCT",
        "NONJUMP_LIBRARY_SIZE",
        "NONJUMP_MEAN_INSERT_SIZE",
        "NONJUMP_STDEV_INSERT_SIZE",
        "CHIMERIC_PAIRS",
        "FRAGMENTS",
        "PCT_JUMPS",
        "PCT_NONJUMPS",
        "PCT_CHIMERAS",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.jump_pairs.into(),
            self.jump_duplicate_pairs.into(),
            self.jump_duplicate_pct.into(),
            self.jump_library_size.into(),
            self.jump_mean_insert_size.into(),
            self.jump_stdev_insert_size.into(),
            self.nonjump_pairs.into(),
            self.nonjump_duplicate_pairs.into(),
            self.nonjump_duplicate_pct.into(),
            self.nonjump_library_size.into(),
            self.nonjump_mean_insert_size.into(),
            self.nonjump_stdev_insert_size.into(),
            self.chimeric_pairs.into(),
            self.fragments.into(),
            self.pct_jumps.into(),
            self.pct_nonjumps.into(),
            self.pct_chimeras.into(),
        ]
    }
}

/// `CollectJumpingLibraryMetrics` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CollectJumpingLibraryMetrics",
    summary: "one row over every input: outward-facing (jump) and inward-facing (non-jump) \
              pairs with their duplicate fraction, library size and insert size mean and SD, \
              chimeric pairs, fragments and the three fractions; the inputs may be several \
              coordinate-sorted BAMs, read in turn.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "A coordinate-sorted BAM or CRAM with duplicates marked; each adds to the list.",
        )
        .required()
        .repeated(),
        Opt::int(
            &["MINIMUM_MAPPING_QUALITY", "MQ"],
            "Skip a pair whose read, or whose MQ tag when present, is below this (default 0).",
        ),
        Opt::int(
            &["TAIL_LIMIT", "T"],
            "Mean and SD stop at the first bin past the mode holding fewer than mode count / \
             TAIL_LIMIT pairs (default 10000).",
        ),
        Opt::int(
            &["CHIMERA_KB_MIN"],
            "An insert above the greater of this and the outward-facing mode is a chimera \
             (default 100000).",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let inputs = inv.all("INPUT");
    if let Some(stream) = inputs.iter().find(|i| is_standard_input(i)) {
        return Err(Error::Usage(format!(
            "INPUT={stream}: CollectJumpingLibraryMetrics reads each input twice (the mode \
             pass, then the counts), so it needs a file, not standard input"
        )));
    }
    let mode = outie_mode(inv, &opts, inputs)?;
    let mut counts = Counts::new(opts, mode);
    for input in inputs {
        let path = Path::new(input);
        counts = inv.read_whole(path, |header| {
            check_coordinate_sorted(path, header)?;
            counts.start_file();
            Ok(counts)
        })?;
    }
    Ok(MetricsOutput::from_rows(&[counts.finish()], None).into())
}

/// `-` (the verb grammar maps it to `/dev/stdin`) or the device itself: a
/// stream the second pass could not read again.
fn is_standard_input(input: &str) -> bool {
    input == "-" || input == "/dev/stdin" || input == "/dev/fd/0"
}

/// Picard's `PicardException` on a header whose sort order is not
/// `coordinate`; the tool has no `ASSUME_SORTED`.
fn check_coordinate_sorted(path: &Path, header: &sam::Header) -> Result<()> {
    let order = pass::sort_order(header);
    if order == "coordinate" {
        return Ok(());
    }
    Err(Error::Metric(format!(
        "{}: the header says the sort order is {order}; CollectJumpingLibraryMetrics needs a \
         coordinate sorted file",
        path.display()
    )))
}

/// `getOutieMode`: the mode of `|TLEN|` over the first `SAMPLE_FOR_MODE /
/// inputs` outward-facing pairs of each input, 0 without any.
fn outie_mode(inv: &Invocation, opts: &Options, inputs: &[String]) -> Result<f64> {
    let per_file = SAMPLE_FOR_MODE / inputs.len().max(1);
    let mut sampler = ModeSampler {
        minimum_mapping_quality: opts.minimum_mapping_quality,
        per_file,
        sampled: 0,
        stop: false,
        histogram: Histogram::new(),
    };
    for input in inputs {
        sampler = inv.read_whole(Path::new(input), |_| {
            sampler.sampled = 0;
            sampler.stop = false;
            Ok(sampler)
        })?;
    }
    Ok(sampler.histogram.mode().map_or(0.0, |m| m as f64))
}

/// What every record of the tool must be: paired, as
/// `SAMRecord.getFirstOfPairFlag` requires.
fn require_paired(record: &bam::Record) -> Result<()> {
    if record.flags().is_segmented() {
        return Ok(());
    }
    Err(Error::MalformedRecord {
        name: record_name(record),
        what: "not paired; CollectJumpingLibraryMetrics reads paired records only",
    })
}

/// The `MQ` tag as `getIntegerAttribute` returns it: `None` without the
/// tag, an error for a value that is not an integer.
fn mate_mapping_quality(record: &bam::Record) -> Result<Option<i64>> {
    match record.data().get(&Tag::MATE_MAPPING_QUALITY) {
        None => Ok(None),
        Some(Err(e)) => Err(e.into()),
        Some(Ok(value)) => value
            .as_int()
            .map(Some)
            .ok_or_else(|| Error::MalformedRecord {
                name: record_name(record),
                what: "MQ tag is not an integer",
            }),
    }
}

/// `getMappingQuality`: 255 for a record without one.
fn mapping_quality(record: &bam::Record) -> i64 {
    record.mapping_quality().map_or(255, |q| i64::from(q.get()))
}

/// The mapping quality test of both passes: the read's own quality and,
/// when the tag is present, the mate's.
fn passes_mapping_quality(record: &bam::Record, minimum: i64) -> Result<bool> {
    let mate_ok = mate_mapping_quality(record)?.is_none_or(|mq| mq >= minimum);
    Ok(mate_ok && mapping_quality(record) >= minimum)
}

fn same_contig(record: &bam::Record) -> Result<bool> {
    let own = record.reference_sequence_id().transpose()?;
    let mate = record.mate_reference_sequence_id().transpose()?;
    Ok(own == mate)
}

fn same_strand(record: &bam::Record) -> bool {
    let flags = record.flags();
    flags.is_reverse_complemented() == flags.is_mate_reverse_complemented()
}

/// True for the first-of-pair record where `doWork` and `getOutieMode`
/// `break`: unmapped, mate unmapped and no reference index.
fn ends_aligned_records(record: &bam::Record) -> bool {
    let flags = record.flags();
    flags.is_unmapped() && flags.is_mate_unmapped() && record.reference_sequence_id().is_none()
}

/// The first pass over each input (`getOutieMode`).
struct ModeSampler {
    minimum_mapping_quality: i64,
    per_file: usize,
    sampled: usize,
    stop: bool,
    histogram: Histogram,
}

impl Accumulator for ModeSampler {
    type Output = ModeSampler;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        require_paired(record)?;
        let flags = record.flags();
        if !flags.is_first_segment() {
            return Ok(());
        }
        if ends_aligned_records(record) {
            self.stop = true;
            return Ok(());
        }
        if flags.is_unmapped() || flags.is_mate_unmapped() {
            return Ok(());
        }
        if passes_mapping_quality(record, self.minimum_mapping_quality)?
            && !same_strand(record)
            && same_contig(record)?
            && orientation(record)? == Orientation::Rf
        {
            self.histogram
                .increment(i64::from(record.template_length()).abs());
            self.sampled += 1;
        }
        Ok(())
    }

    fn done(&self) -> bool {
        self.stop || self.sampled >= self.per_file
    }

    fn finish(self) -> Result<ModeSampler> {
        Ok(self)
    }
}

/// The counts of `doWork`, carried from one input to the next.
struct Counts {
    opts: Options,
    /// `chimeraSizeMinimum`: the greater of the sampled mode and `CHIMERA_KB_MIN`.
    chimera_size_minimum: f64,
    stop: bool,
    fragments: i64,
    innies: i64,
    outies: i64,
    innie_dupes: i64,
    outie_dupes: i64,
    cross_contig: i64,
    super_sized: i64,
    tandem: i64,
    innie_histogram: Histogram,
    outie_histogram: Histogram,
}

impl Counts {
    fn new(opts: Options, mode: f64) -> Self {
        let chimera_size_minimum = mode.max(opts.chimera_kb_min as f64);
        Self {
            opts,
            chimera_size_minimum,
            stop: false,
            fragments: 0,
            innies: 0,
            outies: 0,
            innie_dupes: 0,
            outie_dupes: 0,
            cross_contig: 0,
            super_sized: 0,
            tandem: 0,
            innie_histogram: Histogram::new(),
            outie_histogram: Histogram::new(),
        }
    }

    fn start_file(&mut self) {
        self.stop = false;
    }

    /// The row: counts, library sizes, the trimmed means and SDs, the
    /// fractions.
    fn finish(self) -> Row {
        let Counts {
            opts,
            fragments,
            innies,
            outies,
            innie_dupes,
            outie_dupes,
            cross_contig,
            super_sized,
            tandem,
            mut innie_histogram,
            mut outie_histogram,
            ..
        } = self;
        let pct = |part: i64, whole: i64| {
            if whole != 0 {
                part as f64 / whole as f64
            } else {
                0.0
            }
        };
        let library_size = |pairs: i64, duplicates: i64| {
            if pairs > 0 && duplicates > 0 {
                estimate_library_size(pairs, pairs - duplicates).unwrap_or(0)
            } else {
                0
            }
        };
        outie_histogram.trim_by_tail_limit(opts.tail_limit);
        innie_histogram.trim_by_tail_limit(opts.tail_limit);
        let chimeric_pairs = cross_contig + super_sized + tandem;
        let total = outies + innies + chimeric_pairs;
        Row {
            jump_pairs: outies,
            jump_duplicate_pairs: outie_dupes,
            jump_duplicate_pct: pct(outie_dupes, outies),
            jump_library_size: library_size(outies, outie_dupes),
            jump_mean_insert_size: outie_histogram.mean(),
            jump_stdev_insert_size: outie_histogram.standard_deviation(),
            nonjump_pairs: innies,
            nonjump_duplicate_pairs: innie_dupes,
            nonjump_duplicate_pct: pct(innie_dupes, innies),
            nonjump_library_size: library_size(innies, innie_dupes),
            nonjump_mean_insert_size: innie_histogram.mean(),
            nonjump_stdev_insert_size: innie_histogram.standard_deviation(),
            chimeric_pairs,
            fragments,
            pct_jumps: pct(outies, total),
            pct_nonjumps: pct(innies, total),
            pct_chimeras: pct(chimeric_pairs, total),
        }
    }
}

impl Accumulator for Counts {
    type Output = Counts;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        require_paired(record)?;
        let flags = record.flags();
        if !flags.is_first_segment() {
            return Ok(());
        }
        if flags.is_unmapped() {
            if !flags.is_mate_unmapped() {
                self.fragments += 1;
            } else if record.reference_sequence_id().is_none() {
                self.stop = true;
            }
            return Ok(());
        }
        if flags.is_mate_unmapped() {
            self.fragments += 1;
            return Ok(());
        }
        if !passes_mapping_quality(record, self.opts.minimum_mapping_quality)? {
            return Ok(());
        }
        let insert_size = i64::from(record.template_length()).abs();
        if insert_size as f64 > self.chimera_size_minimum {
            self.super_sized += 1;
        } else if same_strand(record) {
            self.tandem += 1;
        } else if !same_contig(record)? {
            self.cross_contig += 1;
        } else {
            let (histogram, pairs, dupes) = match orientation(record)? {
                Orientation::Rf => (
                    &mut self.outie_histogram,
                    &mut self.outies,
                    &mut self.outie_dupes,
                ),
                Orientation::Fr => (
                    &mut self.innie_histogram,
                    &mut self.innies,
                    &mut self.innie_dupes,
                ),
                // The same-strand test above already took every tandem pair.
                Orientation::Tandem => {
                    return Err(Error::MalformedRecord {
                        name: record_name(record),
                        what: "same-strand pair reached the orientation test",
                    });
                }
            };
            histogram.increment(insert_size);
            *pairs += 1;
            if flags.is_duplicate() {
                *dupes += 1;
            }
        }
        Ok(())
    }

    fn done(&self) -> bool {
        self.stop
    }

    fn finish(self) -> Result<Counts> {
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics_file::format_double;

    fn counts() -> Counts {
        Counts::new(Options::default(), 0.0)
    }

    #[test]
    fn empty_input_prints_picards_row() {
        let row = counts().finish();
        assert_eq!(row.jump_pairs, 0);
        assert_eq!(format_double(row.jump_mean_insert_size), "?");
        assert_eq!(format_double(row.jump_stdev_insert_size), "-0");
        assert_eq!(row.pct_jumps, 0.0);
        assert_eq!(row.jump_library_size, 0);
    }

    #[test]
    fn library_size_needs_pairs_and_duplicates() {
        let mut c = counts();
        c.outies = 10;
        c.outie_dupes = 0;
        c.innies = 10;
        c.innie_dupes = 2;
        for key in [300, 300, 300, 310, 320, 1000] {
            c.innie_histogram.increment(key);
        }
        let row = c.finish();
        assert_eq!(row.jump_library_size, 0);
        assert_eq!(
            row.nonjump_library_size,
            estimate_library_size(10, 8).unwrap()
        );
        assert_eq!(row.nonjump_duplicate_pct, 0.2);
        assert_eq!(row.pct_jumps, 0.5);
        // TAIL_LIMIT 10000: mode 300 (count 3), 301 is missing, so the mean
        // is over the bins up to the mode only.
        assert_eq!(row.nonjump_mean_insert_size, 300.0);
    }

    #[test]
    fn standard_input_is_refused() {
        assert!(is_standard_input("-"));
        assert!(is_standard_input("/dev/stdin"));
        assert!(!is_standard_input("in.bam"));
        let args: Vec<String> = ["I=-", "O=out.txt"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let err = TOOL.execute(&args).expect_err("refused");
        assert!(err.to_string().contains("standard input"), "{err}");
    }

    #[test]
    fn chimera_threshold_is_the_greater_of_mode_and_option() {
        let opts = Options {
            chimera_kb_min: 1000,
            ..Options::default()
        };
        assert_eq!(
            Counts::new(opts.clone(), 2500.0).chimera_size_minimum,
            2500.0
        );
        assert_eq!(Counts::new(opts, 20.0).chimera_size_minimum, 1000.0);
    }
}
