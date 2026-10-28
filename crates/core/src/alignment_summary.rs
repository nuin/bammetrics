//! `CollectAlignmentSummaryMetrics` (docs/01 section 2.4), mirroring Picard
//! 3.5.0 `AlignmentSummaryMetricsCollector` at the `ALL_READS` level.
//!
//! Each primary record is described once (`Read::describe`) and the
//! description is added to its categories: `FIRST_OF_PAIR` or
//! `SECOND_OF_PAIR` plus `PAIR` for a paired read, `UNPAIRED` otherwise.

use std::io::{self, Write};
use std::path::{Path, PathBuf};

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::{bam, sam};

use crate::adapter::{Adapters, DEFAULT_ADAPTERS};
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metrics_file::{HistogramColumn, HistogramTable, MetricsFile, format_double};
use crate::pass::{self, Accumulator};
use crate::record::{Orientation, orientation, record_name};
use crate::reference::Reference;
use crate::sequence::{bases_equal, is_no_call};

/// `MAPPING_QUALITY_THRESHOLD`: MAPQ for a high-quality (HQ) alignment.
const HQ_MAPQ: u8 = 20;
/// `BASE_QUALITY_THRESHOLD`: base quality counted in `PF_HQ_ALIGNED_Q20_BASES`.
const Q20: i8 = 20;
/// `ChimeraUtil.DEFAULT_INSERT_SIZE_LIMIT`.
pub const DEFAULT_MAX_INSERT_SIZE: i32 = 100_000;
/// `ReservedTagConstants.XN`: the Illumina noise-read tag.
const NOISE: Tag = Tag::new(b'X', b'N');

#[derive(Clone, Debug)]
pub struct Options {
    /// `REFERENCE_SEQUENCE`. Without it no base is a mismatch, so the
    /// mismatch and error rates are 0, as Picard prints them.
    pub reference: Option<PathBuf>,
    /// `MAX_INSERT_SIZE`: a mapped pair with a larger |TLEN| is chimeric.
    pub max_insert_size: i32,
    /// `EXPECTED_PAIR_ORIENTATIONS`: a pair in any other orientation is chimeric.
    pub expected_orientations: Vec<Orientation>,
    /// `ADAPTER_SEQUENCE`.
    pub adapter_sequences: Vec<String>,
    /// `COLLECT_ALIGNMENT_INFORMATION`: false counts reads only (and no-calls
    /// for `BAD_CYCLES`); every alignment field stays 0.
    pub collect_alignment_information: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            reference: None,
            max_insert_size: DEFAULT_MAX_INSERT_SIZE,
            expected_orientations: vec![Orientation::Fr],
            adapter_sequences: DEFAULT_ADAPTERS.iter().map(|s| (*s).to_owned()).collect(),
            collect_alignment_information: true,
        }
    }
}

/// `AlignmentSummaryMetrics.Category`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Category {
    FirstOfPair,
    SecondOfPair,
    Pair,
    Unpaired,
}

impl Category {
    pub fn as_str(self) -> &'static str {
        match self {
            Category::FirstOfPair => "FIRST_OF_PAIR",
            Category::SecondOfPair => "SECOND_OF_PAIR",
            Category::Pair => "PAIR",
            Category::Unpaired => "UNPAIRED",
        }
    }
}

/// `AlignmentSummaryMetrics` columns in Picard 3.5.0 file order.
pub const COLUMNS: [&str; 36] = [
    "CATEGORY",
    "TOTAL_READS",
    "PF_READS",
    "PCT_PF_READS",
    "PF_NOISE_READS",
    "PF_READS_ALIGNED",
    "PCT_PF_READS_ALIGNED",
    "PF_ALIGNED_BASES",
    "PF_HQ_ALIGNED_READS",
    "PF_HQ_ALIGNED_BASES",
    "PF_HQ_ALIGNED_Q20_BASES",
    "PF_HQ_MEDIAN_MISMATCHES",
    "PF_MISMATCH_RATE",
    "PF_HQ_ERROR_RATE",
    "PF_INDEL_RATE",
    "MEAN_READ_LENGTH",
    "SD_READ_LENGTH",
    "MEDIAN_READ_LENGTH",
    "MAD_READ_LENGTH",
    "MIN_READ_LENGTH",
    "MAX_READ_LENGTH",
    "MEAN_ALIGNED_READ_LENGTH",
    "READS_ALIGNED_IN_PAIRS",
    "PCT_READS_ALIGNED_IN_PAIRS",
    "PF_READS_IMPROPER_PAIRS",
    "PCT_PF_READS_IMPROPER_PAIRS",
    "BAD_CYCLES",
    "STRAND_BALANCE",
    "PCT_CHIMERAS",
    "PCT_ADAPTER",
    "PCT_SOFTCLIP",
    "PCT_HARDCLIP",
    "AVG_POS_3PRIME_SOFTCLIP_LENGTH",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

/// One `AlignmentSummaryMetrics` row. Java `long` fields are `i64`, Java
/// `double` fields `f64`; a field Picard leaves unset is 0.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub category: Category,
    pub total_reads: i64,
    pub pf_reads: i64,
    pub pct_pf_reads: f64,
    pub pf_noise_reads: i64,
    pub pf_reads_aligned: i64,
    pub pct_pf_reads_aligned: f64,
    pub pf_aligned_bases: i64,
    pub pf_hq_aligned_reads: i64,
    pub pf_hq_aligned_bases: i64,
    pub pf_hq_aligned_q20_bases: i64,
    pub pf_hq_median_mismatches: f64,
    pub pf_mismatch_rate: f64,
    pub pf_hq_error_rate: f64,
    pub pf_indel_rate: f64,
    pub mean_read_length: f64,
    pub sd_read_length: f64,
    pub median_read_length: f64,
    pub mad_read_length: f64,
    pub min_read_length: f64,
    pub max_read_length: f64,
    pub mean_aligned_read_length: f64,
    pub reads_aligned_in_pairs: i64,
    pub pct_reads_aligned_in_pairs: f64,
    pub pf_reads_improper_pairs: i64,
    pub pct_pf_reads_improper_pairs: f64,
    pub bad_cycles: i64,
    pub strand_balance: f64,
    pub pct_chimeras: f64,
    pub pct_adapter: f64,
    pub pct_softclip: f64,
    pub pct_hardclip: f64,
    pub avg_pos_3prime_softclip_length: f64,
}

impl Row {
    /// The row as Picard prints it, in `COLUMNS` order. `SAMPLE`, `LIBRARY`
    /// and `READ_GROUP` are empty at the `ALL_READS` level.
    fn picard_fields(&self) -> Vec<String> {
        let d = |v: f64| format_double(v);
        vec![
            self.category.as_str().to_owned(),
            self.total_reads.to_string(),
            self.pf_reads.to_string(),
            d(self.pct_pf_reads),
            self.pf_noise_reads.to_string(),
            self.pf_reads_aligned.to_string(),
            d(self.pct_pf_reads_aligned),
            self.pf_aligned_bases.to_string(),
            self.pf_hq_aligned_reads.to_string(),
            self.pf_hq_aligned_bases.to_string(),
            self.pf_hq_aligned_q20_bases.to_string(),
            d(self.pf_hq_median_mismatches),
            d(self.pf_mismatch_rate),
            d(self.pf_hq_error_rate),
            d(self.pf_indel_rate),
            d(self.mean_read_length),
            d(self.sd_read_length),
            d(self.median_read_length),
            d(self.mad_read_length),
            d(self.min_read_length),
            d(self.max_read_length),
            d(self.mean_aligned_read_length),
            self.reads_aligned_in_pairs.to_string(),
            d(self.pct_reads_aligned_in_pairs),
            self.pf_reads_improper_pairs.to_string(),
            d(self.pct_pf_reads_improper_pairs),
            self.bad_cycles.to_string(),
            d(self.strand_balance),
            d(self.pct_chimeras),
            d(self.pct_adapter),
            d(self.pct_softclip),
            d(self.pct_hardclip),
            d(self.avg_pos_3prime_softclip_length),
            String::new(),
            String::new(),
            String::new(),
        ]
    }
}

/// The rows Picard writes and the read-length histograms of the `PAIR` and
/// `UNPAIRED` categories (`CollectAlignmentSummaryMetrics.finish`).
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub rows: Vec<Row>,
    /// `PAIRED_TOTAL_LENGTH_COUNT`: read lengths of PF paired reads.
    pub paired_read_lengths: Histogram,
    /// `PAIRED_ALIGNED_LENGTH_COUNT`: bases in M, I, = and X operations.
    pub paired_aligned_lengths: Histogram,
    /// `UNPAIRED_TOTAL_LENGTH_COUNT`.
    pub unpaired_read_lengths: Histogram,
    /// `UNPAIRED_ALIGNED_LENGTH_COUNT`.
    pub unpaired_aligned_lengths: Histogram,
}

/// Write a report as a Picard `AlignmentSummaryMetrics` file. Histograms
/// with no reads are left out (`MetricsFile.printHistogram`).
pub fn write_metrics(
    report: &Report,
    headers: Vec<String>,
    out: &mut impl Write,
) -> io::Result<()> {
    let column = |label: &str, histogram| HistogramColumn {
        label: label.to_owned(),
        histogram,
    };
    MetricsFile {
        headers,
        class: "picard.analysis.AlignmentSummaryMetrics",
        columns: &COLUMNS,
        rows: report.rows.iter().map(Row::picard_fields).collect(),
        histogram: Some(HistogramTable {
            key_label: "READ_LENGTH",
            columns: vec![
                column("PAIRED_TOTAL_LENGTH_COUNT", &report.paired_read_lengths),
                column(
                    "PAIRED_ALIGNED_LENGTH_COUNT",
                    &report.paired_aligned_lengths,
                ),
                column("UNPAIRED_TOTAL_LENGTH_COUNT", &report.unpaired_read_lengths),
                column(
                    "UNPAIRED_ALIGNED_LENGTH_COUNT",
                    &report.unpaired_aligned_lengths,
                ),
            ],
        }),
    }
    .write(out)
}

/// What one primary record adds to a category. `describe` fills it once per
/// record; `Collector::add` applies it to each category the record is in.
#[derive(Default)]
struct Read {
    qc_fail: bool,
    noise: bool,
    length: i64,
    /// `getUnclippedBaseCount`: bases in M, I, = and X operations.
    aligned_length: i64,
    adapter: bool,
    hard_clipped: i64,
    /// Present for a mapped read: soft-clipped bases and 3' soft-clipped bases.
    clips: Option<(i64, i64)>,
    /// Present for a mapped read when alignment information is collected.
    alignment: Option<Alignment>,
    /// Present when the read takes the mismatch branch of `collectQualityData`.
    bases: Option<BaseCounts>,
    /// Keys added to the bad-cycle histogram.
    bad_cycles: Vec<i64>,
}

struct Alignment {
    improper_pair: bool,
    positive_strand: bool,
    aligned_in_pair: bool,
    /// `Some(chimeric)` when the read enters the `PCT_CHIMERAS` denominator.
    chimera: Option<bool>,
}

#[derive(Default)]
struct BaseCounts {
    hq_read: bool,
    aligned: i64,
    hq_aligned: i64,
    hq_q20: i64,
    mismatches: i64,
    hq_mismatches: i64,
    indels: i64,
}

/// Settings and scratch space shared by every record.
struct Describer {
    adapters: Adapters,
    max_insert_size: i64,
    expected: Vec<Orientation>,
    do_ref_metrics: bool,
    bases: Vec<u8>,
    cigar: Vec<(Kind, usize)>,
}

impl Describer {
    fn describe(&mut self, record: &bam::Record, reference: Option<&[u8]>) -> Result<Read> {
        let flags = record.flags();
        // htsjdk `getMappingQuality` reports a missing MAPQ as 255.
        let mapq = record.mapping_quality().map_or(255, |q| q.get());
        let unmapped = flags.is_unmapped();
        let reverse = flags.is_reverse_complemented();
        self.bases.clear();
        self.bases.extend(record.sequence().iter());
        self.cigar.clear();
        for op in record.cigar().iter() {
            let op = op?;
            self.cigar.push((op.kind(), op.len()));
        }
        let mut read = Read {
            qc_fail: flags.is_qc_fail(),
            noise: is_noise(record)?,
            length: self.bases.len() as i64,
            aligned_length: self
                .cigar
                .iter()
                .filter(|(k, _)| {
                    matches!(
                        k,
                        Kind::Match
                            | Kind::Insertion
                            | Kind::SequenceMatch
                            | Kind::SequenceMismatch
                    )
                })
                .map(|&(_, n)| n as i64)
                .sum(),
            adapter: self
                .adapters
                .is_adapter(&self.bases, unmapped, reverse, mapq),
            hard_clipped: self.total(Kind::HardClip),
            ..Read::default()
        };
        if !unmapped {
            read.clips = Some((
                self.total(Kind::SoftClip),
                three_prime_soft_clip(&self.cigar, reverse),
            ));
            if self.do_ref_metrics {
                read.alignment = Some(self.alignment(record, mapq)?);
            }
        }
        if unmapped || read.qc_fail || !self.do_ref_metrics {
            // `collectQualityData` on a read that is not an aligned PF read:
            // every no-call counts toward its machine cycle.
            for (i, &base) in self.bases.iter().enumerate() {
                if is_no_call(base) {
                    read.bad_cycles.push(cycle(reverse, self.bases.len(), i));
                }
            }
        } else {
            read.bases = Some(self.base_counts(record, reference, mapq, &mut read.bad_cycles)?);
        }
        Ok(read)
    }

    /// Summed length of one CIGAR operation (`getTotalCigarOperatorCount`).
    fn total(&self, kind: Kind) -> i64 {
        self.cigar
            .iter()
            .filter(|(k, _)| *k == kind)
            .map(|&(_, n)| n as i64)
            .sum()
    }

    /// The `doRefMetrics` block of `collectReadData` for a mapped PF read.
    fn alignment(&self, record: &bam::Record, mapq: u8) -> Result<Alignment> {
        let flags = record.flags();
        let paired = flags.is_segmented();
        let aligned_in_pair = paired && !flags.is_mate_unmapped();
        let chimera = if aligned_in_pair {
            // Both ends must pass MAPQ 20; a missing MQ tag passes.
            let mate_mapq = mate_mapping_quality(record)?;
            if mate_mapq.is_none_or(|mq| mq >= i64::from(HQ_MAPQ) && mapq >= HQ_MAPQ) {
                Some(self.is_chimeric(record)?)
            } else {
                None
            }
        } else if mapq >= HQ_MAPQ {
            // A fragment, or a pair with one end mapped: chimeric when split
            // (an SA tag).
            Some(has_sa(record))
        } else {
            None
        };
        Ok(Alignment {
            improper_pair: paired && !flags.is_properly_segmented(),
            positive_strand: !flags.is_reverse_complemented(),
            aligned_in_pair,
            chimera,
        })
    }

    /// `ChimeraUtil.isChimeric(rec, maxInsertSize, expectedOrientations)` for
    /// a read whose pair is mapped: |TLEN| above the limit, mates on two
    /// contigs, or an orientation outside the expected set or an SA tag.
    /// The orientation is computed only for mates on one contig.
    fn is_chimeric(&self, record: &bam::Record) -> Result<bool> {
        if i64::from(record.template_length()).abs() > self.max_insert_size {
            return Ok(true);
        }
        let contig = record.reference_sequence_id().transpose()?;
        let mate_contig = record.mate_reference_sequence_id().transpose()?;
        if contig != mate_contig {
            return Ok(true);
        }
        Ok(!self.expected.contains(&orientation(record)?) || has_sa(record))
    }

    /// The mismatch branch of `collectQualityData` for an aligned PF read.
    /// Each alignment block compares read and reference base by base, up to
    /// the end of the contig. Mismatches and no-calls count toward
    /// `getCycle(negativeStrand, readLength, i)` where `i` is the offset in
    /// the *block*, not in the read: Picard passes the block loop variable,
    /// so every block starts again at cycle 1 (or the last cycle when reverse).
    fn base_counts(
        &self,
        record: &bam::Record,
        reference: Option<&[u8]>,
        mapq: u8,
        bad_cycles: &mut Vec<i64>,
    ) -> Result<BaseCounts> {
        let malformed = |what| Error::MalformedRecord {
            name: record_name(record),
            what,
        };
        let reverse = record.flags().is_reverse_complemented();
        let hq = mapq >= HQ_MAPQ;
        let qualities = record.quality_scores();
        let qualities = qualities.as_ref();
        // htsjdk decodes a 0xFF quality string as no qualities.
        let qualities = if qualities.first() == Some(&0xFF) {
            &[][..]
        } else {
            qualities
        };
        let start = record
            .alignment_start()
            .transpose()?
            .ok_or_else(|| malformed("mapped read without an alignment start"))?;
        let mut counts = BaseCounts {
            hq_read: hq,
            ..BaseCounts::default()
        };
        let ref_length = reference.map_or(usize::MAX, <[u8]>::len);
        let (mut read_pos, mut ref_pos) = (0usize, start.get() - 1);
        for &(kind, len) in &self.cigar {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for i in 0..len {
                        if ref_pos + i >= ref_length {
                            break;
                        }
                        let read_base = *self
                            .bases
                            .get(read_pos + i)
                            .ok_or_else(|| malformed("CIGAR longer than the read"))?;
                        let mismatch =
                            reference.is_some_and(|r| !bases_equal(read_base, r[ref_pos + i]));
                        counts.aligned += 1;
                        if hq {
                            counts.hq_aligned += 1;
                            let quality = *qualities
                                .get(read_pos + i)
                                .ok_or_else(|| malformed("aligned read without base qualities"))?;
                            // A Java byte: qualities above 127 are negative.
                            if quality as i8 >= Q20 {
                                counts.hq_q20 += 1;
                            }
                            counts.hq_mismatches += i64::from(mismatch);
                        }
                        counts.mismatches += i64::from(mismatch);
                        if mismatch || is_no_call(read_base) {
                            bad_cycles.push(cycle(reverse, self.bases.len(), i));
                        }
                    }
                    read_pos += len;
                    ref_pos += len;
                }
                Kind::Insertion => {
                    counts.indels += 1;
                    read_pos += len;
                }
                Kind::Deletion => {
                    counts.indels += 1;
                    ref_pos += len;
                }
                Kind::SoftClip => read_pos += len,
                Kind::Skip => ref_pos += len,
                Kind::HardClip | Kind::Pad => {}
            }
        }
        Ok(counts)
    }
}

/// `CoordMath.getCycle`: 1-based machine cycle of read offset `i`.
fn cycle(reverse: bool, read_length: usize, i: usize) -> i64 {
    if reverse {
        read_length as i64 - i as i64
    } else {
        i as i64 + 1
    }
}

/// `get3PrimeSoftClippedBases`: soft clips after the first non-clipping
/// operation, reading the CIGAR from the 3' end's side (reversed for a
/// reverse-strand read). A CIGAR of clips only gives 0.
fn three_prime_soft_clip(cigar: &[(Kind, usize)], reverse: bool) -> i64 {
    let mut seen_aligned = false;
    let mut clipped = 0i64;
    let ops: Box<dyn Iterator<Item = &(Kind, usize)>> = if reverse {
        Box::new(cigar.iter().rev())
    } else {
        Box::new(cigar.iter())
    };
    for &(kind, len) in ops {
        match kind {
            Kind::SoftClip if seen_aligned => clipped += len as i64,
            Kind::SoftClip | Kind::HardClip => {}
            _ => seen_aligned = true,
        }
    }
    clipped
}

/// `isNoiseRead`: an `XN` tag equal to `Integer` 1. Any BAM integer type
/// decodes to a Java `Integer` (or `Long` above `Integer.MAX_VALUE`); a
/// float or string 1 is not noise.
fn is_noise(record: &bam::Record) -> Result<bool> {
    match record.data().get(&NOISE) {
        None => Ok(false),
        Some(value) => Ok(value?.as_int() == Some(1)),
    }
}

fn has_sa(record: &bam::Record) -> bool {
    record.data().get(&Tag::OTHER_ALIGNMENTS).is_some()
}

/// `getIntegerAttribute("MQ")`: `None` when absent, an error when the tag
/// is not an integer.
fn mate_mapping_quality(record: &bam::Record) -> Result<Option<i64>> {
    match record.data().get(&Tag::MATE_MAPPING_QUALITY) {
        None => Ok(None),
        Some(value) => value?
            .as_int()
            .map(Some)
            .ok_or_else(|| Error::MalformedRecord {
                name: record_name(record),
                what: "MQ tag is not an integer",
            }),
    }
}

/// `IndividualAlignmentSummaryMetricsCollector` for one category.
#[derive(Default)]
struct Collector {
    total_reads: i64,
    pf_reads: i64,
    pf_noise_reads: i64,
    pf_reads_aligned: i64,
    pf_reads_improper_pairs: i64,
    reads_aligned_in_pairs: i64,
    positive_strand: i64,
    chimeras: i64,
    chimeras_denominator: i64,
    adapter_reads: i64,
    soft_clipped: i64,
    three_prime_soft_clipped: i64,
    reads_with_three_prime_soft_clips: i64,
    hard_clipped: i64,
    pf_hq_aligned_reads: i64,
    pf_aligned_bases: i64,
    pf_hq_aligned_bases: i64,
    pf_hq_aligned_q20_bases: i64,
    indels: i64,
    read_lengths: Histogram,
    aligned_lengths: Histogram,
    mismatches: Histogram,
    hq_mismatches: Histogram,
    bad_cycles: Histogram,
}

impl Collector {
    /// `collectReadData` then `collectQualityData`.
    fn add(&mut self, read: &Read) {
        self.total_reads += 1;
        if !read.qc_fail {
            self.pf_reads += 1;
            self.pf_noise_reads += i64::from(read.noise);
            self.read_lengths.increment(read.length);
            self.aligned_lengths.increment(read.aligned_length);
            self.adapter_reads += i64::from(read.adapter);
            self.hard_clipped += read.hard_clipped;
            if let Some((soft, three_prime)) = read.clips {
                self.soft_clipped += soft;
                if three_prime > 0 {
                    self.three_prime_soft_clipped += three_prime;
                    self.reads_with_three_prime_soft_clips += 1;
                }
            }
            if let Some(a) = &read.alignment {
                self.pf_reads_aligned += 1;
                self.pf_reads_improper_pairs += i64::from(a.improper_pair);
                self.positive_strand += i64::from(a.positive_strand);
                self.reads_aligned_in_pairs += i64::from(a.aligned_in_pair);
                if let Some(chimeric) = a.chimera {
                    self.chimeras_denominator += 1;
                    self.chimeras += i64::from(chimeric);
                }
            }
        }
        for &c in &read.bad_cycles {
            self.bad_cycles.increment(c);
        }
        if let Some(b) = &read.bases {
            self.pf_hq_aligned_reads += i64::from(b.hq_read);
            self.pf_aligned_bases += b.aligned;
            self.pf_hq_aligned_bases += b.hq_aligned;
            self.pf_hq_aligned_q20_bases += b.hq_q20;
            self.mismatches.increment(b.mismatches);
            self.hq_mismatches.increment(b.hq_mismatches);
            self.indels += b.indels;
        }
    }

    /// `IndividualAlignmentSummaryMetricsCollector.finish`. Ratios go
    /// through `MathUtil.divide` (0 for a zero denominator) except
    /// `PCT_PF_READS` and `PCT_ADAPTER`, which divide directly. With no
    /// bisulfite reads every aligned base is a non-bisulfite base.
    fn finish(&self, category: Category, do_ref_metrics: bool) -> Result<Row> {
        let mut row = Row {
            category,
            total_reads: self.total_reads,
            pf_reads: self.pf_reads,
            pct_pf_reads: 0.0,
            pf_noise_reads: self.pf_noise_reads,
            pf_reads_aligned: self.pf_reads_aligned,
            pct_pf_reads_aligned: 0.0,
            pf_aligned_bases: self.pf_aligned_bases,
            pf_hq_aligned_reads: self.pf_hq_aligned_reads,
            pf_hq_aligned_bases: self.pf_hq_aligned_bases,
            pf_hq_aligned_q20_bases: self.pf_hq_aligned_q20_bases,
            pf_hq_median_mismatches: 0.0,
            pf_mismatch_rate: 0.0,
            pf_hq_error_rate: 0.0,
            pf_indel_rate: 0.0,
            mean_read_length: 0.0,
            sd_read_length: 0.0,
            median_read_length: 0.0,
            mad_read_length: 0.0,
            min_read_length: 0.0,
            max_read_length: 0.0,
            mean_aligned_read_length: 0.0,
            reads_aligned_in_pairs: self.reads_aligned_in_pairs,
            pct_reads_aligned_in_pairs: 0.0,
            pf_reads_improper_pairs: self.pf_reads_improper_pairs,
            pct_pf_reads_improper_pairs: 0.0,
            bad_cycles: 0,
            strand_balance: 0.0,
            pct_chimeras: 0.0,
            pct_adapter: 0.0,
            pct_softclip: 0.0,
            pct_hardclip: 0.0,
            avg_pos_3prime_softclip_length: 0.0,
        };
        if self.pf_reads == 0 {
            if self.total_reads > 0 {
                return Err(Error::Usage(format!(
                    "Input file contains no PF_READS ({} reads in {}, all QC-failed)",
                    self.total_reads,
                    category.as_str()
                )));
            }
            return Ok(row);
        }
        let total = self.total_reads as f64;
        let pf = self.pf_reads as f64;
        let h = &self.read_lengths;
        row.pct_pf_reads = pf / total;
        row.pct_adapter = self.adapter_reads as f64 / pf;
        row.mean_read_length = h.mean();
        row.sd_read_length = h.standard_deviation();
        row.median_read_length = h.median();
        row.mad_read_length = h.median_absolute_deviation();
        row.min_read_length = h.min().unwrap_or_default() as f64;
        row.max_read_length = h.max().unwrap_or_default() as f64;
        // A cycle is bad when no-calls and mismatches there reach 80% of
        // TOTAL_READS (QC-failed reads included in both).
        row.bad_cycles = self
            .bad_cycles
            .keys()
            .filter(|&k| self.bad_cycles.get(k).unwrap_or_default() / total >= 0.8)
            .count() as i64;
        if do_ref_metrics {
            let aligned = self.pf_reads_aligned as f64;
            let total_bases = h.sum();
            row.pct_pf_reads_aligned = divide(aligned, pf);
            row.pct_reads_aligned_in_pairs = divide(self.reads_aligned_in_pairs as f64, aligned);
            row.pct_pf_reads_improper_pairs = divide(self.pf_reads_improper_pairs as f64, aligned);
            row.mean_aligned_read_length = self.aligned_lengths.mean();
            row.strand_balance = divide(self.positive_strand as f64, aligned);
            row.pct_chimeras = divide(self.chimeras as f64, self.chimeras_denominator as f64);
            row.pf_indel_rate = divide(self.indels as f64, self.pf_aligned_bases as f64);
            row.pf_mismatch_rate = divide(self.mismatches.sum(), self.pf_aligned_bases as f64);
            row.pf_hq_error_rate =
                divide(self.hq_mismatches.sum(), self.pf_hq_aligned_bases as f64);
            row.pct_hardclip = divide(self.hard_clipped as f64, total_bases);
            row.pct_softclip = divide(self.soft_clipped as f64, total_bases);
            row.avg_pos_3prime_softclip_length = divide(
                self.three_prime_soft_clipped as f64,
                self.reads_with_three_prime_soft_clips as f64,
            );
            row.pf_hq_median_mismatches = self.hq_mismatches.median();
        }
        Ok(row)
    }
}

/// Picard `MathUtil.divide`: 0 when the denominator is within 1e-6 of 0.
fn divide(numerator: f64, denominator: f64) -> f64 {
    if (0.0 - denominator).abs() > 0.000_001 {
        numerator / denominator
    } else {
        0.0
    }
}

/// `CollectAlignmentSummaryMetrics` as a pass [`Accumulator`]: one
/// collector per category, plus the reference walker when `R` is given.
pub struct AlignmentSummary {
    reference: Option<Reference>,
    describer: Describer,
    first: Collector,
    second: Collector,
    pair: Collector,
    unpaired: Collector,
}

impl AlignmentSummary {
    /// Opens the reference, when there is one, and checks it against the
    /// BAM header.
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        let reference = opts
            .reference
            .as_deref()
            .map(|r| Reference::open(r, header))
            .transpose()?;
        Ok(Self {
            reference,
            describer: Describer {
                adapters: Adapters::new(&opts.adapter_sequences),
                max_insert_size: i64::from(opts.max_insert_size),
                expected: opts.expected_orientations.clone(),
                do_ref_metrics: opts.collect_alignment_information,
                bases: Vec::new(),
                cigar: Vec::new(),
            },
            first: Collector::default(),
            second: Collector::default(),
            pair: Collector::default(),
            unpaired: Collector::default(),
        })
    }
}

impl Accumulator for AlignmentSummary {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        // `SinglePassSamProgram` fetches the contig of every placed record,
        // secondary and unmapped ones included, so the sort check covers them.
        let bases = match (
            &mut self.reference,
            record.reference_sequence_id().transpose()?,
        ) {
            (Some(r), Some(id)) => Some(r.contig(id)?),
            _ => None,
        };
        let flags = record.flags();
        // `AlignmentSummaryMetricsCollector.acceptRecord` drops secondary and
        // supplementary records before any counting.
        if flags.is_secondary() || flags.is_supplementary() {
            return Ok(());
        }
        let read = self.describer.describe(record, bases)?;
        if flags.is_segmented() {
            if flags.is_first_segment() {
                self.first.add(&read);
            } else {
                self.second.add(&read);
            }
            self.pair.add(&read);
        } else {
            self.unpaired.add(&read);
        }
        Ok(())
    }

    fn finish(self) -> Result<Report> {
        let do_ref = self.describer.do_ref_metrics;
        let first_row = self.first.finish(Category::FirstOfPair, do_ref)?;
        let second_row = self.second.finish(Category::SecondOfPair, do_ref)?;
        let mut pair_row = self.pair.finish(Category::Pair, do_ref)?;
        let unpaired_row = self.unpaired.finish(Category::Unpaired, do_ref)?;

        // `GroupAlignmentSummaryMetricsPerUnitMetricCollector.addMetricsToFile`:
        // paired rows only when FIRST_OF_PAIR has reads, with PAIR's BAD_CYCLES
        // the sum of the two ends; UNPAIRED when it has reads or when nothing
        // else is written (an all-zero row for an empty file).
        let mut rows = Vec::with_capacity(4);
        let has_pairs = first_row.total_reads > 0;
        if has_pairs {
            pair_row.bad_cycles = first_row.bad_cycles + second_row.bad_cycles;
            rows.extend([first_row, second_row, pair_row]);
        }
        if unpaired_row.total_reads > 0 || !has_pairs {
            rows.push(unpaired_row);
        }
        Ok(Report {
            rows,
            paired_read_lengths: self.pair.read_lengths,
            paired_aligned_lengths: self.pair.aligned_lengths,
            unpaired_read_lengths: self.unpaired.read_lengths,
            unpaired_aligned_lengths: self.unpaired.aligned_lengths,
        })
    }
}

/// `StringHeader` command line in Picard's shape.
pub fn command_header(opts: &Options, bam: &Path, output: Option<&Path>) -> String {
    let mut command = format!("CollectAlignmentSummaryMetrics INPUT={}", bam.display());
    if let Some(output) = output {
        command.push_str(&format!(" OUTPUT={}", output.display()));
    }
    if let Some(reference) = &opts.reference {
        command.push_str(&format!(" REFERENCE_SEQUENCE={}", reference.display()));
    }
    let orientations: Vec<&str> = opts
        .expected_orientations
        .iter()
        .map(|o| o.as_str())
        .collect();
    command.push_str(&format!(
        " MAX_INSERT_SIZE={} EXPECTED_PAIR_ORIENTATIONS=[{}] ADAPTER_SEQUENCE=[{}] METRIC_ACCUMULATION_LEVEL=[ALL_READS] IS_BISULFITE_SEQUENCED=false COLLECT_ALIGNMENT_INFORMATION={}",
        opts.max_insert_size,
        orientations.join(", "),
        opts.adapter_sequences.join(", "),
        opts.collect_alignment_information
    ));
    command
}

/// Compute the metrics of one BAM with the default BGZF worker count.
pub fn run_over_bam(path: impl AsRef<Path>, opts: &Options) -> Result<Report> {
    pass::run(path.as_ref(), pass::default_threads(), |header| {
        AlignmentSummary::new(opts, header)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn three_prime_soft_clip_follows_strand() {
        let cigar = [
            (Kind::SoftClip, 5),
            (Kind::Match, 90),
            (Kind::SoftClip, 3),
            (Kind::HardClip, 2),
        ];
        assert_eq!(three_prime_soft_clip(&cigar, false), 3);
        assert_eq!(three_prime_soft_clip(&cigar, true), 5);
        assert_eq!(three_prime_soft_clip(&[(Kind::SoftClip, 10)], false), 0);
    }

    #[test]
    fn cycles_count_from_the_machine_start() {
        assert_eq!(cycle(false, 100, 0), 1);
        assert_eq!(cycle(true, 100, 0), 100);
        assert_eq!(cycle(true, 100, 99), 1);
    }

    #[test]
    fn math_util_divide_guards_zero() {
        assert_eq!(divide(3.0, 0.0), 0.0);
        assert_eq!(divide(3.0, 0.000_000_5), 0.0);
        assert_eq!(divide(3.0, 4.0), 0.75);
    }

    #[test]
    fn empty_category_is_an_all_zero_row() {
        let row = Collector::default()
            .finish(Category::Unpaired, true)
            .expect("no reads is not an error");
        assert_eq!(row.total_reads, 0);
        assert_eq!(row.mean_read_length, 0.0);
        let mut out = Vec::new();
        let report = Report {
            rows: vec![row],
            paired_read_lengths: Histogram::new(),
            paired_aligned_lengths: Histogram::new(),
            unpaired_read_lengths: Histogram::new(),
            unpaired_aligned_lengths: Histogram::new(),
        };
        write_metrics(&report, Vec::new(), &mut out).expect("write to Vec");
        let text = String::from_utf8(out).expect("utf8");
        assert!(text.ends_with("\nUNPAIRED\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t\t\t\n\n\n"));
    }

    #[test]
    fn qc_failed_only_category_is_an_error() {
        let mut collector = Collector::default();
        collector.add(&Read {
            qc_fail: true,
            ..Read::default()
        });
        assert!(matches!(
            collector.finish(Category::Unpaired, true),
            Err(Error::Usage(m)) if m.contains("no PF_READS")
        ));
    }

    #[test]
    fn bad_cycle_needs_eighty_percent_of_all_reads() {
        let mut collector = Collector::default();
        for i in 0..5 {
            collector.add(&Read {
                length: 10,
                bad_cycles: if i < 4 { vec![3, 7] } else { vec![3] },
                ..Read::default()
            });
        }
        let row = collector.finish(Category::Unpaired, true).expect("finish");
        // cycle 3: 5/5, cycle 7: 4/5 = 0.8
        assert_eq!(row.bad_cycles, 2);
    }
}
