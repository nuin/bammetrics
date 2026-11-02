//! Per-record description for `CollectAlignmentSummaryMetrics`: what one
//! primary record adds to each category it is in (`collectReadData` and
//! `collectQualityData` of `IndividualAlignmentSummaryMetricsCollector`).

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;

use crate::adapter::Adapters;
use crate::error::{Error, Result};
use crate::record::{Orientation, orientation, record_name};
use crate::sequence::{bases_equal, bisulfite_bases_equal, is_no_call};

/// `MAPPING_QUALITY_THRESHOLD`: MAPQ for a high-quality (HQ) alignment.
const HQ_MAPQ: u8 = 20;
/// `BASE_QUALITY_THRESHOLD`: base quality counted in `PF_HQ_ALIGNED_Q20_BASES`.
const Q20: i8 = 20;
/// `ReservedTagConstants.XN`: the Illumina noise-read tag.
const NOISE: Tag = Tag::new(b'X', b'N');

/// What one primary record adds to a category. `describe` fills it once per
/// record; `Collector::add` applies it to each category the record is in.
#[derive(Default)]
pub(super) struct Read {
    pub(super) qc_fail: bool,
    pub(super) noise: bool,
    pub(super) length: i64,
    /// `getUnclippedBaseCount`: bases in M, I, = and X operations.
    pub(super) aligned_length: i64,
    pub(super) adapter: bool,
    pub(super) hard_clipped: i64,
    /// Present for a mapped read: soft-clipped bases and 3' soft-clipped bases.
    pub(super) clips: Option<(i64, i64)>,
    /// Present for a mapped read when alignment information is collected.
    pub(super) alignment: Option<Alignment>,
    /// Present when the read takes the mismatch branch of `collectQualityData`.
    pub(super) bases: Option<BaseCounts>,
    /// Keys added to the bad-cycle histogram.
    pub(super) bad_cycles: Vec<i64>,
}

pub(super) struct Alignment {
    pub(super) improper_pair: bool,
    pub(super) positive_strand: bool,
    pub(super) aligned_in_pair: bool,
    /// `Some(chimeric)` when the read enters the `PCT_CHIMERAS` denominator.
    pub(super) chimera: Option<bool>,
}

#[derive(Default)]
pub(super) struct BaseCounts {
    pub(super) hq_read: bool,
    pub(super) aligned: i64,
    /// Aligned bases that are not bisulfite conversions: the
    /// `PF_MISMATCH_RATE` denominator.
    pub(super) non_bisulfite: i64,
    pub(super) hq_aligned: i64,
    /// HQ aligned bases that are not bisulfite conversions: the
    /// `PF_HQ_ERROR_RATE` denominator.
    pub(super) hq_non_bisulfite: i64,
    pub(super) hq_q20: i64,
    pub(super) mismatches: i64,
    pub(super) hq_mismatches: i64,
    pub(super) indels: i64,
}

/// Settings and scratch space shared by every record.
pub(super) struct Describer {
    pub(super) adapters: Adapters,
    pub(super) max_insert_size: i64,
    pub(super) expected: Vec<Orientation>,
    pub(super) do_ref_metrics: bool,
    pub(super) bisulfite: bool,
    pub(super) bases: Vec<u8>,
    pub(super) cigar: Vec<(Kind, usize)>,
}

impl Describer {
    pub(super) fn describe(
        &mut self,
        record: &bam::Record,
        reference: Option<&[u8]>,
    ) -> Result<Read> {
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
    ///
    /// With `IS_BISULFITE_SEQUENCED` Picard tests the conversion against the
    /// reference base at the read offset counted from the contig start
    /// (`refBases[readBaseIndex]`), not the aligned reference base. That
    /// index is reproduced here, including its failure when it runs past the
    /// contig end.
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
                        let read_index = read_pos + i;
                        let read_base = *self
                            .bases
                            .get(read_index)
                            .ok_or_else(|| malformed("CIGAR longer than the read"))?;
                        let mut mismatch =
                            reference.is_some_and(|r| !bases_equal(read_base, r[ref_pos + i]));
                        let bisulfite_match = match reference {
                            Some(r) if self.bisulfite => {
                                let shifted = *r.get(read_index).ok_or_else(|| {
                                    malformed("bisulfite check reads past the contig end")
                                })?;
                                bisulfite_bases_equal(reverse, read_base, shifted)
                            }
                            _ => false,
                        };
                        let bisulfite_base = mismatch && bisulfite_match;
                        mismatch &= !bisulfite_match;
                        counts.aligned += 1;
                        counts.non_bisulfite += i64::from(!bisulfite_base);
                        if hq {
                            counts.hq_aligned += 1;
                            counts.hq_non_bisulfite += i64::from(!bisulfite_base);
                            let quality = *qualities
                                .get(read_index)
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
}
