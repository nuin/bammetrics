// bammetrics patch: reachable through `io::writer::raw`, as noodles wrote it.
#![allow(missing_docs)]

mod convert;
mod feature;

use bstr::BString;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::{MappingQuality, data::field::Tag},
        record_buf::data::field::Value,
    },
};

pub use self::feature::Feature;
use crate::record::{Flags, MateFlags};

#[derive(Clone, Debug, Default, PartialEq)]
// bammetrics patch: the fields are public, for `io::writer::raw`.
pub struct Record {
    pub bam_flags: sam::alignment::record::Flags,
    pub cram_flags: Flags,
    pub reference_sequence_id: Option<usize>,
    pub read_length: usize,
    pub alignment_start: Option<Position>,
    pub read_group_id: Option<usize>,
    pub name: Option<BString>,
    pub mate_flags: MateFlags,
    pub mate_reference_sequence_id: Option<usize>,
    pub mate_alignment_start: Option<Position>,
    pub template_length: i32,
    pub mate_distance: Option<usize>,
    pub data: Vec<(Tag, Value)>,
    pub features: Vec<Feature>,
    pub mapping_quality: Option<MappingQuality>,
    pub sequence: Vec<u8>,
    pub quality_scores: Vec<u8>,
}

impl Record {
    pub fn alignment_span(&self) -> usize {
        calculate_alignment_span(self.read_length, &self.features)
    }

    pub fn alignment_end(&self) -> Option<Position> {
        self.alignment_start.and_then(|start| {
            let end = usize::from(start) + self.alignment_span() - 1;
            Position::new(end)
        })
    }
}

fn calculate_alignment_span(read_length: usize, features: &[Feature]) -> usize {
    features
        .iter()
        .fold(read_length, |span, feature| match feature {
            Feature::Insertion { bases, .. } => span - bases.len(),
            Feature::InsertBase { .. } => span - 1,
            Feature::Deletion { len, .. } => span + len,
            Feature::ReferenceSkip { len, .. } => span + len,
            Feature::SoftClip { bases, .. } => span - bases.len(),
            _ => span,
        })
}
