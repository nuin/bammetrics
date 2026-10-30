//! Record filters shared by the metrics. Each field is one Picard behaviour
//! (a `SamRecordFilter`, or an `INCLUDE_*` option); the default accepts
//! every record.

use noodles::bam;

/// Which records a metric counts.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RecordFilter {
    /// Drop records with flag 0x4 (`AlignedFilter`).
    pub skip_unmapped: bool,
    /// Drop flag 0x100 (`NotPrimaryAlignmentFilter`).
    pub skip_secondary: bool,
    /// Drop flag 0x800.
    pub skip_supplementary: bool,
    /// Drop flag 0x400 (`INCLUDE_DUPLICATES=false`).
    pub skip_duplicates: bool,
    /// Drop flag 0x200 (vendor-failed reads).
    pub skip_qc_fail: bool,
    /// Drop records below this mapping quality; 0 keeps all. Unavailable
    /// mapping quality (255) always passes, as in `MappingQualityFilter`.
    pub minimum_mapping_quality: u8,
}

impl RecordFilter {
    pub fn accepts(&self, record: &bam::Record) -> bool {
        let f = record.flags();
        if (self.skip_unmapped && f.is_unmapped())
            || (self.skip_secondary && f.is_secondary())
            || (self.skip_supplementary && f.is_supplementary())
            || (self.skip_duplicates && f.is_duplicate())
            || (self.skip_qc_fail && f.is_qc_fail())
        {
            return false;
        }
        record
            .mapping_quality()
            .is_none_or(|q| u8::from(q) >= self.minimum_mapping_quality)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record::MappingQuality;
    use noodles::sam::alignment::record_buf::RecordBuf;

    /// A record with these flags and mapping quality, round-tripped through BAM bytes.
    fn record(flags: Flags, mapq: Option<u8>) -> bam::Record {
        let mut buf = RecordBuf::default();
        *buf.flags_mut() = flags;
        *buf.mapping_quality_mut() = mapq.and_then(MappingQuality::new);
        let header = noodles::sam::Header::default();
        let mut writer = bam::io::Writer::new(Vec::new());
        writer.write_header(&header).expect("header");
        writer
            .write_alignment_record(&header, &buf)
            .expect("encode");
        let bytes = writer.into_inner().finish().expect("finish");
        let mut reader = bam::io::Reader::new(bytes.as_slice());
        reader.read_header().expect("header");
        let mut record = bam::Record::default();
        reader.read_record(&mut record).expect("decode");
        record
    }

    #[test]
    fn default_accepts_everything() {
        let all = Flags::UNMAPPED | Flags::SECONDARY | Flags::DUPLICATE | Flags::QC_FAIL;
        assert!(RecordFilter::default().accepts(&record(all, Some(0))));
    }

    #[test]
    fn each_flag_filter_drops_its_flag_only() {
        let filter = RecordFilter {
            skip_duplicates: true,
            ..RecordFilter::default()
        };
        assert!(!filter.accepts(&record(Flags::DUPLICATE, Some(60))));
        assert!(filter.accepts(&record(Flags::SECONDARY, Some(60))));
    }

    #[test]
    fn unavailable_mapping_quality_passes() {
        let filter = RecordFilter {
            minimum_mapping_quality: 30,
            ..RecordFilter::default()
        };
        assert!(filter.accepts(&record(Flags::empty(), None)));
        assert!(filter.accepts(&record(Flags::empty(), Some(30))));
        assert!(!filter.accepts(&record(Flags::empty(), Some(29))));
    }
}
