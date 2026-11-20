use super::*;
use noodles::core::Position;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};

const DEFAULT: Thresholds = Thresholds {
    c_quality: 20,
    next_base_quality: 10,
    minimum_read_length: 5,
    max_mismatch_rate: 0.1,
};

/// A mapped record at 1-based `start` with `cigar`, its bases and qualities.
fn record(
    start: usize,
    cigar: &[(Kind, usize)],
    bases: &[u8],
    quals: &[u8],
    reverse: bool,
) -> bam::Record {
    let mut buf = RecordBuf::default();
    *buf.flags_mut() = if reverse {
        Flags::REVERSE_COMPLEMENTED
    } else {
        Flags::empty()
    };
    *buf.reference_sequence_id_mut() = Some(0);
    *buf.alignment_start_mut() = Position::new(start);
    *buf.cigar_mut() = Cigar::from(
        cigar
            .iter()
            .map(|&(k, n)| Op::new(k, n))
            .collect::<Vec<_>>(),
    );
    *buf.sequence_mut() = Sequence::from(bases.to_vec());
    *buf.quality_scores_mut() = QualityScores::from(quals.to_vec());
    let header: sam::Header = "@SQ\tSN:c\tLN:100\n".parse().expect("header");
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

fn observed(record: &bam::Record, contig: &[u8], t: Thresholds) -> Observation {
    let mut o = Observation::default();
    observe(record, contig, t, &mut Scratch::default(), &mut o).expect("observes");
    o
}

#[test]
fn forward_cpgs_count_at_the_zero_based_c_and_t_is_converted() {
    // Reference ACGTTCGA: CpGs at 0-based 1 and 5.
    let contig = b"ACGTTCGA";
    let r = record(1, &[(Kind::Match, 8)], b"ATGTTCGA", &[30; 8], false);
    let o = observed(&r, contig, DEFAULT);
    assert_eq!(o.outcome, Outcome::Counted);
    assert_eq!(o.cpgs, [(1, true), (5, false)]);
    assert_eq!((o.non_cpg, o.non_cpg_converted), (0, 0));
}

#[test]
fn reverse_reads_report_the_forward_c_of_the_site() {
    // Reverse strand: the G of a forward CG is the C read on that strand;
    // a read A there is a conversion.
    let contig = b"ACGTTCGA";
    let r = record(1, &[(Kind::Match, 8)], b"ACATTCGA", &[30; 8], true);
    let o = observed(&r, contig, DEFAULT);
    assert_eq!(o.cpgs, [(5, false), (1, true)]);
}

#[test]
fn non_cpg_c_uses_the_whole_read_qualities_in_file_order() {
    // Reference CCAT reversed is ATGG: no C on the reverse strand.
    // Forward: C at 0 (next C) and C at 1 (next A) are non-CpG.
    let contig = b"CCATGG";
    let r = record(1, &[(Kind::Match, 6)], b"TCATGG", &[30; 6], false);
    let o = observed(&r, contig, DEFAULT);
    assert_eq!((o.non_cpg, o.non_cpg_converted), (2, 1));
    assert!(o.cpgs.is_empty());
    // Picard tests the non-CpG qualities on the record's own array at the
    // block index: with a soft clip of 2, index 0 of the block reads the
    // clipped base's quality.
    let r = record(
        1,
        &[(Kind::SoftClip, 2), (Kind::Match, 6)],
        b"GGTCATGG",
        &[5, 30, 30, 30, 30, 30, 30, 30],
        false,
    );
    let o = observed(&r, contig, DEFAULT);
    assert_eq!((o.non_cpg, o.non_cpg_converted), (1, 0));
}

#[test]
fn short_and_mismatched_reads_are_ignored() {
    let contig = b"ACGTACGTACGT";
    let r = record(1, &[(Kind::Match, 4)], b"ACGT", &[30; 4], false);
    assert_eq!(observed(&r, contig, DEFAULT).outcome, Outcome::Short);
    // Two mismatches in 12 bases: round(1.2) = 1 allowed.
    let r = record(1, &[(Kind::Match, 12)], b"AAGTACGTAAGT", &[30; 12], false);
    assert_eq!(observed(&r, contig, DEFAULT).outcome, Outcome::Mismatched);
    // C to T is no mismatch on the forward strand.
    let r = record(1, &[(Kind::Match, 12)], b"ATGTATGTATGT", &[30; 12], false);
    assert_eq!(observed(&r, contig, DEFAULT).outcome, Outcome::Counted);
}

#[test]
fn low_qualities_drop_the_site_and_the_last_base_starts_nothing() {
    let contig = b"AACG";
    let r = record(1, &[(Kind::Match, 4)], b"AACG", &[30, 30, 19, 30], false);
    assert!(observed(&r, contig, DEFAULT).cpgs.is_empty());
    let r = record(1, &[(Kind::Match, 4)], b"AACG", &[30, 30, 30, 9], false);
    assert!(observed(&r, contig, DEFAULT).cpgs.is_empty());
    let r = record(1, &[(Kind::Match, 3)], b"AAC", &[30; 3], false);
    let o = observed(
        &r,
        b"AACG",
        Thresholds {
            minimum_read_length: 1,
            ..DEFAULT
        },
    );
    assert!(o.cpgs.is_empty() && o.non_cpg == 0);
}

#[test]
fn lower_case_reference_counts_and_n_does_not() {
    let r = record(1, &[(Kind::Match, 5)], b"ACGTA", &[30; 5], false);
    assert_eq!(observed(&r, b"acgta", DEFAULT).cpgs, [(1, false)]);
    let r = record(1, &[(Kind::Match, 5)], b"ANGTA", &[30; 5], false);
    let o = observed(&r, b"ANGTA", DEFAULT);
    assert!(o.cpgs.is_empty() && o.non_cpg == 0);
}

#[test]
fn a_read_past_the_contig_end_fails() {
    let r = record(3, &[(Kind::Match, 5)], b"ACGTA", &[30; 5], false);
    let mut o = Observation::default();
    let error =
        observe(&r, b"ACGTA", DEFAULT, &mut Scratch::default(), &mut o).expect_err("past the end");
    assert!(error.to_string().contains("past the end"), "{error}");
}

#[test]
fn the_prefix_gets_one_dot() {
    assert_eq!(
        file_path(Path::new("out/s"), SUMMARY_EXTENSION),
        Path::new("out/s.rrbs_summary_metrics")
    );
    assert_eq!(
        file_path(Path::new("out/s."), DETAIL_EXTENSION),
        Path::new("out/s.rrbs_detail_metrics")
    );
}

#[test]
fn the_summary_takes_picard_histogram_statistics() {
    let site = |position, total, converted| Site {
        position,
        total,
        converted,
    };
    let contigs = [
        (0, Sites::pack(&mut [site(1, 1, 1), site(5, 4, 2)])),
        (1, Sites::pack(&mut [site(9, 2, 0), site(12, 3, 3)])),
    ];
    let row = summary_row(Counters::default(), &contigs, Labels::default());
    assert_eq!((row.cpg_bases_seen, row.cpg_bases_converted), (10, 6));
    assert_eq!(row.mean_cpg_coverage, 2.5);
    // Median of 1, 2, 3, 4 is 2.5, truncated.
    assert_eq!(row.median_cpg_coverage, 2);
    // Read counts 1, 2, 2, 3 and 7: the median is the third.
    let odd = [(
        0,
        Sites::pack(&mut [
            site(1, 7, 0),
            site(2, 2, 0),
            site(3, 1, 0),
            site(4, 3, 0),
            site(5, 2, 0),
        ]),
    )];
    let row = summary_row(Counters::default(), &odd, Labels::default());
    assert_eq!((row.median_cpg_coverage, row.mean_cpg_coverage), (2, 3.0));
    let empty = summary_row(Counters::default(), &[], Labels::default());
    assert!(empty.mean_cpg_coverage.is_nan());
    assert_eq!(
        (empty.median_cpg_coverage, empty.pct_cpg_bases_converted),
        (0, 0.0)
    );
}

/// An observation with the CpG sites `cpgs`.
fn with_cpgs(cpgs: &[(u32, bool)]) -> Observation {
    Observation {
        cpgs: cpgs.to_vec(),
        ..Observation::default()
    }
}

#[test]
fn a_unit_closes_each_contig_it_leaves() {
    let mut unit = Unit::default();
    unit.add(2, &with_cpgs(&[(7, true), (3, false)]));
    unit.add(2, &with_cpgs(&[(7, false)]));
    unit.add(0, &with_cpgs(&[(4, true)]));
    assert_eq!(unit.closed.len(), 1);
    assert_eq!(unit.open.len(), 1);
    let (counters, contigs) = unit.into_contigs();
    assert_eq!(counters.reads_aligned, 3);
    let site = |position, total, converted| Site {
        position,
        total,
        converted,
    };
    assert_eq!(
        contigs,
        [
            (2, Sites::pack(&mut [site(3, 1, 0), site(7, 2, 1)])),
            (0, Sites::pack(&mut [site(4, 1, 1)]))
        ]
    );
}

#[test]
fn detail_rows_are_made_as_they_are_written() {
    let labels = Labels {
        sample: Some("s".to_owned()),
        library: None,
        read_group: None,
    };
    let site = |position, total, converted| Site {
        position,
        total,
        converted,
    };
    let details = Details {
        names: vec!["b".to_owned(), "a".to_owned()],
        units: vec![(
            labels,
            vec![
                (1, Sites::pack(&mut [site(4, 3, 0)])),
                (0, Sites::pack(&mut [site(2, 4, 1)])),
            ],
        )],
        rows: 2,
    };
    let mut rows = Vec::new();
    details
        .each(&mut |row| {
            rows.push(row.to_vec());
            Ok(())
        })
        .expect("rows");
    let expected = |name: &str, position, total_sites, converted_sites| {
        DetailRow {
            sequence_name: name.to_owned(),
            position,
            total_sites,
            converted_sites,
            labels: Labels {
                sample: Some("s".to_owned()),
                ..Labels::default()
            },
        }
        .values()
    };
    assert_eq!(rows, [expected("a", 4, 3, 0), expected("b", 2, 4, 1)]);
}

#[test]
fn names_sort_as_java_strings() {
    assert_eq!(java_compare("chr10", "chrA"), std::cmp::Ordering::Less);
    assert_eq!(java_compare("B", "a"), std::cmp::Ordering::Less);
}
