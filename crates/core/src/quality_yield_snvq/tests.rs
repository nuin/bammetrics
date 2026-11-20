use super::*;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

/// An unmapped record with `bases`, Phred `qualities` and string tags,
/// round-tripped through BAM bytes.
fn record(flags: u16, bases: &[u8], qualities: &[u8], tags: &[(&[u8; 2], &str)]) -> bam::Record {
    let mut buf = RecordBuf::builder()
        .set_name("r1")
        .set_flags(Flags::from_bits_retain(flags | 0x4))
        .set_sequence(Sequence::from(bases.to_vec()))
        .set_quality_scores(QualityScores::from(qualities.to_vec()))
        .build();
    for (tag, text) in tags {
        buf.data_mut()
            .insert(Tag::new(tag[0], tag[1]), BufValue::from((*text).to_owned()));
    }
    let header = sam::Header::default();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_header(&header).expect("header");
    writer
        .write_alignment_record(&header, &buf)
        .expect("encode");
    let bytes = writer.into_inner().finish().expect("finish");
    let mut reader = bam::io::Reader::new(bytes.as_slice());
    reader.read_header().expect("header");
    let mut out = bam::Record::default();
    reader.read_record(&mut out).expect("decode");
    out
}

const SNVQ: [(&[u8; 2], &str); 4] = [
    (b"qa", "IIIII"),
    (b"qc", "55555"),
    (b"qg", "+++++"),
    (b"qt", "?????"),
];

fn run(opts: &Options, records: &[bam::Record]) -> Result<Report> {
    let mut acc = QualityYieldSnvq::new(opts);
    for r in records {
        acc.accept(r)?;
    }
    acc.finish()
}

#[test]
fn three_observations_per_called_base_and_four_per_n() {
    // Qualities 0, 10, 20, 30, 40 over A, C, G, T, N.
    let r = record(0, b"ACGTN", &[0, 10, 20, 30, 40], &SNVQ);
    let report = run(&Options::default(), &[r]).expect("runs");
    let row = report.row;
    assert_eq!((row.total_reads, row.bases.total), (1, 5));
    assert_eq!((row.bases.q20, row.bases.q30, row.bases.q40), (3, 2, 1));
    // qa 40, qc 20, qg 10, qt 30: A skips qa, C skips qc, G skips qg, T skips qt.
    assert_eq!(row.snvq.total, 16);
    // q40: 4 positions read qa (C, G, T, N); q20: 4 qa + 4 qc + 4 qt = 12.
    assert_eq!((row.snvq.q40, row.snvq.q30, row.snvq.q20), (4, 8, 12));
    assert_eq!(row.snvq.pf, 16);
}

#[test]
fn a_qc_failed_read_counts_in_the_totals_only() {
    let r = record(0x200, b"AAAA", &[40; 4], &SNVQ);
    let row = run(&Options::default(), &[r]).expect("runs").row;
    assert_eq!((row.total_reads, row.pf_reads), (1, 0));
    assert_eq!((row.bases.q40, row.bases.pf_q40), (4, 0));
    assert_eq!((row.snvq.total, row.snvq.pf, row.snvq.pf_q20), (12, 0, 0));
    let values = row.values();
    assert_eq!(values[2], Value::Int(4));
    // PCT_PF_Q20_BASES over no PF base is 0.
    assert_eq!(values[14], Value::Float(0.0));
}

#[test]
fn a_missing_tag_names_the_record_and_the_tag() {
    let r = record(0, b"AC", &[30, 30], &SNVQ[..3]);
    let error = run(&Options::default(), &[r]).expect_err("no qt");
    assert!(
        error.to_string().contains("record r1: no qt tag"),
        "{error}"
    );
    // A secondary record is skipped before its tags are read.
    let secondary = record(0x100, b"AC", &[30, 30], &[]);
    assert!(run(&Options::default(), &[secondary]).is_ok());
}

#[test]
fn a_short_tag_fails_only_where_a_base_reads_past_it() {
    // Every base is A, so qa is never read: a short qa passes.
    let tags = [
        (b"qa", "I"),
        (b"qc", "IIII"),
        (b"qg", "IIII"),
        (b"qt", "IIII"),
    ];
    let r = record(0, b"AAAA", &[30; 4], &tags);
    assert_eq!(
        run(&Options::default(), &[r]).expect("runs").row.snvq.total,
        12
    );
    let r = record(0, b"AACA", &[30; 4], &tags);
    let error = run(&Options::default(), &[r]).expect_err("qa read at index 2");
    assert!(
        error.to_string().contains("read index 2 of the qa tag"),
        "{error}"
    );
}

#[test]
fn qualities_must_match_the_bases() {
    let alternate = Options {
        alternate_quality_attribute: Some("XQ".to_owned()),
        ..Options::default()
    };
    let mut tags = SNVQ.to_vec();
    tags.push((b"XQ", "II"));
    let r = record(0, b"ACG", &[30; 3], &tags);
    let error = run(&alternate, &[r]).expect_err("2 against 3");
    assert!(
        error
            .to_string()
            .contains("quality string length does not match"),
        "{error}"
    );
    tags.pop();
    tags.push((b"XQ", "5I+"));
    let r = record(0, b"ACG", &[0; 3], &tags);
    let row = run(&alternate, &[r]).expect("runs").row;
    assert_eq!((row.bases.q20, row.bases.q30, row.bases.q40), (2, 1, 1));
    let bad_name = Options {
        alternate_quality_attribute: Some("XQZ".to_owned()),
        ..Options::default()
    };
    let r = record(0, b"A", &[30], &SNVQ);
    let error = run(&bad_name, &[r]).expect_err("three letters");
    assert!(error.to_string().contains("length() == 2: XQZ"), "{error}");
}

#[test]
fn the_histogram_keys_qualities_and_read_indexes_on_one_column() {
    let opts = Options {
        include_bq_histogram: true,
        ..Options::default()
    };
    let a = record(0, b"AC", &[0, 40], &SNVQ);
    let b = record(0, b"A", &[20], &SNVQ);
    let report = run(&opts, &[a, b]).expect("runs");
    let section = report.histogram().expect("asked for");
    let labels: Vec<&str> = section.series.iter().map(|s| s.label.as_str()).collect();
    assert_eq!(labels, SERIES);
    // Index 0 holds Q0 and Q20: mean 10; index 1 holds Q40.
    assert_eq!(section.series[2].histogram.get(0), Some(10.0));
    assert_eq!(section.series[2].histogram.get(1), Some(40.0));
    assert_eq!(section.series[0].histogram.get(40), Some(1.0));
    // Six SNVQ observations at index 0 (qc 20, qg 10, qt 30 twice).
    assert_eq!(section.series[3].histogram.get(0), Some(20.0));
    assert!(
        run(&Options::default(), &[])
            .expect("runs")
            .histogram()
            .is_none()
    );
}

#[test]
fn java_bytes_widen_latin1_characters() {
    assert_eq!(java_bytes(b"I5").as_ref(), b"I5");
    assert_eq!(java_bytes(&[b'I', 0xE9]).as_ref(), &[b'I', 0xC3, 0xA9]);
    assert_eq!(phred(b'!'), Some(0));
    assert_eq!(phred(b'~'), Some(93));
    assert_eq!(phred(b' '), None);
    assert_eq!(phred(0xC3), None);
}
