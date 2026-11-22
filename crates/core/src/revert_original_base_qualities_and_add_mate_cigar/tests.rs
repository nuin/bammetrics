use std::path::{Path, PathBuf};

use super::*;
use crate::bam_raw::build;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data")
        .join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-revert-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("scratch");
    dir
}

/// A record at 1-based `start` with `cigar` (`(length, op)`) and MAPQ.
fn record(
    name: &str,
    flags: u16,
    reference: i32,
    start: i32,
    cigar: &[(u32, u8)],
    mapq: u8,
) -> RawRecord {
    let bases = cigar
        .iter()
        .filter(|(_, op)| matches!(op, 0 | 1 | 4 | 7 | 8))
        .map(|(n, _)| *n as usize)
        .sum::<usize>();
    let mut r = build(
        name,
        flags,
        reference,
        start - 1,
        cigar,
        &vec![30; bases],
        b"",
    );
    r.set_mapping_quality(mapq);
    r
}

fn tag(record: &RawRecord, name: &[u8; 2]) -> Option<TagValue> {
    record.tag(*name).expect("tag")
}

#[test]
fn a_mapped_pair_gets_each_others_fields() {
    let mut a = record("p", PAIRED | FIRST, 0, 100, &[(50, 0)], 60);
    let mut b = record(
        "p",
        PAIRED | SECOND | REVERSE,
        0,
        300,
        &[(10, 4), (40, 0)],
        20,
    );
    set_mate_info(&mut a, &mut b, true).expect("mates");
    assert_eq!((a.mate_reference_id(), a.mate_alignment_start()), (0, 300));
    assert_eq!(a.flags() & (MATE_REVERSE | MATE_UNMAPPED), MATE_REVERSE);
    assert_eq!(tag(&a, b"MQ"), Some(TagValue::Int(20)));
    assert_eq!(tag(&a, b"MC"), Some(TagValue::String(b"10S40M".to_vec())));
    assert_eq!(tag(&b, b"MC"), Some(TagValue::String(b"50M".to_vec())));
    // 5' ends 100 and 339: 339 - 100 + 1.
    assert_eq!((a.template_length(), b.template_length()), (240, -240));
}

#[test]
fn a_pair_on_two_contigs_has_tlen_0() {
    let mut a = record("p", PAIRED | FIRST, 0, 100, &[(50, 0)], 60);
    let mut b = record("p", PAIRED | SECOND, 1, 300, &[(50, 0)], 60);
    set_mate_info(&mut a, &mut b, true).expect("mates");
    assert_eq!((a.template_length(), b.template_length()), (0, 0));
    assert_eq!(a.mate_reference_id(), 1);
}

#[test]
fn an_unmapped_mate_is_placed_beside_the_mapped_one() {
    let mut a = record("p", PAIRED | FIRST | REVERSE, 0, 100, &[(50, 0)], 60);
    let mut b = record("p", PAIRED | SECOND | UNMAPPED, -1, 0, &[], 0);
    set_mate_info(&mut a, &mut b, true).expect("mates");
    assert_eq!((b.reference_id(), b.alignment_start()), (0, 100));
    assert_eq!((b.mate_reference_id(), b.mate_alignment_start()), (0, 100));
    assert_ne!(b.flags() & MATE_REVERSE, 0);
    assert_eq!(b.flags() & MATE_UNMAPPED, 0);
    assert_eq!(tag(&b, b"MQ"), Some(TagValue::Int(60)));
    assert_eq!(tag(&b, b"MC"), Some(TagValue::String(b"50M".to_vec())));
    assert_ne!(a.flags() & MATE_UNMAPPED, 0);
    assert_eq!((tag(&a, b"MQ"), tag(&a, b"MC")), (None, None));
    assert_eq!((a.template_length(), b.template_length()), (0, 0));
}

#[test]
fn an_unmapped_pair_is_unplaced() {
    let mut a = record("p", PAIRED | FIRST | UNMAPPED, 0, 100, &[], 0);
    let mut b = record("p", PAIRED | SECOND | UNMAPPED | REVERSE, 0, 100, &[], 0);
    set_mate_info(&mut a, &mut b, true).expect("mates");
    for r in [&a, &b] {
        assert_eq!((r.reference_id(), r.alignment_start()), (-1, 0));
        assert_eq!((r.mate_reference_id(), r.mate_alignment_start()), (-1, 0));
        assert_ne!(r.flags() & MATE_UNMAPPED, 0);
    }
    assert_ne!(a.flags() & MATE_REVERSE, 0);
    assert_eq!(b.flags() & MATE_REVERSE, 0);
}

#[test]
fn a_supplementary_record_takes_its_mates_fields() {
    let mut group = vec![
        record("p", PAIRED | FIRST, 0, 100, &[(50, 0)], 60),
        record("p", PAIRED | SECOND, 0, 400, &[(50, 0)], 30),
        record(
            "p",
            PAIRED | FIRST | SUPPLEMENTARY,
            1,
            900,
            &[(20, 4), (30, 0)],
            10,
        ),
        record("p", PAIRED | SUPPLEMENTARY, 1, 950, &[(30, 0), (20, 4)], 10),
        record("p", PAIRED | SECOND | SECONDARY, 1, 10, &[(50, 0)], 0),
    ];
    mate_group(&mut group).expect("group");
    assert_eq!(group[2].mate_alignment_start(), 400);
    assert_eq!(group[2].template_length(), -group[1].template_length());
    assert_eq!(tag(&group[2], b"MQ"), Some(TagValue::Int(30)));
    assert_eq!(
        group[3].mate_alignment_start(),
        100,
        "neither flag: the first read"
    );
    assert_eq!(
        tag(&group[4], b"MC"),
        None,
        "a secondary record is left alone"
    );
}

#[test]
fn a_name_with_two_first_reads_or_a_missing_mate_fails() {
    let first = record("p", PAIRED | FIRST, 0, 100, &[(50, 0)], 60);
    let second = record("p", PAIRED | SECOND, 0, 300, &[(50, 0)], 60);
    let two = mate_group(&mut [first.clone(), first.clone(), second.clone()]);
    assert!(format!("{:?}", two.err()).contains("first of the pair"));
    let alone = mate_group(&mut [first]);
    assert!(format!("{:?}", alone.err()).contains("Missing second read of pair"));
    let other = mate_group(&mut [second]);
    assert!(format!("{:?}", other.err()).contains("Missing first read of pair"));
    let unpaired = record("u", 0, 0, 10, &[(5, 0)], 60);
    mate_group(&mut [unpaired]).expect("unpaired");
}

#[test]
fn oq_is_read_as_phred_scores_or_fails() {
    let mut r = record("r", 0, 0, 10, &[(3, 0)], 60);
    let mut tags = r.tags().expect("tags");
    tags.set(*b"OQ", TagValue::String(b"!+I".to_vec()));
    r.set_tags(&tags).expect("set");
    assert_eq!(original_qualities(&r).expect("oq"), Some(vec![0, 10, 40]));
    tags.set(*b"OQ", TagValue::String(Vec::new()));
    r.set_tags(&tags).expect("set");
    assert_eq!(original_qualities(&r).expect("empty"), None);
    tags.set(*b"OQ", TagValue::String(b"I I".to_vec()));
    r.set_tags(&tags).expect("set");
    assert!(original_qualities(&r).is_err());
    tags.set(*b"OQ", TagValue::Int(5));
    r.set_tags(&tags).expect("set");
    assert!(original_qualities(&r).is_err());
}

#[test]
fn the_gate_decides_as_can_skip_sam_file() {
    let gate = |examine| Skip {
        examine,
        restore: true,
        seen: 0,
        decision: None,
    };
    let unpaired = record("u", 0, 0, 10, &[(5, 0)], 60);
    let mut with_mc = record("p", PAIRED | FIRST, 0, 10, &[(5, 0)], 60);
    let mut tags = with_mc.tags().expect("tags");
    tags.set(*b"MC", TagValue::String(b"5M".to_vec()));
    with_mc.set_tags(&tags).expect("set");
    let without_mc = record("p", PAIRED | FIRST, 0, 10, &[(5, 0)], 60);

    let mut skip = gate(10);
    assert_eq!(skip.look(&unpaired).expect("look"), None);
    assert_eq!(skip.look(&with_mc).expect("look"), Some(true));
    let mut run = gate(10);
    assert_eq!(run.look(&without_mc).expect("look"), Some(false));
    assert_eq!(run.decision, Some(Decision::FoundNoMc));
    let mut none = gate(1);
    assert_eq!(none.look(&unpaired).expect("look"), None);
    assert_eq!(none.look(&unpaired).expect("look"), Some(false));
    assert_eq!(none.decision, Some(Decision::NoEvidence));
    let mut end = gate(10);
    assert!(end.at_end());
    assert_eq!(
        Decision::NoEvidence.message(7),
        "Found no evidence of OQ or mate with no mate cigar in the first 7 records.  Will \
         continue..."
    );
}

#[test]
fn a_file_in_shape_writes_nothing() {
    let dir = scratch("skip");
    let out = dir.join("out.bam");
    let argv = [
        format!("I={}", data("oxog.bam").display()),
        format!("O={}", out.display()),
    ];
    TOOL.run_cli(&argv).expect("run");
    assert!(!out.exists());
    let argv = [
        argv[0].clone(),
        argv[1].clone(),
        "MAX_RECORDS_TO_EXAMINE=0".to_owned(),
    ];
    TOOL.run_cli(&argv).expect("run");
    let names: Vec<_> = std::fs::read_dir(&dir)
        .expect("dir")
        .map(|e| e.expect("entry").file_name().into_string().expect("name"))
        .collect();
    let mut names = names;
    names.sort();
    assert_eq!(names, ["out.bai", "out.bam", "out.bam.md5"]);
    std::fs::remove_dir_all(&dir).expect("scratch");
}
