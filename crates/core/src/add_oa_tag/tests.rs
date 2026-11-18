//! The `OA` entry of [`super::AddOa`] on synthetic records, with each value
//! Picard 3.5.0 printed for the same SAM lines.

use std::path::PathBuf;

use super::*;
use crate::bam_raw::{Tags, build};
use crate::record_edit::Mode;

fn header() -> sam::Header {
    "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:500\n"
        .parse()
        .expect("header")
}

fn edit(list: Option<Intervals>) -> AddOa {
    let header = header();
    let setup = Setup {
        header: &header,
        reference: None,
        cram: false,
        mode: Mode::Rewrite,
    };
    AddOa::new(&setup, list)
}

fn tags(fields: Vec<([u8; 2], TagValue)>) -> Vec<u8> {
    let mut out = Vec::new();
    Tags::from_unsorted(fields).encode(&mut out).expect("tags");
    out
}

/// The `OA` value after the edit, `None` when the record was left as it was.
fn oa(edit: &mut AddOa, mut record: RawRecord) -> Option<String> {
    let before = record.clone();
    match edit.record(&mut record).expect("edit") {
        Outcome::Same => {
            assert!(record == before);
            None
        }
        _ => match record.tag(OA).expect("tags") {
            Some(TagValue::String(text)) => Some(String::from_utf8(text).expect("ASCII")),
            other => panic!("OA is {other:?}"),
        },
    }
}

const M: u8 = 0;
const I: u8 = 1;
const S: u8 = 4;

#[test]
fn each_entry_is_the_one_picard_writes() {
    let mut edit = edit(None);
    let nm = |v| tags(vec![(NM, v)]);
    let cases: Vec<(RawRecord, &str)> = vec![
        (
            build("nm", 0, 0, 9, &[(10, M)], &[30; 10], &nm(TagValue::Int(2))),
            "c1,10,+,10M,60,2;",
        ),
        (
            build("nonm", 16, 0, 19, &[(4, M), (2, I), (4, M)], &[30; 10], b""),
            "c1,20,-,4M2I4M,60,;",
        ),
        (
            build(
                "oa2",
                256,
                0,
                39,
                &[(5, S), (5, M)],
                &[30; 10],
                &tags(vec![(
                    OA,
                    TagValue::String(b"c2,5,+,10M,60,1;c1,7,-,3S7M,10,;".to_vec()),
                )]),
            ),
            "c2,5,+,10M,60,1;c1,7,-,3S7M,10,;c1,40,+,5S5M,60,;",
        ),
        (
            build(
                "nmf",
                2048,
                0,
                49,
                &[(10, M)],
                &[30; 10],
                &nm(TagValue::Float(1.5)),
            ),
            "c1,50,+,10M,60,1.5;",
        ),
        (
            build(
                "nmz",
                0,
                0,
                59,
                &[(10, M)],
                &[30; 10],
                &nm(TagValue::String(b"abc".to_vec())),
            ),
            "c1,60,+,10M,60,abc;",
        ),
        (
            build(
                "nma",
                0,
                0,
                69,
                &[(10, M)],
                &[30; 10],
                &nm(TagValue::Char(b'x')),
            ),
            "c1,70,+,10M,60,x;",
        ),
        (
            build(
                "oai",
                0,
                0,
                79,
                &[(10, M)],
                &[30; 10],
                &tags(vec![(OA, TagValue::Int(5))]),
            ),
            "5c1,80,+,10M,60,;",
        ),
        (
            build(
                "big",
                0,
                1,
                0,
                &[(10, M)],
                &[30; 10],
                &nm(TagValue::Int(3_000_000_000)),
            ),
            "c2,1,+,10M,60,3000000000;",
        ),
        (
            build("cstar", 0, 0, 94, &[], &[30; 10], b""),
            "c1,95,+,*,60,;",
        ),
        (
            build("placed", 20, 0, 99, &[], &[30; 10], &nm(TagValue::Int(1))),
            "*,0,-,*,255,;",
        ),
        (
            build(
                "unmapped",
                4,
                -1,
                -1,
                &[],
                &[30; 10],
                &tags(vec![(OA, TagValue::String(b"c1,1,+,10M,3,;".to_vec()))]),
            ),
            "c1,1,+,10M,3,;*,0,+,*,255,;",
        ),
    ];
    for (record, expected) in cases {
        assert_eq!(oa(&mut edit, record).as_deref(), Some(expected));
    }
}

#[test]
fn the_tag_goes_to_its_place_in_htsjdk_order() {
    let mut edit = edit(None);
    let fields = tags(vec![
        (*b"XA", TagValue::Int(1)),
        (*b"AA", TagValue::String(b"a".to_vec())),
        (NM, TagValue::Int(300)),
        (*b"zz", TagValue::Int(7)),
    ]);
    let mut record = build("order", 0, 0, 89, &[(10, M)], &[30; 10], &fields);
    edit.record(&mut record).expect("edit");
    let names: Vec<[u8; 2]> = record
        .tags()
        .expect("tags")
        .iter()
        .map(|(name, _)| *name)
        .collect();
    assert_eq!(names, [*b"AA", OA, *b"XA", NM, *b"zz"]);
}

#[test]
fn an_array_value_and_a_comma_fail() {
    let mut edit = edit(None);
    let array = TagValue::Array(b'c', 2, vec![1, 2]);
    for name in [NM, OA] {
        let mut record = build(
            "arr",
            0,
            0,
            4,
            &[(2, M)],
            &[30; 2],
            &tags(vec![(name, array.clone())]),
        );
        let error = edit.record(&mut record).expect_err("array").to_string();
        assert!(error.contains("Java object reference"), "{error}");
    }
    let header: sam::Header = "@SQ\tSN:c,1\tLN:10\n".parse().expect("header");
    let setup = Setup {
        header: &header,
        reference: None,
        cram: false,
        mode: Mode::Rewrite,
    };
    let mut comma = AddOa::new(&setup, None);
    let mut record = build("r3", 0, 0, 4, &[(2, M)], &[30; 2], b"");
    let error = comma.record(&mut record).expect_err("comma").to_string();
    assert_eq!(
        error,
        "Reference name for record r3 contains a comma character."
    );
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-oa-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("scratch");
    dir
}

#[test]
fn intervals_select_the_records_over_them() {
    let dir = scratch("intervals");
    let path = dir.join("iv.interval_list");
    std::fs::write(
        &path,
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c3\tLN:10\n@SQ\tSN:c2\tLN:500\n\
         c1\t15\t15\t+\tt1\nc1\t96\t98\t+\tt2\nc2\t1\t1\t+\tt3\nc3\t1\t5\t+\tt4\n\
         c1\t200\t199\t+\tempty\n",
    )
    .expect("list");
    let mut edit = edit(Some(Intervals::read(&path).expect("read")));
    let cases = [
        (build("over", 0, 0, 9, &[(10, M)], &[30; 10], b""), true),
        (build("next", 0, 0, 15, &[(10, M)], &[30; 10], b""), false),
        (build("ends_at", 0, 0, 5, &[(10, M)], &[30; 10], b""), true),
        (
            build("ins_only", 0, 0, 14, &[(10, I)], &[30; 10], b""),
            false,
        ),
        (build("cstar", 0, 0, 14, &[], &[30; 10], b""), false),
        (build("c2", 0, 1, 0, &[(10, M)], &[30; 10], b""), true),
        (build("placed", 4, 0, 14, &[], &[30; 10], b""), false),
        (build("unmapped", 4, -1, -1, &[], &[30; 10], b""), false),
        (build("empty", 0, 0, 195, &[(10, M)], &[30; 10], b""), false),
    ];
    for (record, tagged) in cases {
        let name = String::from_utf8_lossy(record.name()).into_owned();
        assert_eq!(oa(&mut edit, record).is_some(), tagged, "{name}");
    }
    let bed = dir.join("iv.bed");
    std::fs::write(&bed, "c1\t14\t15\n").expect("bed");
    let error = Intervals::read(&bed).err().expect("a BED").to_string();
    assert!(
        error.contains("Interval list file must contain header."),
        "{error}"
    );
    std::fs::remove_dir_all(&dir).expect("cleanup");
}
