//! The mate pass, the grouping without a sort and the in-place target.

use std::path::{Path, PathBuf};

use super::*;
use crate::bam_raw::{TagValue, build};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data")
        .join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("bammetrics-fix-mate-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("scratch");
    dir
}

/// A pair `name` on reference 0: first at `a`, second at `b` (0-based),
/// ten bases each, the second reverse.
fn pair(name: &str, a: i32, b: i32) -> Vec<RawRecord> {
    vec![
        build(name, 0x41, 0, a, &[(10, 0)], &[30; 10], b""),
        build(name, 0x91, 0, b, &[(10, 0)], &[30; 10], b""),
    ]
}

#[test]
fn each_switch_picks_its_mate_pass() {
    for (mate_cigar, ignore) in [(true, true), (true, false), (false, true), (false, false)] {
        let mut group = pair("p", 100, 300);
        regroup(mate_cigar, ignore)(&mut group).expect("mates");
        assert_eq!(group[0].mate_alignment_start(), 301);
        assert_eq!(group[1].mate_alignment_start(), 101);
        assert_eq!(group[0].template_length(), 210);
        let mc = group[0].tag(*b"MC").expect("tags");
        assert_eq!(
            mc,
            mate_cigar.then(|| TagValue::String(b"10M".to_vec())),
            "{mate_cigar} {ignore}"
        );
        assert_eq!(group[0].tag(*b"MQ").expect("tags"), Some(TagValue::Int(60)));
    }
}

#[test]
fn without_add_mate_cigar_an_old_mc_is_removed() {
    let mut group = pair("p", 100, 300);
    let mut tags = group[0].tags().expect("tags");
    tags.set(*b"MC", TagValue::String(b"5M".to_vec()));
    group[0].set_tags(&tags).expect("tags");
    regroup(false, true)(&mut group).expect("mates");
    assert_eq!(group[0].tag(*b"MC").expect("tags"), None);
}

#[test]
fn a_missing_mate_is_handed_on_or_fails() {
    let mut alone = vec![build("o", 0x41, 0, 5, &[(10, 0)], &[30; 10], b"")];
    regroup(true, true)(&mut alone).expect("ignored");
    assert_eq!(alone[0].mate_alignment_start(), 0);
    assert_eq!(alone[0].tag(*b"MC").expect("tags"), None);
    let error = regroup(true, false)(&mut alone).expect_err("missing mate");
    assert!(
        error.to_string().contains("Missing second read of pair"),
        "{error}"
    );
}

#[test]
fn adjacent_groups_follow_the_records_as_they_come() {
    let mut adjacent = Adjacent {
        group: Vec::new(),
        regroup: regroup(true, true),
    };
    let mut out = Vec::new();
    let mut push = |record: &mut RawRecord| -> Result<()> {
        out.push(std::mem::take(record));
        Ok(())
    };
    // The mates of `a` are apart, so only `b` gets its mate fields.
    let mut records = pair("a", 10, 400);
    records.insert(1, pair("b", 20, 30).remove(0));
    records.insert(2, pair("b", 20, 30).remove(1));
    for record in records {
        adjacent.push(record, &mut push).expect("push");
    }
    adjacent.flush(&mut push).expect("flush");
    let names: Vec<&[u8]> = out.iter().map(RawRecord::name).collect();
    assert_eq!(names, [&b"a"[..], b"b", b"b", b"a"]);
    assert_eq!(out[0].mate_alignment_start(), 0);
    assert_eq!(out[1].mate_alignment_start(), 31);
}

#[test]
fn the_writer_is_presorted_for_queryname_alone() {
    assert_eq!(writer_order(SortOrder::Queryname), WriterOrder::Presorted);
    assert_eq!(writer_order(SortOrder::Unsorted), WriterOrder::Presorted);
    assert_eq!(
        writer_order(SortOrder::Coordinate),
        WriterOrder::Sort(crate::sam_order::Order::Coordinate)
    );
    assert_eq!(writer_order(SortOrder::Unknown), WriterOrder::Unknown);
}

#[test]
fn in_place_takes_a_regular_bam_and_resolves_links() {
    let dir = scratch("in-place");
    let bam = dir.join("x.bam");
    std::fs::copy(data("markdup.bam"), &bam).expect("copy");
    let link = dir.join("link.bam");
    let _ = std::fs::remove_file(&link);
    std::os::unix::fs::symlink(&bam, &link).expect("link");
    assert_eq!(
        in_place_target(&link).expect("target"),
        bam.canonicalize().expect("resolved")
    );
    let error = in_place_target(&data("markdup.sam")).expect_err("SAM text");
    assert!(matches!(error, Error::Usage(_)), "{error}");
    let error = in_place_target(&data("insert_size.cram")).expect_err("CRAM");
    assert!(matches!(error, Error::Usage(_)), "{error}");
    let named_sam = dir.join("bytes.sam");
    std::fs::copy(data("markdup.bam"), &named_sam).expect("copy");
    assert!(in_place_target(&named_sam).is_err());
    std::fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn a_fixed_file_keeps_its_records_and_rewrites_in_place() {
    let dir = scratch("run");
    let bam = dir.join("in.bam");
    std::fs::copy(data("revert.bam"), &bam).expect("copy");
    let inv = TOOL
        .parse(&[format!("INPUT={}", bam.display())])
        .expect("parse");
    let opts = FixOptions::from_invocation(&inv).expect("options");
    assert!(opts.in_place);
    let counters = run(&opts).expect("run");
    assert!(counters.records_out > 0);
    let header =
        crate::pass::read_header_at(&bam, crate::sam_text::Stringency::Silent).expect("header");
    assert_eq!(SortOrder::of(&header), SortOrder::Coordinate);
    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .expect("dir")
        .map(|e| e.expect("entry").file_name())
        .collect();
    assert_eq!(leftovers, [std::ffi::OsString::from("in.bam")]);
    std::fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn create_index_needs_coordinate_order() {
    let dir = scratch("index");
    let out = dir.join("out.bam");
    let inv = TOOL
        .parse(&[
            format!("INPUT={}", data("markdup.bam").display()),
            format!("OUTPUT={}", out.display()),
            "SORT_ORDER=queryname".to_owned(),
            "CREATE_INDEX=true".to_owned(),
        ])
        .expect("parse");
    let error = run(&FixOptions::from_invocation(&inv).expect("options")).expect_err("index");
    assert!(error.to_string().contains("Can't CREATE_INDEX"), "{error}");
    assert!(!out.exists());
    std::fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn two_inputs_need_an_output() {
    let inv = TOOL
        .parse(&[
            format!("INPUT={}", data("markdup.bam").display()),
            format!("INPUT={}", data("markdup.bam").display()),
        ])
        .expect("parse");
    assert!(matches!(
        FixOptions::from_invocation(&inv),
        Err(Error::Usage(_))
    ));
}
