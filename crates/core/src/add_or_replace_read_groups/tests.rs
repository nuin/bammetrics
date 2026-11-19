use std::path::{Path, PathBuf};

use super::*;
use crate::bam_raw::build;
use crate::pass;
use crate::sam_text::Stringency;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data")
        .join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-read-groups-{}-{name}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).expect("scratch");
    dir
}

fn args(extra: &[&str]) -> Vec<String> {
    [
        "I=x.bam",
        "O=y.bam",
        "RGLB=l",
        "RGPL=ILLUMINA",
        "RGPU=u",
        "RGSM=s",
    ]
    .iter()
    .chain(extra)
    .map(|a| (*a).to_owned())
    .collect()
}

fn options(extra: &[&str]) -> Result<ReadGroupOptions> {
    ReadGroupOptions::from_invocation(&TOOL.parse(&args(extra)).expect("parse"))
}

fn usage(extra: &[&str]) -> String {
    match options(extra) {
        Err(Error::Usage(text)) => text,
        other => panic!("{extra:?}: {other:?}"),
    }
}

#[test]
fn the_read_group_holds_htsjdks_fields_in_order() {
    let opts = options(&["RGFO=f", "RGPI=+5", "RGCN=c", "RGDT=2020-01-02"]).expect("options");
    let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:9\n@RG\tID:a\n@RG\tID:b\n"
        .parse()
        .expect("header");
    let edit = ReadGroups { options: opts };
    let out = crate::bam_raw::header_text(&edit.header(&header).expect("edit")).expect("text");
    assert_eq!(
        out,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:9\n@RG\tID:1\tLB:l\tPL:ILLUMINA\tSM:s\tPU:u\t\
         CN:c\tDT:2020-01-02T00:00:00+0000\tPI:5\tFO:f\n"
    );
}

#[test]
fn a_sort_order_sets_so() {
    let edit = ReadGroups {
        options: options(&["SO=queryname"]).expect("options"),
    };
    let header = sam::Header::default();
    let out = crate::bam_raw::header_text(&edit.header(&header).expect("edit")).expect("text");
    assert!(out.starts_with("@HD\tVN:1.6\tSO:queryname\n"), "{out}");
}

#[test]
fn every_bad_value_is_named_in_one_error() {
    let text = usage(&["RGID=a\tb", "RGCN=", "RGKS=é"]);
    assert_eq!(text.lines().count(), 3, "{text}");
    assert!(text.contains("for RGID, 'a\tb'"), "{text}");
    assert!(text.contains("for RGCN, ''"), "{text}");
    assert!(text.contains("for RGKS, 'é'"), "{text}");
}

#[test]
fn required_and_typed_options_fail_as_picard_fails() {
    let missing = ReadGroupOptions::from_invocation(
        &TOOL
            .parse(&["I=x".to_owned(), "O=y".to_owned()])
            .expect("parse"),
    );
    assert!(matches!(missing, Err(Error::Usage(t)) if t == "Argument 'RGLB' is required"));
    let null = ["I=x", "O=y", "RGLB=l", "RGPL=P", "RGPU=u", "RGSM=null"].map(str::to_owned);
    let null = ReadGroupOptions::from_invocation(&TOOL.parse(&null).expect("parse"));
    assert!(
        matches!(null, Err(Error::Usage(t)) if t == "non-null value must be provided for 'RGSM'")
    );
    assert!(usage(&["RGPI=abc"]).contains("Problem constructing Integer"));
    assert!(usage(&["RGDT=garbage"]).contains("Problem constructing Iso8601Date"));
    assert_eq!(
        usage(&["SO=bogus"]),
        "'bogus' is not a valid value for SortOrder."
    );
}

#[test]
fn rgid_null_removes_the_rg_tag_and_names_the_group_null() {
    let opts = options(&["RGID=null"]).expect("options");
    assert_eq!(opts.id, None);
    let mut edit = ReadGroups { options: opts };
    let mut record = build("r", 0, 0, 9, &[(5, 0)], b"IIIII", b"RGZx\0");
    edit.record(&mut record).expect("edit");
    assert_eq!(record.tag(*b"RG").expect("tag"), None);
    let header = edit.header(&sam::Header::default()).expect("header");
    assert!(header.read_groups().contains_key(&b"null"[..]));
}

#[test]
fn each_record_gets_rg_and_htsjdks_qualities() {
    let mut edit = ReadGroups {
        options: options(&["RGID=g"]).expect("options"),
    };
    let mut record = build("r", 0, 0, 9, &[(3, 0)], &[0xFF, 30, 30], b"XIC\x05RGZold\0");
    assert_eq!(edit.record(&mut record).expect("edit"), Outcome::Changed);
    assert_eq!(
        record.tag(*b"RG").expect("tag"),
        Some(TagValue::String(b"g".to_vec()))
    );
    assert_eq!(record.qualities(), &[0xFF, 0xFF, 0xFF]);
}

#[test]
fn a_step_sorts_for_a_sort_order_that_names_a_sort() {
    let argv = |so: &str| vec!["I=x".to_owned(), format!("SORT_ORDER={so}")];
    assert!(sorts(&argv("queryname")));
    assert!(sorts(&argv("duplicate")));
    assert!(!sorts(&argv("unsorted")));
    assert!(!sorts(&argv("unknown")));
    assert!(!sorts(&["I=x".to_owned()]));
}

#[test]
fn the_writer_sorts_only_an_order_that_differs_from_the_inputs() {
    let writer = |order: Option<SortOrder>, input: SortOrder| {
        WriterOrder::Header {
            order,
            presorted_if_same: true,
        }
        .against(input)
    };
    assert_eq!(writer(None, SortOrder::Coordinate), WriterOrder::Presorted);
    assert_eq!(
        writer(Some(SortOrder::Coordinate), SortOrder::Coordinate),
        WriterOrder::Presorted
    );
    assert_eq!(
        writer(Some(SortOrder::Queryname), SortOrder::Coordinate),
        WriterOrder::Sort(Order::Queryname)
    );
    assert_eq!(
        writer(Some(SortOrder::Unsorted), SortOrder::Coordinate),
        WriterOrder::Presorted
    );
    assert_eq!(
        writer(Some(SortOrder::Unknown), SortOrder::Coordinate),
        WriterOrder::Unknown
    );
}

#[test]
fn a_run_sorts_into_queryname_order_and_sets_every_rg() {
    let dir = scratch("run");
    let out = dir.join("out.bam");
    let argv: Vec<String> = [
        format!("I={}", data("markdup.bam").display()),
        format!("O={}", out.display()),
        "RGLB=l".to_owned(),
        "RGPL=P".to_owned(),
        "RGPU=u".to_owned(),
        "RGSM=s".to_owned(),
        "SORT_ORDER=queryname".to_owned(),
        "MAX_RECORDS_IN_RAM=4".to_owned(),
        format!("TMP_DIR={}", dir.display()),
    ]
    .to_vec();
    TOOL.run_cli(&argv).expect("run");
    let names = pass::with_record_stream_at(&out, None, 1, Stringency::Silent, |header, r, _| {
        assert_eq!(SortOrder::of(header), SortOrder::Queryname);
        let mut names = Vec::new();
        let mut record = RawRecord::default();
        while record.read_from(r).map_err(Error::Record)? {
            assert_eq!(
                record.tag(*b"RG").map_err(Error::Record)?,
                Some(TagValue::String(b"1".to_vec()))
            );
            names.push(record.name().to_vec());
        }
        Ok(names)
    })
    .expect("read back");
    assert_eq!(names.len(), 36);
    assert!(names.windows(2).all(|w| w[0] <= w[1]));
    let left: Vec<_> = std::fs::read_dir(&dir)
        .expect("dir")
        .map(|e| e.expect("entry").file_name())
        .collect();
    assert_eq!(left, ["out.bam"], "spilled runs removed");
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn an_unknown_order_over_two_records_writes_nothing() {
    let dir = scratch("unknown");
    let out = dir.join("out.bam");
    let argv: Vec<String> = [
        format!("I={}", data("markdup.bam").display()),
        format!("O={}", out.display()),
        "RGLB=l".to_owned(),
        "RGPL=P".to_owned(),
        "RGPU=u".to_owned(),
        "RGSM=s".to_owned(),
        "SORT_ORDER=unknown".to_owned(),
    ]
    .to_vec();
    let error = TOOL.run_cli(&argv).expect_err("unknown");
    assert!(error.to_string().contains("ClassCastException"), "{error}");
    assert!(!out.exists());
    std::fs::remove_dir_all(&dir).expect("scratch");
}

/// The invocation of a chain step: no OUTPUT.
fn step(extra: &[&str]) -> Invocation {
    let argv: Vec<String> = args(extra)
        .into_iter()
        .filter(|arg| !arg.starts_with("O="))
        .collect();
    TOOL.parse(&argv).expect("parse")
}

/// The step resolves its writer against the stream's header: the stream's
/// own order streams the records with the presorted writer's check, so a
/// stream sorted by samtools keeps its order of tied records; another
/// order sorts.
#[test]
fn a_step_streams_a_stream_already_in_its_order() {
    let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:900\n"
        .parse()
        .expect("header");
    let frame = Frame::new(header, "in.bam".to_owned());
    let temp = TempDir::new(None, "read-groups-test");
    let inv = |extra: &[&str]| step(extra);
    let given: [&[&str]; 3] = [&[], &["SORT_ORDER=coordinate"], &["SORT_ORDER=unsorted"]];
    for extra in given {
        let built = stage(&inv(extra), &frame, &temp).expect("stage");
        assert!(matches!(built, Stage::Transform(_)), "{extra:?}");
    }
    let built = stage(&inv(&["SORT_ORDER=queryname"]), &frame, &temp).expect("stage");
    assert!(matches!(built, Stage::Barrier(_)));
    // Two records at one position keep their order; a record before the
    // last one fails as htsjdk's assertPresorted fails.
    let Ok(Stage::Transform(mut transform)) = stage(&inv(&[]), &frame, &temp) else {
        panic!("a transform");
    };
    let mut batch = vec![
        build("b", 0, 0, 99, &[(5, 0)], b"IIIII", b""),
        build("a", 0, 0, 99, &[(5, 0)], b"IIIII", b""),
    ];
    transform.apply(&mut batch).expect("tied records");
    assert_eq!(batch[0].name(), b"b");
    let mut late = vec![build("c", 0, 0, 9, &[(5, 0)], b"IIIII", b"")];
    let message = transform
        .apply(&mut late)
        .expect_err("out of order")
        .to_string();
    assert!(
        message.contains("Alignments added out of order"),
        "{message}"
    );
}

/// htsjdk's `setHeaderStrict`: a reference or mate reference index outside
/// the dictionary fails the step as it fails the tool.
#[test]
fn a_reference_index_outside_the_dictionary_fails() {
    let header: sam::Header = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c\tLN:900\n"
        .parse()
        .expect("header");
    let frame = Frame::new(header, "in.bam".to_owned());
    let temp = TempDir::new(None, "read-groups-test");
    let inv = step(&[]);
    for (reference, mate, index) in [(77, -1, "77"), (0, 99, "99"), (-5, -1, "-5")] {
        let Ok(Stage::Transform(mut transform)) = stage(&inv, &frame, &temp) else {
            panic!("a transform");
        };
        let mut record = build("far", 0, 0, 9, &[(5, 0)], b"IIIII", b"");
        record.set_reference_ids(reference, mate);
        let message = transform
            .apply(&mut vec![record])
            .expect_err("outside")
            .to_string();
        assert!(
            message.contains(&format!("Reference name for '{index}' not found")),
            "{message}"
        );
    }
}
