//! The filters, the read list reader, the pair logic and the `.reads`
//! names.

use std::path::{Path, PathBuf};

use super::*;
use crate::bam_raw::build;

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-filter-sam-{}-{name}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).expect("scratch");
    dir
}

fn names(records: &VecDeque<RawRecord>) -> Vec<String> {
    records
        .iter()
        .map(|r| format!("{}:{}", String::from_utf8_lossy(r.name()), r.flags()))
        .collect()
}

fn run_filter(rule: Rule, by_pair: bool, records: Vec<RawRecord>) -> Result<Vec<String>> {
    let mut filtering = Filtering::new(rule, by_pair);
    for record in records {
        filtering.take(record)?;
    }
    filtering.end()?;
    Ok(names(&filtering.ready))
}

#[test]
fn a_pair_is_kept_or_dropped_together() {
    let records = vec![
        build("a", 0x41, 0, 1, &[(5, 0)], &[], b""),
        build("a", 0x85, 0, 1, &[], &[], b""),
        build("b", 0x41, 0, 1, &[(5, 0)], &[], b""),
        build("b", 0x81, 0, 1, &[(5, 0)], &[], b""),
        build("c", 0x0, 0, 1, &[(5, 0)], &[], b""),
        build("d", 0x4, -1, -1, &[], &[], b""),
    ];
    let kept = run_filter(Rule::Aligned(true), true, records.clone()).expect("filter");
    assert_eq!(kept, ["b:65", "b:129", "c:0"]);
    let kept = run_filter(Rule::Aligned(false), true, records).expect("filter");
    assert_eq!(kept, ["a:65", "a:133", "d:4"]);
}

#[test]
fn a_second_alone_is_kept_and_a_last_first_is_judged_alone() {
    let records = vec![
        build("a", 0x85, -1, -1, &[], &[], b""),
        build("b", 0x45, -1, -1, &[], &[], b""),
    ];
    let kept = run_filter(Rule::Aligned(true), true, records).expect("filter");
    assert_eq!(kept, ["a:133"]);
}

#[test]
fn a_record_of_both_ends_is_judged_again_after_its_first() {
    // `b` is flagged first and second: after `a` keeps it as its mate it
    // is a first of pair that waits for the record after it.
    let records = vec![
        build("a", 0x41, 0, 1, &[(5, 0)], &[], b""),
        build("a", 0xC1, 0, 1, &[(5, 0)], &[], b""),
        build("a", 0x85, -1, -1, &[], &[], b""),
    ];
    let kept = run_filter(Rule::Aligned(true), true, records).expect("filter");
    assert_eq!(kept, ["a:65"]);
}

#[test]
fn the_pair_filters_want_mates_in_queryname_order() {
    let records = vec![
        build("a", 0x41, 0, 1, &[(5, 0)], &[], b""),
        build("b", 0x81, 0, 1, &[(5, 0)], &[], b""),
    ];
    let error = run_filter(Rule::Aligned(true), true, records).expect_err("no mate");
    assert!(
        error.to_string().contains("Read names do not match"),
        "{error}"
    );
    let records = vec![
        build("b", 0x0, 0, 1, &[(5, 0)], &[], b""),
        build("a", 0x0, 0, 1, &[(5, 0)], &[], b""),
    ];
    let error = run_filter(Rule::Aligned(true), true, records).expect_err("order");
    assert!(error.to_string().contains("should come after"), "{error}");
}

#[test]
fn only_a_string_tag_matches_a_value() {
    let rule = || Rule::Tag {
        tag: java_tag("XT"),
        values: HashSet::from([b"1".to_vec()]),
        include: true,
        tags: HtsjdkTags::default(),
    };
    let records = vec![
        build("z", 0x4, -1, -1, &[], &[], b"XTZ1\0"),
        build("i", 0x4, -1, -1, &[], &[], b"XTc\x01"),
        build("a", 0x4, -1, -1, &[], &[], b"XTA1"),
        build("none", 0x4, -1, -1, &[], &[], b""),
    ];
    assert_eq!(run_filter(rule(), false, records).expect("filter"), ["z:4"]);
    let error = run_filter(
        Rule::Tag {
            tag: java_tag("XTX"),
            values: HashSet::new(),
            include: false,
            tags: HtsjdkTags::default(),
        },
        false,
        vec![build("z", 0x4, -1, -1, &[], &[], b"")],
    )
    .expect_err("tag length");
    assert!(error.to_string().contains("length() == 2"), "{error}");
    assert!(matches!(java_tag("X\u{100}"), JavaTag::NoRecordHasIt));
}

#[test]
fn read_names_take_the_first_word_of_each_line() {
    assert_eq!(first_name(b"r1 extra"), Some(b"r1".to_vec()));
    assert_eq!(first_name(b"r1\tx"), Some(b"r1".to_vec()));
    assert_eq!(first_name(b" r2"), None);
    assert_eq!(first_name(b"  \t"), None);
    assert_eq!(first_name(b"#r5"), Some(b"#r5".to_vec()));
    assert_eq!(
        first_name("\u{e9}t\u{e9}".as_bytes()),
        Some(vec![0xE9, b't', 0xE9])
    );
    assert_eq!(first_name("\u{100}".as_bytes()), None);
    assert_eq!(first_name(b"\xFFbad"), None);
    let dir = scratch("names");
    let list = dir.join("names.txt");
    std::fs::write(&list, b"r1\r\nr2\rr3\n\n  \nr4 x\n").expect("write");
    let mut got = read_names(&list).expect("names").names();
    got.sort();
    assert_eq!(
        got,
        [
            b"r1".to_vec(),
            b"r2".to_vec(),
            b"r3".to_vec(),
            b"r4".to_vec()
        ]
    );
    let gz = dir.join("names.txt.gz");
    let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(&gz).expect("create"));
    writer.write_all(b"g1\ng2\n").expect("gz");
    writer.finish().expect("finish");
    let mut got = read_names(&gz).expect("gz names").names();
    got.sort();
    assert_eq!(got, [b"g1".to_vec(), b"g2".to_vec()]);
    let empty = dir.join("empty.txt");
    std::fs::write(&empty, b"").expect("write");
    assert!(read_names(&empty).is_err());
    std::fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn the_name_set_holds_each_name_once_up_to_its_limit() {
    let mut set = NameSet::new(3);
    for name in [&b"r1"[..], b"r2", b"r1", b"", b"r2"] {
        set.insert(name).expect("under the limit");
    }
    assert!(set.contains(b"r1") && set.contains(b"r2") && set.contains(b""));
    assert!(!set.contains(b"r") && !set.contains(b"r12"));
    let error = set.insert(b"r3").expect_err("a fourth name");
    assert!(error.contains("more than 3 distinct reads"), "{error}");
    set.insert(b"r1").expect("a name held again");
    let mut names = set.names();
    names.sort();
    assert_eq!(names, [b"".to_vec(), b"r1".to_vec(), b"r2".to_vec()]);
}

#[test]
fn a_reads_line_writes_each_byte_from_0x80_up_as_its_character_in_utf8() {
    let dir = scratch("reads-utf8");
    let path = dir.join("in.reads");
    let mut reads = ReadsText::create(&path, &sam::Header::default()).expect("create");
    reads
        .write(&build("\u{e9}1", 0x4, -1, -1, &[], &[], b"XTZ\xC3\xA9\0"))
        .expect("write");
    reads.commit().expect("commit");
    let text = std::fs::read(&path).expect("read");
    assert!(text.starts_with(b"\xC3\x83\xC2\xA91\t4\t"), "{text:?}");
    assert!(text.ends_with(b"XT:Z:\xC3\x83\xC2\xA9\n"), "{text:?}");
    std::fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn reads_files_sit_beside_the_output() {
    assert_eq!(
        reads_path(Path::new("out/x.bam"), Path::new("/data/in.sam")),
        PathBuf::from("out/in.reads")
    );
    assert_eq!(
        reads_path(Path::new("x.bam"), Path::new("x.bam")),
        PathBuf::from("x.reads")
    );
    assert_eq!(
        reads_path(Path::new("o/x.bam"), Path::new(".hidden")),
        PathBuf::from("o/.hidden.reads")
    );
}

#[test]
fn an_option_for_another_filter_is_refused() {
    let parse = |args: &[&str]| {
        let args: Vec<String> = ["INPUT=in.bam", "OUTPUT=out.bam"]
            .iter()
            .chain(args)
            .map(|s| (*s).to_owned())
            .collect();
        FilterOptions::from_invocation(&TOOL.parse(&args).expect("parse"))
    };
    let error = parse(&["FILTER=includeAligned", "TAG=RG"]).expect_err("tag");
    assert!(
        error.to_string().contains("TAG may only be specified"),
        "{error}"
    );
    let error = parse(&["FILTER=includeReadList"]).expect_err("list");
    assert!(
        error
            .to_string()
            .contains("READ_LIST_FILE must be specified"),
        "{error}"
    );
    let error = parse(&["FILTER=includePairedIntervals"]).expect_err("intervals");
    assert!(
        error.to_string().contains("INTERVAL_LIST must be"),
        "{error}"
    );
    let error = parse(&["FILTER=includeJavascript"]).expect_err("js");
    assert!(matches!(error, Error::Usage(_)), "{error}");
    let error = parse(&["FILTER=nope"]).expect_err("unknown");
    assert!(matches!(error, Error::Usage(_)), "{error}");
    let options = parse(&["FILTER=excludeTagValues", "TAG=RG"]).expect("no values");
    assert!(options.tag_values.is_empty());
}

#[test]
fn intervals_keep_both_ends_of_a_pair() {
    let dir = scratch("intervals");
    let list = dir.join("t.interval_list");
    std::fs::write(
        &list,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\nchr1\t101\t200\t+\tt\n",
    )
    .expect("write");
    let mut header = sam::Header::default();
    header.reference_sequences_mut().insert(
        "chr1".into(),
        noodles::sam::header::record::value::Map::<
            noodles::sam::header::record::value::map::ReferenceSequence,
        >::new(std::num::NonZeroUsize::new(100_000).expect("length")),
    );
    let intervals = Intervals::new(&list, &header).expect("list");
    let mut over = build("o", 0x63, 0, 150, &[(10, 0)], &[], b"MCZ10M\0");
    over.set_reference_ids(0, 0);
    over.set_mate_alignment_start(5001);
    assert!(!intervals.out(&over).expect("test"));
    let mut mate_over = build("m", 0xA3, 0, 5000, &[(10, 0)], &[], b"MCZ10M\0");
    mate_over.set_reference_ids(0, 0);
    mate_over.set_mate_alignment_start(195);
    assert!(!intervals.out(&mate_over).expect("test"));
    let mut far = build("f", 0x63, 0, 5000, &[(10, 0)], &[], b"MCZ10M\0");
    far.set_reference_ids(0, 0);
    far.set_mate_alignment_start(6001);
    assert!(intervals.out(&far).expect("test"));
    let mut no_mc = build("n", 0x63, 0, 5000, &[(10, 0)], &[], b"");
    no_mc.set_reference_ids(0, 0);
    no_mc.set_mate_alignment_start(6001);
    assert!(intervals.out(&no_mc).is_err());
    let single = build("s", 0x0, 0, 150, &[(10, 0)], &[], b"");
    assert!(intervals.out(&single).expect("test"));
    std::fs::remove_dir_all(&dir).expect("clean");
}
