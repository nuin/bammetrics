//! The walks, the tally, the header lines and the Picard argument form, on
//! SAM text written here. Picard's own results for these shapes are in
//! `tests/concordance/test_compare_s_a_ms.py`.

use std::path::{Path, PathBuf};

use super::*;

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-compare-sams-{}-{name}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).expect("scratch");
    dir
}

const SQ: &str = "@SQ\tSN:c1\tLN:5000\n@SQ\tSN:c2\tLN:3000\n";
const RG: &str = "@RG\tID:g1\tSM:s\tLB:lib\tPU:u\n";

fn header(order: &str) -> String {
    format!("@HD\tVN:1.6\tSO:{order}\n{SQ}{RG}")
}

/// One record of 10 bases; `pos` 0 is unplaced.
fn rec(name: &str, flag: u16, pos: i32, mapq: u8) -> String {
    let (contig, cigar) = match (pos, flag & 4) {
        (0, _) => ("*", "*"),
        (_, 4) => ("c1", "*"),
        _ => ("c1", "10M"),
    };
    format!(
        "{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:g1\n"
    )
}

fn file(dir: &Path, name: &str, header: &str, records: &[String]) -> String {
    let path = dir.join(name);
    std::fs::write(&path, format!("{header}{}", records.concat())).expect("write");
    path.to_string_lossy().into_owned()
}

/// The comparison of two files with `options`, and the lines it said.
fn run(left: &str, right: &str, options: &[&str]) -> (Result<Comparison>, Vec<String>) {
    let mut args = vec![format!("INPUT={left}"), format!("INPUT={right}")];
    args.push("VALIDATION_STRINGENCY=SILENT".to_owned());
    args.extend(options.iter().map(|o| (*o).to_owned()));
    let inv = TOOL.parse(&args).expect("parse");
    let mut lines = Vec::new();
    let result = compare(&inv, &mut |line| lines.push(line.to_owned()));
    (result, lines)
}

fn counts(left: &str, right: &str, options: &[&str]) -> Counts {
    run(left, right, options).0.expect("compare").row.counts
}

#[test]
fn mapping_quality_keys_order_as_their_texts() {
    let mut pairs: Vec<(u8, u8)> = (0..=255u8)
        .flat_map(|left| (0..=255u8).map(move |right| (left, right)))
        .collect();
    pairs.sort_by_key(|&(left, right)| mq_pair_key(left, right));
    let texts: Vec<String> = pairs
        .iter()
        .map(|&(left, right)| mq_pair_text(mq_pair_key(left, right)))
        .collect();
    let mut sorted = texts.clone();
    sorted.sort();
    assert_eq!(texts, sorted);
    assert_eq!(texts.len(), 65536);
    assert_eq!(mq_pair_text(mq_pair_key(60, 0)), "60,0");
    assert_eq!(mq_pair_text(mq_pair_key(255, 105)), "255,105");
    assert!(mq_pair_key(10, 3) < mq_pair_key(2, 2), "10,3 before 2,2");
}

#[test]
fn positional_arguments_become_inputs() {
    let args = |words: &[&str]| -> Vec<String> { words.iter().map(|w| (*w).to_owned()).collect() };
    let turned = |words: &[&str]| picard_arguments(&args(words)).expect("arguments");
    assert_eq!(
        turned(&["a.bam", "O=x.txt", "b.bam", "COMPARE_MQ=true"]),
        ["INPUT=a.bam", "O=x.txt", "INPUT=b.bam", "COMPARE_MQ=true"]
    );
    assert_eq!(
        turned(&[
            "--OUTPUT",
            "x.txt",
            "a.bam",
            "--LENIENT_DUP",
            "true",
            "b.bam",
            "--COMPARE_MQ"
        ]),
        [
            "--OUTPUT",
            "x.txt",
            "INPUT=a.bam",
            "--LENIENT_DUP",
            "true",
            "INPUT=b.bam",
            "--COMPARE_MQ",
            "true"
        ]
    );
    // A boolean takes the next word only when it is true or false; a
    // flag alone is said in full for the registry.
    assert_eq!(
        turned(&["--COMPARE_MQ", "--QUIET", "a.bam"]),
        ["--COMPARE_MQ", "true", "--QUIET", "true", "INPUT=a.bam"]
    );
    let flag_first = turned(&["--LENIENT_DUP", "a.bam", "b.bam", "--COMPARE_MQ", "FALSE"]);
    assert_eq!(
        flag_first,
        [
            "--LENIENT_DUP",
            "true",
            "INPUT=a.bam",
            "INPUT=b.bam",
            "--COMPARE_MQ",
            "FALSE"
        ]
    );
    let inv = TOOL.parse(&flag_first).expect("parse");
    let options = Options::from_invocation(&inv).expect("options");
    assert_eq!(inv.all("INPUT"), ["a.bam", "b.bam"]);
    assert!(options.lenient_dup && !options.compare_mq);
    for named in [&["I=a.bam"][..], &["--INPUT", "a.bam"], &["INPUT=a.bam"]] {
        let error = picard_arguments(&args(named)).expect_err("INPUT by name");
        assert!(error.to_string().contains("Unrecognized option"), "{error}");
    }
}

#[test]
fn two_inputs_are_needed() {
    for inputs in [1, 3] {
        let args: Vec<String> = (0..inputs).map(|i| format!("INPUT=x{i}.bam")).collect();
        let inv = TOOL.parse(&args).expect("parse");
        let error = Options::from_invocation(&inv).expect_err("two inputs");
        assert!(matches!(error, Error::Usage(_)), "{error}");
    }
    let args = ["INPUT=a", "INPUT=b", "LOW_MQ_THRESHOLD=null"].map(str::to_owned);
    let inv = TOOL.parse(&args).expect("parse");
    assert!(Options::from_invocation(&inv).is_err());
}

#[test]
fn a_coordinate_walk_pairs_records_by_place_then_by_key() {
    let dir = scratch("coordinate");
    let head = header("coordinate");
    let base = [
        rec("a", 0, 100, 60),
        rec("b", 0, 200, 60),
        rec("c", 0, 300, 60),
        rec("u", 4, 0, 0),
    ];
    let left = file(&dir, "left.sam", &head, &base);
    let same = counts(&left, &left, &[]);
    assert_eq!((same.mappings_match, same.unmapped_both), (3, 1));
    assert!(same.all_visited_alignments_equal());

    // b moved, c missing on the right, d only on the right, u placed.
    let other = [
        rec("a", 16, 100, 60),
        rec("b", 0, 250, 60),
        rec("d", 0, 260, 60),
        rec("u", 0, 400, 60),
    ];
    let right = file(&dir, "right.sam", &head, &other);
    let found = counts(&left, &right, &[]);
    assert_eq!(
        found,
        Counts {
            mappings_match: 1,
            mappings_differ: 1,
            unmapped_left: 1,
            missing_left: 1,
            missing_right: 1,
            ..Counts::default()
        }
    );
    let back = counts(&right, &left, &[]);
    assert_eq!(
        (back.unmapped_right, back.missing_left, back.missing_right),
        (1, 1, 1)
    );
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn a_key_given_twice_keeps_its_last_record() {
    let dir = scratch("twice");
    let head = header("coordinate");
    let plain = [
        rec("a", 0, 100, 60),
        rec("b", 0, 200, 60),
        rec("c", 0, 300, 60),
    ];
    let twice = [
        rec("a", 0, 100, 60),
        rec("a", 0, 100, 60),
        rec("b", 0, 200, 60),
    ];
    let plain = file(&dir, "plain.sam", &head, &plain);
    let twice = file(&dir, "twice.sam", &head, &twice);
    let left_twice = counts(&twice, &plain, &[]);
    assert_eq!(
        (
            left_twice.mappings_match,
            left_twice.missing_left,
            left_twice.missing_right
        ),
        (2, 1, 0)
    );
    let right_twice = counts(&plain, &twice, &[]);
    assert_eq!(
        (
            right_twice.mappings_match,
            right_twice.missing_left,
            right_twice.missing_right
        ),
        (2, 1, 1)
    );
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn secondary_and_supplementary_records_are_skipped() {
    let dir = scratch("skipped");
    let head = header("coordinate");
    let left = [
        rec("a", 0, 100, 60),
        rec("a", 256, 150, 60),
        rec("b", 0, 200, 60),
    ];
    let right = [
        rec("a", 0, 100, 60),
        rec("b", 0, 200, 60),
        rec("b", 2048, 250, 60),
    ];
    let left = file(&dir, "left.sam", &head, &left);
    let right = file(&dir, "right.sam", &head, &right);
    let found = run(&left, &right, &[]).0.expect("compare");
    assert_eq!(found.row.counts.mappings_match, 2);
    assert!(found.row.are_equal);
    assert_eq!(found.status(), 0);
    assert_eq!(found.verdict(), "SAM files match.");
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn a_queryname_walk_merges_and_fails_when_the_right_file_ends_first() {
    let dir = scratch("queryname");
    let head = header("queryname");
    let all = [
        rec("a", 0, 100, 60),
        rec("b", 0, 200, 60),
        rec("c", 0, 300, 60),
    ];
    let whole = file(&dir, "whole.sam", &head, &all);
    let no_first = file(&dir, "no_first.sam", &head, &all[1..]);
    let no_last = file(&dir, "no_last.sam", &head, &all[..2]);
    let lacking = counts(&whole, &no_first, &[]);
    assert_eq!((lacking.mappings_match, lacking.missing_right), (2, 1));
    let extra = counts(&no_last, &whole, &[]);
    assert_eq!((extra.mappings_match, extra.missing_left), (2, 1));
    let (failed, _) = run(&whole, &no_last, &[]);
    let error = failed.expect_err("Picard throws here");
    assert!(
        error.to_string().contains("Picard 3.5.0 fails here"),
        "{error}"
    );
    // First and second of a pair are two keys; an unpaired record sorts first.
    let pairs = [
        rec("p", 0, 100, 60),
        rec("p", 65, 100, 60),
        rec("p", 129, 300, 60),
    ];
    let paired = file(&dir, "paired.sam", &head, &pairs);
    assert_eq!(counts(&paired, &paired, &[]).mappings_match, 3);
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn an_unsorted_walk_holds_the_left_file() {
    let dir = scratch("unsorted");
    let head = header("unsorted");
    let all = [
        rec("a", 0, 100, 60),
        rec("b", 0, 200, 60),
        rec("c", 0, 300, 60),
    ];
    let mut reversed = all.clone();
    reversed.reverse();
    let left = file(&dir, "left.sam", &head, &all);
    let right = file(&dir, "right.sam", &head, &reversed);
    let found = run(&left, &right, &[]).0.expect("compare");
    assert_eq!(found.row.counts.mappings_match, 3);
    assert_eq!(found.peak_pending, 3);
    let fewer = file(&dir, "fewer.sam", &head, &all[..1]);
    assert_eq!(counts(&left, &fewer, &[]).missing_right, 2);
    assert_eq!(counts(&fewer, &left, &[]).missing_left, 2);
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn files_in_one_order_hold_no_waiting_record() {
    let dir = scratch("pending");
    let head = header("coordinate");
    let records: Vec<String> = (0..2000)
        .map(|i| rec(&format!("r{i:05}"), 0, 1 + i / 4, 60))
        .collect();
    let left = file(&dir, "left.sam", &head, &records);
    let same = run(&left, &left, &[]).0.expect("compare");
    assert_eq!(same.row.counts.mappings_match, 2000);
    assert_eq!(same.peak_pending, 0);
    // No read in common: every record of both files waits.
    let others: Vec<String> = (0..2000)
        .map(|i| rec(&format!("s{i:05}"), 0, 1 + i / 4, 60))
        .collect();
    let right = file(&dir, "right.sam", &head, &others);
    let disjoint = run(&left, &right, &[]).0.expect("compare");
    assert_eq!(disjoint.peak_pending, 4000);
    assert_eq!(disjoint.row.counts.missing_left, 2000);
    assert_eq!(disjoint.row.counts.missing_right, 2000);
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn leniency_tests_both_qualities_and_the_threshold_is_inclusive() {
    let dir = scratch("lenient");
    let head = header("coordinate");
    let qualities = [(3, 3), (3, 2), (4, 3), (255, 255), (255, 60), (0, 0)];
    let left: Vec<String> = qualities
        .iter()
        .enumerate()
        .map(|(i, q)| rec(&format!("r{i}"), 0, 100 + 100 * i as i32, q.0))
        .collect();
    let right: Vec<String> = qualities
        .iter()
        .enumerate()
        .map(|(i, q)| rec(&format!("r{i}"), 0, 105 + 100 * i as i32, q.1))
        .collect();
    let left = file(&dir, "left.sam", &head, &left);
    let right = file(&dir, "right.sam", &head, &right);
    let matched = |options: &[&str]| counts(&left, &right, options).mappings_match;
    assert_eq!(matched(&[]), 0);
    assert_eq!(matched(&["LENIENT_LOW_MQ_ALIGNMENT=true"]), 3);
    assert_eq!(
        matched(&["LENIENT_LOW_MQ_ALIGNMENT=true", "LOW_MQ_THRESHOLD=2"]),
        1
    );
    assert_eq!(
        matched(&["LENIENT_LOW_MQ_ALIGNMENT=true", "LOW_MQ_THRESHOLD=4"]),
        4
    );
    assert_eq!(matched(&["LENIENT_UNKNOWN_MQ_ALIGNMENT=true"]), 1);
    assert_eq!(matched(&["LOW_MQ_THRESHOLD=60"]), 0);
    let output = run(&left, &right, &["COMPARE_MQ=true"])
        .0
        .expect("compare")
        .output;
    let section = output.main.expect("a file").histogram.expect("a histogram");
    let keys: Vec<String> = section
        .records()
        .iter()
        .map(|row| row[0].1.picard_text())
        .collect();
    assert_eq!(keys, ["0,0", "255,255", "255,60", "3,2", "3,3", "4,3"]);
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn header_differences_are_picards_three_lines() {
    let dir = scratch("header");
    let records = [rec("a", 0, 100, 60)];
    let left = file(&dir, "left.sam", &header("unsorted"), &records);
    let other = "@HD\tVN:1.5\tSO:unsorted\tCR:maker\n@SQ\tSN:c1\tLN:5001\tM5:aa\n\
                 @SQ\tSN:cX\tLN:3000\n@RG\tID:g1\tSM:t\tPU:u\n@PG\tID:p\n";
    let right = file(&dir, "right.sam", other, &records);
    let (result, lines) = run(&left, &right, &[]);
    assert!(!result.expect("compare").row.are_equal);
    let expected = [
        "File format version differs.",
        "File 1: 1.6",
        "File 2: 1.5",
        "File creator differs.",
        "File 1: null",
        "File 2: maker",
        "Length of sequence c1 differs.",
        "File 1: 5000",
        "File 2: 5001",
        "MD5 of sequence c1 differs.",
        "File 1: null",
        "File 2: aa",
        "Name of sequence record 2 differs.",
        "File 1: c2",
        "File 2: cX",
        "Sample for read group g1 differs.",
        "File 1: s",
        "File 2: t",
        "Library for read group g1 differs.",
        "File 1: lib",
        "File 2: null",
        "Number of program records differs.",
        "File 1: 0",
        "File 2: 1",
    ];
    assert_eq!(lines, expected);
    let (lenient, lines) = run(&left, &right, &["LENIENT_HEADER=true"]);
    assert!(!lenient.expect("compare").row.are_equal);
    assert!(
        !lines
            .iter()
            .any(|l| l.starts_with("MD5") || l.contains("program"))
    );
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn sort_orders_decide_whether_records_are_compared() {
    let dir = scratch("orders");
    let records = [rec("a", 0, 100, 60)];
    let coordinate = file(&dir, "c.sam", &header("coordinate"), &records);
    let queryname = file(&dir, "q.sam", &header("queryname"), &records);
    let (result, lines) = run(&coordinate, &queryname, &[]);
    assert_eq!(result.expect("compare").row.counts, Counts::default());
    assert_eq!(
        lines[3..],
        [
            "Sort Order differs.",
            "File 1: coordinate",
            "File 2: queryname",
            "Cannot compare alignments if sort orders differ."
        ]
    );
    let longer = format!("@HD\tVN:1.6\tSO:coordinate\n{SQ}@SQ\tSN:c3\tLN:9\n{RG}");
    let other = file(&dir, "d.sam", &longer, &records);
    let (result, lines) = run(&coordinate, &other, &[]);
    assert_eq!(result.expect("compare").row.counts.mappings_match, 0);
    assert_eq!(
        lines.last().map(String::as_str),
        Some("Cannot compare coordinate-sorted SAM files because sequence dictionaries differ.")
    );
    let unknown = file(&dir, "u.sam", &header("unknown"), &records);
    let error = run(&unknown, &unknown, &[]).0.expect_err("unknown");
    assert_eq!(
        error.to_string(),
        "Unrecognized sort order (unknown) found."
    );
    let no_hd = file(&dir, "n.sam", &format!("{SQ}{RG}"), &records);
    let unsorted = file(&dir, "s.sam", &header("unsorted"), &records);
    let (result, lines) = run(&no_hd, &unsorted, &[]);
    assert_eq!(result.expect("compare").row.counts.mappings_match, 1);
    assert_eq!(
        lines,
        ["Sort order differs.", "File 1: null", "File 2: unsorted"]
    );
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn lenient_dup_forgives_a_swapped_representative() {
    let dir = scratch("dup");
    let head = header("coordinate");
    let set = |marks: [u16; 3]| -> Vec<String> {
        ["a", "b", "c"]
            .iter()
            .zip(marks)
            .map(|(name, mark)| rec(name, mark, 100, 60))
            .collect()
    };
    let left = file(&dir, "left.sam", &head, &set([0, 1024, 1024]));
    let swapped = file(&dir, "swapped.sam", &head, &set([1024, 0, 1024]));
    let unmarked = file(&dir, "unmarked.sam", &head, &set([0, 0, 0]));
    assert_eq!(counts(&left, &swapped, &[]).duplicate_markings_differ, 2);
    let lenient = ["LENIENT_DUP=true"];
    let spilled = ["LENIENT_DUP=true", "MAX_RECORDS_IN_RAM=1"];
    for options in [&lenient[..], &spilled[..]] {
        assert_eq!(
            counts(&left, &swapped, options).duplicate_markings_differ,
            0
        );
        assert_eq!(
            counts(&left, &unmarked, options).duplicate_markings_differ,
            2
        );
    }
    let found = run(&left, &swapped, &lenient).0.expect("compare");
    assert!(found.row.are_equal);
    let error = run(
        &left,
        &swapped,
        &["LENIENT_DUP=true", "MAX_RECORDS_IN_RAM=0"],
    )
    .0;
    assert!(matches!(error, Err(Error::Usage(_))));
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn lenient_dup_reads_mc_only_when_two_pairs_differ() {
    let dir = scratch("mc");
    let head = header("coordinate");
    let pair = |first: u16, second: u16| -> Vec<String> {
        let tail = "ACGTACGTAC\tIIIIIIIIII\tRG:Z:g1\n";
        vec![
            format!("p\t{}\tc1\t100\t60\t10M\t=\t300\t210\t{tail}", 99 | first),
            format!(
                "p\t{}\tc1\t300\t60\t10M\t=\t100\t-210\t{tail}",
                147 | second
            ),
        ]
    };
    let plain = file(&dir, "plain.sam", &head, &pair(0, 0));
    let one = file(&dir, "one.sam", &head, &pair(1024, 0));
    let two = file(&dir, "two.sam", &head, &pair(1024, 1024));
    let lenient = ["LENIENT_DUP=true"];
    assert_eq!(counts(&plain, &one, &lenient).duplicate_markings_differ, 1);
    let error = run(&plain, &two, &lenient).0.expect_err("MC is needed");
    assert!(
        error.to_string().contains("Mate CIGAR (Tag MC) not found"),
        "{error}"
    );
    std::fs::remove_dir_all(&dir).expect("remove");
}

#[test]
fn the_output_may_not_be_an_input() {
    let dir = scratch("output");
    let left = file(
        &dir,
        "left.sam",
        &header("unsorted"),
        &[rec("a", 0, 100, 60)],
    );
    let before = std::fs::read(&left).expect("read");
    let (result, _) = run(&left, &left, &[&format!("OUTPUT={left}")]);
    assert!(result.is_err());
    assert_eq!(std::fs::read(&left).expect("read"), before);
    std::fs::remove_dir_all(&dir).expect("remove");
}

/// A header with one read group whose `DS` is `ds`.
fn described(ds: &str) -> String {
    format!("@HD\tVN:1.6\tSO:coordinate\n{SQ}@RG\tID:g1\tSM:s\tLB:lib\tDS:{ds}\n")
}

/// htsjdk decodes SAM text as UTF-8 and compares strings. The header the
/// records are read against holds one byte per UTF-16 unit, in which
/// `A` and U+0141 are one value, and U+0109 and U+010A a tab and a line
/// end.
#[test]
fn sam_text_header_values_are_compared_as_the_file_spells_them() {
    let dir = scratch("wide-header");
    let records = [rec("a", 0, 100, 60)];
    let pairs = [
        ("caf\u{e9}", "caf\u{e8}"),
        ("A", "\u{141}"),
        ("\u{141}", "\u{241}"),
        ("\u{6f22}\u{5b57}", "\u{6f22}\u{5b50}"),
        ("a\u{109}b", "ab"),
        ("a\u{10a}b", "ab"),
        ("a\u{100}b", "ab"),
        ("x\u{1f600}", "x\u{1f601}"),
    ];
    for (ours, theirs) in pairs {
        let left = file(&dir, "left.sam", &described(ours), &records);
        let right = file(&dir, "right.sam", &described(theirs), &records);
        let (result, lines) = run(&left, &right, &[]);
        let row = result.expect("compare").row;
        assert!(!row.are_equal, "{ours} against {theirs}");
        assert_eq!(row.counts.mappings_match, 1);
        let expected = [
            "DS for read group g1 differs.".to_owned(),
            format!("File 1: {ours}"),
            format!("File 2: {theirs}"),
        ];
        assert_eq!(lines, expected);
        let (same, lines) = run(&left, &left, &[]);
        assert!(
            same.expect("compare").row.are_equal,
            "{ours} against itself"
        );
        assert!(lines.is_empty(), "{lines:?}");
    }
    std::fs::remove_dir_all(&dir).expect("remove");
}

/// Read names of SAM text are htsjdk's strings: `A` and U+0141 are two
/// reads, though a BAM record holds the byte 0x41 for both.
#[test]
fn sam_text_read_names_are_compared_as_the_file_spells_them() {
    let dir = scratch("wide-names");
    let pairs = [
        ("A", "\u{141}"),
        ("caf\u{e9}", "caf\u{e8}"),
        ("\u{6f22}\u{5b57}", "\u{6f22}\u{5b50}"),
        ("x\u{1f600}", "x\u{1f601}"),
    ];
    for order in ["coordinate", "queryname", "unsorted"] {
        let head = header(order);
        for (ours, theirs) in pairs {
            // The queryname merge fails when the right file ends first.
            let (ours, theirs) = if ours.encode_utf16().lt(theirs.encode_utf16()) {
                (ours, theirs)
            } else {
                (theirs, ours)
            };
            // A secondary record first: every record read has an index.
            let records = |name: &str| [rec("s", 256, 50, 60), rec(name, 0, 100, 60)];
            let left = file(&dir, "left.sam", &head, &records(ours));
            let right = file(&dir, "right.sam", &head, &records(theirs));
            let found = counts(&left, &right, &[]);
            assert_eq!(
                (
                    found.missing_left,
                    found.missing_right,
                    found.mappings_match
                ),
                (1, 1, 0),
                "{order}: {ours} against {theirs}"
            );
            let same = counts(&right, &right, &[]);
            assert_eq!((same.mappings_match, same.missing_left), (1, 0), "{theirs}");
        }
    }
    // The queryname merge orders names by UTF-16 unit: U+1F600 (a
    // surrogate pair) before U+FFE5, so the right file does not end first.
    let head = header("queryname");
    let left = file(&dir, "left.sam", &head, &[rec("\u{1f600}", 0, 100, 60)]);
    let right = file(&dir, "right.sam", &head, &[rec("\u{ffe5}", 0, 100, 60)]);
    let found = counts(&left, &right, &[]);
    assert_eq!((found.missing_left, found.missing_right), (1, 1));
    assert!(run(&right, &left, &[]).0.is_err());
    std::fs::remove_dir_all(&dir).expect("remove");
}

/// An unmapped record with a mapped mate and no `MC` is scored by htsjdk
/// only against a record at its place: alone there it counts, with a
/// second one there Picard fails, and so does this.
#[test]
fn lenient_dup_scores_an_unmapped_record_only_on_a_tie() {
    let dir = scratch("unscored");
    let head = header("coordinate");
    let tail = "ACGTACGTAC\tIIIIIIIIII\tRG:Z:g1";
    let pair = |name: &str, at: i32, mark: u16| -> Vec<String> {
        vec![
            format!("{name}\t73\tc1\t{at}\t60\t10M\t=\t{at}\t0\t{tail}\n"),
            format!(
                "{name}\t{}\tc1\t{at}\t0\t*\t=\t{at}\t0\t{tail}\n",
                133 | mark
            ),
        ]
    };
    let two = |second: i32, mark: u16| [pair("u", 100, mark), pair("v", second, mark)].concat();
    let lenient = ["LENIENT_DUP=true"];
    let left = file(&dir, "left.sam", &head, &two(200, 0));
    let right = file(&dir, "right.sam", &head, &two(200, 1024));
    assert_eq!(counts(&left, &right, &lenient).duplicate_markings_differ, 2);
    let left = file(&dir, "left.sam", &head, &two(100, 0));
    let right = file(&dir, "right.sam", &head, &two(100, 1024));
    let error = run(&left, &right, &lenient).0.expect_err("a tie is scored");
    assert!(
        error.to_string().contains("Mate CIGAR (Tag MC) not found"),
        "{error}"
    );
    std::fs::remove_dir_all(&dir).expect("remove");
}
