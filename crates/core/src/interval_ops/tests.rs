use std::path::{Path, PathBuf};

use super::ops::{OverlapTree, compare};
use super::*;

const DICT: &str = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:500\n@SQ\tSN:c3\tLN:0\n";

fn list(records: &str) -> IntervalList {
    IntervalList::parse(&format!("{DICT}{records}")).expect("parses")
}

fn iv(contig: usize, start: i32, end: i32, negative: bool, name: Option<&str>) -> Interval {
    Interval {
        contig,
        start,
        end,
        negative_strand: negative,
        name: name.map(str::to_owned),
    }
}

/// The interval lines of a list, as written.
fn body(list: &IntervalList) -> Vec<String> {
    list.to_text()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(str::to_owned)
        .collect()
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-interval-ops-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

#[test]
fn the_header_is_decoded_and_encoded_as_htsjdk_does() {
    let text = "@CO\tfirst\n@SQ\tSN:c1\tLN:+100\tM5:ab\n@HD\tSO:coordinate\tVN:1.0\tXX:y\n\
                @SQ\tSN:nolength\n@XX\tunknown\n@RG\tID:r\tSM:s\n@PG\tID:7\tPN:p\n@CO\tsecond\n\
                c1\t5\t9\t+\tn\n";
    let parsed = IntervalList::parse(text).expect("parses");
    assert_eq!(
        parsed.to_text(),
        "@HD\tVN:1.6\tSO:coordinate\tXX:y\n@SQ\tSN:c1\tLN:100\tM5:ab\n@RG\tID:r\tSM:s\n\
         @PG\tID:7\tPN:p\n@CO\tfirst\n@CO\tsecond\nc1\t5\t9\t+\tn\n"
    );
    assert_eq!(parsed.header.sort_order(), SortOrder::Coordinate);
    let bogus = IntervalList::parse("@HD\tVN:1.6\tSO:bogus\n").expect("parses");
    assert_eq!(bogus.header.sort_order(), SortOrder::Unknown);
    assert_eq!(bogus.to_text(), "@HD\tVN:1.6\tSO:unknown\n");
    let none = IntervalList::parse("@SQ\tSN:c1\tLN:5\n").expect("parses");
    assert_eq!(none.header.sort_order(), SortOrder::Unsorted);
    assert_eq!(none.to_text(), "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:5\n");
}

#[test]
fn header_errors_are_htsjdks() {
    let message = |text: &str| IntervalList::parse(text).expect_err("fails");
    assert_eq!(
        message("@SQ\tSN:c1\tLN:5\n@SQ\tSN:c1\tLN:5\n"),
        "Cannot add sequence that already exists in SAMSequenceDictionary: c1"
    );
    assert_eq!(message("@SQ\tSN:c1\tLN:x\n"), "For input string: \"x\"");
    assert!(message("@SQ\tSN:*c\tLN:5\n").contains("doesn't match regex"));
    assert!(message("@PG\tID:1\n@PG\tID:1\n").contains("already exists"));
    assert!(message("@RG\tID:1\n@RG\tID:1\n").contains("already exists"));
    assert_eq!(
        message("c1\t1\t2\t+\tn\n"),
        "Interval list file must contain header. "
    );
    assert_eq!(message(""), "Interval list file must contain header. ");
    // A name is cut at its first space.
    let cut = IntervalList::parse("@SQ\tSN:c1 extra\tLN:5\nc1\t1\t2\t+\tn\n").expect("parses");
    assert_eq!(cut.header.sequences[0].name, "c1");
}

#[test]
fn records_follow_the_codec() {
    let records = "c2\t10\t20\t-\tb\n\n  \n@SQ\tSN:late\tLN:9\ncX\t1\t2\t+\tskipped\n\
                   c1\t+5\t6\t+\tn\t\r\nc1\t7\t6\t+\tzero\rc1\t1\t1000\t+\tx y|z\n";
    let parsed = list(records);
    assert_eq!(
        body(&parsed),
        [
            "c2\t10\t20\t-\tb",
            "c1\t5\t6\t+\tn",
            "c1\t7\t6\t+\tzero",
            "c1\t1\t1000\t+\tx y|z"
        ]
    );
    let error =
        |records: &str| IntervalList::parse(&format!("{DICT}{records}")).expect_err("fails");
    assert_eq!(
        error("c1\t1\t5\t+\t\n"),
        "Invalid interval record contains 4 fields: c1\t1\t5\t+\t"
    );
    assert_eq!(
        error("c1\t0\t5\t+\tn\n"),
        "Coordinate less than 1: start value of 0 is less than 1 and thus illegal"
    );
    assert!(error("c1\t9\t7\t+\tn\n").starts_with("Start value of 9 is greater than end + 1"));
    assert_eq!(error("c1\t1\t5\t.\tn\n"), "Invalid strand field: .");
    assert_eq!(
        error("c2\t1\t501\t+\tn\n"),
        "interval with end: 501 extends beyond end of sequence with length: 500"
    );
    assert_eq!(error("c1\tfive\t9\t+\tn\n"), "For input string: \"five\"");
    // A contig of length 0 takes any end.
    assert_eq!(body(&list("c3\t1\t99\t+\tn\n")), ["c3\t1\t99\t+\tn"]);
}

#[test]
fn the_comparator_is_interval_coordinate_comparator() {
    let mut intervals = [
        iv(1, 1, 5, false, Some("a")),
        iv(0, 5, 9, true, Some("a")),
        iv(0, 5, 9, false, Some("b")),
        iv(0, 5, 9, false, None),
        iv(0, 5, 9, false, Some("\u{ff61}")),
        iv(0, 5, 9, false, Some("\u{1f600}")),
        iv(0, 5, 8, false, Some("z")),
    ];
    intervals.sort_by(compare);
    let names: Vec<_> = intervals.iter().map(|i| i.name.as_deref()).collect();
    // A supplementary character sorts by its high surrogate, before U+FF61.
    assert_eq!(
        names,
        [
            Some("z"),
            None,
            Some("b"),
            Some("\u{1f600}"),
            Some("\u{ff61}"),
            Some("a"),
            Some("a")
        ]
    );
}

#[test]
fn merging_follows_the_merger_iterator() {
    let input = vec![
        iv(0, 10, 9, false, Some("zero")),
        iv(0, 11, 20, true, Some("a")),
        iv(0, 21, 30, false, Some("b")),
        iv(0, 25, 26, false, Some("a")),
        iv(0, 32, 31, false, Some("gap")),
        iv(1, 1, 2, false, None),
    ];
    let joined: Vec<_> = merge_intervals(input.clone(), true, true)
        .into_iter()
        .map(|i| (i.contig, i.start, i.end, i.negative_strand, i.name))
        .collect();
    assert_eq!(
        joined,
        [
            (0, 10, 30, false, Some("zero|a|b".to_owned())),
            (0, 32, 31, false, Some("gap".to_owned())),
            (1, 1, 2, false, None),
        ]
    );
    // Without abutting merges, a zero-length interval at 10 still overlaps
    // nothing at 11; 11-20 and 21-30 stay apart.
    let apart = merge_intervals(input.clone(), false, true);
    assert_eq!(apart.len(), 5);
    // Without name joining, the first member's name and the current start.
    let first = merge_intervals(input, true, false);
    assert_eq!(first[0].name.as_deref(), Some("zero"));
    // An unsorted list under `SO:coordinate` is merged as it stands.
    let trusted = IntervalList::parse("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:100\nc1\t50\t60\t+\ta\nc1\t5\t4\t+\tz\nc1\t61\t70\t+\tb\n")
    .expect("parses");
    assert_eq!(
        body(&IntervalList {
            header: trusted.header.clone(),
            intervals: unique_intervals(&trusted, true, true)
        }),
        ["c1\t50\t60\t+\ta", "c1\t5\t4\t+\tz", "c1\t61\t70\t+\tb"]
    );
}

#[test]
fn an_interval_ending_before_it_starts_joins_by_its_unmoved_end() {
    // Picard VcfToIntervalList on c1:100 END=110, c1:120 END=110 and
    // c1:111 END=130 writes `c1 100 130 + a|b|c`: the second record's end
    // 110 lies inside 100-110, though 119-111 (moved one base) does not.
    let spans = |input: Vec<Interval>, concatenate: bool| -> Vec<(i32, i32, Option<String>)> {
        merged(input, true, concatenate)
            .map(|i| (i.start, i.end, i.name))
            .collect()
    };
    let two = vec![
        iv(0, 100, 110, false, Some("a")),
        iv(0, 120, 110, false, Some("b")),
    ];
    assert_eq!(
        spans(two.clone(), true),
        [(100, 110, Some("a|b".to_owned()))]
    );
    assert_eq!(spans(two, false), [(100, 110, Some("a".to_owned()))]);
    let three = vec![
        iv(0, 100, 110, false, Some("a")),
        iv(0, 120, 110, false, Some("b")),
        iv(0, 111, 130, false, Some("c")),
    ];
    assert_eq!(spans(three, true), [(100, 130, Some("a|b|c".to_owned()))]);
    // Without abutting merges only the overlap test runs.
    let apart = vec![
        iv(0, 100, 110, false, Some("a")),
        iv(0, 111, 120, false, Some("b")),
    ];
    assert_eq!(merged(apart, false, true).count(), 2);
}

#[test]
fn the_feature_reader_fails_on_an_empty_line_after_the_header() {
    // Picard IntervalListToBed: "Unable to parse header with error: Index
    // 0 out of bounds for length 0" for each of the first five; a blank
    // line inside or after the records is skipped.
    let dir = scratch("blank");
    let path = dir.join("in.interval_list");
    let read = |text: &str| {
        std::fs::write(&path, text).expect("writes");
        IntervalList::read_features(&path)
    };
    for bad in [
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:100\n\nc1\t1\t5\t+\tn\n",
        "@HD\tVN:1.6\n\n@SQ\tSN:c1\tLN:100\nc1\t1\t5\t+\tn\n",
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:100\n\n",
        "\n@HD\tVN:1.6\n@SQ\tSN:c1\tLN:100\n",
        "\n",
    ] {
        let error = read(bad).expect_err(bad).to_string();
        assert!(error.contains("Index 0 out of bounds"), "{error}");
    }
    for good in [
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:100\nc1\t1\t5\t+\tn\n\nc1\t7\t9\t+\tm\n",
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:100\nc1\t1\t5\t+\tn\n\n",
        "",
    ] {
        assert!(read(good).is_ok(), "{good:?}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn padding_clamps_and_refuses_a_negative_value() {
    let input = list("c1\t3\t998\t-\ta\nc2\t100\t99\t+\tz\nc3\t1\t10\t+\tq\n");
    assert_eq!(
        body(&input.padded(5, 5).expect("pads")),
        ["c1\t1\t1000\t-\ta", "c2\t95\t104\t+\tz", "c3\t1\t0\t+\tq"]
    );
    assert!(matches!(
        input.padded(-1, 0),
        Err(crate::error::Error::Usage(_))
    ));
    // Java `int` arithmetic: the end wraps, then the contig length is the
    // larger value.
    let wrapped = input.padded(0, i32::MAX).expect("pads");
    assert!(wrapped.intervals[0].end < 0);
}

#[test]
fn set_operations_keep_htsjdks_names_and_orders() {
    let a = list("c1\t10\t20\t+\ta1\nc1\t15\t30\t-\ta2\nc1\t15\t30\t+\tsame\nc2\t5\t9\t+\ta3\n");
    let b = list("c1\t18\t40\t-\tb1\nc2\t1\t500\t+\tb2\nc1\t25\t24\t+\tzero\n");
    // `i.intersect(j)` for i of b: b's strand and "b intersection a".
    assert_eq!(
        body(&intersection(&a, &b).expect("intersects")),
        [
            "c1\t18\t30\t-\tb1 intersection a1|b1 intersection a2",
            "c2\t5\t9\t+\tb2 intersection a3"
        ]
    );
    let concatenated = concatenate(vec![a.clone(), b.clone()]).expect("concatenates");
    assert_eq!(concatenated.header.sort_order(), SortOrder::Unsorted);
    assert_eq!(concatenated.len(), 7);
    assert_eq!(
        body(&union(vec![a.clone()]).expect("unites")),
        ["c1\t10\t30\t+\ta1|same|a2", "c2\t5\t9\t+\ta3"]
    );
    assert_eq!(
        body(&subtract(&a, &b).expect("subtracts")),
        ["c1\t10\t17\t+\tinterval-1 intersection a1|interval-1 intersection a2"]
    );
    assert_eq!(
        body(&overlaps(&a, &b).expect("overlaps")),
        [
            "c1\t10\t20\t+\ta1",
            "c1\t15\t30\t-\ta2",
            "c1\t15\t30\t+\tsame",
            "c2\t5\t9\t+\ta3"
        ]
    );
    assert_eq!(
        body(&difference(&a, &b).expect("differs")),
        [
            "c1\t10\t17\t+\tinterval-1 intersection a1|interval-1 intersection a2",
            "c1\t31\t40\t+\tinterval-2 intersection b1",
            "c2\t1\t4\t+\tinterval-3 intersection b2",
            "c2\t10\t500\t+\tinterval-4 intersection b2"
        ]
    );
    let other = IntervalList::parse("@SQ\tSN:c1\tLN:999\n@SQ\tSN:c2\tLN:500\n@SQ\tSN:c3\tLN:0\n")
        .expect("parses");
    assert!(intersection(&a, &other).is_err());
    assert!(concatenate(Vec::new()).is_err());
}

#[test]
fn invert_names_each_gap_and_splits_at_a_zero_length_interval() {
    let input = list("c1\t100\t200\t+\tx\nc1\t300\t299\t+\tz\nc1\t50\t120\t-\ty\n");
    assert_eq!(
        body(&invert(&input)),
        [
            "c1\t1\t49\t+\tinterval-1",
            "c1\t201\t299\t+\tinterval-2",
            "c1\t300\t1000\t+\tinterval-3",
            "c2\t1\t500\t+\tinterval-4"
        ]
    );
    assert_eq!(
        body(&invert(&list(""))),
        ["c1\t1\t1000\t+\tinterval-1", "c2\t1\t500\t+\tinterval-2"]
    );
}

#[test]
fn bands_break_at_multiples() {
    let input = list("c1\t5\t25\t+\tn\nc1\t10\t19\t+\tinside\nc1\t30\t29\t+\tzero\n");
    let mut intervals = input.intervals.clone();
    intervals.push(iv(0, 9, 10, false, None));
    let broken = break_at_band_multiples(&intervals, 10).expect("breaks");
    let shown: Vec<_> = broken
        .iter()
        .map(|i| (i.start, i.end, i.name.clone().unwrap_or_default()))
        .collect();
    assert_eq!(
        shown,
        [
            (5, 9, "n.1".to_owned()),
            (10, 19, "n.2".to_owned()),
            (20, 25, "n.3".to_owned()),
            (10, 19, "inside".to_owned()),
            (30, 29, "zero".to_owned()),
            (9, 9, "null.1".to_owned()),
            (10, 10, "null.2".to_owned()),
        ]
    );
    assert_eq!(
        break_at_band_multiples(&intervals, 0).expect("kept"),
        intervals
    );
    let whole = [iv(0, 1, i32::MAX, false, None)];
    assert!(break_at_band_multiples(&whole, 1).is_err());
}

#[test]
fn the_overlap_tree_finds_what_a_scan_finds() {
    let mut seed: u64 = 7;
    let mut next = |below: i32| {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 33) % below as u64) as i32
    };
    for _ in 0..50 {
        let intervals: Vec<Interval> = (0..60)
            .map(|_| {
                let start = next(500) + 1;
                iv(next(2) as usize, start, start + next(80) - 1, false, None)
            })
            .collect();
        let tree = OverlapTree::new(&intervals, 2);
        for _ in 0..40 {
            let (contig, start) = (next(2) as usize, next(600));
            let end = start + next(50) - 3;
            let mut found = Vec::new();
            tree.for_each_overlap(contig, start, end, |at| found.push(at));
            found.sort_unstable();
            let mut expected: Vec<usize> = Vec::new();
            for (at, i) in intervals.iter().enumerate() {
                let earlier_twin = intervals[..at]
                    .iter()
                    .any(|j| (j.contig, j.start, j.end) == (i.contig, i.start, i.end));
                if i.contig == contig
                    && i.start <= i.end
                    && start <= end
                    && i.start <= end
                    && i.end >= start
                    && !earlier_twin
                {
                    expected.push(at);
                }
            }
            assert_eq!(found, expected);
        }
    }
}

#[test]
fn the_dictionary_check_allows_unknown_lengths_and_compares_md5_as_numbers() {
    let a = Header::decode(["@SQ\tSN:c1\tLN:0\tM5:00FF", "@SQ\tSN:c2\tLN:5"]).expect("decodes");
    let b = Header::decode(["@SQ\tSN:c1\tLN:9\tM5:ff", "@SQ\tSN:c2\tLN:5"]).expect("decodes");
    assert!(a.assert_same_dictionary(&b).is_ok());
    let c = Header::decode(["@SQ\tSN:c1\tLN:9\tM5:fe", "@SQ\tSN:c2\tLN:5"]).expect("decodes");
    assert!(b.assert_same_dictionary(&c).is_err());
    let d = Header::decode(["@SQ\tSN:c1\tLN:9"]).expect("decodes");
    assert_eq!(
        b.assert_same_dictionary(&d).expect_err("sizes"),
        "Sequence dictionaries are not the same size (2, 1)"
    );
}

#[test]
fn programs_and_comments_are_added_as_interval_list_tools_adds_them() {
    let mut header = Header::decode(["@PG\tID:1\tPN:a", "@PG\tID:3\tPN:b"]).expect("decodes");
    header.add_program("cmd", "IntervalListTools");
    header.add_comment("note");
    header.add_comment("@CO\tkept");
    header.add_comment("@COx");
    assert_eq!(
        header.encode(),
        "@HD\tVN:1.6\n@PG\tID:1\tPN:a\n@PG\tID:3\tPN:b\n@PG\tID:2\tCL:cmd\tPN:IntervalListTools\n\
         @CO\tnote\n@CO\tkept\n@CO\t@COx\n"
    );
}

#[test]
fn scatter_by_subdivision_splits_to_the_ideal_weight() {
    let input = list("c1\t201\t300\t+\tb\nc1\t1\t100\t+\ta\n");
    let lists = scatter(&input, ScatterMode::IntervalSubdivision, 3).expect("scatters");
    let shown: Vec<Vec<(i32, i32)>> = lists
        .parts
        .iter()
        .map(|p| p.iter().map(|i| (i.start, i.end)).collect())
        .collect();
    assert_eq!(lists.header, input.uniqued(true).header);
    assert_eq!((lists.base_count(), lists.interval_count()), (200, 4));
    assert_eq!(
        shown,
        [vec![(1, 66)], vec![(67, 100), (201, 232)], vec![(233, 300)]]
    );
    // By count, more lists than intervals: one per interval.
    let by_count = scatter(&input, ScatterMode::IntervalCount, 5).expect("scatters");
    assert_eq!(by_count.len(), 2);
    assert!(scatter(&input, ScatterMode::IntervalCount, 0).is_err());
    assert_eq!(
        scatter_count_for_content(ScatterMode::IntervalSubdivision, &input, 0),
        -1
    );
    assert_eq!(
        scatter_count_for_content(ScatterMode::IntervalSubdivision, &input, 80),
        3
    );
    assert_eq!(
        scatter_path(Path::new("out"), 12_345, 7),
        Path::new("out/temp_12345_of_7/scattered.interval_list")
    );
    for mode in ScatterMode::ALL {
        assert_eq!(ScatterMode::parse(mode.name()), Some(mode));
    }
}

#[test]
fn a_scatter_is_written_as_one_set() {
    let dir = scratch("scatter");
    let input_path = dir.join("in.interval_list");
    let input = list("c1\t1\t100\t+\ta\nc1\t201\t300\t+\tb\n");
    write_list(&input, &input_path, &[]).expect("writes");
    assert_eq!(
        IntervalList::read(&input_path).expect("reads"),
        input,
        "a written list reads back"
    );
    let lists = scatter(&input, ScatterMode::IntervalCount, 2).expect("scatters");
    let out = dir.join("out");
    let written =
        write_scatter(&out, &lists, 2, std::slice::from_ref(&input_path)).expect("writes");
    assert_eq!(
        written,
        [
            out.join("temp_0001_of_2/scattered.interval_list"),
            out.join("temp_0002_of_2/scattered.interval_list")
        ]
    );
    // A file of the set that is an input refuses the whole set before any
    // file or directory is made.
    let clash = dir.join("clash");
    let as_input = scatter_path(&clash, 2, 2);
    std::fs::create_dir_all(as_input.parent().expect("parent")).expect("dir");
    std::fs::write(&as_input, "x").expect("write");
    assert!(write_scatter(&clash, &lists, 2, std::slice::from_ref(&as_input)).is_err());
    assert!(!scatter_path(&clash, 1, 2).exists());
    assert_eq!(std::fs::read_to_string(&as_input).expect("kept"), "x");
    // A list written over its own input is refused.
    assert!(write_list(&input, &input_path, std::slice::from_ref(&input_path)).is_err());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_gzip_name_is_written_compressed_and_read_back() {
    let dir = scratch("gzip");
    let path = dir.join("list.interval_list.gz");
    let input = list("c2\t5\t9\t-\tq\n");
    write_list(&input, &path, &[]).expect("writes");
    let bytes = std::fs::read(&path).expect("reads");
    assert_eq!(&bytes[..2], &[0x1f, 0x8b]);
    assert_eq!(IntervalList::read(&path).expect("inflates"), input);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_vcf_gives_one_interval_per_site() {
    let dir = scratch("vcf");
    let path = dir.join("v.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n##contig=<ID=c1,length=1000,assembly=b37,md5=abc>\n##contig=<ID=c2>\n\
         ##FILTER=<ID=q10,Description=\"x\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         c1\t10\t.\tA\tG\t.\t.\t.\nc1\t20\trs1;rs2\tACG\tA\t.\tPASS\t.\n\
         c1\t30\t.\tA\t<DEL>\t.\tq10\tEND=40\nc2\t5\tid5\tA\t<NON_REF>\t.\t.\tEND=9\n\
         c1\t50\t.\tA\tT\t.\t.\t.\n",
    )
    .expect("write");
    // Picard 3.5.0 VcfToIntervalList on this file.
    assert_eq!(
        IntervalList::from_vcf(&path, false)
            .expect("reads")
            .to_text(),
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\tAS:b37\n@SQ\tSN:c2\tLN:0\n\
         c1\t10\t10\t+\tinterval-1\nc1\t20\t22\t+\trs1;rs2\nc2\t5\t9\t+\tid5\nc1\t50\t50\t+\tinterval-2\n"
    );
    assert_eq!(
        body(&IntervalList::from_vcf(&path, true).expect("reads")),
        [
            "c1\t10\t10\t+\tinterval-1",
            "c1\t20\t22\t+\trs1;rs2",
            "c1\t30\t40\t+\tinterval-2",
            "c2\t5\t9\t+\tid5",
            "c1\t50\t50\t+\tinterval-3"
        ]
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_bed_is_read_as_bed_to_interval_list_reads_it() {
    let dir = scratch("bed");
    let dict = dir.join("ref.dict");
    std::fs::write(&dict, "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\tM5:ab\n").expect("write");
    let bed = dir.join("x.bed");
    std::fs::write(&bed, "c1\t9\t20\tn\t0\t-\nc1\t0\t5\n").expect("write");
    let opts = crate::bed_to_interval_list::Options {
        sequence_dictionary: dict,
        sort: true,
        unique: false,
        drop_missing_contigs: false,
        keep_length_zero_intervals: false,
    };
    let parsed = IntervalList::from_bed(&bed, &opts).expect("reads");
    assert_eq!(
        parsed.to_text(),
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\tM5:ab\nc1\t10\t20\t-\tn\nc1\t1\t5\t+\t.\n"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `qc --panel` and the hs tools read a list this layer writes: the
/// interval list reader of `intervals.rs` takes it against a BAM header
/// with the same `@SQ` lines.
#[test]
fn the_panel_reader_takes_a_written_list() {
    let dir = scratch("panel");
    let path = dir.join("panel.interval_list");
    let mut written = IntervalList::parse(
        "@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:500\nc2\t5\t9\t-\tq\nc1\t1\t10\t+\tp\n",
    )
    .expect("parses")
    .sorted();
    written.header.add_program("cmd", "IntervalListTools");
    written.header.add_comment("made here");
    write_list(&written, &path, &[]).expect("writes");
    let header: noodles::sam::Header = "@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:500\n"
        .parse()
        .expect("header");
    let read = crate::intervals::read_named_files(&[&path], &header).expect("reads");
    let shown: Vec<_> = read
        .iter()
        .map(|i| {
            (
                i.interval.contig,
                i.interval.start,
                i.interval.end,
                i.name.as_str(),
            )
        })
        .collect();
    assert_eq!(shown, [(0, 1, 10, "p"), (1, 5, 9, "q")]);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn reading_stops_at_each_count_cap() {
    let limits = |header_lines, header_fields, intervals| read::Limits {
        header_lines,
        header_fields,
        intervals,
    };
    let text = format!("{DICT}c1\t1\t5\t+\ta\nc1\t7\t9\t+\tb\n");
    let text = text.as_bytes();
    assert!(IntervalList::parse_limited(text, &limits(4, 7, 2), true).is_ok());
    let header = IntervalList::parse_limited(text, &limits(3, 7, 2), true).expect_err("lines");
    assert!(header.contains("more than 3 lines"), "{header}");
    // `VN`, then `SN` and `LN` on each of three `@SQ` lines.
    let fields = IntervalList::parse_limited(text, &limits(4, 6, 2), true).expect_err("fields");
    assert!(fields.contains("more than 6 fields"), "{fields}");
    let intervals =
        IntervalList::parse_limited(text, &limits(4, 7, 1), true).expect_err("intervals");
    assert!(intervals.contains("more than 1 intervals"), "{intervals}");
}

#[test]
fn a_long_bad_line_is_cut_in_the_message() {
    let line = "x".repeat(5000);
    let error = IntervalList::parse(&format!("{DICT}{line}\n")).expect_err("refused");
    assert!(error.len() < 400, "{} bytes", error.len());
    // The record error quotes 200 characters of the line (`error::excerpt`);
    // `cut` bounds any other message at 300 bytes.
    assert!(error.ends_with("x..."), "{error}");
    assert_eq!(read::cut("y".repeat(5000)).len(), 300 + "... (5000 bytes)".len());
}

#[test]
fn bytes_that_are_not_utf8_read_as_replacement_characters() {
    let mut text = DICT.as_bytes().to_vec();
    text.extend_from_slice(b"c1\t1\t5\t+\tn\xff\xfeame\r\nc2\t2\t3\t-\t\xe2\x82\n");
    let parsed = IntervalList::parse_bytes(&text).expect("parses");
    assert_eq!(
        body(&parsed),
        ["c1\t1\t5\t+\tn\u{fffd}\u{fffd}ame", "c2\t2\t3\t-\t\u{fffd}"]
    );
}

#[test]
fn an_intersection_past_its_budget_is_refused() {
    // Three nested intervals on each side: nine pairs.
    let a = list("c1\t1\t100\t+\ta1\nc1\t2\t99\t+\ta2\nc1\t3\t98\t+\ta3\n");
    let b = list("c1\t4\t97\t+\tb1\nc1\t5\t96\t+\tb2\nc1\t6\t95\t+\tb3\n");
    let whole = intersection(&a, &b).expect("within the caps");
    assert_eq!(
        ops::intersection_limited(&a, &b, 9, 1 << 20).expect("nine pairs"),
        whole
    );
    let pairs = ops::intersection_limited(&a, &b, 8, 1 << 20).expect_err("pairs");
    assert!(
        pairs.to_string().contains("more than 8 overlapping pairs"),
        "{pairs}"
    );
    // Each name is `bN intersection aN`, 18 bytes.
    assert!(ops::intersection_limited(&a, &b, 9, 9 * 18).is_ok());
    assert!(ops::intersection_limited(&a, &b, 9, 9 * 18 - 1).is_err());
}

#[test]
fn band_pieces_with_long_names_are_refused_before_they_are_made() {
    let name = "n".repeat(20_000);
    let long = [iv(0, 1, 20_000, false, Some(&name))];
    let error = break_at_band_multiples(&long, 1).expect_err("too many name bytes");
    assert!(error.to_string().contains("bytes of names"), "{error}");
    assert_eq!(
        break_at_band_multiples(&long, 10_000)
            .expect("three pieces")
            .len(),
        3
    );
}

#[test]
fn a_scatter_is_capped_by_the_lists_it_makes() {
    let input = list("c1\t1\t5\t+\ta\nc2\t1\t5\t+\tb\n");
    // Picard writes as many lists as the intervals allow.
    let lists = scatter(&input, ScatterMode::IntervalCount, 100_001).expect("two lists");
    assert_eq!(lists.len(), 2);
    let big = IntervalList::parse("@HD\tVN:1.6\n@SQ\tSN:c1\tLN:200000\nc1\t1\t200000\t+\tw\n")
        .expect("parses");
    let error = scatter(&big, ScatterMode::IntervalSubdivision, 200_000).expect_err("over");
    assert!(
        error.to_string().contains("more than 100000 lists"),
        "{error}"
    );
}

#[test]
fn the_scatter_bytes_are_the_bytes_written() {
    let dir = scratch("scatter-bytes");
    let input = list("c1\t1\t100\t+\ta\nc1\t201\t300\t-\t.\nc2\t9\t8\t+\tlong name\n");
    let lists = scatter(&input, ScatterMode::IntervalCount, 2).expect("scatters");
    let written = write_scatter(&dir.join("out"), &lists, 2, &[]).expect("writes");
    let on_disk: u64 = written
        .iter()
        .map(|p| std::fs::metadata(p).expect("written").len())
        .sum();
    assert_eq!(scatter_bytes(&lists), on_disk);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_vcf_contig_outside_the_name_pattern_is_refused() {
    let dir = scratch("vcf-name");
    let path = dir.join("v.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n##contig=<ID=tig 7,length=999>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\ntig 7\t10\trs1\tA\tC\t.\tPASS\t.\n",
    )
    .expect("write");
    let error = IntervalList::from_vcf(&path, false).expect_err("refused");
    assert!(
        error
            .to_string()
            .contains("Sequence name 'tig 7' doesn't match regex"),
        "{error}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// One gzip member holding `data` in a stored deflate block, as `gzip`
/// writes one without the BGZF extra field.
fn plain_gzip(data: &[u8]) -> Vec<u8> {
    let mut crc = !0u32;
    for &byte in data {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = if crc & 1 == 1 {
                (crc >> 1) ^ 0xedb8_8320
            } else {
                crc >> 1
            };
        }
    }
    let length = u16::try_from(data.len()).expect("one stored block");
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3, 1];
    out.extend_from_slice(&length.to_le_bytes());
    out.extend_from_slice(&(!length).to_le_bytes());
    out.extend_from_slice(data);
    out.extend_from_slice(&(!crc).to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    out
}

#[test]
fn a_plain_gzip_vcf_is_read() {
    let dir = scratch("vcf-gzip");
    let text = "##fileformat=VCFv4.2\n##contig=<ID=c1,length=500>\n\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                c1\t10\trs1\tA\tC\t.\tPASS\t.\nc1\t20\t.\tA\tC\t.\tq\t.\n";
    let path = dir.join("v.vcf.gz");
    std::fs::write(&path, plain_gzip(text.as_bytes())).expect("write");
    assert_eq!(
        body(&IntervalList::from_vcf(&path, true).expect("plain gzip")),
        ["c1\t10\t10\t+\trs1", "c1\t20\t20\t+\tinterval-1"]
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_vcf_header_past_its_line_cap_is_refused() {
    let dir = scratch("vcf-header");
    let path = dir.join("v.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n##contig=<ID=c1,length=500>\n##contig=<ID=c2,length=500>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
    .expect("write");
    assert!(crate::vcf::VcfReader::open_capped_header(&path, 1 << 20, 100, 4).is_ok());
    assert!(crate::vcf::VcfReader::open_capped_header(&path, 1 << 20, 100, 3).is_err());
    let _ = std::fs::remove_dir_all(&dir);
}
