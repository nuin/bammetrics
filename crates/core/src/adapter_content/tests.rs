use std::fmt::Write as _;
use std::fs;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::*;
use crate::sequence::reverse_complement;

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:10000
@RG\tID:a\tSM:s1\tLB:l1\tPU:u1
@RG\tID:b\tSM:s1\tLB:l2\tPU:u2
";

const UNIVERSAL: &str = "AGATCGGAAGAG";
const NEXTERA: &str = "CTGTCTCTTATA";
const LENGTH: usize = 50;

/// A read of `length` bases in sequencing order with no adapter key: the
/// four-base cycle `ACGT`, which holds no run of 12 equal bases and none of
/// the four other keys.
fn background(length: usize) -> String {
    "ACGT".repeat(length / 4 + 1)[..length].to_owned()
}

/// A read with `adapter` at 0-based `start`, the rest background.
fn with_adapter(adapter: &str, start: usize, length: usize) -> String {
    let mut read = background(length).into_bytes();
    read[start..start + adapter.len()].copy_from_slice(adapter.as_bytes());
    String::from_utf8(read).expect("ascii")
}

/// One SAM line of a read given in sequencing order: with flag 0x10 the
/// record stores its reverse complement.
fn line(out: &mut String, name: &str, flag: u16, read: &str, tlen: i64, rg: &str) {
    let stored = if flag & 0x10 != 0 {
        String::from_utf8(reverse_complement(read.as_bytes())).expect("ascii")
    } else {
        read.to_owned()
    };
    let (contig, pos, mapq, cigar) = if flag & 0x4 != 0 {
        ("*", 0, 0, "*".to_owned())
    } else {
        ("chr1", 100, 60, format!("{}M", read.len()))
    };
    let (mate_contig, mate_pos) = if flag & 0x1 != 0 && flag & 0x8 == 0 {
        ("=", 100)
    } else {
        ("*", 0)
    };
    let qual = "I".repeat(read.len());
    writeln!(
        out,
        "{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar}\t{mate_contig}\t{mate_pos}\t{tlen}\t{stored}\t{qual}\tRG:Z:{rg}"
    )
    .expect("write");
}

/// The planted unpaired file: 10 counted reads of 50 bases.
/// * 4 carry the Universal key at cycle 11 (start 10): 2 forward, 1
///   reverse strand, 1 unmapped dimer of read group b.
/// * 2 carry the Nextera key at cycle 21, one on the reverse strand.
/// * 1 is 8 bases long, shorter than the key.
/// * 3 carry nothing.
///
/// Then records that must not count: a secondary and a supplementary copy
/// of a Universal read.
fn unpaired() -> String {
    let mut sam = String::from(HEADER);
    let universal = with_adapter(UNIVERSAL, 10, LENGTH);
    let nextera = with_adapter(NEXTERA, 20, LENGTH);
    line(&mut sam, "u1", 0, &universal, 0, "a");
    line(&mut sam, "u2", 0, &universal, 0, "a");
    line(&mut sam, "u3", 16, &universal, 0, "a");
    line(&mut sam, "dimer", 4, &universal, 0, "b");
    line(&mut sam, "n1", 0, &nextera, 0, "a");
    line(&mut sam, "n2", 16, &nextera, 0, "b");
    line(&mut sam, "short", 0, &background(8), 0, "a");
    for name in ["c1", "c2", "c3"] {
        line(&mut sam, name, 0, &background(LENGTH), 0, "b");
    }
    line(&mut sam, "secondary", 0x100, &universal, 0, "a");
    line(&mut sam, "supplementary", 0x800, &universal, 0, "a");
    sam
}

fn run_report(sam: &str, args: &[&str]) -> Result<Report> {
    static RUNS: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-adapter-content-{}-{}",
        std::process::id(),
        RUNS.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("in.sam");
    fs::write(&path, sam).expect("write");
    let mut all = vec![format!("I={}", path.display())];
    all.extend(args.iter().map(|a| (*a).to_owned()));
    let report = TOOL
        .parse(&all)
        .and_then(|inv| Options::from_invocation(&inv))
        .and_then(|opts| {
            pass::run(Path::new(&path), 1, |header| {
                Ok(AdapterContent::new(&opts, header))
            })
        });
    let _ = fs::remove_dir_all(&dir);
    report
}

fn find<'r>(report: &'r Report, adapter: &str, end: ReadEnd, sample: Option<&str>) -> &'r Row {
    report
        .rows
        .iter()
        .find(|r| r.adapter == adapter && r.read_end == end && r.labels.sample.as_deref() == sample)
        .unwrap_or_else(|| panic!("no row {adapter} {end:?}"))
}

fn series<'r>(report: &'r Report, label: &str) -> &'r Histogram {
    let table = report.table.as_ref().expect("table");
    &table
        .series
        .iter()
        .find(|s| s.label == label)
        .unwrap_or_else(|| panic!("no series {label}"))
        .histogram
}

#[test]
fn planted_unpaired_reads_give_the_planted_percentages() {
    let report = run_report(&unpaired(), &[]).expect("runs");
    // FastQC's six adapters, one read end.
    assert_eq!(report.rows.len(), 6);
    let universal = find(
        &report,
        "Illumina Universal Adapter",
        ReadEnd::Unpaired,
        None,
    );
    assert_eq!(universal.kmer, UNIVERSAL);
    assert_eq!((universal.reads, universal.reads_with_adapter), (10, 4));
    assert_eq!(universal.pct_reads_with_adapter, 40.0);
    assert_eq!(universal.mean_adapter_start, 11.0);
    assert_eq!(universal.pct_at_last_cycle, 40.0);
    assert!(universal.pair_read_through.is_nan());
    let nextera = find(
        &report,
        "Nextera Transposase Sequence",
        ReadEnd::Unpaired,
        None,
    );
    assert_eq!(nextera.reads_with_adapter, 2);
    assert_eq!(nextera.pct_reads_with_adapter, 20.0);
    assert_eq!(nextera.mean_adapter_start, 21.0);
    let poly_a = find(&report, "PolyA", ReadEnd::Unpaired, None);
    assert_eq!(poly_a.reads_with_adapter, 0);
    assert!(poly_a.mean_adapter_start.is_nan());
    assert_eq!(poly_a.pct_reads_with_adapter, 0.0);
}

#[test]
fn the_table_is_cumulative_and_stops_where_a_key_can_last_start() {
    let report = run_report(&unpaired(), &[]).expect("runs");
    let table = report.table.as_ref().expect("table");
    assert_eq!(table.key_label, KEY_LABEL);
    assert_eq!(table.series.len(), 6);
    let universal = series(&report, "All_Reads.UNPAIRED.Illumina Universal Adapter");
    // 50 bases, a 12-base key: the last start is cycle 39.
    assert_eq!(universal.max(), Some(39));
    assert_eq!(universal.min(), Some(1));
    for cycle in 1..=10 {
        assert_eq!(universal.get(cycle), Some(0.0), "cycle {cycle}");
    }
    for cycle in 11..=39 {
        assert_eq!(universal.get(cycle), Some(40.0), "cycle {cycle}");
    }
    let nextera = series(&report, "All_Reads.UNPAIRED.Nextera Transposase Sequence");
    assert_eq!(nextera.get(20), Some(0.0));
    assert_eq!(nextera.get(21), Some(20.0));
}

#[test]
fn levels_split_the_reads_by_read_group() {
    let report = run_report(&unpaired(), &["LEVEL=null", "LEVEL=LIBRARY"]).expect("runs");
    assert_eq!(report.rows.len(), 12);
    let by_library = |library: &str| {
        report
            .rows
            .iter()
            .find(|r| {
                r.adapter == "Illumina Universal Adapter"
                    && r.labels.library.as_deref() == Some(library)
            })
            .expect("row")
    };
    // Library l1 (read group a): 5 reads, 3 Universal. l2: 5 reads, the dimer.
    assert_eq!(
        (by_library("l1").reads, by_library("l1").reads_with_adapter),
        (5, 3)
    );
    assert_eq!(by_library("l1").pct_reads_with_adapter, 60.0);
    assert_eq!(
        (by_library("l2").reads, by_library("l2").reads_with_adapter),
        (5, 1)
    );
    assert_eq!(by_library("l2").pct_reads_with_adapter, 20.0);
    assert_eq!(
        series(&report, "l2.UNPAIRED.Illumina Universal Adapter").get(11),
        Some(20.0)
    );
}

/// 4 pairs of 50-base reads: pair 1 has an insert of 30 (both reads run
/// into the adapter at cycle 31, the first read forward, the second
/// reverse), pair 2 an insert of 300, pair 3 is not proper, pair 4 has a
/// QC-failed and duplicate first read with the Universal key at cycle 1.
fn paired() -> String {
    let mut sam = String::from(HEADER);
    let through = with_adapter(UNIVERSAL, 30, LENGTH);
    line(&mut sam, "p1", 0x1 | 0x2 | 0x40 | 0x20, &through, 30, "a");
    line(&mut sam, "p1", 0x1 | 0x2 | 0x80 | 0x10, &through, -30, "a");
    let plain = background(LENGTH);
    line(&mut sam, "p2", 0x1 | 0x2 | 0x40 | 0x20, &plain, 300, "a");
    line(&mut sam, "p2", 0x1 | 0x2 | 0x80 | 0x10, &plain, -300, "a");
    line(&mut sam, "p3", 0x1 | 0x40, &plain, 20, "a");
    line(&mut sam, "p3", 0x1 | 0x80, &plain, -20, "a");
    let at_start = with_adapter(UNIVERSAL, 0, LENGTH);
    line(
        &mut sam,
        "p4",
        0x1 | 0x40 | 0x200 | 0x400,
        &at_start,
        0,
        "a",
    );
    line(&mut sam, "p4", 0x1 | 0x80, &plain, 0, "a");
    sam
}

#[test]
fn read_ends_are_apart_and_read_through_follows_tlen() {
    let report = run_report(&paired(), &[]).expect("runs");
    assert_eq!(report.rows.len(), 12, "FIRST and SECOND, no UNPAIRED");
    let first = find(&report, "Illumina Universal Adapter", ReadEnd::First, None);
    assert_eq!((first.reads, first.reads_with_adapter), (4, 2));
    assert_eq!(first.mean_adapter_start, 16.0);
    // Pairs 1 and 2 are proper with TLEN set; only pair 1 is under 50.
    assert_eq!(first.pair_read_through, 0.5);
    let second = find(&report, "Illumina Universal Adapter", ReadEnd::Second, None);
    assert_eq!((second.reads, second.reads_with_adapter), (4, 1));
    assert_eq!(second.mean_adapter_start, 31.0);
    assert_eq!(second.pair_read_through, 0.5);
    let curve = series(&report, "All_Reads.FIRST.Illumina Universal Adapter");
    assert_eq!(curve.get(1), Some(25.0));
    assert_eq!(curve.get(30), Some(25.0));
    assert_eq!(curve.get(31), Some(50.0));
}

#[test]
fn qc_failed_and_duplicate_reads_leave_on_request() {
    let first = |args: &[&str]| {
        let report = run_report(&paired(), args).expect("runs");
        let row = find(&report, "Illumina Universal Adapter", ReadEnd::First, None);
        (row.reads, row.reads_with_adapter)
    };
    assert_eq!(first(&["PF_READS_ONLY=true"]), (3, 1));
    assert_eq!(first(&["INCLUDE_DUPLICATES=false"]), (3, 1));
    assert_eq!(
        first(&["INCLUDE_DUPLICATES=true", "PF_READS_ONLY=false"]),
        (4, 2)
    );
}

#[test]
fn own_adapters_replace_the_default_set() {
    let report = run_report(
        &unpaired(),
        &[
            "ADAPTER_SEQUENCE=mine=agatcgga",
            "ADAPTER_SEQUENCE=CTGTCTCT",
            "ADAPTER_KMER=8",
        ],
    )
    .expect("runs");
    let names: Vec<&str> = report.rows.iter().map(|r| r.adapter.as_str()).collect();
    assert_eq!(names, ["mine", "CTGTCTCT"]);
    assert_eq!(report.rows[0].kmer, "AGATCGGA");
    assert_eq!(report.rows[0].reads_with_adapter, 4);
    assert_eq!(report.rows[1].reads_with_adapter, 2);
    // A 50-base read and an 8-base key: the last start is cycle 43.
    assert_eq!(series(&report, "All_Reads.UNPAIRED.mine").max(), Some(43));
}

#[test]
fn a_file_with_no_record_writes_unpaired_rows_and_no_table() {
    let report = run_report(HEADER, &[]).expect("runs");
    assert_eq!(report.rows.len(), 6);
    assert!(
        report
            .rows
            .iter()
            .all(|r| r.read_end == ReadEnd::Unpaired && r.reads == 0)
    );
    assert!(report.rows[0].pct_reads_with_adapter.is_nan());
    assert!(report.table.is_none());
}

#[test]
fn the_table_stops_at_its_cap_and_the_summary_does_not() {
    let mut sam = String::from(HEADER);
    let long = with_adapter(UNIVERSAL, 1_050, 1_200);
    line(&mut sam, "long", 4, &long, 0, "a");
    line(
        &mut sam,
        "early",
        4,
        &with_adapter(UNIVERSAL, 5, 1_200),
        0,
        "a",
    );
    let report = run_report(&sam, &[]).expect("runs");
    let row = find(
        &report,
        "Illumina Universal Adapter",
        ReadEnd::Unpaired,
        None,
    );
    assert_eq!(row.reads_with_adapter, 2);
    assert_eq!(row.pct_reads_with_adapter, 100.0);
    assert_eq!(row.pct_at_last_cycle, 50.0);
    assert_eq!(row.mean_adapter_start, (1_051.0 + 6.0) / 2.0);
    let curve = series(&report, "All_Reads.UNPAIRED.Illumina Universal Adapter");
    assert_eq!(curve.max(), Some(MAX_TABLE_CYCLES as i64));
}

#[test]
fn the_matcher_finds_the_first_exact_occurrence() {
    let adapters = vec![
        parse_adapter("A=AAAA", 4).expect("A"),
        parse_adapter("B=ACGT", 4).expect("B"),
        parse_adapter("C=ACGT", 4).expect("C"),
    ];
    let matcher = Matcher::new(&adapters, 4);
    let mut out = Vec::new();
    matcher.first_starts(b"TTAAAAAAACGTACGT", &mut out);
    assert_eq!(out, [(0, 2), (1, 8), (2, 8)]);
    // An N breaks a run: no AAAA across it.
    matcher.first_starts(b"AANAAACGNT", &mut out);
    assert_eq!(out, []);
    // One mismatch is no match.
    matcher.first_starts(b"AAAT", &mut out);
    assert_eq!(out, []);
    let long = "ACGT".repeat(8);
    let matcher = Matcher::new(&[parse_adapter(&long, 32).expect("32")], 32);
    matcher.first_starts(format!("GG{long}").as_bytes(), &mut out);
    assert_eq!(out, [(0, 2)]);
}

#[test]
fn bad_options_are_usage_errors() {
    let error = |args: &[&str]| {
        let mut all = vec!["I=x.bam".to_owned()];
        all.extend(args.iter().map(|a| (*a).to_owned()));
        let inv = TOOL.parse(&all).expect("parses");
        Options::from_invocation(&inv)
            .expect_err("refused")
            .to_string()
    };
    assert!(error(&["ADAPTER_KMER=0"]).contains("1 to 32"));
    assert!(error(&["ADAPTER_KMER=33"]).contains("1 to 32"));
    assert!(error(&["ADAPTER_KMER=13"]).contains("FastQC's adapters have 12 bases"));
    assert!(
        error(&["ADAPTER_KMER=13", "ADAPTER_SEQUENCE=x=ACGTACGTACGT"])
            .contains("12 bases, ADAPTER_KMER is 13")
    );
    assert!(error(&["ADAPTER_SEQUENCE=x=ACGTNACGTACGT"]).contains("must be A, C, G or T"));
    assert!(error(&["ADAPTER_SEQUENCE==ACGTACGTACGT"]).contains("the name is empty"));
    assert!(
        error(&[
            "ADAPTER_SEQUENCE=x=ACGTACGTACGT",
            "ADAPTER_SEQUENCE=x=TTTTTTTTTTTT"
        ])
        .contains("given twice")
    );
    let many: Vec<String> = (0..=MAX_ADAPTERS)
        .map(|i| format!("ADAPTER_SEQUENCE=a{i}=ACGTACGTACGT"))
        .collect();
    let many: Vec<&str> = many.iter().map(String::as_str).collect();
    assert!(error(&many).contains("at most 64"));
    assert!(error(&["INCLUDE_DUPLICATES=maybe"]).contains("true or false"));
}

#[test]
fn rows_and_columns_agree() {
    let report = run_report(&unpaired(), &[]).expect("runs");
    for row in &report.rows {
        assert_eq!(row.values().len(), Row::COLUMNS.len());
    }
    let output = MetricsOutput::from_rows(&report.rows, report.table);
    let mut text = Vec::new();
    output.write(Vec::new(), &mut text).expect("write");
    let text = String::from_utf8(text).expect("utf8");
    assert!(text.contains("## METRICS CLASS\tbammetrics.qc.AdapterContentMetrics"));
    assert!(text.contains("## HISTOGRAM\tjava.lang.Integer"));
    assert!(text.contains("READ_CYCLE\tAll_Reads.UNPAIRED.Illumina Universal Adapter"));
}

/// BAM's 4-bit packing of `bases`, two to a byte, the first in the high
/// half.
fn pack(bases: &[u8]) -> Vec<u8> {
    const CODES: &[u8] = b"=ACMGRSVTWYHKDBN";
    let code = |b: &u8| CODES.iter().position(|c| c == b).expect("base") as u8;
    bases
        .chunks(2)
        .map(|pair| (code(&pair[0]) << 4) | pair.get(1).map_or(0, code))
        .collect()
}

#[test]
fn packed_records_scan_as_their_bases_in_sequencing_order() {
    let matcher = Matcher::new(
        &Options::from_invocation(&TOOL.parse(&["I=x".to_owned()]).expect("parse"))
            .expect("opts")
            .adapters,
        DEFAULT_KMER,
    );
    let mut seed = 7u64;
    let (mut packed_hits, mut text_hits) = (Vec::new(), Vec::new());
    for round in 0..500 {
        let length = 20 + round % 61;
        let mut read: Vec<u8> = (0..length)
            .map(|_| {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                b"ACGTN"[(seed >> 33) as usize % 5]
            })
            .collect();
        let key = DEFAULT_ADAPTERS[round % 6].1.as_bytes();
        let at = round % (length - 11);
        read[at..at + 12].copy_from_slice(key);
        // The record stores the reverse complement on the reverse strand.
        for reverse in [false, true] {
            let stored = if reverse {
                reverse_complement(&read)
            } else {
                read.clone()
            };
            matcher.scan_record(&pack(&stored), stored.len(), reverse, &mut packed_hits);
            matcher.first_starts(&read, &mut text_hits);
            assert_eq!(packed_hits, text_hits, "round {round} reverse {reverse}");
            assert!(text_hits.iter().any(|&(a, s)| a == round % 6 && s <= at));
        }
    }
}

#[test]
fn a_cut_record_reads_its_missing_bases_as_n() {
    let matcher = Matcher::new(&[parse_adapter("A=AGATCGGAAGAG", 12).expect("A")], 12);
    let read = b"AGATCGGAAGAGAGATCGGAAGAG";
    let packed = pack(read);
    let mut hits = Vec::new();
    // 24 bases need 12 bytes; with 8 the last 8 bases are N.
    matcher.scan_record(&packed[..8], read.len(), false, &mut hits);
    assert_eq!(hits, [(0, 0)]);
    matcher.scan_record(&packed[..8], read.len(), true, &mut hits);
    assert_eq!(hits, []);
    matcher.scan_record(&[], read.len(), true, &mut hits);
    assert_eq!(hits, []);
}
