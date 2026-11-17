use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use super::*;

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:100000
@RG\tID:a\tSM:s1\tLB:l1\tPU:u1
@RG\tID:b\tSM:s1\tLB:l2\tPU:u2
";

/// A directory of this test's own.
fn dir(test: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("bammetrics-telomere-{}-{test}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("dir");
    dir
}

/// A 60-base read: `repeats` copies of TTAGGG, then filler with no G or C.
fn read(repeats: usize) -> String {
    let mut seq = "TTAGGG".repeat(repeats);
    seq.push_str(&"A".repeat(60 - seq.len()));
    seq
}

/// One SAM line; an unmapped read when `flag` has 0x4.
fn line(out: &mut String, name: &str, flag: u16, seq: &str, rg: &str) {
    let (contig, pos, mapq, cigar) = if flag & 4 != 0 {
        ("*", 0, 0, "*".to_owned())
    } else {
        ("chr1", 100, 60, format!("{}M", seq.len()))
    };
    writeln!(
        out,
        "{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t*\tRG:Z:{rg}"
    )
    .expect("write");
}

fn run(sam: &str, args: &[&str]) -> Result<Vec<Bucket>> {
    static RUNS: AtomicUsize = AtomicUsize::new(0);
    let dir = dir(&format!("run{}", RUNS.fetch_add(1, Ordering::Relaxed)));
    let path = dir.join("in.sam");
    fs::write(&path, sam).expect("write");
    let mut all = vec![format!("I={}", path.display())];
    all.extend(args.iter().map(|a| (*a).to_owned()));
    let opts = Options::from_invocation(&TOOL.parse(&all)?)?;
    let buckets = pass::run(Path::new(&path), 1, |header| {
        Ok(Telomere::new(&opts, header))
    });
    let _ = fs::remove_dir_all(&dir);
    buckets
}

fn rows(sam: &str, args: &[&str]) -> Vec<Row> {
    run(sam, args)
        .expect("runs")
        .into_iter()
        .map(|b| b.row)
        .collect()
}

#[test]
fn occurrences_do_not_overlap_and_take_the_larger_strand() {
    assert_eq!(count_occurrences(b"TTAGGGTTAGGG", b"TTAGGG"), 2);
    assert_eq!(count_occurrences(b"TTAGGTTAGGG", b"TTAGGG"), 1);
    assert_eq!(count_occurrences(b"AAAAAAA", b"AAA"), 2, "non-overlapping");
    assert_eq!(count_occurrences(b"AA", b"AAA"), 0);
    assert_eq!(count_occurrences(b"", b"TTAGGG"), 0);
    assert_eq!(count_occurrences(b"ACGT", b""), 0);
    let read = b"CCCTAACCCTAACCCTAATTAGGG";
    assert_eq!(count_repeats(read, b"TTAGGG", b"CCCTAA"), 3);
    assert_eq!(reverse_complement(b"TTAGGG"), b"CCCTAA");
    assert_eq!(reverse_complement(b"ANC="), b"NGNT");
    // An N inside a run breaks one copy.
    assert_eq!(
        count_repeats(b"TTAGGGTTANGGTTAGGG", b"TTAGGG", b"CCCTAA"),
        2
    );
}

/// The one-pass count of repeats up to 8 bases equals the two greedy scans
/// on random reads, self-overlapping repeats included.
#[test]
fn the_one_pass_count_equals_two_scans() {
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
    let mut next = |n: u64| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) % n
    };
    let repeats: [&[u8]; 7] = [
        b"TTAGGG",
        b"AAA",
        b"ACAC",
        b"AAAAAAAA",
        b"TTTAGGGTT",
        b"A",
        b"GG",
    ];
    for _ in 0..4000 {
        let length = next(160) as usize;
        // A small alphabet makes runs and overlaps common.
        let alphabet: &[u8] = if next(2) == 0 { b"ACGTN" } else { b"AG" };
        let read: Vec<u8> = (0..length)
            .map(|_| alphabet[next(alphabet.len() as u64) as usize])
            .collect();
        for repeat in repeats {
            let reverse = reverse_complement(repeat);
            let slow = count_occurrences(&read, repeat).max(count_occurrences(&read, &reverse));
            assert_eq!(
                count_repeats(&read, repeat, &reverse),
                slow,
                "{} in {}",
                String::from_utf8_lossy(repeat),
                String::from_utf8_lossy(&read)
            );
        }
    }
}

#[test]
fn the_gc_bin_is_open_at_48_and_closed_at_52() {
    assert!(!in_gc_bin(48, 100));
    assert!(in_gc_bin(49, 100));
    assert!(in_gc_bin(52, 100));
    assert!(!in_gc_bin(53, 100));
    assert!(!in_gc_bin(0, 0));
    // 150-base reads: 72 is 48 percent, 78 is 52 percent.
    assert!(!in_gc_bin(72, 150));
    assert!(in_gc_bin(73, 150));
    assert!(in_gc_bin(78, 150));
    assert!(!in_gc_bin(79, 150));
}

/// TelSeq's arithmetic: `floor((gc - 0.4f) / 0.02f)` in bins 4 and 5 keeps
/// the fractions the integer test keeps, for every length up to 2,000.
#[test]
fn the_gc_bin_equals_telseq_float_bins() {
    let (lower, size) = (f64::from(0.4f32), f64::from(0.02f32));
    for length in 1..=2000usize {
        for gc in 0..=length {
            let fraction = gc as f64 / length as f64;
            let bin = ((fraction - lower) / size).floor();
            let telseq = bin == 4.0 || bin == 5.0;
            assert_eq!(in_gc_bin(gc, length), telseq, "{gc}/{length}");
        }
    }
}

#[test]
fn the_estimate_is_telseq_formula() {
    let formula = Formula::default();
    let kb = estimate_kb(30, 1000, formula).expect("estimate");
    assert_eq!(kb, 30.0 / 1000.0 * 332_720_800.0 / 1000.0 / 46.0);
    assert_eq!(estimate_kb(30, 0, formula), None);
    let other = Formula {
        telomere_ends: 2,
        genome_length_at_gc: 1000,
    };
    assert_eq!(estimate_kb(1, 1, other), Some(0.5));
}

#[test]
fn the_estimate_is_empty_when_lengths_vary() {
    let mut counts = Counts {
        total_reads: 100,
        telomeric_reads: 10,
        reads_gc_48_52: 20,
        ..Counts::default()
    };
    counts.lengths.increment_by(100, 95.0);
    counts.lengths.increment_by(89, 5.0);
    let five = row(&counts, Formula::default());
    assert_eq!(five.read_length, Some(100));
    assert_eq!(five.pct_off_length, 0.05);
    assert!(five.estimated_length_kb.is_some(), "5 percent off is kept");
    counts.lengths = Histogram::new();
    counts.lengths.increment_by(100, 94.0);
    counts.lengths.increment_by(89, 6.0);
    let six = row(&counts, Formula::default());
    assert_eq!(six.estimated_length_kb, None, "6 percent off is empty");
    counts.lengths = Histogram::new();
    counts.lengths.increment_by(100, 50.0);
    counts.lengths.increment_by(90, 50.0);
    let near = row(&counts, Formula::default());
    assert_eq!(near.read_length, Some(90), "the shorter mode on a tie");
    assert_eq!(near.pct_off_length, 0.5, "100 is 11 percent over 90");
    let empty = row(&Counts::default(), Formula::default());
    assert!(empty.pct_telomeric.is_nan());
    assert_eq!(empty.read_length, None);
    assert_eq!(empty.estimated_length_kb, None);
    assert_eq!(empty.values().len(), Row::COLUMNS.len());
}

/// The planted file: reads at 0, 6, 7 and 20 repeats on both strands,
/// unmapped telomeric reads, an N in a run, and records that must not count.
fn planted() -> String {
    let mut sam = String::from(HEADER);
    let mut n = 0;
    let mut add = |sam: &mut String, flag: u16, seq: &str, rg: &str| {
        n += 1;
        line(sam, &format!("r{n}"), flag, seq, rg);
    };
    for (repeats, copies) in [(0, 4), (6, 3), (7, 2)] {
        for k in 0..copies {
            add(
                &mut sam,
                if k % 2 == 1 { 16 } else { 0 },
                &read(repeats),
                "a",
            );
        }
    }
    // 20 repeats fill a 120-base read; stored reverse complemented on 16.
    let twenty = "TTAGGG".repeat(20);
    add(&mut sam, 0, &twenty, "b");
    add(&mut sam, 16, &"CCCTAA".repeat(20), "b");
    add(&mut sam, 4, &twenty, "b");
    // 7 copies with an N in one of them: 6 left.
    let broken = read(7).replacen("TTAGGG", "TTANGG", 1);
    add(&mut sam, 4, &broken, "b");
    // Must not count: duplicate, secondary, supplementary, QC fail.
    for flag in [0x400, 0x100, 0x800, 0x200] {
        add(&mut sam, flag, &twenty, "a");
    }
    sam
}

#[test]
fn planted_reads_give_planted_counts() {
    let rows = rows(&planted(), &[]);
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(row.total_reads, 13);
    assert_eq!(row.telomeric_reads, 5, "two at 7, three at 20");
    assert_eq!(row.pct_telomeric, 5.0 / 13.0);
    assert_eq!(
        row.reads_gc_48_52, 3,
        "only the 20-repeat reads hold 50 percent GC"
    );
    assert_eq!(row.read_length, Some(60));
    assert_eq!(row.pct_off_length, 3.0 / 13.0);
    assert_eq!(row.estimated_length_kb, None, "3 of 13 reads are 120 bases");
    let duplicates = rows_of(&planted(), &["COUNT_DUPLICATES=true"]);
    assert_eq!(
        (duplicates.total_reads, duplicates.telomeric_reads),
        (14, 6)
    );
    let k20 = rows_of(&planted(), &["K=20"]);
    assert_eq!(k20.telomeric_reads, 3);
    let k1 = rows_of(&planted(), &["TELOMERE_REPEAT_THRESHOLD=1"]);
    assert_eq!(k1.telomeric_reads, 9);
}

fn rows_of(sam: &str, args: &[&str]) -> Row {
    rows(sam, args).remove(0)
}

#[test]
fn the_histogram_counts_reads_by_repeat_count_up_to_the_read_length() {
    let buckets = run(&planted(), &[]).expect("runs");
    let repeats = &buckets[0].repeats;
    assert_eq!(buckets[0].label, "All_Reads.read_count");
    let got: Vec<(i64, f64)> = repeats
        .keys()
        .map(|k| (k, repeats.get(k).unwrap_or(0.0)))
        .collect();
    let mut want: Vec<(i64, f64)> = (0..=10).map(|k| (k, 0.0)).collect();
    want[0].1 = 4.0;
    want[6].1 = 4.0;
    want[7].1 = 2.0;
    want.push((20, 3.0));
    assert_eq!(got, want);
    let mut text = Vec::new();
    output(buckets).write(Vec::new(), &mut text).expect("write");
    let text = String::from_utf8(text).expect("utf8");
    assert!(
        text.contains(
            "## HISTOGRAM\tjava.lang.Integer\nrepeat_count\tAll_Reads.read_count\n0\t4\n1\t0\n"
        ),
        "{text}"
    );
    assert!(
        text.ends_with("7\t2\n8\t0\n9\t0\n10\t0\n20\t3\n\n"),
        "{text}"
    );
}

#[test]
fn every_level_gets_its_rows() {
    let rows = rows(&planted(), &["LEVEL=SAMPLE", "LEVEL=READ_GROUP"]);
    let shape: Vec<(Option<&str>, Option<&str>, i64, i64)> = rows
        .iter()
        .map(|r| {
            (
                r.labels.sample.as_deref(),
                r.labels.read_group.as_deref(),
                r.total_reads,
                r.telomeric_reads,
            )
        })
        .collect();
    assert_eq!(
        shape,
        [
            (None, None, 13, 5),
            (Some("s1"), None, 13, 5),
            (Some("s1"), Some("u1"), 9, 2),
            (Some("s1"), Some("u2"), 4, 3),
        ]
    );
    // Read group u2 holds three 120-base reads and one of 60.
    assert_eq!(rows[3].read_length, Some(120));
    assert_eq!(rows[3].pct_off_length, 0.25);
    assert_eq!(rows[3].reads_gc_48_52, 3);
    // Read group u1: nine 60-base reads, no read in the GC bin.
    assert_eq!(rows[2].estimated_length_kb, None);
}

#[test]
fn a_uniform_file_gets_the_estimate() {
    let mut sam = String::from(HEADER);
    // 100-base reads: 2 telomeric, 4 of exactly 50 percent GC.
    let gc50 = "GC".repeat(25) + &"AT".repeat(25);
    let tel = "TTAGGG".repeat(16) + "AAAA";
    for (i, seq) in [&gc50, &gc50, &gc50, &gc50, &tel, &tel]
        .into_iter()
        .enumerate()
    {
        line(&mut sam, &format!("u{i}"), 4, seq, "a");
    }
    let row = rows_of(&sam, &[]);
    assert_eq!(row.reads_gc_48_52, 4);
    assert_eq!(row.telomeric_reads, 2);
    assert_eq!(
        row.estimated_length_kb,
        Some(2.0 / 4.0 * 332_720_800.0 / 1000.0 / 46.0)
    );
    let ends = rows_of(&sam, &["TELOMERE_ENDS=92", "GENOME_LENGTH_AT_GC=1000"]);
    assert_eq!(
        ends.estimated_length_kb,
        Some(2.0 / 4.0 * 1000.0 / 1000.0 / 92.0)
    );
}

#[test]
fn another_repeat_counts_its_own_motif() {
    let mut sam = String::from(HEADER);
    // Arabidopsis TTTAGGG: 8 copies, and its reverse complement CCCTAAA.
    line(&mut sam, "p1", 4, &"TTTAGGG".repeat(8), "a");
    line(&mut sam, "p2", 4, &"CCCTAAA".repeat(8), "a");
    line(&mut sam, "h1", 4, &"TTAGGG".repeat(8), "a");
    let plant = rows_of(&sam, &["REPEAT=tttaggg"]);
    assert_eq!(plant.telomeric_reads, 2);
    // TTTAGGG holds TTAGGG, so all three reads are telomeric for TTAGGG.
    let human = rows_of(&sam, &[]);
    assert_eq!(human.telomeric_reads, 3);
    let none = rows_of(&sam, &["REPEAT=TTAGGC"]);
    assert_eq!(none.telomeric_reads, 0);
}

/// A 3,000-base read with a one-base repeat fills the histogram to key
/// 1,000 and prints its own count above.
#[test]
fn long_reads_fill_the_histogram_to_the_cap() {
    let mut sam = String::from(HEADER);
    line(&mut sam, "long", 4, &"A".repeat(3000), "a");
    let buckets = run(&sam, &["REPEAT=A"]).expect("runs");
    let repeats = &buckets[0].repeats;
    let keys: Vec<i64> = repeats.keys().collect();
    assert_eq!(keys.len(), 1002);
    assert_eq!(keys[1000], MAX_ZERO_FILLED_KEY);
    assert_eq!(repeats.get(3000), Some(1.0));
}

#[test]
fn options_reject_bad_values() {
    let message = |args: &[&str]| run(&planted(), args).expect_err("refused").to_string();
    assert!(message(&["REPEAT=TTAGGN"]).contains("A, C, G and T"));
    assert!(message(&["REPEAT="]).contains("A, C, G and T"));
    let long = format!("REPEAT={}", "A".repeat(MAX_REPEAT_LENGTH + 1));
    assert!(message(&[&long]).contains("1 to 100 bases"));
    assert!(message(&["K=0"]).contains("TELOMERE_REPEAT_THRESHOLD=0: must be at least 1"));
    assert!(message(&["TELOMERE_ENDS=-1"]).contains("must be at least 1"));
    assert!(message(&["GENOME_LENGTH_AT_GC=0"]).contains("must be at least 1"));
    assert!(message(&["COUNT_DUPLICATES=maybe"]).contains("expected true or false"));
    assert!(message(&["LEVEL=NOPE"]).contains("NOPE"));
}
