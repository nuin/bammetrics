//! `CollectExonCoverageMetrics` on planted reads: a SAM file whose depth at
//! every base is known by construction (reads the Hs filters drop, a
//! deletion, an overlapping pair), and a GTF and a BED12 file whose exons
//! are known. Each cell of each row is computed here from the planted
//! intervals alone. The Hs agreement on real data is
//! `tests/concordance/test_exon_coverage.py`.

use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use bammetrics_core::metric::{MetricsOutput, ToolOutput, Value};
use bammetrics_core::registry;

fn work_dir() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-exon-coverage-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("work dir");
    dir
}

/// One SAM line: unpaired unless `mate` is given as (position, TLEN).
fn read(
    name: &str,
    flag: u16,
    (contig, pos): (&str, i64),
    mapq: u8,
    cigar: &str,
    qual: char,
    mate: Option<(i64, i64)>,
) -> String {
    let bases: i64 = cigar
        .split(|c: char| c.is_ascii_alphabetic())
        .zip(cigar.matches(|c: char| c.is_ascii_alphabetic()))
        .filter(|(_, op)| matches!(*op, "M" | "I" | "S" | "=" | "X"))
        .map(|(n, _)| n.parse::<i64>().expect("length"))
        .sum();
    let length = usize::try_from(bases).expect("length");
    let (rnext, pnext, tlen) = match mate {
        Some((p, t)) => ("=", p, t),
        None => ("*", 0, 0),
    };
    format!(
        "{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{}\t{}\n",
        "A".repeat(length),
        qual.to_string().repeat(length)
    )
}

/// The reads, and the intervals that each add one to the depth.
fn planted() -> (String, Vec<(&'static str, i64, i64)>) {
    let mut sam =
        String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n@SQ\tSN:chr2\tLN:5000\n");
    sam += &read("f1", 0, ("chr1", 101), 60, "100M", 'I', None);
    sam += &read("dup", 1024, ("chr1", 101), 60, "100M", 'I', None);
    sam += &read("lowmq", 0, ("chr1", 101), 10, "100M", 'I', None);
    sam += &read("lowq", 0, ("chr1", 101), 60, "100M", '#', None);
    sam += &read("second", 256, ("chr1", 101), 60, "100M", 'I', None);
    sam += &read("f2", 16, ("chr1", 151), 60, "100M", 'I', None);
    sam += &read("f3", 0, ("chr1", 151), 60, "50M", 'I', None);
    sam += &read("del", 0, ("chr1", 301), 60, "20M10D20M", 'I', None);
    sam += &read("qcfail", 512, ("chr1", 301), 60, "50M", 'I', None);
    sam += &read("p", 99, ("chr1", 501), 60, "50M", 'I', Some((531, 80)));
    sam += &read("p", 147, ("chr1", 531), 60, "50M", 'I', Some((501, -80)));
    for name in ["g1", "g2", "g3"] {
        sam += &read(name, 16, ("chr2", 1001), 60, "100M", 'I', None);
    }
    sam += &read("g4", 0, ("chr2", 1051), 60, "100M", 'I', None);
    let counted = vec![
        ("chr1", 101, 200),
        ("chr1", 151, 250),
        ("chr1", 151, 200),
        ("chr1", 301, 320),
        ("chr1", 331, 350),
        ("chr1", 501, 580),
        ("chr2", 1001, 1100),
        ("chr2", 1001, 1100),
        ("chr2", 1001, 1100),
        ("chr2", 1051, 1150),
    ];
    (sam, counted)
}

fn gtf_line(contig: &str, start: i64, end: i64, strand: &str, attrs: &str) -> String {
    format!("{contig}\tplanted\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
}

/// GA (plus strand, chr1): T1 with a one-base exon on the deletion, T2
/// overlapping T1. GB (minus strand, single exon) on `2`, which the BAM
/// calls `chr2`. GC on `chrUn`, which the BAM lacks.
fn gtf() -> String {
    let ga = |t: &str, n: i64| {
        format!("gene_id \"ENSGA\"; gene_name \"GA\"; transcript_id \"{t}\"; exon_number \"{n}\";")
    };
    let mut text = String::new();
    text += &gtf_line("chr1", 101, 200, "+", &ga("T1", 1));
    text += &gtf_line("chr1", 321, 321, "+", &ga("T1", 2));
    text += &gtf_line("chr1", 151, 260, "+", &ga("T2", 1));
    text += &gtf_line("chr1", 501, 600, "+", &ga("T2", 2));
    text += &gtf_line(
        "2",
        1001,
        1200,
        "-",
        "gene_id \"ENSGB\"; gene_name \"GB\"; transcript_id \"T3\";",
    );
    text += &gtf_line(
        "chrUn",
        1,
        10,
        "+",
        "gene_id \"ENSGC\"; gene_name \"GC\"; transcript_id \"T4\";",
    );
    text
}

/// The expected cells of a set of bases: mean, median, min, max, then the
/// bases at or above 1, 2 and 3 and their fractions.
fn expected(counted: &[(&str, i64, i64)], contig: &str, spans: &[(i64, i64)]) -> Vec<f64> {
    let mut depths: Vec<i64> = spans
        .iter()
        .flat_map(|&(s, e)| s..=e)
        .map(|pos| {
            counted
                .iter()
                .filter(|(c, s, e)| *c == contig && (*s..=*e).contains(&pos))
                .count() as i64
        })
        .collect();
    depths.sort_unstable();
    let n = depths.len();
    let median = if n % 2 == 1 {
        depths[n / 2] as f64
    } else {
        (depths[n / 2 - 1] + depths[n / 2]) as f64 / 2.0
    };
    let mut cells = vec![
        depths.iter().sum::<i64>() as f64 / n as f64,
        median,
        depths[0] as f64,
        depths[n - 1] as f64,
    ];
    let at: Vec<f64> = [1, 2, 3]
        .iter()
        .map(|&m| depths.iter().filter(|&&d| d >= m).count() as f64)
        .collect();
    cells.extend(&at);
    cells.extend(at.iter().map(|a| a / n as f64));
    cells
}

/// (LEVEL, TRANSCRIPT_ID, CHROM, the exon union, EXON_COUNT) of a summary row.
type Summary<'a> = (&'a str, &'a str, &'a str, &'a [(i64, i64)], i64);

fn rows(output: &MetricsOutput) -> Vec<BTreeMap<&str, &Value>> {
    output
        .rows
        .iter()
        .map(|row| output.columns.iter().copied().zip(row).collect())
        .collect()
}

fn depth_cells(row: &BTreeMap<&str, &Value>) -> Vec<f64> {
    [
        "MEAN_DEPTH",
        "MEDIAN_DEPTH",
        "MINIMUM_DEPTH",
        "MAXIMUM_DEPTH",
        "BASES_1X",
        "BASES_2X",
        "BASES_3X",
        "PCT_1X",
        "PCT_2X",
        "PCT_3X",
    ]
    .iter()
    .map(|c| row[c].as_float().unwrap_or_else(|| panic!("{c}")))
    .collect()
}

fn text(row: &BTreeMap<&str, &Value>, column: &str) -> String {
    match row[column] {
        Value::Str(s) => s.clone(),
        Value::Int(i) => i.to_string(),
        Value::Null => String::new(),
        other => panic!("{column}: {other:?}"),
    }
}

fn run(annotation: &str, extra: &[&str]) -> ToolOutput {
    let dir = work_dir();
    let (sam, _) = planted();
    // One file per call: the tests run in parallel in one directory.
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let call = CALLS.fetch_add(1, Ordering::Relaxed);
    let input = dir.join(format!("planted{call}.sam"));
    fs::write(&input, sam).expect("sam");
    let mut args = vec![
        format!("I={}", input.display()),
        format!("ANNOTATION={annotation}"),
        "MIN_DEPTH=1".to_owned(),
        "MIN_DEPTH=2".to_owned(),
        "MIN_DEPTH=3".to_owned(),
        format!(
            "SUMMARY_OUTPUT={}",
            dir.join(format!("summary{call}.txt")).display()
        ),
    ];
    args.extend(extra.iter().map(|s| (*s).to_owned()));
    registry::find("CollectExonCoverageMetrics")
        .expect("registered")
        .execute(&args)
        .unwrap_or_else(|e| panic!("{e}"))
        .1
}

#[test]
fn every_exon_row_is_the_planted_depth() {
    let dir = work_dir();
    let annotation = dir.join("planted.gtf");
    fs::write(&annotation, gtf()).expect("gtf");
    let output = run(&annotation.display().to_string(), &[]);
    let (_, counted) = planted();
    let exons = output.main.expect("exon file");
    let rows = rows(&exons);
    let want: [(&str, &str, &str, &str, i64, i64, &str); 6] = [
        ("GA", "T1", "1", "chr1", 101, 200, "+"),
        ("GA", "T1", "2", "chr1", 321, 321, "+"),
        ("GA", "T2", "1", "chr1", 151, 260, "+"),
        ("GA", "T2", "2", "chr1", 501, 600, "+"),
        ("GB", "T3", "1", "chr2", 1001, 1200, "-"),
        ("GC", "T4", "1", "chrUn", 1, 10, "+"),
    ];
    assert_eq!(rows.len(), want.len());
    for (row, (gene, transcript, number, chrom, start, end, strand)) in rows.iter().zip(want) {
        assert_eq!(text(row, "GENE_NAME"), gene);
        assert_eq!(text(row, "TRANSCRIPT_ID"), transcript);
        assert_eq!(text(row, "EXON_NUMBER"), number);
        assert_eq!(text(row, "CHROM"), chrom);
        assert_eq!(text(row, "START"), start.to_string());
        assert_eq!(text(row, "END"), end.to_string());
        assert_eq!(text(row, "STRAND"), strand);
        assert_eq!(text(row, "LENGTH"), (end - start + 1).to_string());
        if chrom == "chrUn" {
            assert_eq!(
                row["MEAN_DEPTH"],
                &Value::Null,
                "not on the BAM: not counted"
            );
            continue;
        }
        assert_eq!(
            depth_cells(row),
            expected(&counted, chrom, &[(start, end)]),
            "{transcript} exon {number}"
        );
    }
    // The planted values themselves, so a wrong model above cannot pass.
    assert_eq!(depth_cells(&rows[0])[..4], [2.0, 2.0, 1.0, 3.0]);
    assert_eq!(
        depth_cells(&rows[1])[..4],
        [0.0, 0.0, 0.0, 0.0],
        "a deleted base"
    );
    assert_eq!(
        depth_cells(&rows[3])[..4],
        [0.8, 1.0, 0.0, 1.0],
        "the pair's overlap once"
    );
    assert_eq!(depth_cells(&rows[4])[..4], [2.0, 2.0, 0.0, 4.0]);
}

#[test]
fn gene_and_transcript_rows_are_the_union_of_their_exons() {
    let dir = work_dir();
    let annotation = dir.join("summary.gtf");
    fs::write(&annotation, gtf()).expect("gtf");
    let output = run(&annotation.display().to_string(), &[]);
    let (_, counted) = planted();
    let (path, summary) = &output.extra[0];
    assert!(path.display().to_string().contains("summary"));
    let rows = rows(summary);
    let want: [Summary; 7] = [
        ("GENE", "", "chr1", &[(101, 260), (321, 321), (501, 600)], 4),
        ("TRANSCRIPT", "T1", "chr1", &[(101, 200), (321, 321)], 2),
        ("TRANSCRIPT", "T2", "chr1", &[(151, 260), (501, 600)], 2),
        ("GENE", "", "chr2", &[(1001, 1200)], 1),
        ("TRANSCRIPT", "T3", "chr2", &[(1001, 1200)], 1),
        ("GENE", "", "chrUn", &[(1, 10)], 1),
        ("TRANSCRIPT", "T4", "chrUn", &[(1, 10)], 1),
    ];
    assert_eq!(rows.len(), 7);
    for (row, (level, transcript, chrom, spans, count)) in rows.iter().zip(want) {
        assert_eq!(text(row, "LEVEL"), level);
        assert_eq!(text(row, "TRANSCRIPT_ID"), transcript);
        assert_eq!(text(row, "CHROM"), chrom);
        assert_eq!(text(row, "EXON_COUNT"), count.to_string());
        let length: i64 = spans.iter().map(|(s, e)| e - s + 1).sum();
        assert_eq!(text(row, "LENGTH"), length.to_string());
        assert_eq!(text(row, "START"), spans[0].0.to_string());
        assert_eq!(text(row, "END"), spans[spans.len() - 1].1.to_string());
        if chrom != "chrUn" {
            assert_eq!(
                depth_cells(row),
                expected(&counted, chrom, spans),
                "{level} {transcript}"
            );
        }
    }
}

#[test]
fn bed12_blocks_give_the_rows_of_the_same_exons_in_a_gtf() {
    let dir = work_dir();
    let gtf_path = dir.join("same.gtf");
    fs::write(&gtf_path, gtf()).expect("gtf");
    // T1 and T2 as BED12 lines, blocks written out of order with offsets
    // from chromStart; Ensembl contig names, read as chr1.
    let bed = dir.join("same.bed");
    fs::write(
        &bed,
        "1\t100\t321\tT1\t0\t+\t100\t321\t0\t2\t1,100,\t220,0,\n\
         1\t150\t600\tT2\t0\t+\t150\t600\t0\t2\t110,100,\t0,350,\n",
    )
    .expect("bed");
    let from_gtf = run(
        &gtf_path.display().to_string(),
        &["TRANSCRIPT=T1", "TRANSCRIPT=T2"],
    );
    let from_bed = run(&bed.display().to_string(), &[]);
    let gtf_rows = rows(from_gtf.main.as_ref().expect("gtf rows"));
    let bed_rows = rows(from_bed.main.as_ref().expect("bed rows"));
    assert_eq!(gtf_rows.len(), 4);
    assert_eq!(bed_rows.len(), 4);
    for (g, b) in gtf_rows.iter().zip(&bed_rows) {
        for column in [
            "TRANSCRIPT_ID",
            "EXON_NUMBER",
            "CHROM",
            "START",
            "END",
            "STRAND",
            "LENGTH",
        ] {
            assert_eq!(text(g, column), text(b, column), "{column}");
        }
        assert_eq!(depth_cells(g), depth_cells(b));
    }
}

#[test]
fn filters_and_features_choose_the_rows() {
    let dir = work_dir();
    let annotation = dir.join("filters.gtf");
    fs::write(&annotation, gtf()).expect("gtf");
    let path = annotation.display().to_string();
    let count = |extra: &[&str]| run(&path, extra).main.expect("rows").rows.len();
    assert_eq!(count(&["GENE=GB"]), 1);
    assert_eq!(count(&["GENE=ENSGA"]), 4);
    assert_eq!(count(&["TRANSCRIPT=T2"]), 2);
    // No tag and no CDS: the most exon bases, T2 (210) over T1 (101).
    assert_eq!(count(&["CANONICAL=true"]), 4, "T2, T3, T4");
    let list = dir.join("genes.txt");
    fs::write(&list, "GB\nGC\n").expect("list");
    assert_eq!(count(&[&format!("GENE_LIST={}", list.display())]), 2);
    let err = registry::find("CollectExonCoverageMetrics")
        .expect("registered")
        .execute(&[
            format!("I={}/planted0.sam", dir.display()),
            format!("ANNOTATION={path}"),
            "GENE=BRCA1".to_owned(),
        ])
        .expect_err("unknown gene");
    assert!(err.to_string().contains("BRCA1"), "{err}");
}
