use std::path::Path;

use super::*;
use crate::indexed_fasta::testing::{scratch, write_fasta};

/// The FASTA of the Picard 3.5.0 runs below: N runs of 1, 2, 4 and 6
/// bases, lower case, IUPAC codes, `.`, an all-N contig.
const C1: &[u8] = b"ACGTNNACGTnnacgtNNNNNNACGTNNNNRYKMACGTN.ACGNNNN";

fn reference(dir: &Path) -> PathBuf {
    write_fasta(
        dir,
        "a.fa",
        &[("c1", C1), ("c2", b"NNNNNNNNNN"), ("c3", b"acgtACGTNacgt")],
        11,
        "\n",
    )
}

fn run(dir: &Path, args: &[&str]) -> (String, Vec<Vec<Value>>) {
    let output = dir.join("out.interval_list");
    let mut argv = vec![
        format!("R={}", dir.join("a.fa").display()),
        format!("O={}", output.display()),
    ];
    argv.extend(args.iter().map(|a| (*a).to_owned()));
    let (_, result) = TOOL.execute(&argv).expect("runs");
    let rows = result.main.expect("rows").rows;
    let text = std::fs::read_to_string(&output).expect("written");
    let body = text
        .lines()
        .filter(|l| !l.starts_with("@SQ"))
        .collect::<Vec<_>>()
        .join("\n");
    (body, rows)
}

fn lines(runs: &[(&str, i32, i32, &str)]) -> String {
    let mut text = vec!["@HD\tVN:1.6\tSO:coordinate".to_owned()];
    text.extend(
        runs.iter()
            .map(|(c, s, e, n)| format!("{c}\t{s}\t{e}\t+\t{n}mer")),
    );
    text.join("\n")
}

#[test]
fn the_intervals_are_picards() {
    let dir = scratch("scatter-picard");
    reference(&dir);
    let both_1 = [
        ("c1", 1, 4, "ACGT"),
        ("c1", 5, 6, "N"),
        ("c1", 7, 10, "ACGT"),
        ("c1", 11, 12, "N"),
        ("c1", 13, 16, "ACGT"),
        ("c1", 17, 22, "N"),
        ("c1", 23, 26, "ACGT"),
        ("c1", 27, 30, "N"),
        ("c1", 31, 38, "ACGT"),
        ("c1", 39, 40, "N"),
        ("c1", 41, 43, "ACGT"),
        ("c1", 44, 47, "N"),
        ("c2", 1, 10, "N"),
        ("c3", 1, 13, "ACGT"),
    ];
    assert_eq!(run(&dir, &[]).0, lines(&both_1));
    assert_eq!(run(&dir, &["MAX_TO_MERGE=1", "OT=BOTH"]).0, lines(&both_1));
    let mut both_0 = both_1[..13].to_vec();
    both_0.extend([
        ("c3", 1, 8, "ACGT"),
        ("c3", 9, 9, "N"),
        ("c3", 10, 13, "ACGT"),
    ]);
    assert_eq!(run(&dir, &["MAX_TO_MERGE=0"]).0, lines(&both_0));
    assert_eq!(run(&dir, &["MAX_TO_MERGE=-1"]).0, lines(&both_0));
    let both_2 = [
        ("c1", 1, 16, "ACGT"),
        ("c1", 17, 22, "N"),
        ("c1", 23, 26, "ACGT"),
        ("c1", 27, 30, "N"),
        ("c1", 31, 43, "ACGT"),
        ("c1", 44, 47, "N"),
        ("c2", 1, 10, "N"),
        ("c3", 1, 13, "ACGT"),
    ];
    assert_eq!(run(&dir, &["N=2"]).0, lines(&both_2));
    let acgt_4 = [
        ("c1", 1, 16, "ACGT"),
        ("c1", 23, 43, "ACGT"),
        ("c3", 1, 13, "ACGT"),
    ];
    assert_eq!(run(&dir, &["N=4", "OT=ACGT"]).0, lines(&acgt_4));
    let n_1: Vec<_> = both_1.iter().copied().filter(|r| r.3 == "N").collect();
    assert_eq!(run(&dir, &["OUTPUT_TYPE=N"]).0, lines(&n_1));
}

#[test]
fn the_rows_count_both_kinds_per_contig() {
    let dir = scratch("scatter-rows");
    reference(&dir);
    let (_, rows) = run(&dir, &["OT=N"]);
    let int = |v: &Value| match v {
        Value::Int(i) => *i,
        other => panic!("{other:?}"),
    };
    let c1: Vec<i64> = rows[0][1..].iter().map(int).collect();
    assert_eq!(c1, [47, 6, 27, 6, 20]);
    let c2: Vec<i64> = rows[1][1..].iter().map(int).collect();
    assert_eq!(c2, [10, 0, 0, 1, 10]);
}

#[test]
fn the_header_is_the_dict_and_a_failed_run_leaves_no_output() {
    let dir = scratch("scatter-fail");
    let path = reference(&dir);
    let dict = crate::create_sequence_dictionary::default_dictionary(&path).expect("name");
    std::fs::write(
        &dict,
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:c1\tLN:47\tM5:ab\tUR:x\n@SQ\tSN:c2\tLN:10\n\
         @SQ\tSN:c3\tLN:13\n@CO\tdropped\n",
    )
    .expect(".dict");
    let output = dir.join("out.interval_list");
    let argv = |extra: &str| {
        vec![
            format!("R={}", path.display()),
            format!("O={}", output.display()),
            extra.to_owned(),
        ]
    };
    TOOL.execute(&argv("N=1")).expect("runs");
    let text = std::fs::read_to_string(&output).expect("written");
    assert!(text.starts_with(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:47\tM5:ab\tUR:x\n@SQ\tSN:c2\tLN:10\n\
         @SQ\tSN:c3\tLN:13\nc1\t1\t4\t+\tACGTmer\n"
    ));
    std::fs::remove_file(&output).expect("removed");
    for bad in ["OT=acgt", "N=x", "N=2147483648"] {
        assert!(TOOL.execute(&argv(bad)).is_err(), "{bad}");
    }
    std::fs::write(
        crate::indexed_fasta::fai_path(&path),
        "c1\t47\t21\t11\t12\nc2\t0\t90\t11\t12\nc3\t13\t110\t11\t12\n",
    )
    .expect(".fai");
    std::fs::write(
        &dict,
        "@SQ\tSN:c1\tLN:47\n@SQ\tSN:c2\tLN:0\n@SQ\tSN:c3\tLN:13\n",
    )
    .expect(".dict");
    let error = TOOL.execute(&argv("N=1")).expect_err("fails").to_string();
    assert!(error.contains("no base"), "{error}");
    assert!(!output.exists());
    std::fs::remove_file(&dict).expect("no .dict");
    let error = TOOL.execute(&argv("N=1")).expect_err("fails").to_string();
    assert!(error.contains("dictionary"), "{error}");
    assert!(!output.exists());
    let error = TOOL
        .execute(&[
            format!("R={}", path.display()),
            format!("O={}", path.display()),
        ])
        .expect_err("refused")
        .to_string();
    assert!(error.contains("is an input"), "{error}");
}

#[test]
fn the_merge_decides_at_the_head_with_three_in_hand() {
    let mut out = Vec::new();
    let run = |contig, start, end, no_call| Run {
        contig,
        start,
        end,
        no_call,
    };
    let mut merger = Merger::new(3, |r: Run| {
        out.push((r.contig, r.start, r.end, r.no_call));
        Ok(())
    });
    for r in [
        run(0, 1, 5, false),
        run(0, 6, 8, true),
        run(0, 9, 9, false),
        run(0, 10, 10, true),
        run(0, 11, 20, false),
        run(0, 21, 30, true),
        run(1, 1, 2, false),
    ] {
        merger.push(r).expect("push");
    }
    merger.finish().expect("finish");
    assert_eq!(
        out,
        [(0, 1, 20, false), (0, 21, 30, true), (1, 1, 2, false)]
    );
}
