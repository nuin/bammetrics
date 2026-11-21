use std::path::{Path, PathBuf};

use super::*;
use crate::indexed_fasta::testing::{scratch, write_fasta};

const C1: &[u8] = b"ACGTNNACGTnnacgtNNNNNNACGTNNNNRYKMACGTN.ACGNNNN";

fn reference(dir: &Path) -> PathBuf {
    write_fasta(
        dir,
        "a.fa",
        &[("c1", C1), ("c2", b"NNNNNNNNNN"), ("c3", b"acgtACGTNacgt")],
        7,
        "\n",
    )
}

fn run(dir: &Path, extra: &[String]) -> String {
    let mut argv = vec![format!("I={}", dir.join("a.fa").display())];
    argv.extend(extra.iter().cloned());
    let (_, result) = TOOL.execute(&argv).expect("runs");
    result.main.expect("output").text.expect("text")
}

#[test]
fn the_counts_are_picards() {
    let dir = scratch("non-n");
    reference(&dir);
    assert_eq!(run(&dir, &[]), "40\n");
    let list = |name: &str, text: &str| {
        let path = dir.join(name);
        std::fs::write(&path, text).expect("list");
        format!("INTERVALS={}", path.display())
    };
    let header = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:47\n@SQ\tSN:c2\tLN:10\n@SQ\tSN:c3\tLN:13\n";
    let il = list(
        "il.interval_list",
        &format!("{header}c1\t30\t45\t+\tx\nc1\t1\t10\t-\tz\nc3\t5\t13\t+\ty\nc1\t3\t8\t+\to\n"),
    );
    assert_eq!(run(&dir, &[il]), "28\n");
    // Matched by index: the list's first @SQ line stands for c1.
    let swap = list(
        "swap.interval_list",
        "@HD\tVN:1.6\n@SQ\tSN:c3\tLN:13\n@SQ\tSN:c1\tLN:47\nc1\t2\t30\t+\tx\nc3\t5\t13\t+\ty\n",
    );
    assert_eq!(run(&dir, &[swap]), "5\n");
    let other = list(
        "other.interval_list",
        "@SQ\tSN:zz\tLN:47\nzz\t1\t47\t+\tx\nzz\t5\t4\t+\te\n",
    );
    assert_eq!(run(&dir, &[other]), "28\n");
    let empty = list("empty.interval_list", header);
    assert_eq!(run(&dir, &[empty]), "0\n");
}

#[test]
fn a_mask_crosses_chunk_and_line_boundaries() {
    let bases: Vec<u8> = (0..300)
        .map(|i| if i % 3 == 0 { b'n' } else { b'a' })
        .collect();
    let intervals = [
        Interval::new(0, 2, 2),
        Interval::new(0, 10, 150),
        Interval::new(0, 290, 400),
    ];
    let inside = |p: i64| {
        intervals
            .iter()
            .any(|i| i64::from(i.start) <= p && p <= i64::from(i.end))
    };
    let expected = (1..=300i64)
        .filter(|&p| inside(p) && bases[(p - 1) as usize] != b'n')
        .count() as i64;
    for size in [1, 7, 64, 300] {
        let mut next = 0;
        let mut total = 0;
        for (k, slice) in bases.chunks(size).enumerate() {
            total += count_bases(slice, (k * size) as i64 + 1, Some(&intervals), &mut next);
        }
        assert_eq!(total, expected, "{size}");
    }
}

#[test]
fn a_missing_index_or_dictionary_fails_as_picard_does() {
    let dir = scratch("non-n-missing");
    let path = reference(&dir);
    let dict = crate::create_sequence_dictionary::default_dictionary(&path).expect("name");
    std::fs::remove_file(crate::indexed_fasta::fai_path(&path)).expect("no .fai");
    assert!(TOOL.execute(&[format!("I={}", path.display())]).is_err());
    std::fs::write(&dict, "@HD\tVN:1.6\n").expect("empty .dict");
    assert!(TOOL.execute(&[format!("I={}", path.display())]).is_err());
    std::fs::remove_file(&dict).expect("no .dict");
    assert!(TOOL.execute(&[format!("I={}", path.display())]).is_err());
}
