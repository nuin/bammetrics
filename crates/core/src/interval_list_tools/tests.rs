use std::fs;
use std::path::{Path, PathBuf};

use super::*;

const LIST: &str = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:500\n\
                    c2\t10\t60\t+\tb\nc1\t100\t200\t-\ta\nc1\t150\t300\t+\tc\n";

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-interval-list-tools-{name}-{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

fn args(list: &[&str]) -> Vec<String> {
    list.iter().map(|s| (*s).to_owned()).collect()
}

fn run(list: &[String]) -> Result<ToolOutput> {
    let inv = TOOL.parse(list)?;
    let output = (TOOL.run)(&inv)?;
    registry::write_files(&inv, &output)?;
    Ok(output)
}

fn path(dir: &Path, name: &str) -> String {
    dir.join(name).display().to_string()
}

/// The `CL` of Picard 3.5.0's `@PG` line for these argument lists, as
/// IntervalListTools wrote it.
#[test]
fn the_command_line_is_picards() {
    let cases: [(&[&str], &str); 4] = [
        (
            &["I=lists/a.interval_list", "O=p.interval_list", "COMMENT=hi", "QUIET=true"],
            "IntervalListTools INPUT=[lists/a.interval_list] OUTPUT=p.interval_list \
             COMMENT=[hi] QUIET=true    PADDING=0 UNIQUE=false DONT_MERGE_ABUTTING=false \
             SORT=true ACTION=CONCAT SCATTER_COUNT=1 INCLUDE_FILTERED=false \
             BREAK_BANDS_AT_MULTIPLES_OF=0 SUBDIVISION_MODE=INTERVAL_SUBDIVISION INVERT=false \
             OUTPUT_VALUE=NONE VERBOSITY=INFO VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 \
             MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false \
             USE_JDK_DEFLATER=false USE_JDK_INFLATER=false",
        ),
        (
            &[
                "QUIET=TRUE",
                "VERBOSITY=ERROR",
                "O=p2.interval_list",
                "COMMENT=a b",
                "I=lists/a.interval_list",
                "PADDING=007",
                "I=lists/b.interval_list",
                "UNIQUE=false",
                "TMP_DIR=/tmp",
                "COUNT_OUTPUT=null",
                "SCATTER_CONTENT=null",
                "R=null",
                "COMPRESSION_LEVEL=5",
            ],
            "IntervalListTools INPUT=[lists/a.interval_list, lists/b.interval_list] \
             OUTPUT=p2.interval_list PADDING=7 UNIQUE=false COMMENT=[a b] SCATTER_CONTENT=null \
             COUNT_OUTPUT=null TMP_DIR=[/tmp] VERBOSITY=ERROR QUIET=true COMPRESSION_LEVEL=5 \
             REFERENCE_SEQUENCE=null    DONT_MERGE_ABUTTING=false SORT=true ACTION=CONCAT \
             SCATTER_COUNT=1 INCLUDE_FILTERED=false BREAK_BANDS_AT_MULTIPLES_OF=0 \
             SUBDIVISION_MODE=INTERVAL_SUBDIVISION INVERT=false OUTPUT_VALUE=NONE \
             VALIDATION_STRINGENCY=STRICT MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false \
             CREATE_MD5_FILE=false USE_JDK_DEFLATER=false USE_JDK_INFLATER=false",
        ),
        (
            &[
                "--INPUT",
                "lists/a.interval_list",
                "--OUTPUT",
                "p3.interval_list",
                "--COMMENT",
                "x",
                "--QUIET",
                "true",
            ],
            "IntervalListTools --INPUT lists/a.interval_list --OUTPUT p3.interval_list \
             --COMMENT x --QUIET true --PADDING 0 --UNIQUE false --DONT_MERGE_ABUTTING false \
             --SORT true --ACTION CONCAT --SCATTER_COUNT 1 --INCLUDE_FILTERED false \
             --BREAK_BANDS_AT_MULTIPLES_OF 0 --SUBDIVISION_MODE INTERVAL_SUBDIVISION --INVERT \
             false --OUTPUT_VALUE NONE --VERBOSITY INFO --VALIDATION_STRINGENCY STRICT \
             --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false \
             --CREATE_MD5_FILE false --help false --version false --showHidden false \
             --USE_JDK_DEFLATER false --USE_JDK_INFLATER false",
        ),
        (
            &[
                "--INPUT",
                "lists/a.interval_list",
                "--I",
                "lists/b.interval_list",
                "--OUTPUT",
                "p4.interval_list",
                "--COMMENT",
                "x y",
                "--COMMENT",
                "z",
                "--PADDING",
                "007",
                "--SCATTER_CONTENT",
                "null",
                "--TMP_DIR",
                "/tmp",
                "--VERBOSITY",
                "ERROR",
                "--QUIET",
                "true",
            ],
            "IntervalListTools --INPUT lists/a.interval_list --INPUT lists/b.interval_list \
             --OUTPUT p4.interval_list --PADDING 7 --COMMENT x y --COMMENT z  --TMP_DIR /tmp \
             --VERBOSITY ERROR --QUIET true --UNIQUE false --DONT_MERGE_ABUTTING false --SORT \
             true --ACTION CONCAT --SCATTER_COUNT 1 --INCLUDE_FILTERED false \
             --BREAK_BANDS_AT_MULTIPLES_OF 0 --SUBDIVISION_MODE INTERVAL_SUBDIVISION --INVERT \
             false --OUTPUT_VALUE NONE --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 \
             --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --help \
             false --version false --showHidden false --USE_JDK_DEFLATER false \
             --USE_JDK_INFLATER false",
        ),
    ];
    for (given, expected) in cases {
        let inv = TOOL.parse(&args(given)).expect("parses");
        assert_eq!(picard_command_line(&inv), expected, "{given:?}");
    }
}

#[test]
fn only_a_whole_vcf_name_is_a_vcf() {
    for name in ["x/.vcf", ".vcf.gz", "d/.vcf.bgz", ".bcf"] {
        assert!(vcf_name(Path::new(name)).is_some(), "{name}");
    }
    for name in ["a.vcf", "a.vcf.gz", "vcf", "/dev/stdin", "a.interval_list"] {
        assert!(vcf_name(Path::new(name)).is_none(), "{name}");
    }
    let dir = scratch("bcf");
    let bcf = dir.join(".bcf");
    fs::write(&bcf, "x").expect("write");
    let error = run(&args(&[&format!("I={}", bcf.display())])).expect_err("refused");
    assert!(error.to_string().contains("reads no BCF"), "{error}");
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_union_is_written_and_counted() {
    let dir = scratch("union");
    let input = path(&dir, "in.interval_list");
    fs::write(&input, LIST).expect("write");
    let out = path(&dir, "out.interval_list");
    let count = path(&dir, "count.txt");
    let output = run(&args(&[
        &format!("I={input}"),
        "ACTION=UNION",
        &format!("O={out}"),
        "OUTPUT_VALUE=BASES",
        &format!("COUNT_OUTPUT={count}"),
        "COMMENT=note",
    ]))
    .expect("runs");
    let text = fs::read_to_string(&out).expect("written");
    assert!(
        text.ends_with("@CO\tnote\nc1\t100\t300\t-\ta|c\nc2\t10\t60\t+\tb\n"),
        "{text}"
    );
    // The header is the input's: the sort order is set on the sorted copy.
    assert!(text.starts_with("@HD\tVN:1.6\n@SQ"), "{text}");
    assert_eq!(fs::read_to_string(&count).expect("count"), "252\n");
    let rows = output.main.expect("rows").materialized().rows.clone();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0][0], Value::Int(1));
    assert_eq!(rows[0][5], Value::Str("a|c".to_owned()));
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_scatter_writes_its_directory_and_counts_every_list() {
    let dir = scratch("scatter");
    let input = path(&dir, "in.interval_list");
    fs::write(&input, LIST).expect("write");
    let out = dir.join("deep/scatter");
    let output = run(&args(&[
        &format!("I={input}"),
        &format!("O={}", out.display()),
        "SCATTER_COUNT=2",
        "M=INTERVAL_COUNT",
        "OUTPUT_VALUE=INTERVALS",
    ]))
    .expect("runs");
    for n in 1..=2 {
        assert!(
            out.join(format!("temp_000{n}_of_2/scattered.interval_list"))
                .is_file()
        );
    }
    let main = output.main.expect("main");
    assert_eq!(main.text.as_deref(), Some("3\n"));
    let rows = main.materialized().rows.clone();
    assert_eq!(rows.iter().filter(|r| r[0] == Value::Int(2)).count(), 2);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_failing_run_writes_nothing() {
    let dir = scratch("refused");
    let input = path(&dir, "in.interval_list");
    fs::write(&input, LIST).expect("write");
    let out = path(&dir, "out.interval_list");
    let refused = [
        // COUNT_OUTPUT is OUTPUT.
        vec![
            format!("I={input}"),
            format!("O={out}"),
            "OUTPUT_VALUE=BASES".to_owned(),
            format!("COUNT_OUTPUT={out}"),
        ],
        // OUTPUT under a missing directory, as Picard fails.
        vec![format!("I={input}"), format!("O={}", path(&dir, "missing/x"))],
        // COUNT_OUTPUT under a missing directory.
        vec![
            format!("I={input}"),
            format!("O={out}"),
            "OUTPUT_VALUE=BASES".to_owned(),
            format!("COUNT_OUTPUT={}", path(&dir, "missing/count")),
        ],
        // OUTPUT is the input.
        vec![format!("I={input}"), format!("O={input}")],
        // The checks of customCommandLineValidation.
        vec![format!("I={input}"), "ACTION=SUBTRACT".to_owned()],
        vec![format!("I={input}"), format!("SI={input}")],
        vec![format!("I={input}"), "COUNT_OUTPUT=x".to_owned()],
        vec![format!("I={input}"), "SCATTER_COUNT=0".to_owned()],
        vec![format!("I={input}"), "PADDING=-1".to_owned()],
        vec![format!("I={input}"), "ACTION=union".to_owned()],
        vec![format!("I={input}"), "ASSUME_SORTED=true".to_owned()],
    ];
    for list in refused {
        assert!(run(&list).is_err(), "{list:?}");
        assert!(!Path::new(&out).exists(), "{list:?}");
        assert!(!dir.join("missing").exists(), "{list:?}");
    }
    assert_eq!(fs::read_to_string(&input).expect("kept"), LIST);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_vcf_named_dot_vcf_is_read_with_include_filtered() {
    let dir = scratch("vcf");
    let vcf = dir.join(".vcf");
    fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n##contig=<ID=c1,length=1000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         c1\t10\trs1\tA\tC\t.\tPASS\t.\nc1\t20\t.\tAG\tA\t.\tq10\t.\n",
    )
    .expect("write");
    let rows = |filtered: &str| {
        let output = run(&args(&[
            &format!("I={}", vcf.display()),
            &format!("INCLUDE_FILTERED={filtered}"),
        ]))
        .expect("runs");
        output.main.expect("rows").materialized().rows.len()
    };
    assert_eq!(rows("false"), 1);
    assert_eq!(rows("true"), 2);
    let _ = fs::remove_dir_all(&dir);
}
