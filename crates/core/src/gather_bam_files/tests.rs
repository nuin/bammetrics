use std::os::unix::fs::DirBuilderExt;

use super::*;

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-gather-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data")
        .join(name)
}

fn strings(paths: &[&Path]) -> Vec<String> {
    paths.iter().map(|p| p.display().to_string()).collect()
}

#[test]
fn lists_unroll_in_the_order_given_as_htsjdk_unrolls_them() {
    let dir = scratch("unroll");
    let inner = dir.join("inner.txt");
    std::fs::write(&inner, "x.cram\r\ny.sam\r").expect("inner");
    let list = dir.join("inputs.list");
    std::fs::write(
        &list,
        format!("  b.bam \t\n\n{}\n\u{1}c.bam\u{1f}\n", inner.display()),
    )
    .expect("list");
    let (inputs, lists) =
        unroll(&strings(&[Path::new("a.bam"), &list, Path::new("d.bam")])).expect("unroll");
    let names: Vec<String> = inputs.iter().map(|p| p.display().to_string()).collect();
    assert_eq!(
        names,
        ["a.bam", "b.bam", "x.cram", "y.sam", "c.bam", "d.bam"]
    );
    assert_eq!(lists, [list.clone(), inner.clone()]);
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn a_list_that_names_itself_ends_with_an_error() {
    let dir = scratch("loop");
    let list = dir.join("self.list");
    std::fs::write(&list, format!("{}\n{}\n", list.display(), list.display())).expect("list");
    let message = unroll(&strings(&[&list])).expect_err("a loop").to_string();
    assert!(message.contains("does a list name itself"), "{message}");
    let missing = dir.join("missing.list");
    assert!(unroll(&strings(&[&missing])).is_err());
    std::fs::write(&list, b"a.bam\n\xff\n").expect("list");
    let message = unroll(&strings(&[&list]))
        .expect_err("not UTF-8")
        .to_string();
    assert!(message.contains("not UTF-8"), "{message}");
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn a_list_stops_at_the_input_cap() {
    let dir = scratch("cap");
    let list = dir.join("many.list");
    std::fs::write(&list, "a.bam\n".repeat(MAX_INPUTS)).expect("list");
    let (inputs, _) = unroll(&strings(&[&list])).expect("at the cap");
    assert_eq!(inputs.len(), MAX_INPUTS);
    std::fs::write(&list, "a.bam\n".repeat(MAX_INPUTS + 1)).expect("list");
    let message = unroll(&strings(&[&list])).expect_err("over").to_string();
    assert!(message.contains("more than 100000 inputs"), "{message}");
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn the_block_copy_index_is_named_by_the_basename() {
    for (output, index) in [
        ("out.bam", "out.bai"),
        ("dir/out.dat", "dir/out.bai"),
        ("dir/out", "dir/out.bai"),
        ("a.b.bam", "a.b.bai"),
        (".bam", ".bam.bai"),
    ] {
        assert_eq!(
            block_copy_index(Path::new(output)),
            Path::new(index),
            "{output}"
        );
    }
}

/// Run the tool with `args` (Picard syntax).
fn run(args: &[String]) -> Result<ToolOutput> {
    let inv = TOOL.parse(args)?;
    (TOOL.run)(&inv)
}

#[test]
fn a_gather_of_one_bam_is_the_bam() {
    let dir = scratch("one");
    let out = dir.join("out.bam");
    let input = data("markdup.bam");
    run(&[
        format!("I={}", input.display()),
        format!("O={}", out.display()),
    ])
    .expect("run");
    assert_eq!(
        std::fs::read(&out).expect("out"),
        std::fs::read(&input).expect("in")
    );
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn bad_gathers_fail_and_leave_no_file() {
    let dir = scratch("bad");
    let out = dir.join("out.bam");
    let fails = |args: Vec<String>, wanted: &str| {
        let message = run(&args).expect_err(wanted).to_string();
        assert!(message.contains(wanted), "{message}");
        assert!(!out.exists(), "{wanted}: an output was left");
    };
    let markdup = format!("I={}", data("markdup.bam").display());
    let o = format!("O={}", out.display());
    // A file with records under another dictionary.
    let other = dir.join("other.bam");
    let mut header = b"@HD\tVN:1.6\n@SQ\tSN:chrX\tLN:5\n".to_vec();
    header.push(b'\n');
    let mut bytes = b"BAM\x01".to_vec();
    bytes.extend_from_slice(&(header.len() as i32).to_le_bytes());
    bytes.extend_from_slice(&header);
    bytes.extend_from_slice(&1i32.to_le_bytes());
    bytes.extend_from_slice(&5i32.to_le_bytes());
    bytes.extend_from_slice(b"chrX\0");
    bytes.extend_from_slice(&5i32.to_le_bytes());
    // Bytes after the header: the file holds records.
    bytes.extend_from_slice(b"records");
    let mut file = Vec::new();
    bgzf_copy::write_blocks(&mut file, &bytes, 5).expect("blocks");
    file.extend_from_slice(&bgzf_copy::EOF_BLOCK);
    std::fs::write(&other, file).expect("other");
    fails(
        vec![markdup.clone(), format!("I={}", other.display()), o.clone()],
        "its sequence dictionary is not the one of the first INPUT",
    );
    let cut = dir.join("cut.bam");
    let bytes = std::fs::read(data("insert_size.bam")).expect("bam");
    std::fs::write(&cut, &bytes[..bytes.len() / 2]).expect("cut");
    fails(
        vec![markdup.clone(), format!("I={}", cut.display()), o.clone()],
        "valid GZIP block at the end",
    );
    fails(
        vec![
            markdup.clone(),
            format!("I={}", dir.join("none.bam").display()),
            o.clone(),
        ],
        "none.bam",
    );
    fails(
        vec![
            markdup.clone(),
            format!("O={}", dir.join("out.sam").display()),
        ],
        "Picard writes BAM bytes under this name",
    );
    let copy = dir.join("copy.bam");
    std::fs::copy(data("markdup.bam"), &copy).expect("copy");
    let message = run(&[
        format!("I={}", copy.display()),
        markdup.clone(),
        format!("O={}", copy.display()),
    ])
    .expect_err("output is an input")
    .to_string();
    assert!(message.contains("writing it would destroy it"), "{message}");
    assert_eq!(
        std::fs::read(&copy).expect("copy"),
        std::fs::read(data("markdup.bam")).expect("in")
    );
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn a_sam_input_gathers_record_by_record() {
    let dir = scratch("sam");
    let out = dir.join("out.bam");
    let message = run(&[
        format!("I={}", data("markdup.bam").display()),
        format!("I={}", data("markdup.sam").display()),
        format!("O={}", out.display()),
    ])
    .expect_err("the second file starts before the first ends")
    .to_string();
    assert!(
        message.contains("Alignments added out of order"),
        "{message}"
    );
    assert!(!out.exists());
    run(&[
        format!("I={}", data("markdup.sam").display()),
        format!("O={}", out.display()),
    ])
    .expect("run");
    let count = |path: &Path| {
        pass::with_record_stream_at(
            path,
            None,
            1,
            crate::sam_text::Stringency::Silent,
            |_, r, _| {
                let mut record = RawRecord::default();
                let mut n = 0;
                while record.read_from(r).map_err(Error::Record)? {
                    n += 1;
                }
                Ok(n)
            },
        )
        .expect("count")
    };
    assert_eq!(count(&out), count(&data("markdup.bam")));
    std::fs::remove_dir_all(&dir).expect("scratch");
}
