//! The fan-out on synthetic headers: the open-file bound, the rename at
//! the end, nothing left on a failure, Picard's names and sizes.

use std::path::{Path, PathBuf};

use super::*;
use crate::bam_raw::build;

/// A scratch directory of this test, readable by this user only.
fn scratch(name: &str) -> PathBuf {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!("bammetrics-split-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .expect("dir")
        .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// The records of a BAM as bytes after inflation, and its header text.
fn read_bam(path: &Path) -> (String, Vec<Vec<u8>>) {
    let file = std::io::BufReader::new(std::fs::File::open(path).expect("open"));
    let mut reader = crate::gzip::GzipReader::new(file);
    let mut data = Vec::new();
    std::io::Read::read_to_end(&mut reader, &mut data).expect("inflate");
    assert_eq!(&data[..4], b"BAM\x01");
    let text_len = i32::from_le_bytes(data[4..8].try_into().expect("4")) as usize;
    let text = String::from_utf8_lossy(&data[8..8 + text_len]).into_owned();
    let mut at = 8 + text_len;
    let refs = i32::from_le_bytes(data[at..at + 4].try_into().expect("4"));
    at += 4;
    for _ in 0..refs {
        let n = i32::from_le_bytes(data[at..at + 4].try_into().expect("4")) as usize;
        at += 8 + n;
    }
    let mut records = Vec::new();
    while at < data.len() {
        let size = i32::from_le_bytes(data[at..at + 4].try_into().expect("4")) as usize;
        records.push(data[at + 4..at + 4 + size].to_vec());
        at += 4 + size;
    }
    (text, records)
}

fn options() -> WriterOptions {
    WriterOptions {
        compression_level: 1,
        create_index: false,
        create_md5: false,
        threads: 1,
    }
}

fn header(libraries: usize) -> sam::Header {
    let mut text = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c1\tLN:1000\n");
    for n in 0..libraries {
        text.push_str(&format!("@RG\tID:rg{n}\tSM:s\tLB:lib {n}\n"));
    }
    text.push_str("@RG\tID:nolb\tSM:s\n");
    text.parse().expect("header")
}

fn tagged(name: &str, rg: Option<&str>) -> RawRecord {
    let tags = rg.map_or_else(Vec::new, |rg| {
        let mut tags = b"RGZ".to_vec();
        tags.extend_from_slice(rg.as_bytes());
        tags.push(0);
        tags
    });
    build(name, 4, -1, -1, &[], &[], &tags)
}

fn library_split(dir: &Path, header: &sam::Header) -> Splitter {
    library_split_on(dir, header, 1)
}

/// A library split whose blocks `threads` workers compress.
fn library_split_on(dir: &Path, header: &sam::Header, threads: usize) -> Splitter {
    let options = WriterOptions {
        threads,
        ..options()
    };
    let fanout = Fanout::new(options, header, Vec::new(), None).expect("fanout");
    let layout = Layout {
        dir,
        extension: "bam",
        format: Format::Bam,
        input: "in.bam",
        total: 0,
        wide: &[],
    };
    Splitter::new(&SplitKind::Library, header, &layout, fanout).expect("splitter")
}

/// 300 libraries, records round-robin: never more than 64 files open, and
/// every file holds its records in order.
#[test]
fn three_hundred_libraries_keep_at_most_the_bound_open() {
    three_hundred_libraries("bound", 1);
}

/// The same files when four workers compress the blocks.
#[test]
fn workers_compress_the_blocks_of_every_file_in_order() {
    three_hundred_libraries("workers", 4);
    // Records enough for several full blocks per file, two files.
    let dir = scratch("blocks");
    let header = header(2);
    let mut split = library_split_on(&dir, &header, 4);
    for i in 0..40_000 {
        let rg = format!("rg{}", i % 2);
        split
            .push(&tagged(&format!("read{i:06}"), Some(&rg)))
            .expect("push");
    }
    split.finish().expect("finish");
    for n in 0..2 {
        let (_, records) = read_bam(&dir.join(format!("lib_{n}.bam")));
        assert_eq!(records.len(), 20_000);
        for (k, record) in records.iter().enumerate() {
            let name = format!("read{:06}\0", k * 2 + n);
            assert_eq!(
                &record[32..32 + name.len()],
                name.as_bytes(),
                "lib_{n} record {k}"
            );
        }
    }
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

fn three_hundred_libraries(name: &str, threads: usize) {
    let dir = scratch(name);
    let header = header(300);
    let mut split = library_split_on(&dir, &header, threads);
    for round in 0..5 {
        for n in 0..300 {
            let record = tagged(&format!("r{round}_{n}"), Some(&format!("rg{n}")));
            split.push(&record).expect("push");
            assert!(split.fanout.open_files() <= MAX_OPEN_FILES);
        }
    }
    split.push(&tagged("orphan", Some("zz"))).expect("orphan");
    split.push(&tagged("none", None)).expect("none");
    // Nothing is in place before the finish.
    assert!(names(&dir).iter().all(|name| name.starts_with('.')));
    let files = split.finish().expect("finish");
    assert_eq!(files.len(), 301);
    let listed = names(&dir);
    assert_eq!(listed.len(), 301, "no temporary file is left");
    assert!(listed.contains(&"lib_7.bam".to_owned()));
    let (text, records) = read_bam(&dir.join("lib_7.bam"));
    assert!(text.contains("@RG\tID:rg7\tSM:s\tLB:lib 7\n"), "{text}");
    assert!(!text.contains("rg8"), "{text}");
    assert_eq!(records.len(), 5);
    for (round, record) in records.iter().enumerate() {
        let name = format!("r{round}_7\0");
        assert_eq!(&record[32..32 + name.len()], name.as_bytes());
    }
    let (text, records) = read_bam(&dir.join("unknown.bam"));
    assert!(text.contains("ID:nolb") && !text.contains("LB:"), "{text}");
    assert_eq!(records.len(), 2);
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn a_failure_leaves_no_file() {
    let dir = scratch("failure");
    let header: sam::Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:c1\tLN:1000\n\
                                @RG\tID:a\tLB:A\n@RG\tID:b\tLB:B\n"
        .parse()
        .expect("header");
    let mut split = library_split(&dir, &header);
    split.push(&tagged("q2", Some("a"))).expect("first");
    split.push(&tagged("q1", Some("b"))).expect("another file");
    let error = split
        .push(&tagged("q1", Some("a")))
        .expect_err("out of queryname order");
    assert!(error.to_string().contains("out of order"), "{error}");
    drop(split);
    assert!(names(&dir).is_empty(), "{:?}", names(&dir));
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn two_libraries_with_one_safe_name_fail_before_any_record() {
    let dir = scratch("clash");
    let header: sam::Header = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n\
                                @RG\tID:a\tLB:a b\n@RG\tID:b\tLB:a_b\n"
        .parse()
        .expect("header");
    let fanout = Fanout::new(options(), &header, Vec::new(), None).expect("fanout");
    let layout = Layout {
        dir: &dir,
        extension: "bam",
        format: Format::Bam,
        input: "in.bam",
        total: 0,
        wide: &[],
    };
    let error = Splitter::new(&SplitKind::Library, &header, &layout, fanout)
        .err()
        .expect("one file for two libraries");
    assert!(error.to_string().contains("a_b.bam"), "{error}");
    assert!(names(&dir).is_empty());
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

/// SAM text with a library above U+00FF: the file takes the whole name,
/// and two names with the same low bytes are two libraries, as in Picard.
#[test]
fn a_wide_library_of_sam_text_names_its_file_whole() {
    let dir = scratch("wide");
    let text = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n\
                @RG\tID:a\tSM:s\tLB:中文\n@RG\tID:b\tSM:s\tLB:\u{3e2d}\u{6587}\n@RG\tID:c\tSM:s\tLB:café\n";
    let reader =
        crate::sam_text::Reader::new(text.as_bytes(), Stringency::Strict, None).expect("reader");
    let header = reader.header().expect("header");
    assert_eq!(reader.wide_lines().len(), 2);
    let fanout = Fanout::new(options(), &header, Vec::new(), None).expect("fanout");
    let layout = Layout {
        dir: &dir,
        extension: "sam",
        format: Format::Bam,
        input: "in.sam",
        total: 0,
        wide: reader.wide_lines(),
    };
    let split = Splitter::new(&SplitKind::Library, &header, &layout, fanout).expect("splitter");
    split.finish().expect("finish");
    assert_eq!(
        names(&dir),
        ["café.sam", "\u{3e2d}\u{6587}.sam", "中文.sam"]
    );
    assert_eq!(file_name_safe_text(" 中/文 "), "中_文");
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn a_header_without_a_library_is_picards_usage_error() {
    let dir = scratch("nolib");
    let header = header(0);
    let fanout = Fanout::new(options(), &header, Vec::new(), None).expect("fanout");
    let layout = Layout {
        dir: &dir,
        extension: "bam",
        format: Format::Bam,
        input: "/abs/in.bam",
        total: 0,
        wide: &[],
    };
    let error = Splitter::new(&SplitKind::Library, &header, &layout, fanout)
        .err()
        .expect("no library");
    assert!(matches!(error, Error::Usage(_)));
    assert_eq!(
        error.to_string(),
        "No individual libraries are specified in the header of /abs/in.bam"
    );
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn safe_names_follow_make_file_name_safe() {
    for (value, safe) in [
        (&b"Lib A"[..], "Lib_A"),
        (b" sp ", "sp"),
        (b"lib/x", "lib_x"),
        (b"a+b,c-d.e", "a+b,c-d.e"),
        (b"caf\xc3\xa9", "caf\u{c3}\u{a9}"),
        (b"caf\xe9", "caf\u{e9}"),
        (b"\x01x\x7f", "x\u{7f}"),
        (
            b"!\"#$%&'()*/:;<=>?@[]\\^`{|}~",
            "___________________________",
        ),
        (b"", ""),
    ] {
        assert_eq!(file_name_safe(value), safe, "{value:?}");
    }
}

#[test]
fn count_sizes_are_picards_double_arithmetic() {
    let options = |files, reads| CountOptions {
        files,
        reads,
        total: 0,
        prefix: "shard".to_owned(),
    };
    assert_eq!(options(3, 0).sizes(12918), (3, 4306));
    assert_eq!(options(0, 1000).sizes(12918), (13, 994));
    assert_eq!(options(0, 2).sizes(0), (0, 0));
    assert_eq!(options(2, 0).sizes(0), (2, 0));
    assert_eq!(options(0, 2).sizes(i64::MAX), (i32::MAX, i32::MAX));
}

/// Records with one name stay in one file past the size.
#[test]
fn a_count_split_keeps_a_name_together() {
    let dir = scratch("count");
    let header: sam::Header = "@HD\tVN:1.4\n@SQ\tSN:c1\tLN:1000\n"
        .parse()
        .expect("header");
    let fanout = Fanout::new(options(), &header, Vec::new(), None).expect("fanout");
    let layout = Layout {
        dir: &dir,
        extension: "bam",
        format: Format::Bam,
        input: "in.bam",
        total: 7,
        wide: &[],
    };
    let count = CountOptions {
        files: 3,
        reads: 0,
        total: 7,
        prefix: "part".to_owned(),
    };
    let mut split =
        Splitter::new(&SplitKind::Count(count), &header, &layout, fanout).expect("splitter");
    for name in ["a", "a", "a", "b", "c", "c", "d"] {
        split.push(&tagged(name, None)).expect("push");
    }
    split.finish().expect("finish");
    assert_eq!(
        names(&dir),
        ["part_0001.bam", "part_0002.bam", "part_0003.bam"]
    );
    let sizes: Vec<usize> = (1..=3)
        .map(|n| read_bam(&dir.join(format!("part_{n:04}.bam"))).1.len())
        .collect();
    // ceil(7 / 3) = 3 per file: `a a a`, then `b c c` (the second c
    // follows its name), then `d`.
    assert_eq!(sizes, [3, 3, 1]);
    let (text, _) = read_bam(&dir.join("part_0001.bam"));
    assert!(text.starts_with("@HD\tVN:1.6\tSO:unsorted\n"), "{text}");
    std::fs::remove_dir_all(&dir).expect("cleanup");
}
