//! The two runners of [`super`] on a synthetic edit and the fixture files.

use std::path::{Path, PathBuf};

use super::*;
use crate::bam_raw::{TagValue, build};
use crate::stage::Stage;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data")
        .join(name)
}

/// A scratch directory of this test, readable by this user only.
fn scratch(name: &str) -> PathBuf {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-record-edit-{}-{name}",
        std::process::id()
    ));
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

/// Adds `@CO edited`, sets `XN:i:<n>` on the n-th record, drops every
/// third record and writes one record `tail` after the last.
struct Numbering {
    seen: i64,
    mode: Mode,
}

impl RecordEdit for Numbering {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        let mut header = input.clone();
        header.add_comment("edited");
        Ok(header)
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        self.seen += 1;
        if self.seen % 3 == 0 {
            return Ok(Outcome::Dropped);
        }
        let mut tags = record.tags().map_err(Error::Record)?;
        tags.set(*b"XN", TagValue::Int(self.seen));
        record.set_tags(&tags).map_err(Error::Record)?;
        Ok(Outcome::Changed)
    }

    fn finish(&mut self, out: &mut Vec<RawRecord>) -> Result<()> {
        let tail = match self.mode {
            Mode::Chain => "tail-chain",
            _ => "tail",
        };
        out.push(build(tail, 4, -1, -1, &[], &[], b""));
        Ok(())
    }
}

fn numbering(setup: &Setup) -> Result<Box<dyn RecordEdit>> {
    Ok(Box::new(Numbering {
        seen: 0,
        mode: setup.mode,
    }))
}

fn options(input: &Path, output: &Path) -> Options {
    Options {
        input: input.to_path_buf(),
        output: output.to_path_buf(),
        reference: None,
        reads: Vec::new(),
        compression_level: 1,
        create_index: false,
        create_md5: false,
        threads: 1,
        stringency: Stringency::Strict,
        mode: Mode::Rewrite,
    }
}

/// The header and records of a BAM or SAM file.
fn read_back(path: &Path) -> (sam::Header, Vec<RawRecord>) {
    pass::with_record_stream_at(path, None, 1, Stringency::Silent, |header, reader, _| {
        let mut records = Vec::new();
        let mut record = RawRecord::default();
        while record.read_from(reader).map_err(Error::Record)? {
            records.push(record.clone());
        }
        Ok((header.clone(), records))
    })
    .expect("read back")
}

fn count(path: &Path) -> usize {
    read_back(path).1.len()
}

#[test]
fn a_run_edits_the_header_and_each_record_and_adds_the_tail() {
    let dir = scratch("run");
    let input = data("insert_size.bam");
    let n = count(&input);
    for name in ["out.bam", "out.sam"] {
        let out = dir.join(name);
        let counters = run(&options(&input, &out), numbering).expect("run");
        assert_eq!(counters.records_in, n as u64);
        assert_eq!(counters.records_out, (n - n / 3 + 1) as u64);
        assert_eq!(counters.changed, (n - n / 3) as u64);
        let (header, records) = read_back(&out);
        assert_eq!(
            header.comments().last().map(|c| c.to_string()),
            Some("edited".to_owned())
        );
        assert_eq!(records.len(), n - n / 3 + 1);
        assert_eq!(records[2].tag(*b"XN").expect("tag"), Some(TagValue::Int(4)));
        assert_eq!(records.last().expect("tail").name(), b"tail");
    }
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn a_block_copy_keeps_record_bytes_and_names_the_index_after_the_bam() {
    let dir = scratch("copy");
    let input = data("insert_size.bam");
    let out = dir.join("out.bam");
    let opts = Options {
        mode: Mode::BlockCopy,
        create_index: true,
        create_md5: true,
        ..options(&input, &out)
    };
    let keep = |_: &Setup| -> Result<Box<dyn RecordEdit>> { Ok(Box::new(Keep)) };
    run(&opts, keep).expect("run");
    assert_eq!(read_back(&out).1, read_back(&input).1);
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .expect("dir")
        .map(|entry| {
            entry
                .expect("entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    names.sort();
    assert_eq!(names, ["out.bam", "out.bam.bai", "out.bam.md5"]);
    std::fs::remove_dir_all(&dir).expect("scratch");
}

/// Hands every record on as it came.
struct Keep;

impl RecordEdit for Keep {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        Ok(input.clone())
    }

    fn record(&mut self, _record: &mut RawRecord) -> Result<Outcome> {
        Ok(Outcome::Same)
    }
}

fn keep(_: &Setup) -> Result<Box<dyn RecordEdit>> {
    Ok(Box::new(Keep))
}

#[test]
fn a_block_copy_refuses_a_cram_and_a_sam_name_and_leaves_no_file() {
    let dir = scratch("refuse");
    let out = dir.join("out.bam");
    let cram = Options {
        mode: Mode::BlockCopy,
        reference: Some(data("ref_100k.fa")),
        ..options(&data("insert_size.cram"), &out)
    };
    let message = run(&cram, keep).expect_err("a CRAM").to_string();
    assert!(message.contains("a CRAM: this tool copies"), "{message}");
    let sam = Options {
        mode: Mode::BlockCopy,
        ..options(&data("insert_size.bam"), &dir.join("out.sam"))
    };
    let message = run(&sam, keep).expect_err("a .sam name").to_string();
    assert!(message.contains("name the OUTPUT .bam"), "{message}");
    assert_eq!(std::fs::read_dir(&dir).expect("dir").count(), 0);
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn an_output_that_is_an_input_is_refused_before_the_first_read() {
    let dir = scratch("same");
    let input = dir.join("in.bam");
    std::fs::copy(data("insert_size.bam"), &input).expect("copy");
    let before = std::fs::read(&input).expect("read");
    let opts = Options {
        reads: vec![data("markdup.bam")],
        ..options(&input, &input)
    };
    let message = run(&opts, keep).expect_err("same file").to_string();
    assert!(message.contains("writing it would destroy it"), "{message}");
    let header = Options {
        reads: vec![dir.join("in.bai")],
        create_index: true,
        ..options(&data("insert_size.bam"), &input)
    };
    std::fs::write(dir.join("in.bai"), b"x").expect("index");
    let message = run(&header, keep).expect_err("index is read").to_string();
    assert!(message.contains("in.bai"), "{message}");
    assert_eq!(std::fs::read(&input).expect("read"), before);
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn a_cram_output_names_its_wave() {
    let dir = scratch("cram-out");
    let opts = options(&data("insert_size.bam"), &dir.join("out.cram"));
    let message = run(&opts, keep).expect_err("CRAM output").to_string();
    assert!(message.contains("CRAM output is wave"), "{message}");
    std::fs::remove_dir_all(&dir).expect("scratch");
}

#[test]
fn the_stage_edits_a_batch_in_place_and_hands_the_header_down() {
    let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
        .parse()
        .expect("header");
    let frame = Frame::new(header.clone(), "in.bam".to_owned());
    let inv = crate::sort_sam::TOOL
        .parse(&["I=in.bam".to_owned(), "SORT_ORDER=coordinate".to_owned()])
        .expect("parse");
    let Stage::Transform(mut transform) = stage(&inv, &frame, numbering).expect("stage") else {
        panic!("a record edit is a transform");
    };
    let out = transform.header(&header).expect("header");
    assert_eq!(out.comments().len(), 1);
    let mut batch: Vec<RawRecord> = (0..7)
        .map(|i| build(&format!("r{i}"), 0, 0, i, &[], &[], b""))
        .collect();
    transform.apply(&mut batch).expect("apply");
    let names: Vec<&[u8]> = batch.iter().map(RawRecord::name).collect();
    assert_eq!(names, [b"r0", b"r1", b"r3", b"r4", b"r6"]);
    let mut tail = Vec::new();
    transform.finish(&mut tail).expect("finish");
    assert_eq!(tail[0].name(), b"tail-chain");
    let counters = transform.counters();
    assert_eq!(
        (counters.records_in, counters.records_out, counters.changed),
        (7, 6, 5)
    );
}

/// Keeps every record; `rewrites` as given.
struct Pass(bool);

impl RecordEdit for Pass {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        Ok(input.clone())
    }

    fn record(&mut self, _record: &mut RawRecord) -> Result<Outcome> {
        Ok(Outcome::Same)
    }

    fn rewrites(&self) -> bool {
        self.0
    }
}

#[test]
fn a_rewriting_step_gives_each_record_the_writers_bin_and_version() {
    let header: sam::Header = "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
        .parse()
        .expect("header");
    let frame = Frame::new(header.clone(), "in.bam".to_owned());
    let inv = crate::sort_sam::TOOL
        .parse(&["I=in.bam".to_owned(), "SORT_ORDER=coordinate".to_owned()])
        .expect("parse");
    // An unplaced record read with bin 4680: htsjdk's writer gives it 0.
    let mut bytes = build("u", 4, -1, -1, &[], &[], b"").as_bytes().to_vec();
    bytes[10..12].copy_from_slice(&4680u16.to_le_bytes());
    let stale = RawRecord::from_bytes(bytes).expect("record");
    fn rewriting(_: &Setup) -> Result<Box<dyn RecordEdit>> {
        Ok(Box::new(Pass(true)))
    }
    fn copying(_: &Setup) -> Result<Box<dyn RecordEdit>> {
        Ok(Box::new(Pass(false)))
    }
    type Build = fn(&Setup) -> Result<Box<dyn RecordEdit>>;
    for (rewrites, edit, bin) in [(true, rewriting as Build, 0u16), (false, copying, 4680)] {
        let Stage::Transform(mut transform) = stage(&inv, &frame, edit).expect("stage") else {
            panic!("a record edit is a transform");
        };
        let mut batch = vec![stale.clone()];
        transform.apply(&mut batch).expect("apply");
        let got = u16::from_le_bytes([batch[0].as_bytes()[10], batch[0].as_bytes()[11]]);
        assert_eq!(got, bin, "rewrites {rewrites}");
        assert_eq!(transform.counters().changed, 0);
        // htsjdk's writer writes VN:1.6; a block copy keeps the version read.
        let version = transform
            .header(&header)
            .expect("header")
            .header()
            .map(|hd| hd.version());
        let expected = if rewrites { (1, 6) } else { (1, 4) };
        assert_eq!(version, Some(Version::new(expected.0, expected.1)));
    }
}

#[test]
fn a_step_refuses_the_options_of_the_file_it_does_not_write() {
    for given in [
        "CREATE_INDEX=true",
        "COMPRESSION_LEVEL=1",
        "CREATE_MD5_FILE=true",
        "O=x.bam",
    ] {
        let inv = crate::sort_sam::TOOL
            .parse(&[
                "I=in.bam".to_owned(),
                "SORT_ORDER=coordinate".to_owned(),
                given.to_owned(),
            ])
            .expect("parse");
        let message = check_stage(&inv).expect_err(given).to_string();
        assert!(message.contains("belongs to the write step"), "{message}");
    }
}

/// An `SO` value htsjdk's `SortOrder` does not name: an error at `STRICT`,
/// `unknown` in the header written at the other levels, in both modes.
#[test]
fn a_sort_order_htsjdk_does_not_name_is_written_unknown() {
    let header: sam::Header = "@HD\tVN:1.6\tSO:weird\n@SQ\tSN:c1\tLN:100\n"
        .parse()
        .expect("header");
    let path = Path::new("in.bam");
    let message = check_sort_order(&header, Stringency::Strict, path)
        .expect_err("strict")
        .to_string();
    assert!(
        message.ends_with(
            "Found non-conforming header SO tag: weird, exiting because \
             VALIDATION_STRINGENCY=STRICT"
        ),
        "{message}"
    );
    check_sort_order(&header, Stringency::Silent, path).expect("silent");
    for rewrite in [false, true] {
        let text = crate::bam_writer::header_bytes(&with_version(header.clone(), rewrite));
        assert!(text.starts_with(b"@HD\tVN:1.6\tSO:unknown\n"), "{rewrite}");
    }
    let named: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n".parse().expect("header");
    check_sort_order(&named, Stringency::Strict, path).expect("a named order");
}
