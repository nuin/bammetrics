use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use super::*;
use crate::vcf::VcfLines;

fn scratch(name: &str) -> PathBuf {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-vcf-writer-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

const HEADER: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=c1,length=50000>\n##contig=<ID=c2,length=3000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

fn records() -> Vec<String> {
    let mut lines = Vec::new();
    for i in 0..3000 {
        lines.push(format!(
            "c1\t{}\t.\tA\tC\t{}\t.\tDP={i}\tGT\t0/1",
            1 + i * 15,
            i % 97
        ));
    }
    lines.push("c1\t46000\t.\tA\t<DEL>\t.\t.\tEND=46100\tGT\t0/1".to_owned());
    for i in 0..200 {
        lines.push(format!("c2\t{}\t.\tG\tT\t.\t.\t.\tGT\t1/1", 1 + i * 10));
    }
    lines
}

fn header() -> VcfHeader {
    VcfHeader::read(&mut HEADER.as_bytes()).expect("a header")
}

fn write(path: &Path, lines: &[String], options: VcfOptions) -> Result<Vec<PathBuf>> {
    let mut writer = VcfWriter::create(path, &header(), options)?;
    for line in lines {
        writer.write_line(line)?;
    }
    writer.finish()
}

fn read_text(path: &Path) -> String {
    let mut lines = VcfLines::open(path).expect("opens");
    let mut text = lines.header().text();
    let mut line = String::new();
    while lines.next_line(&mut line).expect("reads") {
        text.push_str(&line);
        text.push('\n');
    }
    text
}

#[test]
fn plain_and_compressed_files_hold_the_same_text() {
    let dir = scratch("same-text");
    let lines = records();
    let plain = dir.join("out.vcf");
    let gz = dir.join("out.vcf.gz");
    assert_eq!(
        write(&plain, &lines, VcfOptions::default()).expect("writes"),
        std::slice::from_ref(&plain)
    );
    let options = VcfOptions {
        index: true,
        md5: true,
        threads: 3,
        ..VcfOptions::default()
    };
    let written = write(&gz, &lines, options).expect("writes");
    assert_eq!(written, [gz.clone(), index_path(&gz), md5_path(&gz)]);
    let text = fs::read_to_string(&plain).expect("plain");
    assert!(text.starts_with("##fileformat=VCFv4.2\n##FORMAT=<ID=GT"));
    assert_eq!(read_text(&gz), text);
    let mut bytes = Vec::new();
    fs::File::open(&gz)
        .and_then(|mut f| f.read_to_end(&mut bytes))
        .expect("gz");
    let md5 = fs::read_to_string(md5_path(&gz)).expect("md5");
    assert_eq!(
        md5,
        crate::md5::hex(&{
            let mut m = crate::md5::Md5::new();
            m.update(&bytes);
            m.finish()
        })
    );
    fs::remove_dir_all(dir).expect("clean");
}

/// Every record found by a query of its own span through the index.
#[test]
fn the_index_finds_every_record() {
    use noodles::core::Region;
    let dir = scratch("tabix");
    let gz = dir.join("out.vcf.gz");
    let lines = records();
    let options = VcfOptions {
        index: true,
        ..VcfOptions::default()
    };
    write(&gz, &lines, options).expect("writes");
    let index = noodles::tabix::fs::read(index_path(&gz)).expect("reads the index");
    let mut reader = noodles::vcf::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(&gz)
        .expect("opens");
    let header = reader.read_header().expect("header");
    for (contig, start, end, expected) in [
        ("c1", 1, 50000, 3001),
        ("c1", 16_000, 16_500, 33),
        ("c1", 46_050, 46_060, 1),
        ("c2", 1, 3000, 200),
        ("c2", 1991, 1991, 1),
    ] {
        let region: Region = format!("{contig}:{start}-{end}").parse().expect("region");
        let found = reader
            .query(&header, &region)
            .expect("query")
            .records()
            .count();
        assert_eq!(found, expected, "{contig}:{start}-{end}");
    }
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn unsorted_records_fail_the_index_and_leave_no_file() {
    let dir = scratch("unsorted");
    let gz = dir.join("out.vcf.gz");
    let mut lines = records();
    lines.swap(10, 20);
    let options = VcfOptions {
        index: true,
        ..VcfOptions::default()
    };
    let error = write(&gz, &lines, options).expect_err("fails");
    assert!(error.to_string().contains("out of order"), "{error}");
    let mut lines = records();
    lines.push("c1\t1\t.\tA\tC\t.\t.\t.\tGT\t0".to_owned());
    let error = write(&gz, &lines, options).expect_err("fails");
    assert!(
        error.to_string().contains("Sequence c1 added out of order"),
        "{error}"
    );
    assert_eq!(fs::read_dir(&dir).expect("dir").count(), 0);
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn names_pick_the_kind() {
    for (name, kind) in [
        ("a.vcf", VcfFormat::Plain),
        ("a.vcf.gz", VcfFormat::Bgzf),
        ("a.VCF.GZ", VcfFormat::Bgzf),
        ("a.vcf.bgz", VcfFormat::Bgzf),
    ] {
        assert_eq!(VcfFormat::of(Path::new(name)).expect(name), kind, "{name}");
    }
    assert!(VcfFormat::of(Path::new("a.bcf")).is_err());
    assert!(VcfFormat::of(Path::new("a.txt")).is_err());
    let index = VcfOptions {
        index: true,
        ..VcfOptions::default()
    };
    assert!(index.check(Path::new("a.vcf")).is_err());
    assert!(index.check(Path::new("a.vcf.gz")).is_ok());
}

#[test]
fn a_bad_record_fails_and_leaves_no_file() {
    let dir = scratch("bad");
    let out = dir.join("out.vcf");
    let lines = vec!["c1\t1\t.\tA\tC\t.\tnofilter\t.\tGT\t0".to_owned()];
    let error = write(&out, &lines, VcfOptions::default()).expect_err("fails");
    assert!(error.to_string().contains("FILTER"), "{error}");
    let mut writer = VcfWriter::create(&out, &header(), VcfOptions::default()).expect("creates");
    writer
        .write_line("c1\t1\t.\tA\tC\t.\t.\tDP=1\tGT\t0")
        .expect("writes");
    drop(writer);
    assert_eq!(fs::read_dir(&dir).expect("dir").count(), 0);
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn text_is_written_as_iso_8859_1() {
    let dir = scratch("latin1");
    let out = dir.join("out.vcf");
    // An ID with U+00E9, U+4E2D, U+1F600 and the U+FFFD of a byte that was
    // not UTF-8: Picard writes E9 and a `?` for each of the others.
    let lines = vec!["c1\t1\tr\u{e9}\u{4e2d}\u{1f600}\u{fffd}\tA\tC\t.\t.\t.\tGT\t0".to_owned()];
    write(&out, &lines, VcfOptions::default()).expect("writes");
    let bytes = fs::read(&out).expect("vcf");
    let record = bytes
        .split(|&b| b == b'\n')
        .find(|line| line.starts_with(b"c1"))
        .expect("record");
    assert_eq!(record, b"c1\t1\tr\xe9???\tA\tC\t.\t.\t.\tGT\t0");
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn the_index_is_built_from_the_written_columns() {
    let line = b"c1\t10\t.\tACG\tA\t.\t.\tDP=1;XEND=9\tGT\t0";
    assert_eq!(tabix::record_span(line), Ok((&b"c1"[..], 10, 12)));
    let line = b"c1\t10\t.\tA\t<DEL>\t.\t.\tDP=1;END=40\t\t";
    assert_eq!(tabix::record_span(line), Ok((&b"c1"[..], 10, 40)));
    assert!(tabix::record_span(b"c1\t10\t.").is_err());
}

#[test]
fn a_record_past_the_tbi_windows_fails_the_index_and_leaves_no_file() {
    let dir = scratch("tbi-limit");
    let out = dir.join("out.vcf.gz");
    let options = VcfOptions {
        index: true,
        ..VcfOptions::default()
    };
    for line in [
        "c1\t1\t.\tA\t<DEL>\t.\t.\tEND=600000000\tGT\t0",
        "c2\t550000000\t.\tA\tC\t.\t.\t.\tGT\t0",
    ] {
        let error = write(&out, &[line.to_owned()], options).expect_err("fails");
        assert!(error.to_string().contains(".tbi"), "{error}");
        assert_eq!(fs::read_dir(&dir).expect("dir").count(), 0);
    }
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn only_names_htsjdk_types_are_written() {
    for (name, format) in [
        ("x.vcf", Some(VcfFormat::Plain)),
        ("x.vcf.gz", Some(VcfFormat::Bgzf)),
        ("x.Vcf.Gz", Some(VcfFormat::Bgzf)),
        ("x.gz", Some(VcfFormat::Bgzf)),
        ("x.vcf.bgzf", Some(VcfFormat::Bgzf)),
        ("x.gzip", Some(VcfFormat::Bgzf)),
        ("x.VCF", None),
        ("x.txt", None),
        ("x", None),
        ("-", None),
        ("x.vcf.", None),
        ("x.vcf.gz.tbi", None),
        ("x.vcf.idx", None),
        ("x.bam", None),
        ("x.bcf", None),
    ] {
        assert_eq!(VcfFormat::of(Path::new(name)).ok(), format, "{name}");
    }
}
