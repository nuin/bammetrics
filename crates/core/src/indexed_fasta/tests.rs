use super::testing::{scratch, write_fasta};
use super::*;

fn indexed(path: &Path) -> IndexedFasta {
    match open(path).expect("opens") {
        Opened::Indexed(fasta) => fasta,
        Opened::Unindexed { .. } => panic!("{} is not indexed", path.display()),
    }
}

fn bases(fasta: &mut IndexedFasta, index: usize) -> Vec<u8> {
    let mut all = Vec::new();
    fasta
        .walk(index, |run| all.extend_from_slice(run))
        .expect("walks");
    all
}

#[test]
fn a_walk_gives_the_bases_at_the_index_geometry() {
    let dir = scratch("walk");
    let long: Vec<u8> = (0..CHUNK + 1234).map(|i| b"ACGTN"[i % 5]).collect();
    for (name, width, end) in [
        ("lf.fa", 7, "\n"),
        ("crlf.fasta", 60, "\r\n"),
        ("one.fna", 5000, "\n"),
    ] {
        let path = write_fasta(
            &dir,
            name,
            &[("c1", b"ACGTNNacgtRYKM.-x"), ("c2", &long), ("c3", b"A")],
            width,
            end,
        );
        let mut fasta = indexed(&path);
        assert_eq!(fasta.entries().len(), 3);
        assert_eq!(fasta.dictionary().expect("dict").sequences.len(), 3);
        assert_eq!(bases(&mut fasta, 0), b"ACGTNNacgtRYKM.-x", "{name}");
        assert_eq!(bases(&mut fasta, 1), long, "{name}");
        assert_eq!(bases(&mut fasta, 2), b"A", "{name}");
        assert_eq!(bases(&mut fasta, 0), b"ACGTNNacgtRYKM.-x", "{name} again");
    }
}

#[test]
fn a_bgzf_fasta_is_read_through_its_gzi() {
    let dir = scratch("bgzf");
    let plain = write_fasta(
        &dir,
        "p.fa",
        &[("c1", b"ACGTNNNNacgt"), ("c2", b"NNA")],
        5,
        "\n",
    );
    let text = std::fs::read(&plain).expect("FASTA");
    let gz = dir.join("p.fa.gz");
    let mut writer = noodles::bgzf::io::Writer::new(File::create(&gz).expect("create"));
    std::io::Write::write_all(&mut writer, &text).expect("write");
    writer.finish().expect("finish");
    std::fs::copy(fai_path(&plain), fai_path(&gz)).expect(".fai");
    // The .dict of p.fa.gz is p.dict, the one p.fa has.
    assert!(matches!(
        open(&gz).expect("opens"),
        Opened::Unindexed { .. }
    ));
    std::fs::write(with_suffix(&gz, ".gzi"), 0u64.to_le_bytes()).expect(".gzi");
    let mut fasta = indexed(&gz);
    assert_eq!(bases(&mut fasta, 0), b"ACGTNNNNacgt");
    assert_eq!(bases(&mut fasta, 1), b"NNA");
}

#[test]
fn plain_gzip_and_unsupported_names_are_refused() {
    let dir = scratch("refused");
    let plain = write_fasta(&dir, "g.fa", &[("c1", b"ACGT")], 60, "\n");
    let gz = dir.join("g.fa.gz");
    std::fs::write(&gz, stored_gzip(&std::fs::read(&plain).expect("FASTA"))).expect("gzip");
    std::fs::copy(fai_path(&plain), fai_path(&gz)).expect(".fai");
    let error = open(&gz).err().expect("refused").to_string();
    assert!(error.contains("BGZF"), "{error}");
    let other = dir.join("g.seq");
    std::fs::copy(&plain, &other).expect("copy");
    let error = open(&other).err().expect("refused").to_string();
    assert!(
        error.contains("not a supported reference file type"),
        "{error}"
    );
}

/// A gzip member holding `bytes` in one stored block: gzip, not BGZF.
fn stored_gzip(bytes: &[u8]) -> Vec<u8> {
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    out.push(1);
    let length = bytes.len() as u16;
    out.extend_from_slice(&length.to_le_bytes());
    out.extend_from_slice(&(!length).to_le_bytes());
    out.extend_from_slice(bytes);
    out.extend_from_slice(&crc32(bytes).to_le_bytes());
    out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
    out
}

fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &b in bytes {
        crc ^= u32::from(b);
        for _ in 0..8 {
            crc = if crc & 1 != 0 {
                (crc >> 1) ^ 0xedb8_8320
            } else {
                crc >> 1
            };
        }
    }
    !crc
}

#[test]
fn the_dictionary_must_match_the_index() {
    let dir = scratch("sanity");
    let path = write_fasta(&dir, "a.fa", &[("c1", b"ACGT"), ("c2", b"NN")], 60, "\n");
    let dict = default_dictionary(&path).expect("name");
    for (text, wanted) in [
        ("@SQ\tSN:c1\tLN:4\n", "different numbers of contigs"),
        (
            "@SQ\tSN:c2\tLN:2\n@SQ\tSN:c1\tLN:4\n",
            "Mismatch between sequence dictionary",
        ),
        (
            "@SQ\tSN:c1\tLN:5\n@SQ\tSN:c2\tLN:2\n",
            "Index length does not match",
        ),
    ] {
        std::fs::write(&dict, text).expect(".dict");
        let error = open(&path).err().expect("refused").to_string();
        assert!(error.contains(wanted), "{text:?}: {error}");
    }
    std::fs::remove_file(&dict).expect("no .dict");
    let appended = with_suffix(&path, ".dict");
    std::fs::write(&appended, "@SQ\tSN:c1\tLN:4\n@SQ\tSN:c2\tLN:2\n").expect("a.fa.dict");
    assert_eq!(
        dictionary_path(&path).expect("name"),
        Some(appended.clone())
    );
    assert!(indexed(&path).dictionary().is_some());
    std::fs::remove_file(&appended).expect("no .dict");
    assert!(indexed(&path).dictionary().is_none());
    std::fs::remove_file(fai_path(&path)).expect("no .fai");
    assert!(matches!(
        open(&path).expect("opens"),
        Opened::Unindexed { dictionary: None }
    ));
}

#[test]
fn index_lines_are_read_as_htsjdk_reads_them() {
    let dir = scratch("fai");
    let fai = dir.join("x.fa.fai");
    std::fs::write(&fai, "c1 extra\t10\t4\t5\t6\r\nc2\t3\t20\t60\t61\n\n  \n").expect(".fai");
    let entries = read_index(&fai).expect("reads");
    assert_eq!(entries[0].name, "c1");
    assert_eq!(entries[0].length, 10);
    assert_eq!(entries[1].line_width, 61);
    for bad in [
        "c1\t10\t4\t5\n",
        "c1\t-1\t4\t5\t6\n",
        "c1\t1\t4\t5\t99999999999\n",
        "c1\t1\t4\t5\t6\nc1\t1\t4\t5\t6\n",
    ] {
        std::fs::write(&fai, bad).expect(".fai");
        assert!(read_index(&fai).is_err(), "{bad:?}");
    }
}

/// The six shapes of the 0.29.0 review, each as Picard 3.5.0 reads it:
/// htsjdk finds `(.+)\t+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)` in each line of a
/// `Scanner`.
#[test]
fn index_lines_follow_the_scanner_pattern() {
    let dir = scratch("fai-scanner");
    let fai = dir.join("x.fa.fai");
    let read = |text: &str| {
        std::fs::write(&fai, text).expect(".fai");
        read_index(&fai)
    };
    for good in [
        "c1\t10\t4\t5\t6\t\nc2\t3\t20\t60\t61\t\n",
        "c1\t10\t4\t5\t6  \nc2\t3\t20\t60\t61  \n",
        "c1\t\t10\t\t4\t\t5\t\t6\nc2\t\t3\t\t20\t\t60\t\t61\n",
        "c1\t10\t4\t5\t6x\nc2\t3\t20\t60\t61x\n",
        "c1\t10\t4\t5\t6\rc2\t3\t20\t60\t61\r",
        "c1\t10 4 5 6\nc2\t3\t20\t60\t61x",
    ] {
        let entries = read(good).unwrap_or_else(|e| panic!("{good:?}: {e}"));
        assert_eq!(entries.len(), 2, "{good:?}");
        assert_eq!(entries[0].length, 10, "{good:?}");
        assert_eq!(entries[1].line_width, 61, "{good:?}");
    }
    let longer = read("c1\t10\t4\t5\t6\t7\n").expect("a sixth number");
    assert_eq!((longer[0].length, longer[0].line_width), (4, 7));
    for (bad, message) in [
        ("c1\t10\t4\t5\t6\nc2\t3\t20\t60\t61", "No line found"),
        ("\nc1\t10\t4\t5\t6\n", "Found invalid line in index file:"),
        (
            "c1\t10\t4\t5\t6\n \nc2\t3\t20\t60\t61\n",
            "invalid line in index file: \n",
        ),
        ("c1 10\t4\t5\t6\n", "invalid line"),
    ] {
        let error = read(bad).expect_err(bad).to_string();
        assert!(error.contains(message.trim_end()), "{bad:?}: {error}");
    }
}

#[test]
fn a_short_file_and_a_zero_line_length_fail() {
    let dir = scratch("short");
    let path = write_fasta(&dir, "s.fa", &[("c1", b"ACGTACGTAC")], 4, "\n");
    let text = std::fs::read(&path).expect("FASTA");
    std::fs::write(&path, &text[..text.len() - 4]).expect("cut");
    let mut fasta = indexed(&path);
    let error = fasta.walk(0, |_| {}).expect_err("fails").to_string();
    assert!(error.contains("ends before"), "{error}");
    std::fs::write(fai_path(&path), "c1\t0\t4\t0\t0\n").expect(".fai");
    std::fs::write(
        default_dictionary(&path).expect("name"),
        "@SQ\tSN:c1\tLN:0\n",
    )
    .expect(".dict");
    let mut fasta = indexed(&path);
    let error = fasta.walk(0, |_| {}).expect_err("fails").to_string();
    assert!(error.contains("/ by zero"), "{error}");
}

/// A named pipe as the FASTA or the `.dict` is refused before it is
/// opened: an open with no writer would block.
#[cfg(unix)]
#[test]
fn a_named_pipe_is_refused_before_it_is_opened() {
    let mkfifo = |path: &Path| {
        let made = std::process::Command::new("mkfifo")
            .arg(path)
            .status()
            .expect("mkfifo runs");
        assert!(made.success());
    };
    let dir = scratch("fifo");
    let path = write_fasta(&dir, "p.fa", &[("c1", b"ACGT")], 4, "\n");
    let dict = default_dictionary(&path).expect("name");
    std::fs::remove_file(&dict).expect("no .dict");
    mkfifo(&dict);
    let error = open(&path).err().expect("a pipe .dict").to_string();
    assert!(error.contains("not a regular file"), "{error}");
    std::fs::remove_file(&dict).expect("no pipe");
    let fifo = dir.join("q.fa");
    mkfifo(&fifo);
    std::fs::copy(fai_path(&path), fai_path(&fifo)).expect(".fai");
    let error = open(&fifo).err().expect("a pipe FASTA").to_string();
    assert!(error.contains("not a regular file"), "{error}");
}
