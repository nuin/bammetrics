use std::os::unix::fs::DirBuilderExt;

use super::*;

/// A directory of its own under the system temporary directory.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-bgzf-copy-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

/// The inflated BAM header: magic, text, one reference `c1` of 1000 bases.
fn header(text: &[u8]) -> Vec<u8> {
    let mut bytes = b"BAM\x01".to_vec();
    bytes.extend_from_slice(&(text.len() as i32).to_le_bytes());
    bytes.extend_from_slice(text);
    bytes.extend_from_slice(&1i32.to_le_bytes());
    bytes.extend_from_slice(&3i32.to_le_bytes());
    bytes.extend_from_slice(b"c1\0");
    bytes.extend_from_slice(&1000i32.to_le_bytes());
    bytes
}

/// A file of the given blocks, each written as one stored block, then the
/// end-of-file block when `eof`.
fn file_of(dir: &Path, name: &str, blocks: &[&[u8]], eof: bool) -> PathBuf {
    let mut bytes = Vec::new();
    for block in blocks {
        if block.is_empty() {
            bytes.extend_from_slice(&EOF_BLOCK);
        } else {
            write_blocks(&mut bytes, block, 0).expect("block");
        }
    }
    if eof {
        bytes.extend_from_slice(&EOF_BLOCK);
    }
    let path = dir.join(name);
    std::fs::write(&path, bytes).expect("write");
    path
}

/// The raw blocks of `bytes`.
fn raw_blocks(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut reader = bytes;
    let mut raw = Vec::new();
    let mut out = Vec::new();
    while read_raw_block(&mut reader, &mut raw).expect("block") {
        out.push(raw.clone());
    }
    out
}

fn copy(path: &Path, skip_header: bool, skip_terminator: bool, level: u8) -> Vec<u8> {
    let source = Source::open(path).expect("source");
    let mut out = Vec::new();
    source
        .copy(&mut out, skip_header, skip_terminator, level)
        .expect("copy");
    out
}

#[test]
fn the_end_of_a_file_is_told_as_htsjdk_tells_it() {
    let dir = scratch("termination");
    let check = |bytes: &[u8]| {
        let path = dir.join("t");
        std::fs::write(&path, bytes).expect("write");
        let mut file = File::open(&path).expect("open");
        termination(&mut file, bytes.len() as u64).expect("read")
    };
    let mut block = Vec::new();
    write_blocks(&mut block, b"some bytes", 0).expect("block");
    assert_eq!(check(&EOF_BLOCK), Termination::Terminator);
    assert_eq!(check(&block), Termination::HealthyLastBlock);
    assert_eq!(check(&block[..block.len() - 1]), Termination::Defective);
    assert_eq!(check(&EOF_BLOCK[..27]), Termination::Defective);
    assert_eq!(check(&[b'x'; 100]), Termination::Defective);
    let mut two = block.clone();
    two.extend_from_slice(&block);
    assert_eq!(check(&two), Termination::HealthyLastBlock);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_first_record_follows_the_header_as_get_file_pointer_gives_it() {
    let dir = scratch("start");
    let head = header(b"@HD\tVN:1.6\n");
    // The header alone in its block: the next block at offset 0.
    let alone = file_of(&dir, "alone.bam", &[&head, b"records"], true);
    let start = Source::open(&alone).expect("open").start;
    let first_block = 18 + 5 + head.len() + 8;
    assert_eq!(start.first_record, (first_block as u64, 0));
    assert_eq!(start.references, vec![(b"c1".to_vec(), 1000)]);
    assert!(start.has_records);
    let empty = file_of(&dir, "empty.bam", &[&head], true);
    assert!(!Source::open(&empty).expect("open").start.has_records);
    // The header and records in one block.
    let mut shared = head.clone();
    shared.extend_from_slice(b"records");
    let shared = file_of(&dir, "shared.bam", &[&shared], true);
    let start = Source::open(&shared).expect("open").start;
    assert_eq!(start.first_record, (0, head.len()));
    // The header over two blocks.
    let (a, b) = head.split_at(10);
    let split = file_of(&dir, "split.bam", &[a, b, b"r"], true);
    let start = Source::open(&split).expect("open").start;
    assert_eq!(start.first_record.1, 0);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_block_of_the_first_record_is_deflated_again_and_the_rest_copied() {
    let dir = scratch("copy");
    let head = header(b"@HD\tVN:1.6\n");
    let path = file_of(&dir, "a.bam", &[&head, b"first", b"second", b"third"], true);
    let bytes = std::fs::read(&path).expect("read");
    let input = raw_blocks(&bytes);
    // Level 0: the deflated block is the stored block htsjdk writes.
    let out = copy(&path, true, true, 0);
    assert_eq!(raw_blocks(&out), input[1..4].to_vec());
    let out = copy(&path, true, false, 0);
    assert_eq!(raw_blocks(&out), input[1..].to_vec());
    // Level 5: the same bytes inflated, the later blocks as they are.
    let out = copy(&path, true, true, 5);
    let blocks = raw_blocks(&out);
    let mut first = Vec::new();
    inflate(&blocks[0], &mut first).expect("inflate");
    assert_eq!(first, b"first");
    assert_eq!(blocks[1..], input[2..4]);
    // No header skip: every block, the end-of-file block left out.
    assert_eq!(copy(&path, false, true, 5), bytes[..bytes.len() - 28]);
    assert_eq!(copy(&path, false, false, 5), bytes);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn records_sharing_the_header_block_are_deflated_from_their_offset() {
    let dir = scratch("shared");
    let mut shared = header(b"@HD\tVN:1.6\n");
    shared.extend_from_slice(b"records");
    let path = file_of(&dir, "a.bam", &[&shared, b"more"], true);
    let input = raw_blocks(&std::fs::read(&path).expect("read"));
    let blocks = raw_blocks(&copy(&path, true, true, 0));
    let mut first = Vec::new();
    inflate(&blocks[0], &mut first).expect("inflate");
    assert_eq!(first, b"records");
    assert_eq!(blocks[1], input[1]);
    assert_eq!(blocks.len(), 2);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_header_only_bam_gives_nothing_after_its_header() {
    let dir = scratch("header-only");
    let path = file_of(&dir, "a.bam", &[&header(b"")], true);
    assert!(copy(&path, true, true, 5).is_empty());
    assert!(copy(&path, true, false, 5).is_empty());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn an_empty_block_after_the_header_is_passed_over() {
    let dir = scratch("empty-block");
    let head = header(b"");
    let path = file_of(&dir, "a.bam", &[&head, b"", b"data", b"tail"], true);
    let input = raw_blocks(&std::fs::read(&path).expect("read"));
    let blocks = raw_blocks(&copy(&path, true, true, 0));
    assert_eq!(blocks, vec![input[2].clone(), input[3].clone()]);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn bad_files_fail_before_and_during_the_copy() {
    let dir = scratch("bad");
    let head = header(b"");
    let fails = |bytes: &[u8], wanted: &str| {
        let path = dir.join("bad.bam");
        std::fs::write(&path, bytes).expect("write");
        let error = Source::open(&path)
            .and_then(|source| source.copy(&mut io::sink(), true, true, 5))
            .expect_err(wanted)
            .to_string();
        assert!(error.contains(wanted), "{error}");
    };
    let good = std::fs::read(file_of(&dir, "good.bam", &[&head, b"r"], true)).expect("read");
    fails(&good[..good.len() - 30], "valid GZIP block at the end");
    let mut not_bam = Vec::new();
    write_blocks(&mut not_bam, b"SAM\x01 text", 0).expect("block");
    not_bam.extend_from_slice(&EOF_BLOCK);
    fails(&not_bam, "Invalid BAM file header");
    // A text length larger than the file: the read ends, no buffer is sized.
    let mut huge = b"BAM\x01".to_vec();
    huge.extend_from_slice(&i32::MAX.to_le_bytes());
    let mut file = Vec::new();
    write_blocks(&mut file, &huge, 0).expect("block");
    file.extend_from_slice(&EOF_BLOCK);
    fails(&file, "Premature EOF");
    let mut negative = b"BAM\x01".to_vec();
    negative.extend_from_slice(&(-1i32).to_le_bytes());
    let mut file = Vec::new();
    write_blocks(&mut file, &negative, 0).expect("block");
    file.extend_from_slice(&EOF_BLOCK);
    fails(&file, "claims -1 bytes");
    // A byte that starts no block between two blocks.
    let mut garbage = Vec::new();
    write_blocks(&mut garbage, &head, 0).expect("block");
    write_blocks(&mut garbage, b"records", 0).expect("block");
    garbage.extend_from_slice(b"junk");
    write_blocks(&mut garbage, b"more", 0).expect("block");
    garbage.extend_from_slice(&EOF_BLOCK);
    fails(&garbage, "no BGZF block starts here");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn blocks_hold_at_most_htsjdks_block_size() {
    let data: Vec<u8> = (0..200_000u32).map(|i| (i * 7 % 251) as u8).collect();
    for level in [0, 1, 5, 9] {
        let mut out = Vec::new();
        let count = write_blocks(&mut out, &data, level).expect("blocks");
        assert_eq!(count, 4, "level {level}");
        let mut inflated = Vec::new();
        let mut all = Vec::new();
        for block in raw_blocks(&out) {
            inflate(&block, &mut inflated).expect("inflate");
            assert!(inflated.len() <= BLOCK_DATA_BYTES);
            all.extend_from_slice(&inflated);
        }
        assert_eq!(all, data);
    }
}
