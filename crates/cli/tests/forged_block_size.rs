//! A BAM record whose `block_size` claims far more bytes than the file
//! holds, on the shared one-pass reader of `crates/core/src/pass.rs`. Two
//! BAMs under 300 bytes, one record each: `block_size` 0x7fffff00 and
//! -5 (0xfffffffb), 100 bytes behind it. Before `bam_raw::BoundedRecords`,
//! noodles sized its record buffer from the field: 2,052 MiB and 4,100 MiB
//! peak RSS in every tool below. Each run must exit 1 with `truncated BAM`
//! under 64 MiB. A CRAM container header claiming 0x7fffff00 bytes (the
//! `bammetrics patch` of `vendor/noodles-cram` `read_container`) must exit 1
//! with `truncated CRAM` under the same limit; it cost 1,578 MiB before.
//! SAM text needs no such test: its records are built from parsed lines.

use std::fs::{self, DirBuilder};
use std::io::Read;
use std::os::unix::fs::DirBuilderExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

const TOOLS: &[&str] = &[
    "adapter-content",
    "telomere",
    "duplicate-by-tile",
    "insert-size",
    "quality-yield",
    "sex-check",
];

const LIMIT: i64 = 64 << 20;

/// The 28-byte BGZF end-of-file block of SAM specification 4.1.2.
const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, b'B', b'C', 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

fn crc32(data: &[u8]) -> u32 {
    let mut crc = !0u32;
    for byte in data {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            crc = if crc & 1 == 1 {
                (crc >> 1) ^ 0xedb8_8320
            } else {
                crc >> 1
            };
        }
    }
    !crc
}

/// `data` as one BGZF block holding one stored (uncompressed) deflate block.
fn bgzf(data: &[u8]) -> Vec<u8> {
    let len = u16::try_from(data.len()).expect("one block");
    let total = u16::try_from(18 + 5 + data.len() + 8).expect("one block");
    let mut block = vec![
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, b'B', b'C', 0x02, 0,
    ];
    block.extend((total - 1).to_le_bytes());
    block.push(0x01);
    block.extend(len.to_le_bytes());
    block.extend((!len).to_le_bytes());
    block.extend(data);
    block.extend(crc32(data).to_le_bytes());
    block.extend(u32::from(len).to_le_bytes());
    block
}

/// A coordinate-sorted BAM with one reference and one record whose
/// `block_size` is `size`, followed by 100 bytes and the end-of-file block.
fn forged(size: u32) -> Vec<u8> {
    let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n";
    let mut data = b"BAM\x01".to_vec();
    data.extend(u32::try_from(text.len()).unwrap().to_le_bytes());
    data.extend(text);
    data.extend(1u32.to_le_bytes());
    data.extend(5u32.to_le_bytes());
    data.extend(b"chr1\0");
    data.extend(100_000u32.to_le_bytes());
    data.extend(size.to_le_bytes());
    data.extend([0u8; 100]);
    let mut bam = bgzf(&data);
    bam.extend(EOF_BLOCK);
    bam
}

/// A new directory per call. The two tests run on two threads of one
/// process and may read the same clock value, so a counter tells them apart.
fn fresh_dir() -> PathBuf {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let call = CALLS.fetch_add(1, Ordering::Relaxed);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-forged-block-size-{}-{nanos}-{call}",
        std::process::id()
    ));
    DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("make the temp dir");
    dir
}

/// Exit code, standard error and peak RSS in bytes of one run.
#[expect(
    clippy::zombie_processes,
    reason = "libc::wait4 reaps the child to get its rusage"
)]
fn measured(args: &[&str]) -> (i32, String, i64) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run bammetrics");
    let mut stderr = String::new();
    child
        .stderr
        .take()
        .expect("stderr")
        .read_to_string(&mut stderr)
        .expect("read stderr");
    let pid = libc::pid_t::try_from(child.id()).expect("pid");
    let mut status = 0;
    // SAFETY: an all-zero `rusage` is a valid value of the C struct.
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    // SAFETY: `pid` is this process's child, not yet waited for, and both
    // pointers are to live locals.
    let waited = unsafe { libc::wait4(pid, &mut status, 0, &mut usage) };
    assert_eq!(waited, pid, "wait4");
    assert!(libc::WIFEXITED(status), "{args:?} did not exit: {status}");
    // ru_maxrss is bytes on macOS and kilobytes on Linux.
    let peak = if cfg!(target_os = "macos") {
        usage.ru_maxrss
    } else {
        usage.ru_maxrss * 1024
    };
    (libc::WEXITSTATUS(status), stderr, peak)
}

#[test]
fn a_block_size_past_the_file_costs_the_bytes_present() {
    let dir = fresh_dir();
    for (name, size) in [
        ("past", 0x7fff_ff00u32),
        ("negative", (-5i32).cast_unsigned()),
    ] {
        let bam = dir.join(format!("{name}.bam"));
        fs::write(&bam, forged(size)).expect("write the BAM");
        for tool in TOOLS {
            let out = dir.join(format!("{name}.{tool}.txt"));
            let args = [*tool, path(&bam), "-o", path(&out)];
            let (code, stderr, peak) = measured(&args);
            assert_eq!(code, 1, "{name} {tool}: {stderr}");
            assert!(stderr.contains("truncated BAM"), "{name} {tool}: {stderr}");
            assert!(peak < LIMIT, "{name} {tool}: peak RSS {peak} bytes");
            assert!(!out.exists(), "{name} {tool} wrote {}", out.display());
        }
    }
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

fn path(p: &Path) -> &str {
    p.to_str().expect("UTF-8 path")
}

/// The offset after the ITF8 at `at`.
fn itf8(bytes: &[u8], at: usize) -> usize {
    at + 1 + bytes[at].leading_ones().min(4) as usize
}

/// The offset after the LTF8 at `at`.
fn ltf8(bytes: &[u8], at: usize) -> usize {
    at + 1 + bytes[at].leading_ones().min(8) as usize
}

/// The offset after the CRAM 3 container header at `at`: the length, four
/// ITF8s, two LTF8s, the block count, the landmarks and the CRC32.
fn container_header_end(bytes: &[u8], at: usize) -> usize {
    let mut i = at + 4;
    for _ in 0..4 {
        i = itf8(bytes, i);
    }
    i = ltf8(bytes, ltf8(bytes, i));
    i = itf8(bytes, i);
    let landmarks = bytes[i];
    assert!(landmarks < 0x80, "a one-byte landmark count");
    i = itf8(bytes, i);
    for _ in 0..landmarks {
        i = itf8(bytes, i);
    }
    i + 4
}

/// `tests/data/insert_size.cram` with the length of its first data
/// container set to 0x7fffff00, the header CRC32 made to match, and 200
/// bytes of the container kept.
fn forged_cram(root: &Path) -> Vec<u8> {
    let cram = fs::read(root.join("tests/data/insert_size.cram")).expect("read the CRAM");
    // The 26-byte file definition, then the header container.
    let header = 26;
    let length = i32::from_le_bytes(cram[header..header + 4].try_into().unwrap());
    let data = container_header_end(&cram, header) + usize::try_from(length).unwrap();
    let mut forged = cram[..data + 200].to_vec();
    forged[data..data + 4].copy_from_slice(&0x7fff_ff00i32.to_le_bytes());
    let end = container_header_end(&forged, data);
    let crc = crc32(&forged[data..end - 4]);
    forged[end - 4..end].copy_from_slice(&crc.to_le_bytes());
    forged
}

#[test]
fn a_cram_container_past_the_file_costs_the_bytes_present() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let dir = fresh_dir();
    let cram = dir.join("past.cram");
    fs::write(&cram, forged_cram(&root)).expect("write the CRAM");
    let reference = root.join("tests/data/ref_100k.fa");
    for tool in ["insert-size", "quality-yield"] {
        let out = dir.join(format!("{tool}.txt"));
        let args = [tool, path(&cram), "-R", path(&reference), "-o", path(&out)];
        let (code, stderr, peak) = measured(&args);
        assert_eq!(code, 1, "{tool}: {stderr}");
        assert!(stderr.contains("truncated CRAM"), "{tool}: {stderr}");
        assert!(peak < LIMIT, "{tool}: peak RSS {peak} bytes");
    }
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}
