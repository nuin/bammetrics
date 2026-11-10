//! Deflater spike (docs/concordance.md, "BAM bytes"): can a Rust deflater
//! write the BGZF blocks Picard wrote?
//!
//!     deflater-spike PICARD.bam [MORE.bam ...]
//!
//! Every BGZF block of each file is inflated, deflated again by each
//! candidate at each of its levels, and the result compared with the bytes
//! the block holds. One line per file, candidate and level goes to standard
//! output, tab separated:
//!
//!     file  encoder  level  equal_blocks  data_blocks  bytes  file_bytes
//!
//! `equal_blocks` counts the blocks whose deflate stream came out byte for
//! byte, `bytes` is the candidate's total and `file_bytes` the total of the
//! streams in the file. A file is reproducible by a candidate only when
//! `equal_blocks` equals `data_blocks` at one level.

use std::env;
use std::fs;
use std::ops::RangeInclusive;
use std::process::ExitCode;

/// The fixed part of a BGZF block before the deflate stream: the gzip
/// header with one `BC` extra field.
const HEADER: usize = 18;
/// CRC32 and ISIZE.
const TRAILER: usize = 8;

/// One BGZF block: the deflate stream it holds and the bytes it inflates to.
struct Block {
    stream: Vec<u8>,
    data: Vec<u8>,
}

/// A candidate deflater: raw deflate of `data` at `level`.
type Encode = fn(&[u8], i32) -> Result<Vec<u8>, String>;

const CANDIDATES: &[(&str, RangeInclusive<i32>, Encode)] = &[
    ("zlib-rs", 1..=9, zlib_rs_deflate),
    ("miniz_oxide", 1..=10, miniz_deflate),
    ("libdeflate", 1..=12, libdeflate_deflate),
];

fn zlib_rs_deflate(data: &[u8], level: i32) -> Result<Vec<u8>, String> {
    let mut out = vec![0; zlib_rs::compress_bound(data.len())];
    let mut encoder = zlib_rs::Deflate::new(level, false, 15);
    match encoder.compress(data, &mut out, zlib_rs::DeflateFlush::Finish) {
        Ok(zlib_rs::Status::StreamEnd) => {
            out.truncate(encoder.total_out() as usize);
            Ok(out)
        }
        other => Err(format!("zlib-rs level {level}: {other:?}")),
    }
}

fn miniz_deflate(data: &[u8], level: i32) -> Result<Vec<u8>, String> {
    Ok(miniz_oxide::deflate::compress_to_vec(data, level as u8))
}

fn libdeflate_deflate(data: &[u8], level: i32) -> Result<Vec<u8>, String> {
    let level = libdeflater::CompressionLvl::new(level)
        .map_err(|_| format!("libdeflate has no level {level}"))?;
    let mut encoder = libdeflater::Compressor::new(level);
    let mut out = vec![0; encoder.deflate_compress_bound(data.len())];
    let len = encoder
        .deflate_compress(data, &mut out)
        .map_err(|e| format!("libdeflate: {e}"))?;
    out.truncate(len);
    Ok(out)
}

fn inflate(stream: &[u8], size: usize) -> Result<Vec<u8>, String> {
    let mut data = vec![0; size];
    let mut decoder = zlib_rs::Inflate::new(false, 15);
    match decoder.decompress(stream, &mut data, zlib_rs::InflateFlush::Finish) {
        Ok(zlib_rs::Status::StreamEnd) if decoder.total_out() as usize == size => Ok(data),
        other => Err(format!("inflate: {other:?}")),
    }
}

/// The blocks of a BGZF file that hold data (the empty end-of-file block
/// is left out).
fn blocks(bytes: &[u8]) -> Result<Vec<Block>, String> {
    let mut found = Vec::new();
    let mut at = 0;
    while at < bytes.len() {
        let head = bytes
            .get(at..at + HEADER)
            .ok_or_else(|| format!("offset {at}: cut block header"))?;
        if head[..4] != [0x1f, 0x8b, 8, 4] || head[10..14] != [6, 0, b'B', b'C'] {
            return Err(format!("offset {at}: not a BGZF block"));
        }
        let total = usize::from(u16::from_le_bytes([head[16], head[17]])) + 1;
        let block = bytes
            .get(at..at + total)
            .filter(|_| total >= HEADER + TRAILER)
            .ok_or_else(|| format!("offset {at}: cut block"))?;
        let size = u32::from_le_bytes([
            block[total - 4],
            block[total - 3],
            block[total - 2],
            block[total - 1],
        ]) as usize;
        if size > 0 {
            let stream = block[HEADER..total - TRAILER].to_vec();
            let data = inflate(&stream, size).map_err(|e| format!("offset {at}: {e}"))?;
            found.push(Block { stream, data });
        }
        at += total;
    }
    Ok(found)
}

fn report(path: &str) -> Result<(), String> {
    let bytes = fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    let blocks = blocks(&bytes).map_err(|e| format!("{path}: {e}"))?;
    let file_bytes: usize = blocks.iter().map(|b| b.stream.len()).sum();
    for (name, levels, encode) in CANDIDATES {
        for level in levels.clone() {
            let mut equal = 0;
            let mut total = 0;
            for block in &blocks {
                let stream = encode(&block.data, level)?;
                total += stream.len();
                equal += usize::from(stream == block.stream);
            }
            println!(
                "{path}\t{name}\t{level}\t{equal}\t{}\t{total}\t{file_bytes}",
                blocks.len()
            );
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    let paths: Vec<String> = env::args().skip(1).collect();
    if paths.is_empty() {
        eprintln!("usage: deflater-spike PICARD.bam [MORE.bam ...]");
        return ExitCode::from(2);
    }
    for path in &paths {
        if let Err(error) = report(path) {
            eprintln!("deflater-spike: {error}");
            return ExitCode::FAILURE;
        }
    }
    ExitCode::SUCCESS
}
