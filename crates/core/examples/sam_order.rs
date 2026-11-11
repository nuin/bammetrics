//! Sort the records of a SAM, BAM or CRAM in an htsjdk order
//! ([`bammetrics_core::sam_order`]) and compare the result, record by
//! record as SAM text, with a file that holds the expected order:
//!
//! ```text
//! cargo run --release -p bammetrics-core --example sam_order -- ORDER INPUT EXPECTED
//! ```
//!
//! `ORDER` is `coordinate`, `queryname`, `duplicate` or `unsorted`.
//! `tests/concordance/test_sam_order.py` runs it on a shuffled corpus
//! against the output of Picard `SortSam`. The sort is stable and in
//! memory, as htsjdk's is for a file under `MAX_RECORDS_IN_RAM` records.
//!
//! Exit 0 when every record matches, 1 at the first difference, 2 when a
//! file cannot be read or a record has no key (htsjdk throws there too).

use std::cmp::Ordering;
use std::path::Path;
use std::process::ExitCode;

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_order::{CoordinateKey, DuplicateKey, Libraries, Order, QuerynameKey};
use bammetrics_core::sam_text::{self, Stringency};
use bammetrics_core::{Error, Result};
use noodles::sam;

/// The header and every record of `path`, read without validation.
fn read(path: &Path) -> Result<(sam::Header, Vec<RawRecord>)> {
    pass::with_record_stream_at(path, None, 1, Stringency::Silent, |header, reader, _| {
        let mut records = Vec::new();
        let mut record = RawRecord::default();
        while record.read_from(reader).map_err(Error::Record)? {
            records.push(record.clone());
        }
        Ok((header.clone(), records))
    })
}

/// The positions of `records` in the order of their keys, ties in input
/// order.
fn sorted<'a, K>(
    records: &'a [RawRecord],
    key: impl Fn(&'a RawRecord) -> Result<K>,
    compare: impl Fn(&K, &K) -> Ordering,
) -> Result<Vec<usize>> {
    let keys = records.iter().map(key).collect::<Result<Vec<K>>>()?;
    let mut order: Vec<usize> = (0..records.len()).collect();
    // Java's sort, as htsjdk's: an order that is not total never panics.
    bammetrics_core::timsort::sort_by(&mut order, |&a, &b| compare(&keys[a], &keys[b]))
        .map_err(|violation| Error::Metric(violation.to_string()))?;
    Ok(order)
}

/// The SAM line of each record, without its line end.
fn lines(header: &sam::Header, records: &[RawRecord]) -> Result<Vec<Vec<u8>>> {
    let references: Vec<Vec<u8>> = header
        .reference_sequences()
        .keys()
        .map(|name| name.to_vec())
        .collect();
    records
        .iter()
        .map(|record| {
            let mut line = Vec::new();
            sam_text::write_record(&mut line, record, &references).map_err(Error::Record)?;
            line.pop();
            Ok(line)
        })
        .collect()
}

fn run(order: Order, input: &Path, expected: &Path) -> Result<bool> {
    let (header, records) = read(input)?;
    let libraries = Libraries::new(&header);
    let positions = match order {
        Order::Unsorted => (0..records.len()).collect(),
        Order::Coordinate => sorted(
            &records,
            |record| Ok(CoordinateKey::of(record)),
            CoordinateKey::compare,
        )?,
        Order::Queryname => sorted(&records, QuerynameKey::of, QuerynameKey::compare)?,
        Order::Duplicate => sorted(
            &records,
            |record| DuplicateKey::of(record, &libraries),
            DuplicateKey::compare,
        )?,
    };
    let ours = lines(&header, &records)?;
    let (expected_header, expected_records) = read(expected)?;
    let theirs = lines(&expected_header, &expected_records)?;
    if ours.len() != theirs.len() {
        println!("{} records sorted, {} expected", ours.len(), theirs.len());
        return Ok(false);
    }
    for (at, (&position, expected)) in positions.iter().zip(&theirs).enumerate() {
        if &ours[position] != expected {
            println!(
                "record {} differs:\n  sorted:   {}\n  expected: {}",
                at + 1,
                String::from_utf8_lossy(&ours[position]),
                String::from_utf8_lossy(expected)
            );
            return Ok(false);
        }
    }
    println!("{} records in {} order", ours.len(), order.name());
    Ok(true)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [order, input, expected] = args.as_slice() else {
        eprintln!("usage: sam_order ORDER INPUT EXPECTED");
        return ExitCode::from(2);
    };
    let Some(order) = Order::parse(order) else {
        eprintln!("sam_order: {order}: expected coordinate, queryname, duplicate or unsorted");
        return ExitCode::from(2);
    };
    match run(order, Path::new(input), Path::new(expected)) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("sam_order: {e}");
            ExitCode::from(2)
        }
    }
}
