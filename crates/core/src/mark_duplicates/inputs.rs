//! `AbstractMarkDuplicatesCommandLineProgram.openInputs`: every `INPUT` as
//! one header and one record stream.
//!
//! One input is read as it is. Several are merged as htsjdk 5.0.0 merges
//! them (`SamFileHeaderMerger` without dictionary merging,
//! `MergingSamRecordIterator`), through [`crate::header_merge`]. Without
//! `ASSUME_SORT_ORDER` every input must declare the first input's sort
//! order.

use std::io::Read;
use std::path::{Path, PathBuf};

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::header_merge::{HeaderMerger, Merger};
pub use crate::header_merge::{HeaderTag, SortOrder, base36, set_header_tag};
use crate::pass;
use crate::sam_text::Stringency;

/// One opened input: its header and its record stream.
type Stream<'a> = (sam::Header, &'a mut (dyn Read + 'static), bool);

/// What `openInputs` gives: the header and the records in order.
pub struct Inputs<'a> {
    pub header: sam::Header,
    /// Some input is a CRAM: htsjdk hands its records over with sorted tags
    /// and rebuilt `MD` and `NM` ([`super::md_nm`]).
    pub cram: bool,
    source: Source<'a>,
}

enum Source<'a> {
    One(&'a mut (dyn Read + 'static)),
    Merged(Box<Merger<'a>>),
}

impl Inputs<'_> {
    /// The next record into `record`; false after the last.
    pub fn next(&mut self, record: &mut RawRecord) -> Result<bool> {
        match &mut self.source {
            Source::One(reader) => record.read_from(*reader).map_err(Error::Record),
            Source::Merged(merger) => merger.next(record),
        }
    }
}

/// Open every input and hand `f` the merged view. `assume` is
/// `ASSUME_SORT_ORDER` (with `ASSUME_SORTED=true` as `coordinate`): it is
/// written into the first header, as Picard does. SAM text is read at
/// `stringency`.
pub fn open<T>(
    paths: &[PathBuf],
    reference: Option<&Path>,
    threads: usize,
    stringency: Stringency,
    assume: Option<SortOrder>,
    f: impl FnOnce(Inputs<'_>) -> Result<T>,
) -> Result<T> {
    let mut f = Some(f);
    let mut call = |streams: Vec<Stream<'_>>| -> Result<T> {
        let f = f.take().expect("called once");
        f(combine(streams, assume)?)
    };
    let decoding = Decoding {
        reference,
        threads,
        stringency,
    };
    open_nested(paths, decoding, Vec::new(), &mut call)
}

/// How every input is decoded.
#[derive(Clone, Copy)]
struct Decoding<'a> {
    reference: Option<&'a Path>,
    threads: usize,
    stringency: Stringency,
}

fn open_nested<T>(
    paths: &[PathBuf],
    decoding: Decoding<'_>,
    opened: Vec<Stream<'_>>,
    f: &mut dyn FnMut(Vec<Stream<'_>>) -> Result<T>,
) -> Result<T> {
    match paths.split_first() {
        None => f(opened),
        Some((first, rest)) => {
            let Decoding {
                reference,
                threads,
                stringency,
            } = decoding;
            pass::with_record_stream_at(
                first,
                reference,
                threads,
                stringency,
                |header, reader, cram| {
                    let mut opened = opened;
                    opened.push((header.clone(), reader, cram));
                    open_nested(rest, decoding, opened, f)
                },
            )
        }
    }
}

fn combine(mut streams: Vec<Stream<'_>>, assume: Option<SortOrder>) -> Result<Inputs<'_>> {
    if let Some(order) = assume {
        set_header_tag(&mut streams[0].0, HeaderTag::SortOrder, order.name());
    }
    if streams.len() == 1 {
        let (header, reader, cram) = streams.pop().expect("one input");
        return Ok(Inputs {
            header,
            cram,
            source: Source::One(reader),
        });
    }
    let cram = streams.iter().any(|s| s.2);
    let (headers, readers): (Vec<sam::Header>, Vec<_>) =
        streams.into_iter().map(|(h, r, _)| (h, r)).unzip();
    let merged = HeaderMerger::merge(&headers)?;
    let order = SortOrder::of(&merged.header);
    if assume.is_none() && order != SortOrder::Unsorted {
        for header in &headers {
            if SortOrder::of(header) != order {
                return Err(Error::Metric(
                    "Files are not compatible with sort order".to_owned(),
                ));
            }
        }
    }
    let header = merged.header.clone();
    let merger = Merger::new(readers, merged, order)?;
    Ok(Inputs {
        header,
        cram,
        source: Source::Merged(Box::new(merger)),
    })
}

#[cfg(test)]
mod tests {
    use crate::sam_order::{coordinate_order, queryname_order};
    use std::cmp::Ordering;

    #[test]
    fn coordinate_order_puts_unplaced_last_and_reverse_after_forward() {
        let a = crate::bam_raw::build("a", 0, 0, 10, &[], &[], b"");
        let b = crate::bam_raw::build("b", 0x10, 0, 10, &[], &[], b"");
        let unplaced = crate::bam_raw::build("c", 0x4, -1, -1, &[], &[], b"");
        assert_eq!(coordinate_order(&a, &b), Ordering::Less);
        assert_eq!(coordinate_order(&unplaced, &a), Ordering::Greater);
        let first = crate::bam_raw::build("x", 0x41, 0, 1, &[], &[], b"");
        let second = crate::bam_raw::build("x", 0x81, 0, 1, &[], &[], b"");
        assert_eq!(queryname_order(&first, &second), Ordering::Less);
    }
}
