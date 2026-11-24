//! `SamFormatConverter`, mirroring Picard 3.5.0
//! (`picard.sam.SamFormatConverter`): every record of `INPUT` (SAM text,
//! BAM, or a CRAM with `REFERENCE_SEQUENCE`) written unchanged to `OUTPUT`,
//! in the format its name gives.
//!
//! Picard opens `INPUT` with `SamReaderFactory` and writes through
//! `SAMFileWriterFactory.makeWriter(header, true, OUTPUT, REFERENCE_SEQUENCE)`,
//! a presorted writer. The name decides the format as htsjdk decides it
//! ([`WriteFormat::of`]): SAM text for a name ending in `.sam`, a CRAM for
//! `.cram`, both lower case, and a BAM for any other name (`out.BAM`,
//! `out.sam.gz`, `out.txt`; htsjdk logs `Unknown file extension, assuming
//! BAM format`). This is [`record_edit::run`] in [`Mode::Rewrite`] with an
//! edit that changes nothing:
//! * the header is the input's as htsjdk's writer writes it (`@HD VN:1.6`,
//!   `SO:unsorted` added to an `@HD` without one); no `@PG` line is added;
//! * each record keeps its fields and tags, with the `bin` htsjdk's codec
//!   computes; a record of a CRAM gets `MD` and `NM` as htsjdk decodes it;
//! * the writer checks the header's order (`assertPresorted`).
//!
//! `CREATE_INDEX` on a header whose `SO` is not `coordinate` fails before
//! any record is written, with Picard's sentence. Picard checks it after it
//! has made `OUTPUT`, so it leaves an empty BAM or a SAM header behind;
//! bammetrics leaves none. A BAM gets `out.bai`, SAM text gets no index
//! (htsjdk writes none and says nothing), and a CRAM refuses the option
//! (`cram_writer::index_refusal`). `MAX_RECORDS_IN_RAM` changes nothing: a
//! presorted writer sorts nothing.
//!
//! Not a chain step: the `write` step of a chain writes SAM text, a BAM or
//! a CRAM by its name, which is all this tool does.

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::header_merge::SortOrder;
use crate::metric::ToolOutput;
use crate::record_edit::{self, Outcome, RecordEdit, Setup};
use crate::registry::{Invocation, Opt, Tool};

pub const TOOL: Tool = Tool {
    name: "SamFormatConverter",
    summary: "no rows; OUTPUT gets every record of INPUT unchanged, as SAM text for a name \
              ending in .sam, a CRAM for one ending in .cram, else a BAM.",
    options: &[
        Opt::path(
            &["OUTPUT", "O"],
            "The file to write: SAM text for a name ending in .sam, a CRAM (with \
             REFERENCE_SEQUENCE) for one ending in .cram, a BAM for any other name.",
        )
        .written()
        .required(),
        record_edit::COMPRESSION_OPTION,
        Opt::bool(
            &["CREATE_INDEX"],
            "Write a .bai beside a BAM OUTPUT; an input whose header is not coordinate sorted \
             fails, as in Picard (default false).",
        ),
        record_edit::MD5_OPTION,
        Opt::ignored_int(
            &["MAX_RECORDS_IN_RAM"],
            "records are written in input order",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// The edit: the header and every record as they are. `create_index`
/// holds Picard's check of the header's order.
struct Identity {
    create_index: bool,
}

impl RecordEdit for Identity {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        if self.create_index && SortOrder::of(input) != SortOrder::Coordinate {
            return Err(Error::Metric(
                "Can't CREATE_INDEX unless sort order is coordinate".to_owned(),
            ));
        }
        Ok(input.clone())
    }

    fn record(&mut self, _record: &mut RawRecord) -> Result<Outcome> {
        Ok(Outcome::Same)
    }

    fn rewrites(&self) -> bool {
        true
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = record_edit::Options::from_invocation(inv)?;
    let create_index = opts.create_index;
    record_edit::run(&opts, |_: &Setup| {
        Ok(Box::new(Identity { create_index }) as Box<dyn RecordEdit>)
    })?;
    Ok(ToolOutput::NONE)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(text: &str) -> sam::Header {
        text.parse().expect("header")
    }

    #[test]
    fn the_header_is_kept() {
        let input = header("@HD\tVN:1.5\tSO:queryname\n@SQ\tSN:chr1\tLN:10\n@CO\tc\n");
        let edit = Identity {
            create_index: false,
        };
        assert_eq!(edit.header(&input).expect("kept"), input);
    }

    #[test]
    fn an_index_needs_a_coordinate_header() {
        let edit = Identity { create_index: true };
        for text in [
            "@HD\tVN:1.6\tSO:queryname\n",
            "@HD\tVN:1.6\n",
            "@HD\tVN:1.6\tSO:unknown\n",
            "@SQ\tSN:chr1\tLN:10\n",
        ] {
            let error = edit.header(&header(text)).expect_err(text);
            assert_eq!(
                error.to_string(),
                "Can't CREATE_INDEX unless sort order is coordinate"
            );
        }
        assert!(
            edit.header(&header("@HD\tVN:1.6\tSO:coordinate\n"))
                .is_ok()
        );
    }

    #[test]
    fn records_are_kept() {
        let mut edit = Identity {
            create_index: false,
        };
        let mut record = RawRecord::default();
        assert_eq!(edit.record(&mut record).expect("kept"), Outcome::Same);
        assert!(edit.rewrites());
    }
}
