//! `AddOrReplaceReadGroups`, mirroring Picard 3.5.0
//! (`picard.sam.AddOrReplaceReadGroups`): every `@RG` line of the input is
//! replaced by the one read group the options give, and every record's
//! `RG` tag is set to `RGID`.
//!
//! The `@RG` line holds `ID`, then `LB`, `PL`, `SM`, `PU`, then the
//! optional `CN`, `DS`, `DT`, `PI`, `PG`, `PM`, `KS`, `FO` that are given,
//! in that order (the order `doWork` sets them). Each text value must match
//! `^[ -~]+$` (printable ASCII and space, at least one character); every
//! value that does not is named in one error, as Picard's
//! `customCommandLineValidation` names them. `RGPI` is read by
//! `Integer.parseInt` and written in decimal. `RGDT` is read by htsjdk's
//! `DateParser` and written in UTC ([`date`]). `RGID=null` is Java's
//! `null`: Picard writes `ID:null` and removes each record's `RG` tag, and
//! so does bammetrics. No `@PG` line is added.
//!
//! The output header is the input's with `SO` set to `SORT_ORDER` when it
//! is given. Picard writes through a writer made with `presorted` true
//! when the output order equals the input's (the order checked, the
//! records as they come) and false otherwise: then `queryname`,
//! `coordinate` and `duplicate` are sorted as `SortSam` sorts them
//! (`MAX_RECORDS_IN_RAM` per run, spilled to `TMP_DIR`), `unsorted` is
//! written as read, and `unknown` fails on a second record
//! ([`record_edit::WriterOrder`]). Every record is encoded again, as
//! htsjdk encodes a record whose tag changed: tags in
//! [`crate::bam_raw::Tags`] order and widths, the `bin` htsjdk computes,
//! qualities whose first byte is `0xFF` written as all `0xFF`.
//!
//! In a chain the step is a streaming transform that leaves `rg`, or with
//! a `SORT_ORDER` of `queryname`, `coordinate` or `duplicate` a barrier
//! that sorts ([`record_edit::sorting_stage`]); `SORT_ORDER=unknown` is a
//! usage error there.

use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::read_group::tag as rg;
use noodles::sam::header::record::value::map::tag::Other;

use crate::bam_raw::{RawRecord, TagValue};
use crate::error::{Error, Result};
use crate::header_merge::{HeaderTag, SortOrder, set_header_tag};
use crate::metric::ToolOutput;
use crate::record_edit::{self, Outcome, RecordEdit, Reorder, Setup, WriterOrder};
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_order::Order;
use crate::spill::TempDir;
use crate::stage::{Frame, Stage};

pub mod date;

pub const TOOL: Tool = Tool {
    name: "AddOrReplaceReadGroups",
    summary: "no rows; OUTPUT gets every record of the input with its RG tag set to RGID under \
              the input header whose @RG lines are replaced by the one read group the options \
              give, sorted when SORT_ORDER differs from the input's: a BAM, or SAM text for a \
              name ending in .sam.",
    options: &[
        Opt::text(&["RGID", "ID"], "Read group ID (default 1)."),
        Opt::text(&["RGLB", "LB"], "Read group library.").required(),
        Opt::text(
            &["RGPL", "PL"],
            "Read group platform (ILLUMINA, SOLID, ...).",
        )
        .required(),
        Opt::text(&["RGPU", "PU"], "Read group platform unit (a run barcode).").required(),
        Opt::text(&["RGSM", "SM"], "Read group sample name.").required(),
        Opt::text(&["RGCN", "CN"], "Read group sequencing center name."),
        Opt::text(&["RGDS", "DS"], "Read group description."),
        Opt::text(
            &["RGDT", "DT"],
            "Read group run date, ISO 8601 (2020-01-02, 2020-01-02T03:04:05+05:30), written \
             in UTC.",
        ),
        Opt::int(&["RGPI", "PI"], "Read group predicted insert size."),
        Opt::text(&["RGPG", "PG"], "Read group program group."),
        Opt::text(&["RGPM", "PM"], "Read group platform model."),
        Opt::text(&["RGKS", "KS"], "Read group key sequence."),
        Opt::text(&["RGFO", "FO"], "Read group flow order."),
        Opt::text(
            &["SORT_ORDER", "SO"],
            "The order of OUTPUT: unsorted, queryname, coordinate, duplicate or unknown \
             (default: the input's). An order that differs from the input's is sorted.",
        ),
        record_edit::OUTPUT_OPTION,
        record_edit::COMPRESSION_OPTION,
        record_edit::INDEX_OPTION,
        record_edit::MD5_OPTION,
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records held in memory before a sorted run spills to TMP_DIR, 1 or more, when \
             the records are sorted (default 500000).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the spilled runs of a sort (default: the system temporary \
             directory).",
        )
        .repeated(),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// htsjdk's `SAMReadGroupRecord` fields after `ID`, in the order `doWork`
/// sets them, each with its option.
const FIELDS: [(&str, Other<rg::Standard>); 12] = [
    ("RGLB", rg::LIBRARY),
    ("RGPL", rg::PLATFORM),
    ("RGSM", rg::SAMPLE),
    ("RGPU", rg::PLATFORM_UNIT),
    ("RGCN", rg::SEQUENCING_CENTER),
    ("RGDS", rg::DESCRIPTION),
    ("RGDT", rg::PRODUCED_AT),
    ("RGPI", rg::PREDICTED_MEDIAN_INSERT_SIZE),
    ("RGPG", rg::PROGRAM),
    ("RGPM", rg::PLATFORM_MODEL),
    ("RGKS", rg::KEY_SEQUENCE),
    ("RGFO", rg::FLOW_ORDER),
];

/// The options `customCommandLineValidation` checks, in its order.
const CHECKED: [&str; 11] = [
    "RGID", "RGLB", "RGPL", "RGPU", "RGSM", "RGCN", "RGDS", "RGKS", "RGFO", "RGPG", "RGPM",
];

/// The required options, as barclay names them.
const REQUIRED: [&str; 4] = ["RGLB", "RGPL", "RGPU", "RGSM"];

/// The read group and order the options give.
#[derive(Clone, Debug)]
struct ReadGroupOptions {
    /// `None` for `RGID=null`.
    id: Option<String>,
    /// `(tag, value)` after `ID`, in htsjdk's order.
    fields: Vec<(Other<rg::Standard>, String)>,
    /// `SORT_ORDER`, `None` when not given.
    order: Option<SortOrder>,
}

/// A text option: `None` when absent or `null`.
fn text(inv: &Invocation, name: &str) -> Option<String> {
    inv.get(name)
        .filter(|value| *value != "null")
        .map(str::to_owned)
}

/// `^[ -~]+$`.
fn printable(value: &str) -> bool {
    !value.is_empty() && value.chars().all(|c| (' '..='~').contains(&c))
}

impl ReadGroupOptions {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        for name in REQUIRED {
            if text(inv, name).is_none() {
                let shown = if inv.get(name).is_some() {
                    format!("non-null value must be provided for '{name}'")
                } else {
                    format!("Argument '{name}' is required")
                };
                return Err(Error::Usage(shown));
            }
        }
        let order = inv
            .get("SORT_ORDER")
            .filter(|value| *value != "null")
            .map(|value| {
                SortOrder::parse(value).ok_or_else(|| {
                    Error::Usage(format!("'{value}' is not a valid value for SortOrder."))
                })
            })
            .transpose()?;
        let bad: Vec<String> = CHECKED
            .iter()
            .filter_map(|name| {
                let value = inv.get(name).filter(|v| *v != "null")?;
                (!printable(value)).then(|| {
                    format!(
                        "The values of tags in a SAM header must adhere to the regular \
                         expression '^[ -~]+$',but the value provided for {name}, '{value}', \
                         doesn't."
                    )
                })
            })
            .collect();
        if !bad.is_empty() {
            return Err(Error::Usage(bad.join("\n")));
        }
        let mut fields = Vec::new();
        for (name, tag) in FIELDS {
            let Some(value) = text(inv, name) else {
                continue;
            };
            let value = match name {
                "RGDT" => date::parse(&value).map(date::format_utc).map_err(|why| {
                    Error::Usage(format!(
                        "RGDT={value}: Problem constructing Iso8601Date from the string \
                         '{value}' ({why})"
                    ))
                })?,
                "RGPI" => crate::java::parse_integer(&value, 32)
                    .map(|n| n.to_string())
                    .map_err(|_| {
                        Error::Usage(format!(
                            "RGPI={value}: Problem constructing Integer from the string \
                             '{value}'"
                        ))
                    })?,
                _ => value,
            };
            fields.push((tag, value));
        }
        let id = match inv.get("RGID") {
            None => Some("1".to_owned()),
            Some("null") => None,
            Some(id) => Some(id.to_owned()),
        };
        Ok(Self { id, fields, order })
    }

    /// The sort of a chain step: a `SORT_ORDER` that names a sort.
    fn chain_sort(&self) -> Option<Order> {
        match self.order? {
            SortOrder::Queryname => Some(Order::Queryname),
            SortOrder::Coordinate => Some(Order::Coordinate),
            SortOrder::Duplicate => Some(Order::Duplicate),
            SortOrder::Unsorted | SortOrder::Unknown => None,
        }
    }
}

/// htsjdk `BAMRecord.getBaseQualities` reads a first byte `0xFF` as no
/// qualities, which its codec writes back as `0xFF` for every base.
pub(crate) fn htsjdk_qualities(record: &mut RawRecord) {
    let qualities = record.qualities_mut();
    if qualities.first() == Some(&0xFF) {
        qualities.fill(0xFF);
    }
}

/// The edit: the read group in the header, `RG` on each record.
struct ReadGroups {
    options: ReadGroupOptions,
}

impl RecordEdit for ReadGroups {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        let mut header = input.clone();
        let mut map = Map::<ReadGroup>::default();
        for (tag, value) in &self.options.fields {
            map.other_fields_mut().insert(*tag, value.as_str().into());
        }
        let id = self.options.id.as_deref().unwrap_or("null");
        header.read_groups_mut().clear();
        header.read_groups_mut().insert(id.into(), map);
        if let Some(order) = self.options.order {
            set_header_tag(&mut header, HeaderTag::SortOrder, order.name());
        }
        Ok(header)
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        let mut tags = record.tags().map_err(Error::Record)?;
        match &self.options.id {
            Some(id) => tags.set(*b"RG", TagValue::String(id.clone().into_bytes())),
            None => tags.remove(*b"RG"),
        }
        record.set_tags(&tags).map_err(Error::Record)?;
        htsjdk_qualities(record);
        Ok(Outcome::Changed)
    }

    fn rewrites(&self) -> bool {
        true
    }
}

fn build(options: ReadGroupOptions) -> impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>> {
    move |_| Ok(Box::new(ReadGroups { options }))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let options = ReadGroupOptions::from_invocation(inv)?;
    let mut opts = record_edit::Options::from_invocation(inv)?;
    let (limits, tmp_dir) = record_edit::sort_limits(inv)?;
    opts.reorder = Reorder {
        regroup: None,
        writer: WriterOrder::Header {
            order: options.order,
            presorted_if_same: true,
        },
        limits,
        tmp_dir,
    };
    record_edit::run(&opts, build(options))?;
    Ok(ToolOutput::NONE)
}

/// The usage errors of the step, asked before the first read of a chain:
/// an option of the written file, a read group value Picard refuses, an
/// unknown order.
pub fn check_stage(inv: &Invocation) -> Result<()> {
    record_edit::check_stage(inv)?;
    let options = ReadGroupOptions::from_invocation(inv)?;
    if options.order == Some(SortOrder::Unknown) {
        return Err(Error::Usage(
            "SORT_ORDER=unknown in a chain: Picard sorts unknown with no comparator and fails \
             on any input whose order is not unknown"
                .to_owned(),
        ));
    }
    if options.chain_sort().is_some() {
        record_edit::stage_limits(inv, 1)?;
    }
    Ok(())
}

/// True for a step that sorts: a `SORT_ORDER` of `queryname`,
/// `coordinate` or `duplicate` makes the step a barrier.
pub fn sorts(argv: &[String]) -> bool {
    argv.iter().rev().find_map(|arg| {
        let (key, value) = arg.split_once('=')?;
        ["SORT_ORDER", "SO"]
            .contains(&key)
            .then_some(matches!(value, "queryname" | "coordinate" | "duplicate"))
    }) == Some(true)
}

/// `AddOrReplaceReadGroups` as a step of a chain. The writer is resolved
/// against the stream's header as the tool resolves it against its
/// input's: a `SORT_ORDER` that is the stream's order (or none) streams the
/// records with the presorted writer's order check, so a stream sorted by
/// another tool keeps its order of tied records; any other order sorts.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let options = ReadGroupOptions::from_invocation(inv)?;
    let writer = WriterOrder::Header {
        order: options.order,
        presorted_if_same: true,
    }
    .against(SortOrder::of(&frame.header));
    let Stage::Transform(transform) = record_edit::stage(inv, frame, build(options))? else {
        return Err(Error::Metric(
            "AddOrReplaceReadGroups: the edit of a step is a transform".to_owned(),
        ));
    };
    if writer == WriterOrder::Presorted {
        return record_edit::presorted_stage(frame, transform, TOOL.name);
    }
    let reorder = Reorder {
        writer,
        limits: record_edit::stage_limits(inv, 1)?,
        ..Reorder::default()
    };
    record_edit::sorting_stage(inv, frame, temp, &reorder, transform)
}

#[cfg(test)]
mod tests;
