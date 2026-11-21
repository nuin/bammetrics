//! `PositionBasedDownsampleSam`, mirroring Picard 3.5.0
//! (`picard.sam.PositionBasedDownsampleSam`): the records whose read name
//! places them inside (or, for `FRACTION` over 0.5, outside) a circle of
//! area `FRACTION` on each tile of the flow cell, written in input order
//! through a presorted writer.
//!
//! The location of a read is tile, x and y from its name through
//! `READ_NAME_REGEX` ([`crate::optical::ReadNameParser`]), starting from
//! Picard's `PhysicalLocationInt` of -1, -1, -1: a name without a location
//! counts on tile -1 at (-1, -1), and a regex that matches but whose x or
//! y is no `int` leaves the fields after the tile as they were. Neither is
//! an error, as in Picard.
//!
//! Two reads of `INPUT`, so standard input and pipes are refused. The
//! first finds each tile's extent: `Coord` starts at 0 for the minimum and
//! the maximum of x and y (so a minimum above 0 stays 0, as Picard has it),
//! then each side widens by `(max - min) / count` in Java `int` division.
//! The second keeps a record when `d = r(u)^2 + r(v)^2` is at most `f / pi`,
//! or above it for `FRACTION` over 0.5: `u` is `(x - minX) / (maxX - minX)`
//! less `offset`, `v` the same for y, `r(t)` is `t - Math.round(t)`, `f` is
//! `FRACTION` (`1 - FRACTION` above 0.5) and `offset` is
//! `sqrt(f / pi - f^2 / 4)`. `STOP_AFTER` (a Java `Long`, `null` by
//! default) ends both reads after that many records; 0 or less reads none.
//! `REMOVE_DUPLICATE_INFORMATION` (default true) clears the duplicate flag
//! of each record written.
//!
//! An input whose header has an `@PG` with `PN:PositionBasedDownsampleSam`
//! fails unless `ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS`. The output
//! header gets `@PG` `ID:PositionBasedDownsampleSam` (made unique as
//! htsjdk's `PgIdGenerator` makes it) with `PN`, `CL` and `VN`.

use std::collections::HashMap;

use noodles::sam;
use noodles::sam::header::record::value::map::program::tag;

use crate::bam_raw::RawRecord;
use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::metric::ToolOutput;
use crate::optical::{Location, ReadNameParser};
use crate::pass;
use crate::record_edit::{self, Mode, Out, Restore};
use crate::registry::{Invocation, Opt, Tool};
use crate::sensitivity::math_round;

pub const TOOL: Tool = Tool {
    name: "PositionBasedDownsampleSam",
    summary: "no rows; OUTPUT gets the records whose read name places them in a circle of \
              area FRACTION on each tile of the flow cell (outside one of area 1 - FRACTION \
              above 0.5), duplicate flags cleared, with a @PG line: a BAM, or SAM text for a \
              name ending in .sam.",
    options: &[
        record_edit::OUTPUT_OPTION,
        Opt::float(
            &["FRACTION", "F"],
            "The approximate fraction of reads to keep, 0 to 1.",
        )
        .required(),
        Opt::bool(
            &["REMOVE_DUPLICATE_INFORMATION"],
            "Clear the duplicate flag of each record written (default true).",
        ),
        Opt::text(
            &["READ_NAME_REGEX"],
            "How tile, x and y are read from a read name: Picard's default (5 or 7 fields \
             split on :), a regex with three groups, or null (no location: every read on \
             tile -1).",
        ),
        Opt::text(
            &["STOP_AFTER"],
            "Read only the first N records in each read (default null, all); 0 or less \
             writes none.",
        ),
        Opt::bool(
            &["ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS"],
            "Run on an input this tool downsampled before (default false: an error).",
        ),
        record_edit::COMPRESSION_OPTION,
        record_edit::INDEX_OPTION,
        record_edit::MD5_OPTION,
        record_edit::RECORDS_IN_RAM_OPTION,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// `PositionBasedDownsampleSam.PG_PROGRAM_NAME`.
pub const PROGRAM: &str = "PositionBasedDownsampleSam";
/// The duplicate flag, which `REMOVE_DUPLICATE_INFORMATION` clears.
const DUPLICATE: u16 = 0x400;

/// The options of one run.
#[derive(Clone, Debug)]
pub struct Options {
    pub write: record_edit::Options,
    pub fraction: f64,
    pub remove_duplicate_information: bool,
    pub read_name_regex: Option<String>,
    /// `STOP_AFTER`; `None` for `null`.
    pub stop_after: Option<i64>,
    pub allow_multiple: bool,
    pub command: String,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let text = inv
            .get("FRACTION")
            .ok_or_else(|| Error::Usage("FRACTION is required".to_owned()))?;
        let fraction = crate::java::parse_double(text)
            .map_err(|e| Error::Usage(format!("FRACTION={text}: {e}")))?;
        // `customCommandLineValidation`, NaN passing as in Java.
        if fraction < 0.0 || fraction > 1.0 {
            return Err(Error::Usage(format!(
                "FRACTION must be a value between 0 and 1, found: {}",
                crate::java::double_to_string(fraction)
            )));
        }
        let stop_after = match inv.get("STOP_AFTER") {
            None | Some("null") => None,
            Some(text) => Some(
                crate::java::parse_integer(text, 64)
                    .map_err(|e| Error::Usage(format!("STOP_AFTER={text}: {e}")))?,
            ),
        };
        Ok(Self {
            write: record_edit::Options::from_invocation(inv)?,
            fraction,
            remove_duplicate_information: inv.flag("REMOVE_DUPLICATE_INFORMATION", true)?,
            read_name_regex: inv.get("READ_NAME_REGEX").map(str::to_owned),
            stop_after,
            allow_multiple: inv.flag("ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS", false)?,
            command: inv.command.clone(),
        })
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    run(&opts)?;
    Ok(ToolOutput::NONE)
}

/// The extent of one tile, `PositionBasedDownsampleSam.Coord`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Coord {
    min_x: i32,
    min_y: i32,
    max_x: i32,
    max_y: i32,
    count: i32,
}

impl Coord {
    fn add(&mut self, at: &Location) {
        self.max_x = self.max_x.max(at.x);
        self.min_x = self.min_x.min(at.x);
        self.max_y = self.max_y.max(at.y);
        self.min_y = self.min_y.min(at.y);
        self.count = self.count.wrapping_add(1);
    }

    /// The end of `fillTileMinMaxCoord`: each side widened by the span over
    /// the count, in Java `int` arithmetic.
    fn widen(&mut self) {
        let span_x = self.max_x.wrapping_sub(self.min_x);
        let span_y = self.max_y.wrapping_sub(self.min_y);
        // Every tile counted one read or more; Java's `int` division wraps.
        let step = |span: i32| span.checked_div(self.count).unwrap_or(0);
        let (step_x, step_y) = (step(span_x), step(span_y));
        self.max_x = self.max_x.wrapping_add(step_x);
        self.min_x = self.min_x.wrapping_sub(step_x);
        self.max_y = self.max_y.wrapping_add(step_y);
        self.min_y = self.min_y.wrapping_sub(step_y);
    }
}

/// `PositionBasedDownsampleSam.CircleSelector`.
#[derive(Clone, Copy, Debug)]
struct CircleSelector {
    radius_squared: f64,
    offset: f64,
    positive: bool,
}

impl CircleSelector {
    fn new(fraction: f64) -> Self {
        let (area, positive) = if fraction > 0.5 {
            (1.0 - fraction, false)
        } else {
            (fraction, true)
        };
        let radius_squared = area / std::f64::consts::PI;
        Self {
            radius_squared,
            offset: (radius_squared - area * area / 4.0).sqrt(),
            positive,
        }
    }

    /// `roundedPart`: the distance to the nearest integer, signed.
    fn rounded_part(value: f64) -> f64 {
        value - math_round(value) as f64
    }

    fn select(&self, at: &Location, tile: &Coord) -> bool {
        let axis = |value: i32, min: i32, max: i32| {
            let place = f64::from(value.wrapping_sub(min)) / f64::from(max.wrapping_sub(min));
            let part = Self::rounded_part(place - self.offset);
            part * part
        };
        let distance = axis(at.x, tile.min_x, tile.max_x) + axis(at.y, tile.min_y, tile.max_y);
        (distance > self.radius_squared) ^ self.positive
    }
}

/// The location `getSamRecordLocation` gives a record: Picard's
/// `PhysicalLocationInt` of -1, -1, -1, filled as far as the name allows.
fn location(parser: &ReadNameParser, record: &RawRecord) -> Result<Location> {
    let mut at = Location::default();
    let name = String::from_utf8_lossy(record.name());
    parser.read_location(&name, &mut at)?;
    Ok(at)
}

/// One read of `INPUT` under `STOP_AFTER`: `start` makes a state from the
/// header and whether the records come from a CRAM, `each` takes every
/// record with it. Returns the state.
fn each_record<S>(
    opts: &Options,
    start: impl FnOnce(&sam::Header, bool) -> Result<S>,
    mut each: impl FnMut(&mut S, &mut RawRecord) -> Result<()>,
) -> Result<S> {
    let write = &opts.write;
    pass::with_record_stream_kind(
        &write.input,
        write.reference.as_deref(),
        write.threads,
        write.stringency,
        |header, reader, kind| {
            record_edit::check_sort_order(header, write.stringency, &write.input)?;
            let cram = kind == Kind::Cram;
            let mut state = start(header, cram)?;
            let mut record = RawRecord::default();
            let mut read = 0i64;
            while opts.stop_after.is_none_or(|stop| read < stop)
                && record
                    .read_from(reader)
                    .map_err(|e| pass::read_error(&write.input, cram, e))?
            {
                read += 1;
                each(&mut state, &mut record)?;
            }
            Ok(state)
        },
    )
}

/// `checkProgramRecords`: an earlier run of this tool on the input fails
/// the run unless `ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS`.
/// Picard logs a warning per earlier `@PG` line; the allowed run prints
/// one line with the first ID and the count.
fn check_programs(header: &sam::Header, allow: bool) -> Result<()> {
    let mut earlier = header.programs().as_ref().iter().filter(|(_, map)| {
        map.other_fields()
            .get(&tag::NAME)
            .is_some_and(|name| name.as_slice() == PROGRAM.as_bytes())
    });
    let Some((first, _)) = earlier.next() else {
        return Ok(());
    };
    let message = format!(
        "Found previous Program Record that indicates that this file has been downsampled \
         already with this program. Operation not supported! Previous PG: {}",
        String::from_utf8_lossy(first)
    );
    if !allow {
        return Err(Error::Metric(message));
    }
    let more = earlier.count();
    if more == 0 {
        eprintln!("bammetrics: {message}");
    } else {
        eprintln!("bammetrics: {message} (and {more} more {PROGRAM} @PG lines)");
    }
    Ok(())
}

/// The tool on its own: refuse a stream, check the header, find each tile's
/// extent in a first read, write the records kept in a second.
pub fn run(opts: &Options) -> Result<()> {
    let input = opts.write.input.to_string_lossy();
    if crate::jumping_library::is_stream(&input) {
        return Err(Error::Usage(format!(
            "INPUT={input}: PositionBasedDownsampleSam reads its input twice (the tile extents, \
             then the records), so it needs a file, not standard input or a pipe"
        )));
    }
    record_edit::check_output(&opts.write)?;
    let header = pass::read_header_at(&opts.write.input, opts.write.stringency)?;
    check_programs(&header, opts.allow_multiple)?;
    let parser = ReadNameParser::new(opts.read_name_regex.as_deref())?;
    let mut tiles = each_record(
        opts,
        |_, _| Ok(HashMap::<i16, Coord>::new()),
        |tiles, record| {
            let at = location(&parser, record)?;
            tiles.entry(at.tile).or_default().add(&at);
            Ok(())
        },
    )?;
    tiles.values_mut().for_each(Coord::widen);
    write_kept(opts, &parser, &mut tiles)
}

/// The file of the second read and what each record written needs.
struct Writing {
    out: Out,
    header: sam::Header,
    restore: Restore,
    references: usize,
}

/// The second read: each record the circle selects, written.
fn write_kept(
    opts: &Options,
    parser: &ReadNameParser,
    tiles: &mut HashMap<i16, Coord>,
) -> Result<()> {
    let selector = CircleSelector::new(opts.fraction);
    let writing = each_record(
        opts,
        |header, cram| {
            let mut out_header = header.clone();
            crate::downsample_sam::add_program(&mut out_header, PROGRAM, &opts.command);
            let out_header = record_edit::with_version(out_header, true);
            Ok(Writing {
                restore: Restore::new(
                    cram,
                    Mode::Rewrite,
                    opts.write.reference.as_deref(),
                    header,
                )?,
                out: Out::create(&opts.write, &out_header)?,
                header: out_header,
                references: header.reference_sequences().len(),
            })
        },
        |writing, record| {
            let at = location(parser, record)?;
            let tile = tiles.entry(at.tile).or_default();
            if !selector.select(&at, tile) {
                return Ok(());
            }
            if opts.remove_duplicate_information {
                record.set_flags(record.flags() & !DUPLICATE);
            }
            record_edit::check_reference_ids(record, writing.references)?;
            writing.restore.apply(record)?;
            writing.out.write(record, Mode::Rewrite)
        },
    )?;
    writing.out.finish(&writing.header, &opts.write)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_tile_widens_by_its_span_over_its_count() {
        let mut tile = Coord::default();
        for (x, y) in [(1000, 2000), (5000, 9000), (3000, 4000)] {
            tile.add(&Location {
                read_group: -1,
                tile: 1101,
                x,
                y,
            });
        }
        assert_eq!(
            (tile.min_x, tile.max_x, tile.min_y, tile.max_y),
            (0, 5000, 0, 9000)
        );
        tile.widen();
        assert_eq!(
            (tile.min_x, tile.max_x, tile.min_y, tile.max_y),
            (-1666, 6666, -3000, 12000)
        );
    }

    #[test]
    fn the_circle_keeps_its_middle_or_its_outside() {
        let tile = Coord {
            min_x: 0,
            min_y: 0,
            max_x: 100,
            max_y: 100,
            count: 4,
        };
        let at = |x, y| Location {
            read_group: -1,
            tile: 1,
            x,
            y,
        };
        let small = CircleSelector::new(0.1);
        let large = CircleSelector::new(0.9);
        // The circle sits at (offset, offset) on the unit torus; 0.9 keeps
        // what 0.1 drops.
        let centre = (small.offset * 100.0).round() as i32;
        for (x, y, inside) in [(centre, centre, true), (centre + 50, centre + 50, false)] {
            assert_eq!(small.select(&at(x, y), &tile), inside);
            assert_eq!(large.select(&at(x, y), &tile), !inside);
        }
    }
}
