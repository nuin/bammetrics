//! `CollectDuplicateByTileMetrics` (docs/06 section 7): duplicate rates per
//! flowcell lane and tile from the read names and the duplicate flags
//! already in the file. Not a Picard tool.
//!
//! One pass. The lane and tile come from the read name through
//! `optical.rs` `ReadNameParser`, as `MarkDuplicates` reads them:
//! the default 5-field (`instrument:lane:tile:x:y`) and 7-field
//! (`instrument:run:flowcell:lane:tile:x:y`) Illumina names, or a
//! `READ_NAME_REGEX` with tile, x and y groups (then `FLOWCELL` and `LANE`
//! are empty). The name is read one char per byte, as `MarkDuplicates`
//! reads it. The tile is Picard's 16-bit value. A name that gives no
//! location, or whose flowcell or lane holds an ASCII control character,
//! counts on the `UNPARSED` row.
//!
//! Records are counted as `CollectDuplicateMetrics` counts them, so the
//! `ALL` row reconciles with its totals: an unmapped record in `UNMAPPED`,
//! then secondary and supplementary records skipped, then a mapped record
//! unpaired or with an unmapped mate in `UNPAIRED_READS`. A pair with both
//! mates mapped is counted once, on its first-of-pair record, in `PAIRS`.
//! The duplicate flag adds to `UNPAIRED_DUPLICATES` or `PAIR_DUPLICATES`
//! by the same rule. `OPTICAL_DUPLICATES` counts duplicate pairs whose
//! first-of-pair record carries `DT:Z:SQ`, the tag `MarkDuplicates` writes
//! under `TAGGING_POLICY=All` or `OpticalOnly`; without any `DT` tag and
//! without a `MarkDuplicates` `@PG` line naming one of those policies the
//! column is empty. Rows: each lane's tiles in tile order then its `LANE`
//! row, lanes by flowcell then lane, then `UNPARSED`, then `ALL`.
//!
//! `TILE_OUTLIER` is `Y` for a tile whose `PCT_DUPLICATION` is above
//! `median + max(OUTLIER_SD x 1.4826 x MAD, 0.005)` of the lane's tiles
//! (`MIN_OUTLIER_EXCESS`), with at least
//! 3 tiles that examined a read; tiles that examined none are left out.

use std::borrow::Cow;
use std::collections::HashMap;

use noodles::sam::alignment::record::data::field::{Tag, Value as FieldValue};
use noodles::sam::header::record::value::map::program::tag;
use noodles::{bam, sam};

use crate::error::{Error, Result};
use crate::mark_duplicates::latin1;
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::optical::ReadNameParser;
use crate::pass::{self, Accumulator};
use crate::registry::{Invocation, Opt, Tool};

/// The metrics class string. Not a Picard class: readers key on the last
/// name.
pub const CLASS: &str = "bammetrics.qc.DuplicateByTileMetrics";

pub const DEFAULT_OUTLIER_SD: f64 = 3.0;

/// `1 / Phi^-1(3/4)`: the MAD of a normal sample times this is its
/// standard deviation.
pub const MAD_TO_SD: f64 = 1.4826;

/// The least `PCT_DUPLICATION` a tile must stand above its lane's median
/// to be an outlier, so a lane of near-equal tiles (MAD 0) flags none.
pub const MIN_OUTLIER_EXCESS: f64 = 0.005;

/// Tiles of one lane that must have examined a read before any is judged.
pub const MIN_TILES_FOR_OUTLIERS: usize = 3;

/// Lane and tile pairs one run may hold. A NovaSeq S4 flowcell has 3,744
/// tiles; more than this means the names are not Illumina names.
pub const MAX_TILES: usize = 100_000;

const UNMAPPED: u16 = 0x4;
const PAIRED: u16 = 0x1;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST: u16 = 0x40;
const SECONDARY_OR_SUPPLEMENTARY: u16 = 0x100 | 0x800;
const DUPLICATE: u16 = 0x400;
const DT: Tag = Tag::new(b'D', b'T');

/// What a row sums.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unit {
    Tile,
    Lane,
    Unparsed,
    All,
}

impl Unit {
    pub fn name(self) -> &'static str {
        match self {
            Unit::Tile => "TILE",
            Unit::Lane => "LANE",
            Unit::Unparsed => "UNPARSED",
            Unit::All => "ALL",
        }
    }
}

/// The smallest and largest x and y seen.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Extent {
    pub min_x: i32,
    pub max_x: i32,
    pub min_y: i32,
    pub max_y: i32,
}

impl Extent {
    fn point(x: i32, y: i32) -> Self {
        Self {
            min_x: x,
            max_x: x,
            min_y: y,
            max_y: y,
        }
    }

    fn merge(self, other: Self) -> Self {
        Self {
            min_x: self.min_x.min(other.min_x),
            max_x: self.max_x.max(other.max_x),
            min_y: self.min_y.min(other.min_y),
            max_y: self.max_y.max(other.max_y),
        }
    }
}

/// The counts of one row.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub reads: i64,
    pub unmapped: i64,
    pub unpaired_reads: i64,
    pub pairs: i64,
    pub unpaired_duplicates: i64,
    pub pair_duplicates: i64,
    pub optical_duplicates: i64,
    pub extent: Option<Extent>,
}

impl Counts {
    fn add(&mut self, other: &Counts) {
        self.reads += other.reads;
        self.unmapped += other.unmapped;
        self.unpaired_reads += other.unpaired_reads;
        self.pairs += other.pairs;
        self.unpaired_duplicates += other.unpaired_duplicates;
        self.pair_duplicates += other.pair_duplicates;
        self.optical_duplicates += other.optical_duplicates;
        self.extent = match (self.extent, other.extent) {
            (Some(a), Some(b)) => Some(a.merge(b)),
            (a, b) => a.or(b),
        };
    }

    fn see(&mut self, x: i32, y: i32) {
        let point = Extent::point(x, y);
        self.extent = Some(self.extent.map_or(point, |e| e.merge(point)));
    }

    /// Duplicate reads: each duplicate pair counts its two reads.
    pub fn duplicates(&self) -> i64 {
        self.unpaired_duplicates + 2 * self.pair_duplicates
    }

    /// Reads examined for duplication, as `CollectDuplicateMetrics`
    /// divides by them.
    pub fn examined(&self) -> i64 {
        self.unpaired_reads + 2 * self.pairs
    }

    /// `PERCENT_DUPLICATION` of these counts: 0 without an examined read.
    pub fn pct_duplication(&self) -> f64 {
        match self.examined() {
            0 => 0.0,
            examined => self.duplicates() as f64 / examined as f64,
        }
    }

    /// Optical duplicate pairs over examined pairs: 0 without a pair.
    pub fn pct_optical(&self) -> f64 {
        match self.pairs {
            0 => 0.0,
            pairs => self.optical_duplicates as f64 / pairs as f64,
        }
    }
}

/// One row of the metrics file.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub unit: Unit,
    pub flowcell: Option<String>,
    pub lane: Option<String>,
    pub tile: Option<i64>,
    pub counts: Counts,
    /// False when the file carries no optical duplicate marks.
    pub optical_known: bool,
    pub outlier_limit: Option<f64>,
    pub tile_outlier: Option<bool>,
    pub notes: Option<String>,
}

impl MetricRow for Row {
    const CLASS: &'static str = CLASS;
    const COLUMNS: &'static [&'static str] = &[
        "UNIT",
        "FLOWCELL",
        "LANE",
        "TILE",
        "READS",
        "UNMAPPED",
        "UNPAIRED_READS",
        "PAIRS",
        "UNPAIRED_DUPLICATES",
        "PAIR_DUPLICATES",
        "DUPLICATES",
        "OPTICAL_DUPLICATES",
        "PCT_DUPLICATION",
        "PCT_OPTICAL",
        "OUTLIER_LIMIT",
        "TILE_OUTLIER",
        "MIN_X",
        "MAX_X",
        "MIN_Y",
        "MAX_Y",
        "NOTES",
    ];

    fn values(&self) -> Vec<Value> {
        let c = &self.counts;
        let extent = |f: fn(&Extent) -> i32| c.extent.as_ref().map(|e| i64::from(f(e)));
        vec![
            self.unit.name().into(),
            self.flowcell.clone().into(),
            self.lane.clone().into(),
            self.tile.into(),
            c.reads.into(),
            c.unmapped.into(),
            c.unpaired_reads.into(),
            c.pairs.into(),
            c.unpaired_duplicates.into(),
            c.pair_duplicates.into(),
            c.duplicates().into(),
            self.optical_known.then_some(c.optical_duplicates).into(),
            c.pct_duplication().into(),
            self.optical_known.then(|| c.pct_optical()).into(),
            self.outlier_limit.into(),
            self.tile_outlier.map(|y| if y { "Y" } else { "N" }).into(),
            extent(|e| e.min_x).into(),
            extent(|e| e.max_x).into(),
            extent(|e| e.min_y).into(),
            extent(|e| e.max_y).into(),
            self.notes.clone().into(),
        ]
    }
}

/// The options of one run.
#[derive(Clone, Debug)]
pub struct Options {
    pub parser: ReadNameParser,
    pub outlier_sd: f64,
}

impl Options {
    pub(crate) fn from_invocation(inv: &Invocation) -> Result<Self> {
        let outlier_sd: f64 = inv.parse("OUTLIER_SD", DEFAULT_OUTLIER_SD)?;
        if !outlier_sd.is_finite() || outlier_sd < 0.0 {
            return Err(Error::Usage(format!(
                "OUTLIER_SD={outlier_sd}: a finite number of at least 0"
            )));
        }
        Ok(Self {
            parser: ReadNameParser::new(inv.get("READ_NAME_REGEX"))?,
            outlier_sd,
        })
    }
}

/// The byte ranges of a default-form name's flowcell and lane fields, and
/// the key of its lane: the bytes from the flowcell through the lane (7
/// fields, always holding a `:`) or the lane alone (5 fields, never
/// holding one). `None` for another field count.
fn lane_fields(name: &[u8]) -> Option<(&[u8], &[u8], &[u8])> {
    let mut colons = [0usize; 6];
    let mut count = 0usize;
    for (i, &b) in name.iter().enumerate() {
        if b == b':' {
            if count == colons.len() {
                return None;
            }
            colons[count] = i;
            count += 1;
        }
    }
    match count {
        6 => {
            let flowcell = &name[colons[1] + 1..colons[2]];
            let lane = &name[colons[2] + 1..colons[3]];
            Some((&name[colons[1] + 1..colons[3]], flowcell, lane))
        }
        4 => {
            let lane = &name[colons[0] + 1..colons[1]];
            Some((lane, &[], lane))
        }
        _ => None,
    }
}

/// The tiles of one flowcell lane.
#[derive(Clone, Debug)]
struct Lane {
    flowcell: String,
    lane: String,
    tiles: HashMap<i16, Counts>,
}

/// What the header says about duplicate marking.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Marking {
    /// An `@PG` line whose ID, PN or CL names a duplicate marker.
    pub marker: bool,
    /// Such a line's CL sets `TAGGING_POLICY` to `All` or `OpticalOnly`.
    pub tagging: bool,
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// True when a lower-case command line sets `TAGGING_POLICY` (or
/// `--tagging-policy`) to `all` or `opticalonly`.
fn tags_optical(command_line: &[u8]) -> bool {
    let mut rest = command_line;
    while let Some(at) = rest
        .windows(14)
        .position(|w| w == b"tagging_policy" || w == b"tagging-policy")
    {
        rest = &rest[at + 14..];
        let value = rest.trim_ascii_start();
        let value = value.strip_prefix(b"=").unwrap_or(value).trim_ascii_start();
        let end = value
            .iter()
            .position(|b| !b.is_ascii_alphabetic())
            .unwrap_or(value.len());
        if matches!(&value[..end], b"all" | b"opticalonly") {
            return true;
        }
    }
    false
}

/// The `@PG` lines of `header` read for duplicate marking.
pub fn marking(header: &sam::Header) -> Marking {
    let mut marking = Marking::default();
    for (id, map) in header.programs().as_ref() {
        let fields = map.other_fields();
        let lower = |bytes: &[u8]| bytes.to_ascii_lowercase();
        let id = lower(id.as_ref());
        let name = fields.get(&tag::NAME).map(|v| lower(v.as_ref()));
        let command_line = fields.get(&tag::COMMAND_LINE).map(|v| lower(v.as_ref()));
        let names_marker = [Some(&id), name.as_ref(), command_line.as_ref()]
            .into_iter()
            .flatten()
            .any(|text| contains(text, b"markdup"));
        if names_marker {
            marking.marker = true;
            if command_line.as_deref().is_some_and(tags_optical) {
                marking.tagging = true;
            }
        }
    }
    marking
}

/// The rows and the notes of one run.
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub rows: Vec<Row>,
    pub notes: Vec<String>,
}

/// `CollectDuplicateByTileMetrics` as a pass [`Accumulator`].
pub struct DuplicateByTile {
    parser: ReadNameParser,
    outlier_sd: f64,
    marking: Marking,
    lanes: HashMap<Vec<u8>, Lane>,
    tiles: usize,
    unparsed: Counts,
    /// Records sent to `UNPARSED` for a control character in the
    /// flowcell or lane field of a name that parses.
    control_names: i64,
    any_duplicate: bool,
    any_dt: bool,
}

impl DuplicateByTile {
    pub fn new(opts: &Options, header: &sam::Header) -> Self {
        Self {
            parser: opts.parser.clone(),
            outlier_sd: opts.outlier_sd,
            marking: marking(header),
            lanes: HashMap::new(),
            tiles: 0,
            unparsed: Counts::default(),
            control_names: 0,
            any_duplicate: false,
            any_dt: false,
        }
    }

    /// The counts of the row `name` belongs to; `MAX_TILES` bounds the
    /// rows a file can make.
    fn counts_of(&mut self, name: &[u8]) -> Result<&mut Counts> {
        // One char per byte, as `MarkDuplicates` reads the name.
        let text = match std::str::from_utf8(name) {
            Ok(text) if text.is_ascii() => Cow::Borrowed(text),
            _ => Cow::Owned(latin1(name)),
        };
        let Some((tile, x, y)) = self.parser.location(&text) else {
            return Ok(&mut self.unparsed);
        };
        let (key, flowcell, lane) = match self.parser {
            ReadNameParser::Default => lane_fields(name).unwrap_or((&[], &[], &[])),
            _ => (&[][..], &[][..], &[][..]),
        };
        // A tab or line end in FLOWCELL or LANE would break the row.
        if key.iter().any(u8::is_ascii_control) {
            self.control_names += 1;
            return Ok(&mut self.unparsed);
        }
        if !self.lanes.contains_key(key) {
            self.lanes.insert(
                key.to_vec(),
                Lane {
                    flowcell: latin1(flowcell),
                    lane: latin1(lane),
                    tiles: HashMap::new(),
                },
            );
        }
        let lane = self.lanes.get_mut(key).expect("inserted above");
        if !lane.tiles.contains_key(&tile) {
            if self.tiles == MAX_TILES {
                return Err(Error::Metric(format!(
                    "more than {MAX_TILES} lane and tile pairs in the read names; \
                     READ_NAME_REGEX does not fit these names"
                )));
            }
            self.tiles += 1;
        }
        let counts = lane.tiles.entry(tile).or_default();
        counts.see(x, y);
        Ok(counts)
    }

    fn report(self) -> Report {
        let mut lanes: Vec<Lane> = self.lanes.into_values().collect();
        lanes.sort_by(|a, b| lane_order(a).cmp(&lane_order(b)));
        let optical_known = self.any_dt || self.marking.tagging;
        let mut rows = Vec::new();
        let mut all = Counts::default();
        for lane in lanes {
            let mut tiles: Vec<(i16, Counts)> = lane.tiles.into_iter().collect();
            tiles.sort_unstable_by_key(|(tile, _)| *tile);
            let judged: Vec<f64> = tiles
                .iter()
                .filter(|(_, c)| c.examined() > 0)
                .map(|(_, c)| c.pct_duplication())
                .collect();
            let limit = outlier_limit(&judged, self.outlier_sd);
            let mut total = Counts::default();
            for (tile, counts) in &tiles {
                total.add(counts);
                let outlier = limit
                    .filter(|_| counts.examined() > 0)
                    .map(|limit| counts.pct_duplication() > limit);
                rows.push(Row {
                    unit: Unit::Tile,
                    flowcell: Some(lane.flowcell.clone()),
                    lane: Some(lane.lane.clone()),
                    tile: Some(i64::from(*tile)),
                    counts: *counts,
                    optical_known,
                    outlier_limit: limit,
                    tile_outlier: outlier,
                    notes: None,
                });
            }
            all.add(&total);
            rows.push(Row {
                unit: Unit::Lane,
                flowcell: Some(lane.flowcell),
                lane: Some(lane.lane),
                tile: None,
                counts: total,
                optical_known,
                outlier_limit: limit,
                tile_outlier: None,
                notes: None,
            });
        }
        all.add(&self.unparsed);
        rows.push(Row {
            unit: Unit::Unparsed,
            flowcell: None,
            lane: None,
            tile: None,
            counts: self.unparsed,
            optical_known,
            outlier_limit: None,
            tile_outlier: None,
            notes: None,
        });
        let mut notes = Vec::new();
        if !self.any_duplicate && !self.marking.marker {
            notes.push(
                "no record carries the duplicate flag and no @PG line names a duplicate \
                 marker: duplicates were probably never marked; run MarkDuplicates first"
                    .to_owned(),
            );
        }
        if !optical_known {
            notes.push(
                "OPTICAL_DUPLICATES is empty: no duplicate carries a DT tag and no MarkDuplicates \
                 @PG line sets TAGGING_POLICY; run MarkDuplicates with TAGGING_POLICY=All (or \
                 OpticalOnly) to fill it"
                    .to_owned(),
            );
        }
        let unread = self.unparsed.reads - self.control_names;
        if unread > 0 {
            let parser = match self.parser {
                ReadNameParser::Default => "the 5- and 7-field Illumina name forms",
                ReadNameParser::Regex(_) => "READ_NAME_REGEX",
                ReadNameParser::None => "READ_NAME_REGEX=null",
            };
            notes.push(format!(
                "{unread} reads have names without a lane and tile under {parser}; they are \
                 on the UNPARSED row"
            ));
        }
        if self.control_names > 0 {
            notes.push(format!(
                "{} reads have a control character (tab, line end) in the flowcell or lane of \
                 their name; they are on the UNPARSED row",
                self.control_names
            ));
        }
        rows.push(Row {
            unit: Unit::All,
            flowcell: None,
            lane: None,
            tile: None,
            counts: all,
            optical_known,
            outlier_limit: None,
            tile_outlier: None,
            notes: (!notes.is_empty()).then(|| notes.join("; ")),
        });
        Report { rows, notes }
    }
}

/// Lanes by flowcell, then by number, then lanes that are no number by
/// text.
fn lane_order(lane: &Lane) -> (&str, u8, u64, &str) {
    match lane.lane.parse::<u64>() {
        Ok(number) => (&lane.flowcell, 0, number, ""),
        Err(_) => (&lane.flowcell, 1, 0, &lane.lane),
    }
}

/// The median of `values`, which must not be empty: the mean of the two
/// middle values for an even count.
fn median(values: &mut [f64]) -> f64 {
    values.sort_unstable_by(f64::total_cmp);
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    }
}

/// `median + max(outlier_sd x 1.4826 x MAD, MIN_OUTLIER_EXCESS)` of the
/// tiles' `PCT_DUPLICATION`, or `None` with fewer than
/// [`MIN_TILES_FOR_OUTLIERS`] tiles.
pub fn outlier_limit(pcts: &[f64], outlier_sd: f64) -> Option<f64> {
    if pcts.len() < MIN_TILES_FOR_OUTLIERS {
        return None;
    }
    let mut values = pcts.to_vec();
    let center = median(&mut values);
    let mut deviations: Vec<f64> = pcts.iter().map(|p| (p - center).abs()).collect();
    let mad = median(&mut deviations);
    Some(center + (outlier_sd * MAD_TO_SD * mad).max(MIN_OUTLIER_EXCESS))
}

impl Accumulator for DuplicateByTile {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = u16::from(record.flags());
        let unmapped = flags & UNMAPPED != 0;
        if !unmapped && flags & SECONDARY_OR_SUPPLEMENTARY != 0 {
            return Ok(());
        }
        let duplicate = !unmapped && flags & DUPLICATE != 0;
        let mut optical = false;
        if duplicate {
            self.any_duplicate = true;
            match record.data().get(&DT) {
                Some(Ok(FieldValue::String(value))) => {
                    self.any_dt = true;
                    optical = AsRef::<[u8]>::as_ref(value) == b"SQ";
                }
                Some(Ok(_)) => self.any_dt = true,
                Some(Err(e)) => return Err(e.into()),
                None => {}
            }
        }
        let name: &[u8] = record.name().map_or(&[], |n| n.as_ref());
        let counts = self.counts_of(name)?;
        counts.reads += 1;
        if unmapped {
            counts.unmapped += 1;
        } else if flags & PAIRED == 0 || flags & MATE_UNMAPPED != 0 {
            counts.unpaired_reads += 1;
            if duplicate {
                counts.unpaired_duplicates += 1;
            }
        } else if flags & FIRST != 0 {
            counts.pairs += 1;
            if duplicate {
                counts.pair_duplicates += 1;
                if optical {
                    counts.optical_duplicates += 1;
                }
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<Report> {
        Ok(self.report())
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectDuplicateByTileMetrics",
    summary: "one DuplicateByTileMetrics row per flowcell lane and tile from the read names, \
              one per lane, one for names without a tile and one for all: reads, pairs and \
              duplicates from the duplicate flags in the file, optical duplicates from \
              MarkDuplicates' DT:Z:SQ tags, and TILE_OUTLIER for a tile whose duplication is \
              far above its lane's median. Not a Picard tool (docs/06).",
    options: &[
        Opt::text(
            &["READ_NAME_REGEX"],
            "Regex with tile, x and y groups, as MarkDuplicates takes it; null for none \
             (default: 5- and 7-field Illumina names, which also give the lane).",
        ),
        Opt::float(
            &["OUTLIER_SD"],
            "A tile is an outlier above its lane's median PCT_DUPLICATION plus this many \
             robust standard deviations, 1.4826 x MAD, and at least 0.005 above it \
             (default 3).",
        ),
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn pass::Collector>> {
    let opts = Options::from_invocation(inv)?;
    let quiet = inv.quiet;
    Ok(pass::boxed(
        DuplicateByTile::new(&opts, header),
        move |report| {
            for note in report.notes.iter().filter(|_| !quiet) {
                eprintln!("bammetrics: {}: {note}", TOOL.name);
            }
            Ok(MetricsOutput::from_rows(&report.rows, None).into())
        },
    ))
}

#[cfg(test)]
mod tests;
