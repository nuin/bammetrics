//! `CollectSexCheckMetrics` (docs/06 section 2): reads per base on the X
//! contig, the Y contig and the autosomes, the two ratios to the autosomes,
//! and the sex they imply, against a declared sex. Not a Picard tool.
//!
//! One pass, one `ALL_READS` row. A read counts when it is primary (not
//! secondary, not supplementary), mapped, not a duplicate, not QC-failed,
//! and at or above `MINIMUM_MAPPING_QUALITY` (a missing MAPQ, 255, passes).
//! With `TARGET_INTERVALS` the lengths are the merged target bases per
//! contig set and a read counts only when its alignment overlaps a target.

use std::path::PathBuf;

use noodles::{bam, sam};

use crate::blocks::alignment_end;
use crate::error::{Error, Result};
use crate::intervals::{self, OverlapIndex};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::pass::{self, Accumulator};
use crate::record::position;
use crate::registry::{Invocation, Opt, Tool};

pub const DEFAULT_MINIMUM_MAPPING_QUALITY: u8 = 20;

/// The metrics class string. Not a Picard class: the last name is what
/// readers key on.
pub const CLASS: &str = "bammetrics.qc.SexCheckMetrics";

/// `FEMALE` when `X_RATIO` is at least this and `Y_RATIO` under [`FEMALE_MAX_Y`].
pub const FEMALE_MIN_X: f64 = 0.75;
pub const FEMALE_MAX_Y: f64 = 0.05;
/// `MALE` when `X_RATIO` is within these and `Y_RATIO` at least [`MALE_MIN_Y`].
pub const MALE_MIN_X: f64 = 0.3;
pub const MALE_MAX_X: f64 = 0.75;
pub const MALE_MIN_Y: f64 = 0.1;

/// A declared or inferred sex.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sex {
    Female,
    Male,
}

impl Sex {
    pub fn name(self) -> &'static str {
        match self {
            Sex::Female => "FEMALE",
            Sex::Male => "MALE",
        }
    }

    /// `MALE` or `FEMALE`, any case.
    pub fn parse(text: &str) -> Option<Self> {
        match text.to_ascii_uppercase().as_str() {
            "MALE" => Some(Sex::Male),
            "FEMALE" => Some(Sex::Female),
            _ => None,
        }
    }
}

/// The thresholds of docs/06 section 2.3 on the two ratios. `None` is
/// `UNKNOWN`, with the reason.
pub fn infer(x_ratio: f64, y_ratio: f64) -> std::result::Result<Sex, String> {
    if !x_ratio.is_finite() || !y_ratio.is_finite() {
        return Err("a ratio is undefined".to_owned());
    }
    if x_ratio >= FEMALE_MIN_X && y_ratio < FEMALE_MAX_Y {
        return Ok(Sex::Female);
    }
    if (MALE_MIN_X..=MALE_MAX_X).contains(&x_ratio) && y_ratio >= MALE_MIN_Y {
        return Ok(Sex::Male);
    }
    Err(format!(
        "X_RATIO {} and Y_RATIO {} fit neither pattern",
        crate::metrics_file::format_double(x_ratio),
        crate::metrics_file::format_double(y_ratio)
    ))
}

/// The one row of the metrics file.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub autosome_reads: i64,
    pub x_reads: i64,
    pub y_reads: i64,
    pub autosome_length: i64,
    pub x_length: i64,
    pub y_length: i64,
    pub x_ratio: f64,
    pub y_ratio: f64,
    pub inferred_sex: String,
    pub declared_sex: Option<String>,
    pub sex_match: Option<String>,
    pub reason: Option<String>,
}

impl MetricRow for Row {
    const CLASS: &'static str = CLASS;
    const COLUMNS: &'static [&'static str] = &[
        "AUTOSOME_READS",
        "X_READS",
        "Y_READS",
        "AUTOSOME_LENGTH",
        "X_LENGTH",
        "Y_LENGTH",
        "X_RATIO",
        "Y_RATIO",
        "INFERRED_SEX",
        "DECLARED_SEX",
        "SEX_MATCH",
        "REASON",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.autosome_reads.into(),
            self.x_reads.into(),
            self.y_reads.into(),
            self.autosome_length.into(),
            self.x_length.into(),
            self.y_length.into(),
            self.x_ratio.into(),
            self.y_ratio.into(),
            self.inferred_sex.clone().into(),
            self.declared_sex.clone().into(),
            self.sex_match.clone().into(),
            self.reason.clone().into(),
        ]
    }
}

/// The counts of one contig set.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Set {
    reads: i64,
    length: i64,
}

impl Set {
    fn per_base(self) -> f64 {
        self.reads as f64 / self.length as f64
    }
}

/// The row from the three sets and the declaration.
pub fn row(autosomes: (i64, i64), x: (i64, i64), y: (i64, i64), declared: Option<Sex>) -> Row {
    let [autosomes, x, y] = [autosomes, x, y].map(|(reads, length)| Set { reads, length });
    let mut reasons: Vec<String> = Vec::new();
    for (set, what) in [(autosomes, "autosome"), (x, "X"), (y, "Y")] {
        if set.length == 0 {
            reasons.push(format!("no {what} bases"));
        }
    }
    if autosomes.length > 0 && autosomes.reads == 0 {
        reasons.push("no autosome reads".to_owned());
    }
    let x_ratio = x.per_base() / autosomes.per_base();
    let y_ratio = y.per_base() / autosomes.per_base();
    let inferred = match infer(x_ratio, y_ratio) {
        Ok(sex) => Ok(sex),
        Err(why) => {
            if reasons.is_empty() {
                reasons.push(why);
            }
            Err(())
        }
    };
    let sex_match = match (inferred, declared) {
        (_, None) => None,
        (Ok(i), Some(d)) if i == d => Some("Y".to_owned()),
        (Ok(_), Some(_)) => Some("N".to_owned()),
        (Err(()), Some(d)) => {
            reasons.push(format!(
                "{} declared, the inferred sex is UNKNOWN",
                d.name()
            ));
            None
        }
    };
    Row {
        autosome_reads: autosomes.reads,
        x_reads: x.reads,
        y_reads: y.reads,
        autosome_length: autosomes.length,
        x_length: x.length,
        y_length: y.length,
        x_ratio,
        y_ratio,
        inferred_sex: inferred.map_or("UNKNOWN", Sex::name).to_owned(),
        declared_sex: declared.map(|d| d.name().to_owned()),
        sex_match,
        reason: (!reasons.is_empty()).then(|| reasons.join("; ")),
    }
}

/// The options of one run.
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub expected_sex: Option<Sex>,
    pub minimum_mapping_quality: u8,
    pub target_intervals: Vec<PathBuf>,
    pub x_contig: Option<String>,
    pub y_contig: Option<String>,
    pub autosomes: Vec<String>,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let expected_sex = match inv.get("EXPECTED_SEX") {
            None => None,
            Some(text) => Some(Sex::parse(text).ok_or_else(|| {
                Error::Usage(format!("EXPECTED_SEX={text}: expected MALE or FEMALE"))
            })?),
        };
        let mapq: i64 = inv.parse(
            "MINIMUM_MAPPING_QUALITY",
            i64::from(DEFAULT_MINIMUM_MAPPING_QUALITY),
        )?;
        Ok(Self {
            expected_sex,
            minimum_mapping_quality: u8::try_from(mapq)
                .map_err(|_| Error::Usage(format!("MINIMUM_MAPPING_QUALITY={mapq}: 0 to 255")))?,
            target_intervals: inv
                .all("TARGET_INTERVALS")
                .iter()
                .map(PathBuf::from)
                .collect(),
            x_contig: inv.get("X_CONTIG").map(str::to_owned),
            y_contig: inv.get("Y_CONTIG").map(str::to_owned),
            autosomes: inv.all("AUTOSOMES").to_vec(),
        })
    }
}

/// Which set a contig belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Group {
    Autosome,
    X,
    Y,
}

/// The default names: `chrX` or `X`, `chrY` or `Y`, `chr1`..`chr22` or
/// `1`..`22`.
fn default_group(name: &str) -> Option<Group> {
    let bare = name.strip_prefix("chr").unwrap_or(name);
    match bare {
        "X" => Some(Group::X),
        "Y" => Some(Group::Y),
        _ => bare
            .parse::<u8>()
            .ok()
            .filter(|n| (1..=22).contains(n) && bare == n.to_string())
            .map(|_| Group::Autosome),
    }
}

/// The group of each contig of the header, by index.
fn groups(opts: &Options, header: &sam::Header) -> Vec<Option<Group>> {
    header
        .reference_sequences()
        .keys()
        .map(|name| {
            let name = name.to_string();
            if opts.x_contig.as_deref() == Some(name.as_str()) {
                Some(Group::X)
            } else if opts.y_contig.as_deref() == Some(name.as_str()) {
                Some(Group::Y)
            } else if opts.autosomes.contains(&name) {
                Some(Group::Autosome)
            } else {
                match default_group(&name) {
                    Some(Group::X) if opts.x_contig.is_none() => Some(Group::X),
                    Some(Group::Y) if opts.y_contig.is_none() => Some(Group::Y),
                    Some(Group::Autosome) if opts.autosomes.is_empty() => Some(Group::Autosome),
                    _ => None,
                }
            }
        })
        .collect()
}

/// `CollectSexCheckMetrics` as a pass [`Accumulator`].
pub struct SexCheck {
    groups: Vec<Option<Group>>,
    sets: [Set; 3],
    targets: Option<OverlapIndex>,
    minimum_mapping_quality: u8,
    expected_sex: Option<Sex>,
}

impl SexCheck {
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        let groups = groups(opts, header);
        let mut sets = [Set::default(); 3];
        let targets = if opts.target_intervals.is_empty() {
            for (index, contig) in header.reference_sequences().values().enumerate() {
                if let Some(group) = groups[index] {
                    sets[group as usize].length += contig.length().get() as i64;
                }
            }
            None
        } else {
            let unique =
                intervals::uniqued(&intervals::read_files(&opts.target_intervals, header)?);
            for interval in &unique {
                if let Some(group) = groups[interval.contig] {
                    sets[group as usize].length += interval.length();
                }
            }
            Some(OverlapIndex::new(unique, groups.len(), 0))
        };
        Ok(Self {
            groups,
            sets,
            targets,
            minimum_mapping_quality: opts.minimum_mapping_quality,
            expected_sex: opts.expected_sex,
        })
    }

    fn counted(&self, record: &bam::Record) -> Result<Option<Group>> {
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            return Ok(None);
        }
        if record
            .mapping_quality()
            .is_some_and(|q| u8::from(q) < self.minimum_mapping_quality)
        {
            return Ok(None);
        }
        let Some(contig) = record.reference_sequence_id().transpose()? else {
            return Ok(None);
        };
        let Some(group) = self.groups.get(contig).copied().flatten() else {
            return Ok(None);
        };
        if let Some(targets) = &self.targets {
            let start = position(record, record.alignment_start(), "missing alignment start")?;
            let end = alignment_end(record)?.unwrap_or(start);
            if targets.overlaps(contig, start, end).is_empty() {
                return Ok(None);
            }
        }
        Ok(Some(group))
    }

    /// The reads and bases of a set as `(reads, length)`.
    fn set(&self, group: Group) -> (i64, i64) {
        let set = self.sets[group as usize];
        (set.reads, set.length)
    }
}

impl Accumulator for SexCheck {
    type Output = Row;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if let Some(group) = self.counted(record)? {
            self.sets[group as usize].reads += 1;
        }
        Ok(())
    }

    fn finish(self) -> Result<Row> {
        Ok(row(
            self.set(Group::Autosome),
            self.set(Group::X),
            self.set(Group::Y),
            self.expected_sex,
        ))
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectSexCheckMetrics",
    summary: "one row of the 12 SexCheckMetrics columns: reads per base on the X contig, the Y \
              contig and the autosomes, their ratios, the sex they imply (FEMALE, MALE or \
              UNKNOWN) and whether it matches EXPECTED_SEX. Not a Picard tool (docs/06).",
    options: &[
        Opt::text(
            &["EXPECTED_SEX"],
            "MALE or FEMALE, the declared sex; fills DECLARED_SEX and SEX_MATCH.",
        ),
        Opt::int(
            &["MINIMUM_MAPPING_QUALITY", "MQ"],
            "Reads below this MAPQ are not counted (default 20).",
        ),
        Opt::path(
            &["TARGET_INTERVALS", "TI"],
            "Capture targets, an interval list or BED; repeat for several. Lengths are then the \
             target bases per contig and only reads over a target count.",
        )
        .repeated(),
        Opt::text(&["X_CONTIG"], "The X contig name (default chrX or X)."),
        Opt::text(&["Y_CONTIG"], "The Y contig name (default chrY or Y)."),
        Opt::text(
            &["AUTOSOMES"],
            "An autosome contig name; repeat for each (default chr1 to chr22 or 1 to 22).",
        )
        .repeated(),
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
    Ok(pass::boxed(SexCheck::new(&opts, header)?, |row| {
        Ok(MetricsOutput::from_rows(&[row], None).into())
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thresholds_give_female_male_and_unknown() {
        assert_eq!(infer(1.0, 0.01), Ok(Sex::Female));
        assert_eq!(infer(0.75, 0.049), Ok(Sex::Female));
        assert_eq!(infer(0.5, 0.3), Ok(Sex::Male));
        assert_eq!(infer(0.3, 0.1), Ok(Sex::Male));
        assert_eq!(infer(0.75, 0.1), Ok(Sex::Male));
        assert!(infer(0.75, 0.07).is_err(), "Y between the two bands");
        assert!(infer(0.6, 0.02).is_err(), "X low with no Y");
        assert!(infer(1.0, 0.2).is_err(), "X high with Y");
        assert!(infer(0.2, 0.5).is_err(), "X under the male band");
        assert_eq!(
            infer(0.6, 0.02).unwrap_err(),
            "X_RATIO 0.6 and Y_RATIO 0.02 fit neither pattern"
        );
        assert_eq!(infer(f64::NAN, 0.0).unwrap_err(), "a ratio is undefined");
    }

    #[test]
    fn rows_carry_the_match_and_the_reason() {
        // 1 read per 10 bases on the autosomes, half that on X, a fifth on Y.
        let male = row((1000, 10_000), (50, 1000), (20, 1000), Some(Sex::Male));
        assert_eq!(male.inferred_sex, "MALE");
        assert_eq!(male.x_ratio, 0.5);
        assert!((male.y_ratio - 0.2).abs() < 1e-12);
        assert_eq!(male.declared_sex.as_deref(), Some("MALE"));
        assert_eq!(male.sex_match.as_deref(), Some("Y"));
        assert_eq!(male.reason, None);
        let mismatch = row((1000, 10_000), (100, 1000), (0, 1000), Some(Sex::Male));
        assert_eq!(mismatch.inferred_sex, "FEMALE");
        assert_eq!(mismatch.sex_match.as_deref(), Some("N"));
        let undeclared = row((1000, 10_000), (100, 1000), (0, 1000), None);
        assert_eq!(undeclared.declared_sex, None);
        assert_eq!(undeclared.sex_match, None);
        let unknown = row((1000, 10_000), (60, 1000), (2, 1000), Some(Sex::Female));
        assert_eq!(unknown.inferred_sex, "UNKNOWN");
        assert_eq!(unknown.sex_match, None, "UNKNOWN is not a mismatch");
        assert_eq!(
            unknown.reason.as_deref(),
            Some(
                "X_RATIO 0.6 and Y_RATIO 0.02 fit neither pattern; FEMALE declared, the inferred sex is UNKNOWN"
            )
        );
        let no_y = row((1000, 10_000), (100, 1000), (0, 0), None);
        assert_eq!(no_y.inferred_sex, "UNKNOWN");
        assert!(no_y.y_ratio.is_nan());
        assert_eq!(no_y.reason.as_deref(), Some("no Y bases"));
        let empty = row((0, 10_000), (0, 1000), (0, 1000), None);
        assert_eq!(empty.reason.as_deref(), Some("no autosome reads"));
        assert_eq!(empty.values().len(), Row::COLUMNS.len());
        assert_eq!(empty.values()[8], Value::from("UNKNOWN"));
        assert_eq!(empty.values()[10], Value::Null);
    }

    #[test]
    fn contig_names_are_recognised_exactly() {
        assert_eq!(default_group("chrX"), Some(Group::X));
        assert_eq!(default_group("X"), Some(Group::X));
        assert_eq!(default_group("chrY"), Some(Group::Y));
        assert_eq!(default_group("Y"), Some(Group::Y));
        assert_eq!(default_group("chr1"), Some(Group::Autosome));
        assert_eq!(default_group("22"), Some(Group::Autosome));
        for other in [
            "chr23",
            "0",
            "chrM",
            "MT",
            "chr01",
            "x",
            "chrX_random",
            "chrUn_gl000220",
        ] {
            assert_eq!(default_group(other), None, "{other}");
        }
    }

    #[test]
    fn options_parse_and_reject_a_bad_sex() {
        let strings =
            |args: &[&str]| -> Vec<String> { args.iter().map(|s| (*s).to_owned()).collect() };
        let opts = Options::from_invocation(
            &TOOL
                .parse(&strings(&[
                    "I=x.bam",
                    "EXPECTED_SEX=female",
                    "MQ=0",
                    "X_CONTIG=23",
                    "AUTOSOMES=1",
                    "AUTOSOMES=2",
                ]))
                .expect("parses"),
        )
        .expect("options");
        assert_eq!(opts.expected_sex, Some(Sex::Female));
        assert_eq!(opts.minimum_mapping_quality, 0);
        assert_eq!(opts.x_contig.as_deref(), Some("23"));
        assert_eq!(opts.autosomes, ["1", "2"]);
        let err = TOOL
            .parse(&strings(&["I=x.bam", "EXPECTED_SEX=unknown"]))
            .and_then(|inv| Options::from_invocation(&inv))
            .expect_err("bad sex");
        assert!(err.to_string().contains("MALE or FEMALE"), "{err}");
        let err = TOOL
            .parse(&strings(&["I=x.bam", "MQ=300"]))
            .and_then(|inv| Options::from_invocation(&inv))
            .expect_err("bad mapq");
        assert!(err.to_string().contains("0 to 255"), "{err}");
    }
}
