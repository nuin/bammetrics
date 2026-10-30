//! `METRIC_ACCUMULATION_LEVEL`: one accumulator per sample, library or read
//! group, the way htsjdk's `MultiLevelCollector` keeps them.
//!
//! A tool keeps its per-unit state in [`Levels<T>`] and calls
//! [`Levels::route`] for each record; the closure sees every bucket the
//! record belongs to. [`Levels::into_buckets`] returns the buckets in the
//! order Picard writes their rows.
//!
//! Picard 3.5.0 behaviour (checked on a 12-read-group BAM):
//! * Row order is `ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`, whatever
//!   order the options came in. Within a level, buckets follow the order their
//!   key first appears among the `@RG` header lines, not alphabetical order.
//! * `READ_GROUP` keys on the platform unit (`PU`), not on the `ID`.
//! * A record with no `RG` tag, or one the header does not list, falls in the
//!   `unknown` bucket of every level but `ALL_READS`, even for a BAM with no
//!   `@RG` lines.
//! * A read group without `LB` (or `PU`) falls in an `unknown` bucket that
//!   comes after every named bucket of that level. At `LIBRARY` level its
//!   `SAMPLE` and `LIBRARY` columns both print `unknown`; at `READ_GROUP` level
//!   a missing `PU` prints `unknown` in all three columns, while a missing `LB`
//!   with a `PU` prints the sample and an empty library.
//! * `METRIC_ACCUMULATION_LEVEL` is a list option with default `ALL_READS`:
//!   each value is added to it, and `null` clears it.
//! * Histogram columns are named `<key>.<orientation>_count`, with
//!   `All_Reads` as the `ALL_READS` key.

use std::collections::HashMap;

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::header::record::value::map::read_group::tag;

use crate::error::{Error, Result};

/// `MetricAccumulationLevel`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Level {
    AllReads,
    Sample,
    Library,
    ReadGroup,
}

impl Level {
    pub fn parse(text: &str) -> Result<Self> {
        match text {
            "ALL_READS" => Ok(Level::AllReads),
            "SAMPLE" => Ok(Level::Sample),
            "LIBRARY" => Ok(Level::Library),
            "READ_GROUP" => Ok(Level::ReadGroup),
            _ => Err(Error::Usage(format!(
                "METRIC_ACCUMULATION_LEVEL {text:?}: expected ALL_READS, SAMPLE, LIBRARY or READ_GROUP"
            ))),
        }
    }

    /// The level's spelling in the command-line header.
    pub fn name(self) -> &'static str {
        match self {
            Level::AllReads => "ALL_READS",
            Level::Sample => "SAMPLE",
            Level::Library => "LIBRARY",
            Level::ReadGroup => "READ_GROUP",
        }
    }
}

/// The level list of a `METRIC_ACCUMULATION_LEVEL` option: `ALL_READS` plus
/// each value given, with `null` clearing what came before it.
pub fn parse_levels<S: AsRef<str>>(values: &[S]) -> Result<Vec<Level>> {
    let mut levels = vec![Level::AllReads];
    for value in values {
        if value.as_ref() == "null" {
            levels.clear();
        } else {
            let level = Level::parse(value.as_ref())?;
            if !levels.contains(&level) {
                levels.push(level);
            }
        }
    }
    Ok(levels)
}

/// The identity fields of one `@RG` line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadGroupInfo {
    pub id: String,
    pub sample: Option<String>,
    pub library: Option<String>,
    pub platform_unit: Option<String>,
}

/// The `@RG` lines of a header, in header order.
pub fn read_groups(header: &sam::Header) -> Vec<ReadGroupInfo> {
    header
        .read_groups()
        .iter()
        .map(|(id, map)| {
            let field = |t| map.other_fields().get(&t).map(|v| v.to_string());
            ReadGroupInfo {
                id: id.to_string(),
                sample: field(tag::SAMPLE),
                library: field(tag::LIBRARY),
                platform_unit: field(tag::PLATFORM_UNIT),
            }
        })
        .collect()
}

/// The `SAMPLE`, `LIBRARY` and `READ_GROUP` columns of a bucket; `None`
/// prints as an empty field.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Labels {
    pub sample: Option<String>,
    pub library: Option<String>,
    pub read_group: Option<String>,
}

/// One accumulator and what identifies it.
#[derive(Clone, Debug)]
pub struct Bucket<T> {
    pub level: Level,
    /// Prefix of the bucket's histogram columns (`All_Reads`, the sample,
    /// library or platform unit name, or `unknown`).
    pub key: String,
    pub labels: Labels,
    pub state: T,
}

const UNKNOWN: &str = "unknown";

/// The accumulators of every requested level.
pub struct Levels<T> {
    buckets: Vec<Bucket<T>>,
    /// Bucket indices that take every record (the `ALL_READS` bucket).
    all: Vec<usize>,
    /// Read group `ID` to the non-`ALL_READS` buckets it feeds.
    by_read_group: HashMap<Vec<u8>, Vec<usize>>,
    /// The `unknown` bucket of each non-`ALL_READS` level, for records
    /// whose read group is missing or not in the header.
    unknown: Vec<usize>,
}

impl<T> Levels<T> {
    /// Build one bucket per unit of each level in `levels`; `new_state` makes
    /// the empty accumulator for a bucket.
    pub fn new(header: &sam::Header, levels: &[Level], mut new_state: impl FnMut() -> T) -> Self {
        let groups = read_groups(header);
        let mut buckets = Vec::new();
        let mut all = Vec::new();
        let mut by_read_group: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
        let mut unknown = Vec::new();
        for level in [
            Level::AllReads,
            Level::Sample,
            Level::Library,
            Level::ReadGroup,
        ] {
            if !levels.contains(&level) {
                continue;
            }
            if level == Level::AllReads {
                all.push(buckets.len());
                buckets.push(Bucket {
                    level,
                    key: "All_Reads".to_owned(),
                    labels: Labels::default(),
                    state: new_state(),
                });
                continue;
            }
            // Named keys in header order, then the `unknown` bucket.
            let mut order: Vec<Option<&str>> = Vec::new();
            for group in &groups {
                let key = unit_key(level, group);
                if !order.contains(&key) {
                    order.push(key);
                }
            }
            if !order.contains(&None) {
                order.push(None);
            }
            order.sort_by_key(Option::is_none);
            let first = buckets.len();
            unknown.push(first + order.len() - 1);
            for key in &order {
                let representative = groups.iter().find(|g| unit_key(level, g) == *key);
                buckets.push(Bucket {
                    level,
                    key: key.unwrap_or(UNKNOWN).to_owned(),
                    labels: labels_for(level, *key, representative),
                    state: new_state(),
                });
            }
            for group in &groups {
                let key = unit_key(level, group);
                let at = order.iter().position(|k| *k == key).unwrap_or_default();
                by_read_group
                    .entry(group.id.clone().into_bytes())
                    .or_default()
                    .push(first + at);
            }
        }
        Self {
            buckets,
            all,
            by_read_group,
            unknown,
        }
    }

    /// Run `f` on every bucket this record belongs to: the `ALL_READS`
    /// bucket, and the sample, library and read group buckets of its `RG` tag.
    /// A record without a read group in the header reaches the `unknown`
    /// bucket of each other level, as in Picard.
    pub fn route(&mut self, record: &bam::Record, mut f: impl FnMut(&mut T)) {
        for &i in &self.all {
            f(&mut self.buckets[i].state);
        }
        let targets = match record.data().get(&Tag::READ_GROUP) {
            Some(Ok(Value::String(id))) => self.by_read_group.get(AsRef::<[u8]>::as_ref(id)),
            _ => None,
        }
        .unwrap_or(&self.unknown);
        for &i in targets {
            f(&mut self.buckets[i].state);
        }
    }

    /// Buckets in the order Picard writes their rows.
    pub fn into_buckets(self) -> Vec<Bucket<T>> {
        self.buckets
    }
}

/// The key a read group has at `level`: its sample, library or platform
/// unit, `None` when the `@RG` line lacks the field.
fn unit_key(level: Level, group: &ReadGroupInfo) -> Option<&str> {
    match level {
        Level::Sample => group.sample.as_deref(),
        Level::Library => group.library.as_deref(),
        _ => group.platform_unit.as_deref(),
    }
}

fn labels_for(level: Level, key: Option<&str>, group: Option<&ReadGroupInfo>) -> Labels {
    let text = |s: &str| Some(s.to_owned());
    match (level, key, group) {
        (Level::Sample, Some(sample), _) => Labels {
            sample: text(sample),
            ..Labels::default()
        },
        (Level::Sample, None, _) => Labels {
            sample: text(UNKNOWN),
            ..Labels::default()
        },
        (Level::Library, Some(library), Some(g)) => Labels {
            sample: g.sample.clone(),
            library: text(library),
            read_group: None,
        },
        (Level::Library, _, _) => Labels {
            sample: text(UNKNOWN),
            library: text(UNKNOWN),
            read_group: None,
        },
        (Level::ReadGroup, Some(unit), Some(g)) => Labels {
            sample: g.sample.clone(),
            library: g.library.clone(),
            read_group: text(unit),
        },
        _ => Labels {
            sample: text(UNKNOWN),
            library: text(UNKNOWN),
            read_group: text(UNKNOWN),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(text: &str) -> sam::Header {
        text.parse().expect("header")
    }

    const HEADER: &str = "@HD\tVN:1.6\n@SQ\tSN:c\tLN:100\n\
        @RG\tID:b1\tSM:S2\tLB:LibB\tPU:pu-b1\n\
        @RG\tID:a1\tSM:S1\tLB:LibA\tPU:pu-a1\n\
        @RG\tID:a2\tSM:S1\tLB:LibA\n\
        @RG\tID:c1\tSM:S1\tPU:pu-c1\n";

    #[test]
    fn levels_option_is_a_list_with_default_all_reads() {
        let none: [&str; 0] = [];
        assert_eq!(parse_levels(&none).unwrap(), [Level::AllReads]);
        assert_eq!(
            parse_levels(&["SAMPLE", "READ_GROUP"]).unwrap(),
            [Level::AllReads, Level::Sample, Level::ReadGroup]
        );
        assert_eq!(
            parse_levels(&["null", "LIBRARY"]).unwrap(),
            [Level::Library]
        );
        assert!(parse_levels(&["BOGUS"]).is_err());
    }

    #[test]
    fn buckets_follow_header_order_with_unknown_last() {
        let h = header(HEADER);
        let levels = Levels::new(
            &h,
            &[
                Level::ReadGroup,
                Level::Library,
                Level::Sample,
                Level::AllReads,
            ],
            || 0u32,
        );
        let buckets = levels.into_buckets();
        let keys: Vec<(Level, &str)> = buckets.iter().map(|b| (b.level, b.key.as_str())).collect();
        assert_eq!(
            keys,
            [
                (Level::AllReads, "All_Reads"),
                (Level::Sample, "S2"),
                (Level::Sample, "S1"),
                (Level::Sample, "unknown"),
                (Level::Library, "LibB"),
                (Level::Library, "LibA"),
                (Level::Library, "unknown"),
                (Level::ReadGroup, "pu-b1"),
                (Level::ReadGroup, "pu-a1"),
                (Level::ReadGroup, "pu-c1"),
                (Level::ReadGroup, "unknown"),
            ]
        );
    }

    #[test]
    fn unknown_labels_match_picard() {
        let h = header(HEADER);
        let buckets = Levels::new(&h, &[Level::Library, Level::ReadGroup], || ()).into_buckets();
        let unknown_library = &buckets[2].labels;
        assert_eq!(unknown_library.sample.as_deref(), Some("unknown"));
        assert_eq!(unknown_library.library.as_deref(), Some("unknown"));
        // pu-c1 has a PU but no LB: sample kept, library empty.
        let no_library = &buckets[5].labels;
        assert_eq!(no_library.sample.as_deref(), Some("S1"));
        assert_eq!(no_library.library, None);
        assert_eq!(no_library.read_group.as_deref(), Some("pu-c1"));
        let unknown_unit = &buckets[6].labels;
        assert_eq!(unknown_unit.read_group.as_deref(), Some("unknown"));
        assert_eq!(unknown_unit.sample.as_deref(), Some("unknown"));
    }
}
