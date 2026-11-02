//! `METRIC_ACCUMULATION_LEVEL`: one accumulator per sample, library or read
//! group, the way Picard's `MultiLevelCollector` keeps them.
//!
//! A tool keeps its per-unit state in [`Levels<T>`] and calls
//! [`Levels::route`] for each record; the closure sees every bucket the
//! record belongs to. [`Levels::into_buckets`] returns the buckets in the
//! order Picard writes their rows.
//!
//! Picard 3.5.0 behaviour (`MultiLevelCollector.Distributor`, checked on a
//! 12-read-group BAM with `CollectInsertSizeMetrics` and
//! `CollectAlignmentSummaryMetrics`):
//! * Row order is `ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`, whatever
//!   order the options came in. Within a level, buckets follow the order their
//!   key first appears among the `@RG` header lines, not alphabetical order.
//! * `READ_GROUP` keys on the platform unit (`PU`), not on the `ID`.
//! * The buckets are made from the header before any record is read. The
//!   first `@RG` line without the level's field (`SM`, `LB` or `PU`) makes a
//!   bucket too, at its header position, labelled with that line's fields.
//!   No record ever reaches it. A tool that writes a row for an empty bucket
//!   (`CollectAlignmentSummaryMetrics` writes an all-zero `UNPAIRED` row)
//!   shows it; `CollectInsertSizeMetrics` writes nothing for it.
//! * A record with no `RG` tag, one the header does not list, or one whose
//!   `@RG` line lacks the level's field goes to the `unknown` bucket. That
//!   bucket is made when the first such record arrives, after every header
//!   bucket, and does not exist when no record needs it. Its labels are
//!   `unknown` up to the level: `SAMPLE` prints `unknown`; `LIBRARY` prints
//!   `unknown` twice; `READ_GROUP` prints `unknown` three times.
//! * A header bucket keeps the `@RG` line's fields up to its level: a missing
//!   `LB` with a `PU` prints the sample and an empty library.
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

impl Labels {
    /// The histogram column prefix of a collector: the read group, else the
    /// library, else the sample, else `All_Reads`.
    fn key(&self) -> String {
        self.read_group
            .as_deref()
            .or(self.library.as_deref())
            .or(self.sample.as_deref())
            .unwrap_or("All_Reads")
            .to_owned()
    }
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

/// Where a read group's records go at one level.
#[derive(Clone, Copy, Debug)]
enum Route {
    /// The header bucket of the read group's key.
    Bucket(usize),
    /// The read group lacks the level's field: the `unknown` bucket.
    Unknown,
}

/// The buckets of one `SAMPLE`, `LIBRARY` or `READ_GROUP` level (Picard's
/// `Distributor`).
struct Unit<T> {
    level: Level,
    /// Bucket keys in creation order; `None` for read groups without the
    /// level's field.
    keys: Vec<Option<String>>,
    buckets: Vec<Bucket<T>>,
    /// Index of the bucket keyed `unknown`, once it exists.
    unknown: Option<usize>,
}

impl<T: Clone> Unit<T> {
    /// One bucket per distinct key among the `@RG` lines, in header order,
    /// the `None` key included.
    fn new(level: Level, groups: &[ReadGroupInfo], empty: &T) -> Self {
        let mut unit = Self {
            level,
            keys: Vec::new(),
            buckets: Vec::new(),
            unknown: None,
        };
        for group in groups {
            let key = unit_key(level, group).map(str::to_owned);
            if !unit.keys.contains(&key) {
                if key.as_deref() == Some(UNKNOWN) {
                    unit.unknown = Some(unit.buckets.len());
                }
                unit.push(key, labels_for(level, group), empty.clone());
            }
        }
        unit
    }

    fn push(&mut self, key: Option<String>, labels: Labels, state: T) {
        self.keys.push(key);
        self.buckets.push(Bucket {
            level: self.level,
            key: labels.key(),
            labels,
            state,
        });
    }

    fn route(&self, group: &ReadGroupInfo) -> Route {
        match unit_key(self.level, group) {
            None | Some(UNKNOWN) => Route::Unknown,
            Some(key) => Route::Bucket(
                self.keys
                    .iter()
                    .position(|k| k.as_deref() == Some(key))
                    .unwrap_or_default(),
            ),
        }
    }

    /// The `unknown` bucket's state, made on first use after every header
    /// bucket.
    fn unknown_state(&mut self, empty: &T) -> &mut T {
        let index = match self.unknown {
            Some(i) => i,
            None => {
                self.push(
                    Some(UNKNOWN.to_owned()),
                    unknown_labels(self.level),
                    empty.clone(),
                );
                let i = self.buckets.len() - 1;
                self.unknown = Some(i);
                i
            }
        };
        &mut self.buckets[index].state
    }
}

/// The accumulators of every requested level.
pub struct Levels<T> {
    all: Option<Bucket<T>>,
    units: Vec<Unit<T>>,
    /// Read group `ID` to its route at each entry of `units`.
    by_read_group: HashMap<Vec<u8>, Vec<Route>>,
    /// The empty accumulator, copied into each `unknown` bucket.
    empty: T,
}

impl<T: Clone> Levels<T> {
    /// One bucket per unit of each level in `levels`, each a copy of `empty`.
    pub fn new(header: &sam::Header, levels: &[Level], empty: T) -> Self {
        let groups = read_groups(header);
        let all = levels.contains(&Level::AllReads).then(|| Bucket {
            level: Level::AllReads,
            key: "All_Reads".to_owned(),
            labels: Labels::default(),
            state: empty.clone(),
        });
        let units: Vec<Unit<T>> = [Level::Sample, Level::Library, Level::ReadGroup]
            .into_iter()
            .filter(|level| levels.contains(level))
            .map(|level| Unit::new(level, &groups, &empty))
            .collect();
        let by_read_group = groups
            .iter()
            .map(|g| {
                let routes = units.iter().map(|u| u.route(g)).collect();
                (g.id.clone().into_bytes(), routes)
            })
            .collect();
        Self {
            all,
            units,
            by_read_group,
            empty,
        }
    }

    /// Run `f` on every bucket this record belongs to: the `ALL_READS`
    /// bucket, and the sample, library and read group buckets of its `RG` tag.
    pub fn route(&mut self, record: &bam::Record, mut f: impl FnMut(&mut T)) {
        if let Some(all) = &mut self.all {
            f(&mut all.state);
        }
        if self.units.is_empty() {
            return;
        }
        let routes = match record.data().get(&Tag::READ_GROUP) {
            Some(Ok(Value::String(id))) => self.by_read_group.get(AsRef::<[u8]>::as_ref(id)),
            _ => None,
        };
        for (i, unit) in self.units.iter_mut().enumerate() {
            match routes.map_or(Route::Unknown, |r| r[i]) {
                Route::Bucket(b) => f(&mut unit.buckets[b].state),
                Route::Unknown => f(unit.unknown_state(&self.empty)),
            }
        }
    }

    /// Buckets in the order Picard writes their rows.
    pub fn into_buckets(self) -> Vec<Bucket<T>> {
        self.all
            .into_iter()
            .chain(self.units.into_iter().flat_map(|u| u.buckets))
            .collect()
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

/// Labels of a header bucket: the group's fields up to the level.
fn labels_for(level: Level, group: &ReadGroupInfo) -> Labels {
    Labels {
        sample: group.sample.clone(),
        library: group.library.clone().filter(|_| level != Level::Sample),
        read_group: group
            .platform_unit
            .clone()
            .filter(|_| level == Level::ReadGroup),
    }
}

/// Labels of an `unknown` bucket: `unknown` up to the level.
fn unknown_labels(level: Level) -> Labels {
    let text = Some(UNKNOWN.to_owned());
    Labels {
        sample: text.clone(),
        library: text.clone().filter(|_| level != Level::Sample),
        read_group: text.filter(|_| level == Level::ReadGroup),
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

    /// An unmapped record with read group `id`, or none, round-tripped
    /// through BAM bytes.
    fn record(id: Option<&str>) -> bam::Record {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::RecordBuf;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;

        let mut buf = RecordBuf::default();
        if let Some(id) = id {
            buf.data_mut()
                .insert(Tag::READ_GROUP, BufValue::from(id.to_owned()));
        }
        let header = sam::Header::default();
        let mut writer = bam::io::Writer::new(Vec::new());
        writer.write_header(&header).expect("header");
        writer
            .write_alignment_record(&header, &buf)
            .expect("encode");
        let bytes = writer.into_inner().finish().expect("finish");
        let mut reader = bam::io::Reader::new(bytes.as_slice());
        reader.read_header().expect("header");
        let mut out = bam::Record::default();
        reader.read_record(&mut out).expect("decode");
        out
    }

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

    fn keys(buckets: &[Bucket<u32>]) -> Vec<(Level, &str, u32)> {
        buckets
            .iter()
            .map(|b| (b.level, b.key.as_str(), b.state))
            .collect()
    }

    #[test]
    fn header_buckets_in_header_order_and_unknown_only_when_used() {
        let all = [
            Level::ReadGroup,
            Level::Library,
            Level::Sample,
            Level::AllReads,
        ];
        let h = header(HEADER);
        let mut levels = Levels::new(&h, &all, 0u32);
        levels.route(&record(Some("a1")), |n| *n += 1);
        let buckets = levels.into_buckets();
        // c1 has no LB: its LIBRARY bucket takes c1's sample as its key and
        // gets nothing. a2 has no PU: likewise at READ_GROUP.
        assert_eq!(
            keys(&buckets),
            [
                (Level::AllReads, "All_Reads", 1),
                (Level::Sample, "S2", 0),
                (Level::Sample, "S1", 1),
                (Level::Library, "LibB", 0),
                (Level::Library, "LibA", 1),
                (Level::Library, "S1", 0),
                (Level::ReadGroup, "pu-b1", 0),
                (Level::ReadGroup, "pu-a1", 1),
                (Level::ReadGroup, "LibA", 0),
                (Level::ReadGroup, "pu-c1", 0),
            ]
        );
    }

    #[test]
    fn records_without_the_field_make_an_unknown_bucket_last() {
        let h = header(HEADER);
        let mut levels = Levels::new(&h, &[Level::Library, Level::ReadGroup], 0u32);
        for id in [Some("c1"), None, Some("zz"), Some("a2")] {
            levels.route(&record(id), |n| *n += 1);
        }
        let buckets = levels.into_buckets();
        let unknown_library = &buckets[3];
        assert_eq!(
            (unknown_library.key.as_str(), unknown_library.state),
            ("unknown", 3)
        );
        assert_eq!(unknown_library.labels.sample.as_deref(), Some("unknown"));
        assert_eq!(unknown_library.labels.library.as_deref(), Some("unknown"));
        assert_eq!(unknown_library.labels.read_group, None);
        // pu-c1 has a PU but no LB: sample kept, library empty.
        let no_library = &buckets[7].labels;
        assert_eq!(no_library.sample.as_deref(), Some("S1"));
        assert_eq!(no_library.library, None);
        assert_eq!(no_library.read_group.as_deref(), Some("pu-c1"));
        let unknown_unit = &buckets[8];
        assert_eq!(unknown_unit.state, 3);
        assert_eq!(unknown_unit.labels.read_group.as_deref(), Some("unknown"));
        assert_eq!(unknown_unit.labels.library.as_deref(), Some("unknown"));
        // The empty a2 bucket keeps a2's sample and library.
        let empty_unit = &buckets[6];
        assert_eq!(empty_unit.state, 0);
        assert_eq!(empty_unit.labels.library.as_deref(), Some("LibA"));
        assert_eq!(empty_unit.labels.read_group, None);
    }
}
