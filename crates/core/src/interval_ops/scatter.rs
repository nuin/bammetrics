//! Picard 3.5.0's `IntervalListScatter` and its five `IntervalListScatterer`s,
//! and the directory `IntervalListTools.writeScatterIntervals` writes.
//!
//! The list is first `sorted()` (`uniqued()` for `INTERVAL_SUBDIVISION`).
//! Each list is filled from a queue of the intervals: `takeSome` decides
//! whether the next interval joins the list, is split (subdivision only),
//! or ends it; the last list asked for takes every interval left. A list
//! may come out empty, and fewer lists than asked come out when the
//! intervals run out first.

use std::collections::VecDeque;
use std::path::{Path, PathBuf};

use super::header::Header;
use super::write::{check_target, stage};
use super::{Interval, IntervalList};
use crate::cleanup::{AtomicFile, MadeDirs};
use crate::error::{Error, Result};

/// The lists a scatter makes at most, so the files it writes at most: a
/// scatter that would make more is refused, whatever `SCATTER_COUNT` asks.
pub const MAX_SCATTER_FILES: i64 = 100_000;

/// The bytes a scatter writes at most, counted as text before any is
/// written: every list carries the whole header, so the GRCh38 dictionary
/// (723 KB) in 1,000 lists is 723 MB.
pub const MAX_SCATTER_BYTES: u64 = 32 << 30;

/// `IntervalListScatterMode`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScatterMode {
    IntervalSubdivision,
    BalancingWithoutIntervalSubdivision,
    BalancingWithoutIntervalSubdivisionWithOverflow,
    IntervalCount,
    IntervalCountWithDistributedRemainder,
}

impl ScatterMode {
    pub const ALL: [ScatterMode; 5] = [
        Self::IntervalSubdivision,
        Self::BalancingWithoutIntervalSubdivision,
        Self::BalancingWithoutIntervalSubdivisionWithOverflow,
        Self::IntervalCount,
        Self::IntervalCountWithDistributedRemainder,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Self::IntervalSubdivision => "INTERVAL_SUBDIVISION",
            Self::BalancingWithoutIntervalSubdivision => "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
            Self::BalancingWithoutIntervalSubdivisionWithOverflow => {
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
            }
            Self::IntervalCount => "INTERVAL_COUNT",
            Self::IntervalCountWithDistributedRemainder => {
                "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER"
            }
        }
    }

    /// The Picard value, case as written.
    pub fn parse(text: &str) -> Option<ScatterMode> {
        Self::ALL.into_iter().find(|mode| mode.name() == text)
    }

    fn by_count(self) -> bool {
        matches!(
            self,
            Self::IntervalCount | Self::IntervalCountWithDistributedRemainder
        )
    }

    /// `intervalWeight`: the length, or 1 by count.
    fn interval_weight(self, interval: &Interval) -> i64 {
        if self.by_count() {
            1
        } else {
            i64::from(interval.length())
        }
    }

    /// `deduceIdealSplitWeight`, an `int`.
    fn ideal_split_weight(self, list: &IntervalList, count: i32) -> i32 {
        if self.by_count() {
            return (list.len() as i64).div_euclid(i64::from(count)).max(1) as i32;
        }
        let base = 1.max(list.unique_base_count().div_euclid(i64::from(count)) as i32);
        match self {
            Self::BalancingWithoutIntervalSubdivision
            | Self::BalancingWithoutIntervalSubdivisionWithOverflow => {
                let longest = list
                    .intervals
                    .iter()
                    .map(Interval::length)
                    .max()
                    .unwrap_or(1);
                longest.max(base)
            }
            _ => base,
        }
    }

    /// `takeSome`: the part of `interval` the list takes and the part left
    /// in the queue.
    fn take_some(
        self,
        interval: Interval,
        ideal: i64,
        current: i64,
        projected: f64,
    ) -> (Option<Interval>, Option<Interval>) {
        match self {
            Self::IntervalSubdivision => {
                let amount = ideal - current;
                if amount >= i64::from(interval.length()) {
                    (Some(interval), None)
                } else if amount == 0 {
                    (None, Some(interval))
                } else {
                    let cut = interval.start.wrapping_add(amount as i32);
                    let left = Interval {
                        end: cut.wrapping_sub(1),
                        ..interval.clone()
                    };
                    let right = Interval {
                        start: cut,
                        ..interval
                    };
                    (Some(left), Some(right))
                }
            }
            Self::BalancingWithoutIntervalSubdivision
            | Self::BalancingWithoutIntervalSubdivisionWithOverflow => {
                let size = current + self.interval_weight(&interval);
                let include = size <= ideal
                    || (self == Self::BalancingWithoutIntervalSubdivisionWithOverflow
                        && (ideal as f64) < projected);
                if include {
                    (Some(interval), None)
                } else {
                    (None, Some(interval))
                }
            }
            Self::IntervalCount => {
                if ideal - current > 0 {
                    (Some(interval), None)
                } else {
                    (None, Some(interval))
                }
            }
            Self::IntervalCountWithDistributedRemainder => {
                if projected > current as f64 {
                    (Some(interval), None)
                } else {
                    (None, Some(interval))
                }
            }
        }
    }
}

/// `listWeight`: the base count, or the number of intervals by count.
pub fn list_weight(mode: ScatterMode, list: &IntervalList) -> i64 {
    if mode.by_count() {
        list.len() as i64
    } else {
        list.base_count()
    }
}

/// `Math.round(double)` of Java 8 and later.
fn java_round(value: f64) -> i64 {
    let bits = value.to_bits() as i64;
    let biased_exponent = (bits & 0x7ff0_0000_0000_0000) >> 52;
    let shift = (52 - 1 + 1023) - biased_exponent;
    if shift & -64 == 0 {
        let mut r = (bits & 0x000f_ffff_ffff_ffff) | 0x0010_0000_0000_0000;
        if bits < 0 {
            r = -r;
        }
        ((r >> shift) + 1) >> 1
    } else {
        // NaN gives 0 and an infinity or a large value the nearest bound,
        // as Java's `(long)` cast.
        value as i64
    }
}

/// `IntervalListTools`'s `SCATTER_COUNT` from `SCATTER_CONTENT`:
/// `(int) Math.round(listWeight / content)`. A content of 0 gives -1 (the
/// list is then written whole); a count above [`MAX_SCATTER_FILES`] is
/// refused by [`scatter`].
pub fn scatter_count_for_content(mode: ScatterMode, list: &IntervalList, content: i32) -> i32 {
    java_round(list_weight(mode, list) as f64 / f64::from(content)) as i32
}

/// The lists of a scatter: the header every list carries, held once, and
/// each list's intervals. A header of 3,000 `@SQ` lines cloned into
/// [`MAX_SCATTER_FILES`] lists would hold gigabytes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Scattered {
    pub header: Header,
    pub parts: Vec<Vec<Interval>>,
}

impl Scattered {
    /// The number of lists.
    pub fn len(&self) -> usize {
        self.parts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }

    /// The summed base counts of the lists (`ScatterSummary.baseCount`).
    pub fn base_count(&self) -> i64 {
        self.parts
            .iter()
            .flatten()
            .map(|i| i64::from(i.length()))
            .sum()
    }

    /// The summed interval counts of the lists.
    pub fn interval_count(&self) -> i64 {
        self.parts.iter().map(|p| p.len() as i64).sum()
    }
}

/// `IntervalListScatter(scatterer, list, count)` iterated to its end. A
/// scatter that makes more than [`MAX_SCATTER_FILES`] lists is refused
/// when the next list starts.
pub fn scatter(list: &IntervalList, mode: ScatterMode, count: i32) -> Result<Scattered> {
    if count < 1 {
        return Err(Error::Usage(format!(
            "Scatter count should be a positive integer, but it was{count}."
        )));
    }
    let list = match mode {
        ScatterMode::IntervalSubdivision => list.uniqued(true),
        _ => list.sorted(),
    };
    let ideal = i64::from(mode.ideal_split_weight(&list, count));
    let mut weight_remaining = list_weight(mode, &list);
    let IntervalList { header, intervals } = list;
    let mut queue = VecDeque::from(intervals);
    let mut out = Vec::new();
    let mut returned: i64 = 0;
    while !queue.is_empty() {
        if out.len() as i64 == MAX_SCATTER_FILES {
            return Err(Error::Usage(format!(
                "SCATTER_COUNT={count} makes more than {MAX_SCATTER_FILES} lists; at most \
                 {MAX_SCATTER_FILES} files are written"
            )));
        }
        returned += 1;
        let mut part = Vec::new();
        let mut current: i64 = 0;
        let mut ended = false;
        while returned < i64::from(count) {
            let Some(interval) = queue.pop_front() else {
                break;
            };
            let projected =
                (weight_remaining - current) as f64 / (i64::from(count) - returned) as f64;
            let (take, rest) = mode.take_some(interval, ideal, current, projected);
            if let Some(rest) = rest {
                queue.push_front(rest);
            }
            let Some(take) = take else {
                weight_remaining -= current;
                ended = true;
                break;
            };
            current += mode.interval_weight(&take);
            part.push(take);
        }
        if !ended {
            part.extend(queue.drain(..));
        }
        out.push(part);
    }
    Ok(Scattered { header, parts: out })
}

/// `createSubDirectoryAndGetScatterFile`: `dir/temp_NNNN_of_COUNT/scattered.interval_list`,
/// `index` counted from 1 and written with at least four digits.
pub fn scatter_path(dir: &Path, index: usize, count: i32) -> PathBuf {
    dir.join(format!("temp_{index:04}_of_{count}"))
        .join("scattered.interval_list")
}

/// Write the lists of [`scatter`] under `dir` at [`scatter_path`], every
/// file under a temporary name first and all renamed into place as one
/// set; a failure leaves none of them and removes the directories it made.
/// No file may be one of `inputs`.
pub fn write_scatter(
    dir: &Path,
    lists: &Scattered,
    count: i32,
    inputs: &[PathBuf],
) -> Result<Vec<PathBuf>> {
    if lists.len() as i64 > MAX_SCATTER_FILES {
        return Err(Error::Usage(format!(
            "{} lists to write; at most {MAX_SCATTER_FILES} are written",
            lists.len()
        )));
    }
    let bytes = scatter_bytes(lists);
    if bytes > MAX_SCATTER_BYTES {
        return Err(Error::Usage(format!(
            "the {} lists would write {bytes} bytes (the header in each); at most \
             {MAX_SCATTER_BYTES} are written",
            lists.len()
        )));
    }
    let paths: Vec<PathBuf> = (1..=lists.len())
        .map(|index| scatter_path(dir, index, count))
        .collect();
    for path in &paths {
        check_target(path, inputs)?;
    }
    let mut dirs = MadeDirs::default();
    let mut staged = Vec::with_capacity(lists.len());
    for (part, path) in lists.parts.iter().zip(&paths) {
        dirs.create_parent(path).map_err(|source| Error::Output {
            path: path.clone(),
            source,
        })?;
        let file = AtomicFile::new(path);
        stage(&lists.header, part, &file)?;
        staged.push(file);
    }
    let written = crate::split_sam::commit_all(staged)?;
    dirs.keep();
    Ok(written)
}

/// The bytes of text the lists of `lists` hold: the header once per list
/// and one line per interval.
pub fn scatter_bytes(lists: &Scattered) -> u64 {
    let header = lists.header.encode().len() as u64;
    let lines: u64 = lists
        .parts
        .iter()
        .flatten()
        .map(|i| {
            let fields = lists.header.sequences[i.contig].name.len()
                + decimal_length(i.start)
                + decimal_length(i.end)
                + i.name.as_ref().map_or(1, String::len);
            // Four tabs, the strand and the line feed.
            fields as u64 + 6
        })
        .sum();
    header * lists.len() as u64 + lines
}

/// The characters of `value` in decimal.
fn decimal_length(value: i32) -> usize {
    let digits = value
        .unsigned_abs()
        .checked_ilog10()
        .map_or(1, |d| d as usize + 1);
    digits + usize::from(value < 0)
}
