//! The planner of docs/07 section 2.4: what each barrier of a chain holds
//! in memory and writes to the temporary directory, estimated before the
//! first read, and the free-space refusal.
//!
//! Each barrier reads its input once and then hands the stream on:
//! * `sort-sam` holds its records ([`Next::Records`]). It keeps at most
//!   `--sort-memory` bytes and spills the rest as sorted runs.
//! * `mark-duplicates` reads its input twice. It replays the source when it
//!   reads the source itself. Behind another transform it writes a spool of
//!   its input and replays that. A pipe source in front of it is copied to
//!   the temporary directory first.
//!
//! The record count comes from the source's size. Every constant below was
//! measured (docs/chaining.md section 11); none is a guess.
//!
//! [`Next::Records`]: bammetrics_core::stage::Next::Records

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use bammetrics_core::mark_duplicates::{self, SortOrder};
use bammetrics_core::sam_sort::ITEM_BYTES;
use bammetrics_core::sam_text::Stringency;
use bammetrics_core::stage::WriteFormat;
use bammetrics_core::{cram, pass};

use super::prepare::{Chain, Planned, Role};
use super::{Spec, VERB, usage};
use crate::error::VerbError;

/// Compressed bytes per record, for the record count of a file that is
/// not read: 10,339,495,604 bytes over 150,538,907 records (the HG002
/// exome BAM) and 15,797,182,294 bytes over 768.6 million records (the
/// NA12878 genome CRAM), both in docs/benchmarks.md.
const BAM_BYTES_PER_RECORD: u64 = 69;
const CRAM_BYTES_PER_RECORD: u64 = 21;
/// Bytes of SAM text per record: 2,289,237,081 bytes over 5,270,660
/// records (chromosomes 20 and 21 of that exome, as `view-sam` prints
/// them). Also what `write x.sam` writes per record.
const SAM_BYTES_PER_RECORD: u64 = 434;
/// Bytes of one record in memory: 2,951,869,900 inflated bytes over
/// 8,445,233 records of chromosomes 20 to 22 of that exome (349.5).
const RECORD_BYTES: u64 = 350;
/// Resident memory per byte a sort counts against `--sort-memory`, in
/// quarters: 4,390 MB resident for the 3,526 MB a sort of those 8,445,233
/// records held in memory (1.25).
const RESIDENT_QUARTERS: u64 = 5;
/// Resident memory of a chain beside its barriers: 425 MB resident at
/// `--sort-memory 256M`, less 1.25 times the 268 MB the sort held.
const BASE_MEMORY: u64 = 100 << 20;
/// Bytes per record a sorted run or a spool takes at BGZF level 3, by the
/// order of the records in it, over those 8,445,233 records: queryname
/// runs 871,834,880 and 877,285,505, a queryname spool 879,914,157 (104.2);
/// coordinate runs 616,683,279 from coordinate input and 639,427,369 from
/// queryname input (75.7). The larger of each.
const SPILL_BYTES_QUERYNAME: u64 = 105;
const SPILL_BYTES_COORDINATE: u64 = 76;
/// BAM bytes per record `write` makes at level 5: 856,813,333 (queryname)
/// and 584,011,111 (coordinate) over those records.
const BAM_OUT_QUERYNAME: u64 = 102;
const BAM_OUT_COORDINATE: u64 = 70;
/// Peak bytes in `TMP_DIR` per record of a `MarkDuplicates` run that
/// spills: 46 on that BAM (6.5 GiB), 59 on that CRAM (42.1 GiB). The
/// larger one, so the estimate errs high.
const MARKDUP_SPILL_PER_RECORD: u64 = 59;
/// Resident bytes of `MarkDuplicates` per record it holds (at most
/// `MAX_RECORDS_IN_RAM`) beside the base: 646 MiB with 4,000,000 held on
/// the 8,445,233 records, 800 MiB on the 150,538,907 of the exome.
const MARKDUP_MEMORY_PER_RECORD: u64 = 200;
/// Runs one merge reads at once (`spill.rs`): more cost a further level.
const MERGE_FAN_IN: u64 = 64;

/// `--sort-memory` when none is given: 2 GiB (docs/07 section 2.1).
pub const DEFAULT_SORT_MEMORY: u64 = 2 << 30;

/// A byte count as `--sort-memory` takes it: digits, then an optional K, M,
/// G or T (1024 based) with an optional `B` or `iB`, any case. `None` for
/// zero, an overflow, or any other text.
pub fn parse_bytes(text: &str) -> Option<u64> {
    let text = text.trim();
    let digits = text
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(text.len());
    let (number, unit) = text.split_at(digits);
    let number: u64 = number.parse().ok()?;
    let unit = unit.trim().to_ascii_uppercase();
    let unit = unit
        .strip_suffix("IB")
        .or_else(|| unit.strip_suffix('B'))
        .unwrap_or(&unit);
    let shift = match unit {
        "" => 0,
        "K" => 10,
        "M" => 20,
        "G" => 30,
        "T" => 40,
        _ => return None,
    };
    number.checked_mul(1 << shift).filter(|bytes| *bytes > 0)
}

/// `--sort-memory`, parsed.
pub(super) fn sort_memory(spec: &Spec) -> Result<u64, VerbError> {
    let Some(text) = &spec.sort_memory else {
        return Ok(DEFAULT_SORT_MEMORY);
    };
    parse_bytes(text).ok_or_else(|| {
        usage(format!("--sort-memory {text:?}: expected bytes above 0"))
            .from("--sort-memory")
            .hint("a number of bytes, or one with K, M, G or T: 512M, 2G, 2GiB")
    })
}

/// `--validation-stringency`, parsed.
pub(super) fn stringency(spec: &Spec) -> Result<Stringency, VerbError> {
    let Some(text) = &spec.stringency else {
        return Ok(Stringency::Strict);
    };
    Stringency::parse(text).ok_or_else(|| {
        usage(format!(
            "--validation-stringency {text}: expected STRICT, LENIENT or SILENT"
        ))
        .from("--validation-stringency")
    })
}

/// What the plan knows about the source without reading a record.
pub struct Source {
    /// `BAM`, `CRAM` or `SAM`.
    pub kind: &'static str,
    /// The size of a regular file.
    pub bytes: Option<u64>,
    /// The header's sort order, for a regular file.
    pub order: Option<SortOrder>,
}

impl Source {
    /// The source of `chain`: its header is read, no record. A pipe is
    /// not opened.
    pub fn of(chain: &Chain) -> Result<Self, VerbError> {
        if !chain.replayable() {
            return Ok(Self {
                kind: "BAM, CRAM or SAM",
                bytes: None,
                order: None,
            });
        }
        let header = pass::read_header_at(&chain.source_path, chain.stringency)
            .map_err(|e| VerbError::input(VERB, e.to_string()).from("BAM"))?;
        Ok(Self {
            kind: if cram::is_cram(&chain.source_path) {
                "CRAM"
            } else if cram::is_sam(&chain.source_path) {
                "SAM"
            } else {
                "BAM"
            },
            bytes: std::fs::metadata(&chain.source_path).ok().map(|m| m.len()),
            order: Some(SortOrder::of(&header)),
        })
    }

    /// The records of the source, estimated from its size.
    pub fn records(&self) -> Option<u64> {
        let per_record = match self.kind {
            "CRAM" => CRAM_BYTES_PER_RECORD,
            "SAM" => SAM_BYTES_PER_RECORD,
            _ => BAM_BYTES_PER_RECORD,
        };
        self.bytes.map(|bytes| bytes / per_record)
    }
}

/// How a barrier hands its stream on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum How {
    /// It holds the records, in memory or in sorted runs (`sort-sam`).
    Holds,
    /// It reads the source file a second time.
    Replays,
    /// It reads a copy of a piped source, made before the first read.
    SpoolsSource,
    /// It writes its input to a spool and reads that a second time.
    Spools,
}

/// One barrier, estimated.
pub struct Barrier {
    /// The step's number.
    pub number: usize,
    pub verb: &'static str,
    pub how: How,
    /// The order its records take (a sort), or the order of its input.
    pub order: Option<SortOrder>,
    /// Records it holds before it spills: `MAX_RECORDS_IN_RAM`, or `None`
    /// for a sort bounded by `--sort-memory` alone.
    pub records_in_ram: Option<u64>,
    /// Bytes of sorted runs or spilled read ends in the temporary
    /// directory; `None` when the record count is unknown.
    pub spill: Option<u64>,
    /// Bytes of its spool (or of the copy of a piped source).
    pub spool: Option<u64>,
    /// Resident memory while it reads its input, the base included.
    pub memory: Option<u64>,
}

/// The estimate of a chain, before the first read.
pub struct Estimate {
    pub source: Source,
    pub records: Option<u64>,
    pub barriers: Vec<Barrier>,
    /// Every file written by `write`, with its estimated bytes.
    pub outputs: Vec<(PathBuf, u64)>,
}

impl Estimate {
    /// Bytes in the temporary directory at most: every spill and spool,
    /// as if they were all held at once. `None` when one is unknown.
    pub fn temp(&self) -> Option<u64> {
        self.barriers
            .iter()
            .map(|b| Some(b.spill?.saturating_add(b.spool?)))
            .sum()
    }

    /// The largest resident memory over the barriers; the base without
    /// one. `None` when the record count is unknown.
    pub fn memory(&self) -> Option<u64> {
        self.barriers
            .iter()
            .map(|b| b.memory)
            .try_fold(BASE_MEMORY, |most, m| m.map(|m| most.max(m)))
    }
}

/// `NAME=VALUE` of a step's Picard arguments, the last one given.
fn picard_value<'a>(step: &'a Planned, name: &str) -> Option<&'a str> {
    step.argv.iter().rev().find_map(|arg| {
        arg.strip_prefix(name)
            .and_then(|rest| rest.strip_prefix('='))
    })
}

fn spill_bytes(order: Option<SortOrder>) -> u64 {
    match order {
        Some(SortOrder::Coordinate) => SPILL_BYTES_COORDINATE,
        _ => SPILL_BYTES_QUERYNAME,
    }
}

/// Resident bytes for `held` bytes a sort counts.
fn resident(held: u64) -> u64 {
    held.saturating_mul(RESIDENT_QUARTERS) / 4
}

/// The estimate of `chain`: the header of a file source is read, no
/// record.
pub fn estimate(chain: &Chain) -> Result<Estimate, VerbError> {
    let source = Source::of(chain)?;
    let records = source.records();
    let mut order = source.order;
    let mut first_stage = true;
    // The memory of the barrier before, beside the base. A barrier's peak
    // adds it: a sort that did not spill still holds its records while the
    // next barrier reads them, and the allocator keeps what the barrier
    // before freed. Measured: `sort-sam` (1 GiB) + `mark-duplicates` +
    // `sort-sam` peaked at 1,771,077,632 bytes on 8,445,233 records, more
    // than any one barrier alone.
    let mut held_before: Option<u64> = Some(0);
    let mut barriers = Vec::new();
    let mut outputs = Vec::new();
    for step in &chain.steps {
        match (step.role, step.tool.map(|tool| tool.name)) {
            (Role::Transform, Some("SortSam")) => {
                let sorted = picard_value(step, "SORT_ORDER").and_then(SortOrder::parse);
                let records_in_ram =
                    picard_value(step, "MAX_RECORDS_IN_RAM").and_then(|n| n.parse::<u64>().ok());
                let per_record = RECORD_BYTES + ITEM_BYTES as u64;
                let budget = records.map(|n| n.saturating_mul(per_record));
                let cap = records_in_ram
                    .map_or(chain.sort_memory, |n| n.saturating_mul(per_record))
                    .min(chain.sort_memory);
                let held = budget.map(|b| b.min(cap));
                let spills = budget.map(|b| b > cap);
                let spill = records.zip(spills).map(|(n, spills)| {
                    if !spills {
                        return 0;
                    }
                    let runs = budget.unwrap_or(0).div_ceil(cap.max(1));
                    // One write per level of merges before the last.
                    let mut levels = 1;
                    let mut left = runs;
                    while left > MERGE_FAN_IN {
                        left = left.div_ceil(MERGE_FAN_IN);
                        levels += 1;
                    }
                    n.saturating_mul(spill_bytes(sorted)).saturating_mul(levels)
                });
                let memory = held
                    .zip(held_before)
                    .map(|(held, before)| BASE_MEMORY + resident(held) + before);
                held_before = held.map(resident);
                barriers.push(Barrier {
                    number: step.number,
                    verb: step.verb,
                    how: How::Holds,
                    order: sorted,
                    records_in_ram,
                    spill,
                    spool: Some(0),
                    memory,
                });
                order = sorted;
                first_stage = false;
            }
            (Role::Transform, _) => {
                let how = match (first_stage, chain.replayable()) {
                    (true, true) => How::Replays,
                    (true, false) => How::SpoolsSource,
                    (false, _) => How::Spools,
                };
                let spool = match how {
                    How::Replays | How::Holds => Some(0),
                    How::SpoolsSource => source.bytes,
                    How::Spools => records.map(|n| n.saturating_mul(spill_bytes(order))),
                };
                let records_in_ram = picard_value(step, "MAX_RECORDS_IN_RAM")
                    .and_then(|n| n.parse::<u64>().ok())
                    .unwrap_or(mark_duplicates::DEFAULT_RECORDS_IN_RAM as u64);
                let spill = records.map(|n| {
                    if n > records_in_ram {
                        n.saturating_mul(MARKDUP_SPILL_PER_RECORD)
                    } else {
                        0
                    }
                });
                let own = records.map(|n| n.min(records_in_ram) * MARKDUP_MEMORY_PER_RECORD);
                let memory = own
                    .zip(held_before)
                    .map(|(own, before)| BASE_MEMORY + own + before);
                held_before = own;
                barriers.push(Barrier {
                    number: step.number,
                    verb: step.verb,
                    how,
                    order,
                    records_in_ram: Some(records_in_ram),
                    spill,
                    spool,
                    memory,
                });
                first_stage = false;
            }
            (Role::Write, _) if !step.stdout => {
                let Some(options) = &step.write else { continue };
                let per_record = match (options.format, order) {
                    (WriteFormat::Sam, _) => SAM_BYTES_PER_RECORD,
                    (WriteFormat::Bam, Some(SortOrder::Coordinate)) => BAM_OUT_COORDINATE,
                    (WriteFormat::Bam, _) => BAM_OUT_QUERYNAME,
                };
                if let Some(n) = records {
                    outputs.push((options.path.clone(), n.saturating_mul(per_record)));
                }
            }
            _ => {}
        }
    }
    Ok(Estimate {
        source,
        records,
        barriers,
        outputs,
    })
}

/// The free bytes of the file system that holds `dir`, for this user.
#[cfg(unix)]
fn free_bytes(dir: &Path) -> Option<u64> {
    use std::os::unix::ffi::OsStrExt;
    let path = std::ffi::CString::new(dir.as_os_str().as_bytes()).ok()?;
    // SAFETY: `path` is a NUL-terminated string that outlives the call, and
    // `stat` is a plain C struct the call fills.
    let mut stat: libc::statvfs = unsafe { std::mem::zeroed() };
    let status = unsafe { libc::statvfs(path.as_ptr(), &mut stat) };
    (status == 0).then(|| u64::from(stat.f_bavail).saturating_mul(stat.f_frsize))
}

#[cfg(not(unix))]
fn free_bytes(_dir: &Path) -> Option<u64> {
    None
}

/// The nearest directory at or above `path` that exists.
fn existing_dir(path: &Path) -> PathBuf {
    let mut dir = if path.as_os_str().is_empty() {
        PathBuf::from(".")
    } else {
        path.to_path_buf()
    };
    while !dir.is_dir() {
        match dir.parent() {
            Some(parent) if !parent.as_os_str().is_empty() => dir = parent.to_path_buf(),
            _ => return PathBuf::from("."),
        }
    }
    dir
}

/// What the temporary directory holds, in a [`Need`].
const SPILLS: &str = "spills and spools";

/// A file system the chain writes to: a directory on it, the bytes the
/// chain needs there, and what they are for.
pub struct Need {
    pub dir: PathBuf,
    pub bytes: u64,
    pub what: Vec<String>,
    pub free: Option<u64>,
}

/// Where the chain writes, by file system (device), with the bytes it
/// needs on each. The temporary directory takes every spill and spool,
/// each output directory its `write` files.
pub fn needs(chain: &Chain, estimate: &Estimate) -> Vec<Need> {
    let temp = chain.tmp_dir.clone().unwrap_or_else(std::env::temp_dir);
    let mut wanted: Vec<(PathBuf, u64, String)> = Vec::new();
    if let Some(bytes) = estimate.temp().filter(|bytes| *bytes > 0) {
        wanted.push((existing_dir(&temp), bytes, SPILLS.to_owned()));
    }
    for (path, bytes) in &estimate.outputs {
        let dir = existing_dir(path.parent().unwrap_or(Path::new(".")));
        wanted.push((dir, *bytes, path.display().to_string()));
    }
    let mut by_device: HashMap<u64, usize> = HashMap::new();
    let mut needs: Vec<Need> = Vec::new();
    for (dir, bytes, what) in wanted {
        let device = device_of(&dir);
        match by_device.get(&device) {
            Some(&index) => {
                needs[index].bytes = needs[index].bytes.saturating_add(bytes);
                needs[index].what.push(what);
            }
            None => {
                by_device.insert(device, needs.len());
                let free = chain.free_space.or_else(|| free_bytes(&dir));
                needs.push(Need {
                    dir,
                    bytes,
                    what: vec![what],
                    free,
                });
            }
        }
    }
    needs
}

#[cfg(unix)]
fn device_of(dir: &Path) -> u64 {
    use std::os::unix::fs::MetadataExt;
    std::fs::metadata(dir).map_or(0, |m| m.dev())
}

#[cfg(not(unix))]
fn device_of(_dir: &Path) -> u64 {
    0
}

/// The needs that exceed their free space.
pub fn short(needs: &[Need]) -> impl Iterator<Item = &Need> {
    needs
        .iter()
        .filter(|need| need.free.is_some_and(|free| need.bytes > free))
}

/// Refuse the chain, before the first read, when a file system it writes
/// to has less free space than the estimate needs there; `--force-space`
/// runs it anyway. Exit 1.
pub fn check_space(chain: &Chain) -> Result<(), VerbError> {
    if chain.force_space {
        return Ok(());
    }
    let estimate = estimate(chain)?;
    let needs = needs(chain, &estimate);
    let Some(need) = short(&needs).next() else {
        return Ok(());
    };
    let temp = need.what.iter().any(|what| what == SPILLS);
    let error = VerbError::input(
        VERB,
        format!(
            "not enough free space in {}: the chain needs about {} there ({}) and {} are \
             free",
            need.dir.display(),
            super::describe::human(need.bytes),
            need.what.join(", "),
            super::describe::human(need.free.unwrap_or(0))
        ),
    );
    let error = if temp {
        error.hint("--tmp-dir DIR puts spills and spools on another file system")
    } else {
        error
    };
    Err(error.hint("--force-space runs the chain anyway"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The `MarkDuplicates` spill estimate for the two inputs it was
    /// measured on: within a third of the 6.5 GiB and 42.1 GiB peaks of
    /// docs/benchmarks.md.
    #[test]
    fn the_markdup_spill_estimate_is_near_the_measured_peaks() {
        let estimate = |bytes: u64, per_record: u64| bytes / per_record * MARKDUP_SPILL_PER_RECORD;
        let exome = estimate(10_339_495_604, BAM_BYTES_PER_RECORD) as f64;
        let genome = estimate(15_797_182_294, CRAM_BYTES_PER_RECORD) as f64;
        let gib = 1024.0 * 1024.0 * 1024.0;
        assert!((exome / gib - 6.5).abs() / 6.5 < 0.34, "{}", exome / gib);
        assert!(
            (genome / gib - 42.1).abs() / 42.1 < 0.34,
            "{}",
            genome / gib
        );
    }

    /// The sort estimates against the runs they came from: 8,445,233
    /// records, 573,136,551 bytes of BAM.
    #[test]
    fn the_sort_estimates_are_the_measured_runs() {
        let records = 573_136_551 / BAM_BYTES_PER_RECORD;
        assert!((8_300_000..8_500_000).contains(&records));
        let queryname = records * SPILL_BYTES_QUERYNAME;
        assert!(queryname >= 871_834_880 * 98 / 100, "{queryname}");
        // At 256 MiB the resident peak was 425 MB.
        let resident = BASE_MEMORY + resident(256 << 20);
        assert!(
            resident.abs_diff(425_099_264) < 425_099_264 / 10,
            "{resident}"
        );
    }

    #[test]
    fn the_nearest_existing_directory_takes_a_missing_one() {
        let here = std::env::temp_dir();
        assert_eq!(existing_dir(&here.join("no/such/dir")), here);
        assert_eq!(existing_dir(Path::new("")), PathBuf::from("."));
        assert!(free_bytes(&here).is_some_and(|free| free > 0));
    }

    #[test]
    fn sizes_take_a_binary_suffix() {
        assert_eq!(parse_bytes("1000"), Some(1000));
        assert_eq!(parse_bytes("2G"), Some(2 << 30));
        assert_eq!(parse_bytes("2GiB"), Some(2 << 30));
        assert_eq!(parse_bytes("512m"), Some(512 << 20));
        assert_eq!(parse_bytes("64KB"), Some(64 << 10));
        assert_eq!(parse_bytes("1T"), Some(1 << 40));
        for bad in ["", "0", "G", "2X", "-1", "1.5G", "99999999999T"] {
            assert_eq!(parse_bytes(bad), None, "{bad}");
        }
    }
}
