//! The plan of a chain as text: `chain ... --dry-run` and `plan` with a
//! chain's arguments. It prints the segments, which transform replays the
//! source or spools it, the files each step writes and the estimated
//! temporary disk. Nothing runs: the source's header is read to tell a BAM
//! from a file that is none, and no record.

use std::fmt::Write as _;

use bammetrics_core::mark_duplicates::{self, SortOrder};
use bammetrics_core::{cram, pass};

use super::VERB;
use super::prepare::{Chain, Planned, PrefixOrigin, Role};
use crate::error::VerbError;

/// Compressed bytes per record, for the record count of a file that is
/// not read: 10,339,495,604 bytes over 150,538,907 records (the HG002
/// exome BAM) and 15,797,182,294 bytes over 768.6 million records (the
/// NA12878 genome CRAM), both in docs/benchmarks.md.
const BAM_BYTES_PER_RECORD: u64 = 69;
const CRAM_BYTES_PER_RECORD: u64 = 21;
/// Bytes of SAM text per record: 2,289,237,081 bytes over 5,270,660
/// records (chromosomes 20 and 21 of that exome, as `view-sam` prints
/// them).
const SAM_BYTES_PER_RECORD: u64 = 434;
/// Peak bytes in `TMP_DIR` per record of a `MarkDuplicates` run that
/// spills: 46 on that BAM (6.5 GiB), 59 on that CRAM (42.1 GiB). The
/// larger one, so the estimate errs high.
const SPILL_BYTES_PER_RECORD: u64 = 59;

/// `1,158`: a byte count with thousands separators.
fn thousands(n: u64) -> String {
    let digits = n.to_string();
    let mut text = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, digit) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            text.push(',');
        }
        text.push(digit);
    }
    text
}

/// A size in the unit that keeps it under 1024: `6.5 GiB`.
fn human(bytes: u64) -> String {
    let units = ["bytes", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit + 1 < units.len() {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{bytes} bytes")
    } else {
        format!("{value:.1} {}", units[unit])
    }
}

/// What the plan knows about the source without reading a record.
struct Source {
    /// `BAM`, `CRAM` or `SAM`.
    kind: &'static str,
    /// The size of a regular file.
    bytes: Option<u64>,
    /// The header's sort order, for a regular file.
    order: Option<SortOrder>,
}

impl Source {
    fn of(chain: &Chain) -> Result<Self, VerbError> {
        if !chain.replayable() {
            return Ok(Self {
                kind: "BAM, CRAM or SAM",
                bytes: None,
                order: None,
            });
        }
        let header = pass::read_header(&chain.source_path)
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
    fn records(&self) -> Option<u64> {
        let per_record = match self.kind {
            "CRAM" => CRAM_BYTES_PER_RECORD,
            "SAM" => SAM_BYTES_PER_RECORD,
            _ => BAM_BYTES_PER_RECORD,
        };
        self.bytes.map(|bytes| bytes / per_record)
    }

    fn line(&self, chain: &Chain) -> String {
        match (self.bytes, self.order) {
            (Some(bytes), Some(order)) => format!(
                "{} ({}, {} bytes, sort order {}; a file, it can be read again)",
                chain.source,
                self.kind,
                thousands(bytes),
                order.name()
            ),
            _ => format!(
                "{} (a pipe: it is read once)",
                if chain.source == "-" {
                    "standard input"
                } else {
                    &chain.source
                }
            ),
        }
    }
}

/// `MAX_RECORDS_IN_RAM` of a transform step: its own value, or the
/// tool's default.
fn records_in_ram(step: &Planned) -> u64 {
    step.argv
        .iter()
        .rev()
        .find_map(|arg| arg.strip_prefix("MAX_RECORDS_IN_RAM="))
        .and_then(|value| value.parse().ok())
        .unwrap_or(mark_duplicates::DEFAULT_RECORDS_IN_RAM as u64)
}

/// The files of a step, as the plan lists them.
fn files(step: &Planned) -> String {
    // `write -` plans /dev/stdout, which reads as standard output below.
    let mut names: Vec<String> = step
        .files
        .iter()
        .filter(|_| !(step.stdout && step.role == Role::Write))
        .map(|file| file.path.display().to_string())
        .collect();
    if let Some(prefix) = &step.late_prefix {
        names.push(format!(
            "{}.* (named when the step finishes)",
            prefix.display()
        ));
    }
    if step.stdout {
        names.push("standard output".to_owned());
    }
    if names.is_empty() {
        "no file".to_owned()
    } else {
        names.join(", ")
    }
}

/// What a step does, in the column after its verb.
fn role(step: &Planned, first_read: bool) -> String {
    match (step.role, &step.leaves) {
        (Role::Transform, Some(_)) if first_read => "transform, first read".to_owned(),
        (Role::Transform, Some(stream)) => format!("stream {stream}"),
        (Role::Write, _) => format!("BAM of {}", step.reads),
        _ => format!("tap on {}", step.reads),
    }
}

/// One step line: `  step 3  insert-size  tap on marked  FILES`.
fn step_line(out: &mut String, widths: (usize, usize), step: &Planned, first_read: bool) {
    let what = if first_read {
        "takes every record, then finds the duplicate sets".to_owned()
    } else {
        files(step)
    };
    let _ = writeln!(
        out,
        "  step {:<2} {:<verb$}  {:<role$}  {what}",
        step.number,
        step.verb,
        role(step, first_read),
        verb = widths.0,
        role = widths.1,
    );
}

/// The header of segment `number`: what it reads.
fn segment_line(chain: &Chain, number: usize, taps_only: bool) -> String {
    let source = if chain.source == "-" {
        "standard input"
    } else {
        &chain.source
    };
    let transform = chain
        .steps
        .iter()
        .filter(|step| step.role == Role::Transform)
        .nth(number.saturating_sub(2));
    match (number, transform) {
        (1, _) if taps_only && !chain.stage_executor => {
            format!("segment 1  one read of {source}, shared by its taps as qc shares it")
        }
        (1, _) if chain.replayable() || chain.segments() == 1 => {
            format!("segment 1  one read of {source}")
        }
        (1, _) => format!(
            "segment 1  one read of {source}, copied to the temporary directory for the \
             second read (spool)"
        ),
        (_, Some(step)) if chain.replayable() => format!(
            "segment {number}  a second read of {source} through step {} (replay: nothing is \
             spooled)",
            step.number
        ),
        (_, Some(step)) => format!(
            "segment {number}  a read of the spooled copy through step {}",
            step.number
        ),
        (_, None) => format!("segment {number}"),
    }
}

/// The temporary disk a run needs, estimated.
fn temp_line(chain: &Chain, source: &Source) -> String {
    let Some(transform) = chain.steps.iter().find(|s| s.role == Role::Transform) else {
        return "temp disk  none: no step spills".to_owned();
    };
    let limit = records_in_ram(transform);
    let place = match &chain.tmp_dir {
        Some(dir) => dir.display().to_string(),
        None => "the system temporary directory".to_owned(),
    };
    let Some(records) = source.records() else {
        return format!(
            "temp disk  in {place}: the whole input once (the spool; a pipe has no size to \
             read), plus about {SPILL_BYTES_PER_RECORD} bytes per record when {} holds more \
             than {} records (MAX_RECORDS_IN_RAM)",
            transform.verb,
            thousands(limit)
        );
    };
    if records <= limit {
        return format!(
            "temp disk  none expected: about {} records, and {} keeps {} in memory before \
             it spills (MAX_RECORDS_IN_RAM)",
            thousands(records),
            transform.verb,
            thousands(limit)
        );
    }
    format!(
        "temp disk  about {} in {place}: {} spills its read ends, about \
         {SPILL_BYTES_PER_RECORD} bytes for each of about {} records (estimated from the \
         file size)",
        human(records * SPILL_BYTES_PER_RECORD),
        transform.verb,
        thousands(records)
    )
}

/// The plan of `chain` as text.
pub fn text(chain: &Chain) -> Result<String, VerbError> {
    let source = Source::of(chain)?;
    let mut out = String::new();
    let named = if chain.source == "-" {
        "standard input"
    } else {
        &chain.source
    };
    let _ = writeln!(out, "chain plan for {named}");
    let _ = writeln!(out, "  source    {}", source.line(chain));
    let origin = match chain.prefix_origin {
        PrefixOrigin::Option => "from -o",
        PrefixOrigin::Directory => "from -o, a directory, and the source name",
        PrefixOrigin::SourceName => "from the source name, no -o",
    };
    if chain.prefix_used {
        let _ = writeln!(out, "  prefix    {} ({origin})", chain.prefix.display());
    } else {
        let _ = writeln!(out, "  prefix    not used: every step names its files");
    }
    if let Some(reference) = &chain.reference {
        let _ = writeln!(out, "  reference {}", reference.display());
    }
    let taps_only = chain.steps.iter().all(|step| step.role == Role::Tap);
    let widths = (
        chain.steps.iter().map(|s| s.verb.len()).max().unwrap_or(0),
        chain
            .steps
            .iter()
            .map(|s| role(s, false).len().max(role(s, true).len()))
            .max()
            .unwrap_or(0),
    );
    for segment in 1..=chain.segments() {
        let _ = writeln!(out, "{}", segment_line(chain, segment, taps_only));
        for step in &chain.steps {
            let emits = step.role == Role::Transform && step.segment + 1 == segment;
            if step.segment == segment {
                step_line(&mut out, widths, step, step.role == Role::Transform);
            } else if emits {
                step_line(&mut out, widths, step, false);
            }
        }
    }
    if let Some(report) = &chain.report {
        let _ = writeln!(out, "also       {} (every metrics table)", report.display());
    }
    let _ = writeln!(out, "{}", temp_line(chain, &source));
    out.push_str("nothing ran: this is a plan\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_read_in_the_unit_that_fits() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(1158), "1,158");
        assert_eq!(thousands(10_339_495_604), "10,339,495,604");
        assert_eq!(human(512), "512 bytes");
        assert_eq!(human(6_962_982_912), "6.5 GiB");
    }

    /// The estimate for the two inputs it was measured on: within a third
    /// of the 6.5 GiB and 42.1 GiB peaks of docs/benchmarks.md.
    #[test]
    fn the_spill_estimate_is_near_the_measured_peaks() {
        let estimate = |bytes: u64, per_record: u64| bytes / per_record * SPILL_BYTES_PER_RECORD;
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
}
