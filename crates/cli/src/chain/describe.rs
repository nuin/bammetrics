//! The plan of a chain as text: `chain ... --dry-run` and `plan` with a
//! chain's arguments. It prints the segments and what each reads, the
//! files each step writes, then one line per barrier from the planner
//! ([`super::planner`]): whether it holds its records, replays the source
//! or spools its input, the bytes it spills and the memory it takes; then
//! the temporary disk, the output disk, the peak memory, and a line when a
//! file system has less free space than the chain needs. Nothing runs: the
//! source's header is read to tell a BAM from a file that is none, and no
//! record.

use std::fmt::Write as _;

use super::planner::{self, Barrier, Estimate, How, Source};
use super::prepare::{Chain, Planned, PrefixOrigin, Role};
use crate::error::VerbError;

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
pub(super) fn human(bytes: u64) -> String {
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

fn source_line(source: &Source, chain: &Chain) -> String {
    match (source.bytes, source.order) {
        (Some(bytes), Some(order)) => format!(
            "{} ({}, {} bytes, sort order {}; a file, it can be read again)",
            chain.source,
            source.kind,
            thousands(bytes),
            order.name()
        ),
        _ => format!("{} (a pipe: it is read once)", shown_source(chain)),
    }
}

fn shown_source(chain: &Chain) -> &str {
    if chain.source == "-" {
        "standard input"
    } else {
        &chain.source
    }
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
    match (step.role, &step.leaves, &step.write) {
        (Role::Transform, Some(_), _) if first_read => "transform, first read".to_owned(),
        (Role::Transform, Some(stream), _) => format!("stream {stream}"),
        (Role::Write, _, Some(options)) => {
            let kind = match options.format {
                bammetrics_core::stage::WriteFormat::Bam => "BAM",
                bammetrics_core::stage::WriteFormat::Sam => "SAM text",
            };
            format!("{kind} of {}", step.reads)
        }
        _ => format!("tap on {}", step.reads),
    }
}

/// What a transform does in its first read.
fn first_read(step: &Planned) -> String {
    let order = step
        .argv
        .iter()
        .rev()
        .find_map(|arg| arg.strip_prefix("SORT_ORDER="));
    match order {
        Some(order) => format!("takes every record and sorts them by {order}"),
        None => "takes every record, then finds the duplicate sets".to_owned(),
    }
}

/// One step line: `  step 3  insert-size  tap on marked  FILES`.
fn step_line(out: &mut String, widths: (usize, usize), step: &Planned, first: bool) {
    let what = if first { first_read(step) } else { files(step) };
    let _ = writeln!(
        out,
        "  step {:<2} {:<verb$}  {:<role$}  {what}",
        step.number,
        step.verb,
        role(step, first),
        verb = widths.0,
        role = widths.1,
    );
}

/// The header of segment `number`: what it reads.
fn segment_line(chain: &Chain, estimate: &Estimate, number: usize, taps_only: bool) -> String {
    let source = shown_source(chain);
    if number == 1 {
        let copied = estimate
            .barriers
            .first()
            .is_some_and(|b| b.how == How::SpoolsSource);
        return if taps_only && !chain.stage_executor {
            format!("segment 1  one read of {source}, shared by its taps as qc shares it")
        } else if copied {
            format!(
                "segment 1  one read of {source}, copied to the temporary directory for the \
                 second read (spool)"
            )
        } else {
            format!("segment 1  one read of {source}")
        };
    }
    let Some(barrier) = estimate.barriers.get(number - 2) else {
        return format!("segment {number}");
    };
    let step = barrier.number;
    match barrier.how {
        How::Holds => format!(
            "segment {number}  the records step {step} holds, sorted by {}",
            barrier.order.map_or("its order", |order| order.name())
        ),
        How::Replays => format!(
            "segment {number}  a second read of {source} through step {step} (replay: nothing \
             is spooled)"
        ),
        How::SpoolsSource => {
            format!("segment {number}  a read of the spooled copy through step {step}")
        }
        How::Spools => format!(
            "segment {number}  a read of the spool of step {step}'s input through step {step} \
             (spool: its input is no file)"
        ),
    }
}

/// `about 1.2 GiB`, or what is unknown.
fn about(bytes: Option<u64>) -> String {
    match bytes {
        Some(bytes) => format!("about {}", human(bytes)),
        None => "unknown (a pipe has no size to read)".to_owned(),
    }
}

/// One line per barrier: how it hands its stream on, what it spills, its
/// memory.
fn barrier_line(chain: &Chain, records: Option<u64>, barrier: &Barrier) -> String {
    let how = match barrier.how {
        How::Holds => "holds its records".to_owned(),
        How::Replays => "replays the source".to_owned(),
        How::SpoolsSource => "replays a copy of the piped source".to_owned(),
        How::Spools => format!("spools its input, {}", about(barrier.spool)),
    };
    let bound = match (barrier.how, barrier.records_in_ram) {
        (How::Holds, None) => format!("--sort-memory {}", human(chain.sort_memory)),
        (How::Holds, Some(n)) => format!(
            "--sort-memory {} and MAX_RECORDS_IN_RAM {}",
            human(chain.sort_memory),
            thousands(n)
        ),
        (_, n) => format!("MAX_RECORDS_IN_RAM {}", thousands(n.unwrap_or(0))),
    };
    let spill = match (barrier.spill, records) {
        (Some(0), Some(n)) => format!("about {} records fit under {bound}", thousands(n)),
        (Some(bytes), _) => format!("spills about {} ({bound})", human(bytes)),
        (None, _) => format!("spills past {bound}"),
    };
    let memory = match barrier.memory {
        Some(bytes) => format!("memory about {}", human(bytes)),
        None => "memory unknown".to_owned(),
    };
    format!(
        "barrier    step {} {}: {how}; {spill}; {memory}",
        barrier.number, barrier.verb
    )
}

/// The lines after the steps: barriers, disk, memory, free space.
fn budget_lines(out: &mut String, chain: &Chain, estimate: &Estimate) {
    for barrier in &estimate.barriers {
        let _ = writeln!(out, "{}", barrier_line(chain, estimate.records, barrier));
    }
    let place = match &chain.tmp_dir {
        Some(dir) => dir.display().to_string(),
        None => "the system temporary directory".to_owned(),
    };
    let temp = match (estimate.barriers.is_empty(), estimate.temp()) {
        (true, _) => "temp disk  none: no step spills".to_owned(),
        (false, Some(0)) => "temp disk  none expected".to_owned(),
        (false, bytes) => format!("temp disk  {} in {place}", about(bytes)),
    };
    let _ = writeln!(out, "{temp}");
    if !estimate.outputs.is_empty() {
        let bytes = estimate.outputs.iter().map(|(_, b)| *b).sum();
        let _ = writeln!(out, "output     about {} written by write", human(bytes));
    }
    let _ = writeln!(out, "memory     {} at the peak", about(estimate.memory()));
    let needs = planner::needs(chain, estimate);
    for need in planner::short(&needs) {
        let _ = writeln!(
            out,
            "space      {} has {} free and the chain needs about {} there: it refuses to run \
             without --force-space",
            need.dir.display(),
            human(need.free.unwrap_or(0)),
            human(need.bytes)
        );
    }
}

/// The plan of `chain` as text.
pub fn text(chain: &Chain) -> Result<String, VerbError> {
    let estimate = planner::estimate(chain)?;
    let mut out = String::new();
    let _ = writeln!(out, "chain plan for {}", shown_source(chain));
    let _ = writeln!(out, "  source    {}", source_line(&estimate.source, chain));
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
        let _ = writeln!(
            out,
            "{}",
            segment_line(chain, &estimate, segment, taps_only)
        );
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
    budget_lines(&mut out, chain, &estimate);
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
}
