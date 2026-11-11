//! The files of a chain, each counted once whatever its spelling (a link,
//! a hard link and `..` name the file they lead to): every planned output
//! claimed by its step, and every file the chain reads. Two steps that
//! write one file, a step that writes a file the chain reads, and two
//! writers of standard output are usage errors, found before the first
//! read.

use std::path::{Path, PathBuf};

use bammetrics_core::stage::Outputs;

use super::prepare::Planned;
use super::{VERB, usage};
use crate::error::VerbError;
use crate::qc::with_suffix;

/// The files a chain reads beside those its steps name.
pub(super) struct Read<'a> {
    pub source: &'a Path,
    pub reference: Option<&'a Path>,
    /// The plan file of `--plan`.
    pub plan_file: Option<&'a Path>,
}

/// Every planned file is claimed by its step: two steps that name one
/// file are a usage error.
pub(super) fn claim_outputs(
    steps: &[Planned],
    report: Option<&Path>,
) -> Result<Outputs, VerbError> {
    let mut outputs = Outputs::default();
    for step in steps {
        for file in &step.files {
            outputs.claim(&step.label, &file.path).map_err(|e| {
                VerbError::core(VERB, &e)
                    .hint("a step's own -o names another file (bammetrics chain --help)")
            })?;
        }
    }
    if let Some(report) = report {
        outputs
            .claim("--also markdown", report)
            .map_err(|e| VerbError::core(VERB, &e))?;
    }
    Ok(outputs)
}

/// The files a chain reads: the source, the reference with its index, the
/// plan file, and every file a step reads by name. The same file through a
/// link or `..` counts as one.
pub(super) fn inputs_of(read: &Read, steps: &[Planned]) -> Outputs {
    let mut inputs = Outputs::default();
    let fai = read
        .reference
        .map(|reference| with_suffix(reference, ".fai"));
    let files = std::iter::once(read.source)
        .chain(read.reference)
        .chain(fai.as_deref())
        .chain(read.plan_file)
        .chain(
            steps
                .iter()
                .flat_map(|step| &step.inputs)
                .map(PathBuf::as_path),
        );
    for path in files {
        // One file read by several steps is claimed once; the error of a
        // second claim says nothing new.
        let _ = inputs.claim("an input", path);
    }
    inputs
}

/// A planned file that is a file the chain reads would destroy it: a usage
/// error.
pub(super) fn refuse_overwriting_inputs(
    inputs: &Outputs,
    steps: &[Planned],
    report: Option<&Path>,
) -> Result<(), VerbError> {
    let also = "--also markdown".to_owned();
    let planned = steps
        .iter()
        .flat_map(|step| step.files.iter().map(move |file| (&step.label, &file.path)))
        .map(|(label, path)| (label, path.as_path()))
        .chain(report.map(|report| (&also, report)));
    for (label, path) in planned {
        if inputs.owner(path).is_some() {
            return Err(usage(format!(
                "{label}: {} is an input; writing it would destroy it",
                path.display()
            )));
        }
    }
    Ok(())
}

/// At most one step writes standard output.
pub(super) fn check_one_stdout(steps: &[Planned]) -> Result<(), VerbError> {
    let mut writers = steps.iter().filter(|step| step.stdout);
    match (writers.next(), writers.next()) {
        (Some(first), Some(second)) => Err(usage(format!(
            "{} and {} both write standard output; it takes one file",
            first.label, second.label
        ))),
        _ => Ok(()),
    }
}
