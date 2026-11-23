//! Kept products (docs/chaining.md section 23): `--keep` on a transform
//! writes the stream it leaves as `PREFIX.<stream>.bam`, with its `.bai`
//! when the stream is coordinate sorted. The file is the BAM `write` writes
//! at that position, so `chain --resume` can read it in place of the
//! source and rerun only the steps after the transform.

use std::path::Path;

use bammetrics_core::bam_writer;
use bammetrics_core::stage::{Products, Step, WriteOptions};

use super::prepare::{Planned, PlannedFile};
use super::run::Written;
use crate::qc::with_suffix;

/// The planned kind of a kept BAM.
pub const KIND: &str = "kept";
/// The planned kind of its `.bai`, written for a coordinate-sorted stream.
pub const INDEX_KIND: &str = "kept-index";

/// A transform whose stream the run writes as a BAM.
#[derive(Clone, Debug)]
pub struct Kept {
    /// The transform's index in `Chain::steps`.
    pub index: usize,
    /// The stream it leaves.
    pub stream: String,
    pub options: WriteOptions,
}

impl Kept {
    pub fn path(&self) -> &Path {
        &self.options.path
    }
}

/// The kept stream of transform `step` (at `index`) under `prefix`: its
/// BAM and `.bai` join the step's planned files, so the claims, the
/// input check and the writability check see them before the first read.
pub(super) fn plan(prefix: &Path, index: usize, step: &mut Planned) -> Kept {
    let stream = step.leaves.clone().unwrap_or_default();
    let path = with_suffix(prefix, &format!(".{stream}.bam"));
    let mut options = WriteOptions::new(&path);
    options.create_index = true;
    step.files.push(PlannedFile {
        path: path.clone(),
        kind: KIND,
    });
    step.files.push(PlannedFile {
        path: bam_writer::index_path(&path),
        kind: INDEX_KIND,
    });
    Kept {
        index,
        stream,
        options,
    }
}

/// The step of the stage layer that writes the kept stream: a `write` at
/// the position right after the transform.
pub(super) fn step(kept: &Kept, threads: usize, label: &str) -> Step {
    let mut step = Step::write(kept.options.clone(), threads);
    step.label = format!("{label} --keep");
    step
}

/// The files the keep write renamed into place. The `.bai` of a stream
/// that is not coordinate sorted is not written, and that is the rule of
/// `--keep`, not news: the write's note about it is dropped.
pub(super) fn written(products: &Products) -> Vec<Written> {
    products
        .files
        .iter()
        .map(|path| Written {
            path: path.clone(),
            kind: if path.extension().is_some_and(|e| e == "bai") {
                "index"
            } else {
                KIND
            },
        })
        .collect()
}

/// The canonical argument that a step keeps its stream, for the manifest:
/// keeping changes what the step writes, so it changes the step's key.
pub(super) fn argv(kept: &Kept) -> String {
    format!("KEEP={}", kept.path().display())
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn the_index_of_a_kept_stream_is_named_as_htsjdk_names_it() {
        let kept = Kept {
            index: 1,
            stream: "marked".to_owned(),
            options: WriteOptions::new(Path::new("out/s.marked.bam")),
        };
        assert_eq!(argv(&kept), "KEEP=out/s.marked.bam");
        let products = Products {
            files: vec![
                PathBuf::from("out/s.marked.bam"),
                PathBuf::from("out/s.marked.bai"),
            ],
            ..Products::default()
        };
        let kinds: Vec<&str> = written(&products).iter().map(|w| w.kind).collect();
        assert_eq!(kinds, ["kept", "index"]);
    }
}
