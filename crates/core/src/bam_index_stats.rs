//! `BamIndexStats`, mirroring Picard 3.5.0 `BamIndexStats.doWork` with
//! htsjdk's `BAMIndexMetaData.printIndexStats` and `getIndexStats`.
//!
//! Picard opens `INPUT` as a BAM, finds its index beside it as
//! `SamFiles.lookForIndex` does (`x.bai`, `x.csi`, `x.bam.bai`, `x.bam.csi`
//! for `x.bam`, then the same names beside the file a symbolic link points
//! to) and reads the metadata pseudo-bin of every reference: the mapped and
//! unmapped record counts `samtools index` stores, cast to a Java `int`, and
//! the trailing count of records without a coordinate. It prints one line
//! per `@SQ` line of the header, `NAME length=\tLENGTH\tAligned= N\tUnaligned= N`,
//! then `NoCoordinateCount= N`. A reference the index lacks metadata for
//! prints `0` and `0`; an index with no reference counts as one empty
//! reference; a header with more references than the index fails
//! (`ArrayIndexOutOfBoundsException` in Picard). Without an index Picard
//! throws `No index for bam file INPUT`.
//!
//! bammetrics prints the same text to standard output (Picard has no
//! `OUTPUT`). A CRAM, or a BAM on a pipe, has no index to read: its
//! records are counted instead, by reference and unmapped flag, which is
//! what `samtools index` would have stored for the same records. Picard
//! fails on a CRAM (KNOWN_DIFFERENCES.md).

use std::path::{Path, PathBuf};

use noodles::csi::BinningIndex;
use noodles::{bam, csi, sam};

use crate::cram;
use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::pass::{self, Accumulator};
use crate::registry::{Invocation, Opt, Tool};

/// The columns of the rows returned to Python and printed by `--format`:
/// one row per reference, then a `*` row whose `unaligned` is the
/// `NoCoordinateCount`, as `samtools idxstats` prints it.
pub const COLUMNS: [&str; 4] = ["reference", "length", "aligned", "unaligned"];

/// The name of the last row, the records without a coordinate.
pub const NO_COORDINATE: &str = "*";

/// The `TypedDict` of one row in `bammetrics.pyi` (written by hand there,
/// one key per [`COLUMNS`] entry), the return type of `bam_index_stats`.
pub const PYTHON_ROW: &str = "_BamIndexRow";

pub const TOOL: Tool = Tool {
    name: "BamIndexStats",
    summary: "one row per reference of the BAM header (reference, length, aligned, unaligned) \
              read from the .bai or .csi beside INPUT, then a * row whose unaligned count is \
              the NoCoordinateCount; the command prints the text Picard prints.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The BAM with its .bai or .csi beside it; a CRAM is read whole with R.",
        )
        .required(),
        Opt::rejected(&["OUTPUT", "O"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

/// The statistics of one reference of the header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    pub name: String,
    pub length: i64,
    pub aligned: i64,
    pub unaligned: i64,
}

/// What Picard prints: one [`Row`] per header reference and the count of
/// records without a coordinate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Stats {
    pub rows: Vec<Row>,
    pub no_coordinate: i64,
}

impl Stats {
    /// The text `printIndexStats` writes to standard output.
    pub fn text(&self) -> String {
        let mut text = String::new();
        for row in &self.rows {
            text.push_str(&format!(
                "{} length=\t{}\tAligned= {}\tUnaligned= {}\n",
                row.name, row.length, row.aligned, row.unaligned
            ));
        }
        text.push_str(&format!("NoCoordinateCount= {}\n", self.no_coordinate));
        text
    }

    /// The rows of [`COLUMNS`]: one per reference, then the `*` row.
    pub fn values(&self) -> Vec<Vec<Value>> {
        let mut rows: Vec<Vec<Value>> = self
            .rows
            .iter()
            .map(|r| {
                vec![
                    Value::Str(r.name.clone()),
                    Value::Int(r.length),
                    Value::Int(r.aligned),
                    Value::Int(r.unaligned),
                ]
            })
            .collect();
        rows.push(vec![
            Value::Str(NO_COORDINATE.to_owned()),
            Value::Int(0),
            Value::Int(0),
            Value::Int(self.no_coordinate),
        ]);
        rows
    }
}

/// The `(int)` cast `BAMIndexMetaData` applies to the chunk offsets that
/// hold the counts: a count past 2^31 wraps as it does in Picard.
fn java_int(count: u64) -> i64 {
    i64::from(count as u32 as i32)
}

/// `SamFiles.lookForIndex`: the index beside `bam` under the names htsjdk
/// tries, in its order, then (`findIndex`) the same names beside the file a
/// symbolic link resolves to.
pub fn find_index(bam: &Path) -> Option<PathBuf> {
    look_for_index(bam).or_else(|| {
        let real = bam.canonicalize().ok()?;
        (real != bam).then(|| look_for_index(&real)).flatten()
    })
}

fn look_for_index(bam: &Path) -> Option<PathBuf> {
    let name = bam.file_name()?.to_str()?;
    let mut candidates = Vec::new();
    if let Some(stem) = name.strip_suffix(".bam") {
        candidates.push(format!("{stem}.bai"));
        candidates.push(format!("{stem}.csi"));
    }
    candidates.push(format!("{name}.bai"));
    candidates.push(format!("{name}.csi"));
    candidates
        .into_iter()
        .map(|candidate| bam.with_file_name(candidate))
        .find(|path| path.is_file())
}

/// `getIndexStats` then `printIndexStats`: the header's references against
/// the index's metadata.
fn from_index(bam: &Path, index_path: &Path) -> Result<Stats> {
    let header = pass::read_header(bam)?;
    let input = |source| Error::Input {
        path: index_path.to_path_buf(),
        source,
    };
    let index: Box<dyn BinningIndex> = if index_path.extension().is_some_and(|e| e == "csi") {
        Box::new(csi::fs::read(index_path).map_err(input)?)
    } else {
        Box::new(bam::bai::fs::read(index_path).map_err(input)?)
    };
    let mut counts: Vec<(u64, u64)> = index
        .reference_sequences()
        .map(|reference| {
            reference.metadata().map_or((0, 0), |m| {
                (m.mapped_record_count(), m.unmapped_record_count())
            })
        })
        .collect();
    if counts.is_empty() {
        // `result[0] = new BAMIndexMetaData()`: one empty reference.
        counts.push((0, 0));
    }
    let no_coordinate = index.unplaced_unmapped_record_count().unwrap_or(0);
    let references = header.reference_sequences();
    if references.len() > counts.len() {
        return Err(Error::Metric(format!(
            "{}: the index holds {} reference sequences, the header {} (Picard: Index {} out \
             of bounds for length {})",
            index_path.display(),
            index.reference_sequences().count(),
            references.len(),
            counts.len(),
            counts.len()
        )));
    }
    let rows = references
        .iter()
        .zip(&counts)
        .map(|((name, contig), (aligned, unaligned))| Row {
            name: name.to_string(),
            length: contig.length().get() as i64,
            aligned: java_int(*aligned),
            unaligned: java_int(*unaligned),
        })
        .collect();
    Ok(Stats {
        rows,
        no_coordinate: no_coordinate as i64,
    })
}

/// The counts `samtools index` would store, taken from the records: a
/// record without a reference has no coordinate, one with a reference is
/// aligned unless its unmapped flag is set.
struct Counter {
    rows: Vec<Row>,
    no_coordinate: u64,
    aligned: Vec<u64>,
    unaligned: Vec<u64>,
}

impl Counter {
    fn new(header: &sam::Header) -> Self {
        let rows: Vec<Row> = header
            .reference_sequences()
            .iter()
            .map(|(name, contig)| Row {
                name: name.to_string(),
                length: contig.length().get() as i64,
                aligned: 0,
                unaligned: 0,
            })
            .collect();
        let n = rows.len();
        Self {
            rows,
            no_coordinate: 0,
            aligned: vec![0; n],
            unaligned: vec![0; n],
        }
    }
}

impl Accumulator for Counter {
    type Output = Stats;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let Some(id) = record.reference_sequence_id().transpose()? else {
            self.no_coordinate += 1;
            return Ok(());
        };
        let counts = if record.flags().is_unmapped() {
            &mut self.unaligned
        } else {
            &mut self.aligned
        };
        match counts.get_mut(id) {
            Some(count) => *count += 1,
            None => {
                return Err(Error::MalformedRecord {
                    name: String::from_utf8_lossy(record.name().unwrap_or_default()).into_owned(),
                    what: "reference sequence id outside the header",
                });
            }
        }
        Ok(())
    }

    fn finish(mut self) -> Result<Stats> {
        for ((row, aligned), unaligned) in
            self.rows.iter_mut().zip(&self.aligned).zip(&self.unaligned)
        {
            row.aligned = java_int(*aligned);
            row.unaligned = java_int(*unaligned);
        }
        Ok(Stats {
            rows: self.rows,
            no_coordinate: self.no_coordinate as i64,
        })
    }
}

/// The statistics of `INPUT`: from the index beside a BAM, from the records
/// of a CRAM or of a file that is no regular file (a pipe).
pub fn stats(inv: &Invocation) -> Result<Stats> {
    let input = inv.input.as_path();
    let name = input.to_string_lossy();
    if name.ends_with(".bai") || name.ends_with(".csi") {
        return Err(Error::Metric(
            "INPUT should be the BAM file name, not its index file".to_owned(),
        ));
    }
    match find_index(input) {
        Some(index) if !cram::is_cram(input) => from_index(input, &index),
        _ if cram::may_be_cram(input) => inv.read_whole(input, |header| Ok(Counter::new(header))),
        _ => Err(Error::Metric(format!(
            "No index for bam file {}",
            input.display()
        ))),
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let stats = stats(inv)?;
    Ok(MetricsOutput::text(&COLUMNS, stats.values(), stats.text()).into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn text_is_picards_line_per_reference_then_the_no_coordinate_count() {
        let stats = Stats {
            rows: vec![
                Row {
                    name: "chr1".to_owned(),
                    length: 100_000,
                    aligned: 10,
                    unaligned: 0,
                },
                Row {
                    name: "chr2".to_owned(),
                    length: 50_000,
                    aligned: 0,
                    unaligned: 3,
                },
            ],
            no_coordinate: 7,
        };
        assert_eq!(
            stats.text(),
            "chr1 length=\t100000\tAligned= 10\tUnaligned= 0\n\
             chr2 length=\t50000\tAligned= 0\tUnaligned= 3\n\
             NoCoordinateCount= 7\n"
        );
        let rows = stats.values();
        assert_eq!(rows.len(), 3);
        assert_eq!(
            rows[2],
            vec![
                Value::Str("*".to_owned()),
                Value::Int(0),
                Value::Int(0),
                Value::Int(7)
            ]
        );
    }

    #[test]
    fn counts_wrap_as_a_java_int() {
        assert_eq!(java_int(10), 10);
        assert_eq!(java_int(1 << 31), i64::from(i32::MIN));
        assert_eq!(java_int((1 << 32) + 5), 5);
    }

    /// An empty file at `path`, created new with mode 0600; an existing
    /// path (a planted link) fails the test.
    fn touch(path: &Path) {
        use std::os::unix::fs::OpenOptionsExt;
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(path)
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    }

    #[test]
    fn index_names_follow_htsjdk_order() {
        use std::os::unix::fs::DirBuilderExt;
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let dir =
            std::env::temp_dir().join(format!("bammetrics-bai-{}-{nanos}", std::process::id()));
        std::fs::DirBuilder::new()
            .mode(0o700)
            .create(&dir)
            .expect("a new temp dir");
        let bam = dir.join("x.bam");
        touch(&bam);
        assert_eq!(find_index(&bam), None);
        touch(&dir.join("x.bam.csi"));
        assert_eq!(find_index(&bam), Some(dir.join("x.bam.csi")));
        touch(&dir.join("x.bam.bai"));
        assert_eq!(find_index(&bam), Some(dir.join("x.bam.bai")));
        touch(&dir.join("x.csi"));
        assert_eq!(find_index(&bam), Some(dir.join("x.csi")));
        touch(&dir.join("x.bai"));
        assert_eq!(find_index(&bam), Some(dir.join("x.bai")));
        let other = dir.join("y.sam");
        touch(&other);
        touch(&dir.join("y.sam.bai"));
        assert_eq!(find_index(&other), Some(dir.join("y.sam.bai")));
        std::fs::remove_dir_all(&dir).expect("remove the temp dir");
    }
}
