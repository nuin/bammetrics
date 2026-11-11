//! `BuildBamIndex`, mirroring Picard 3.5.0 `picard.sam.BuildBamIndex` with
//! htsjdk 5.0.0 `BAMIndexer.createIndex`.
//!
//! Picard opens `INPUT` as a BAM, requires `SO:coordinate` in its header
//! (`Input bam file must be sorted by coordinate` for `queryname`,
//! `unsorted`, `unknown` and a header without the tag), reads every record
//! and writes the `.bai` to `OUTPUT`. [`crate::bai`] holds the format and
//! what htsjdk writes into it; the bytes are Picard's for the same BAM
//! bytes. A reference of length 0 (`@SQ LN:0`) is indexed as Picard indexes
//! it: [`crate::bai::IndexHeader`] reads the header, which noodles' own
//! reader refuses.
//!
//! Without `OUTPUT`, Picard 3.5.0 writes into the current directory, under
//! the file name of `INPUT`: `x.bai` for `dir/x.bam`, `NAME.bai` for any
//! other name. Its help says the index goes beside the BAM; every input
//! takes the branch its source keeps for a URL. The Picard face does what
//! Picard does ([`with_default_output`]). The `build-bam-index` verb and the
//! Python function put the index beside the BAM ([`default_output`]), where
//! htsjdk and samtools look for it.
//!
//! What bammetrics does differently (KNOWN_DIFFERENCES.md):
//!
//! - Records are not validated. Picard's default `VALIDATION_STRINGENCY`
//!   refuses a record that starts past its reference, has a reference and
//!   no position, or is mapped without a CIGAR; bammetrics indexes it as
//!   Picard does under `VALIDATION_STRINGENCY=SILENT`, and fails only where
//!   the indexer itself fails (a reference out of order, an end past 2^29,
//!   a start in a window past the reference).
//! - A failed run writes nothing: the index is written under a temporary
//!   name and renamed into place. Picard leaves an empty `OUTPUT`.
//! - An empty BGZF block in the middle of the file does not end the read.
//! - `OUTPUT` may not be `INPUT`, and a CRAM is refused by name.
//! - The verb and the Python function need an output name for a pipe: a
//!   pipe has no name to put the index beside.
//!
//! The Python function returns the rows `BamIndexStats` prints for the
//! index written: one per reference, then the `*` row.

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::bai::{BamSource, Counts};
use crate::bam_index_stats::{COLUMNS, Row, Stats};
use crate::cleanup::AtomicFile;
use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput};
use crate::registry::{self, Invocation, Opt, Tool};

pub const TOOL: Tool = Tool {
    name: "BuildBamIndex",
    summary: "writes the .bai of a coordinate-sorted BAM, byte for byte the index Picard \
              writes for the same file: to OUTPUT, or beside INPUT as x.bai for x.bam (the \
              BuildBamIndex command without OUTPUT writes that name into the current \
              directory, as Picard does). Returns one row per reference of the header \
              (reference, length, aligned, unaligned), then a * row whose unaligned count is \
              the records without a coordinate.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The BAM to index; its header must say SO:coordinate.",
        )
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The index file (default: beside INPUT, with .bai in place of .bam, or INPUT.bai \
             for another name).",
        )
        .written(),
        Opt::ignored(
            &["REFERENCE_SEQUENCE", "R"],
            "a BAM is indexed without its reference",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

/// The index beside `input`, named as Picard names it: `.bai` in place of
/// a final `.bam`, else appended.
pub fn default_output(input: &Path) -> PathBuf {
    let mut name = input.as_os_str().to_owned();
    if let Some(stem) = input.to_str().and_then(|text| text.strip_suffix(".bam")) {
        name = stem.into();
    }
    name.push(".bai");
    PathBuf::from(name)
}

/// Where Picard 3.5.0 writes the index of `input` without `OUTPUT`:
/// [`default_output`] of its file name, in the current directory.
pub fn picard_default_output(input: &Path) -> PathBuf {
    let name = input.file_name().unwrap_or(input.as_os_str());
    default_output(Path::new(name))
}

/// The Picard face's arguments with `OUTPUT` added when none is given:
/// [`picard_default_output`] of `INPUT`.
pub fn with_default_output(args: &[String]) -> Result<Vec<String>> {
    let inv = TOOL.parse(args)?;
    let mut args = args.to_vec();
    if inv.get("OUTPUT").is_none() {
        let output = picard_default_output(&inv.input);
        args.push(format!("OUTPUT={}", output.display()));
    }
    Ok(args)
}

/// `OUTPUT`, or the default name beside `INPUT`. A device or a pipe has no
/// name to put an index beside.
fn output_path(inv: &Invocation) -> Result<PathBuf> {
    if let Some(output) = inv.get("OUTPUT") {
        return Ok(PathBuf::from(output));
    }
    let regular = std::fs::metadata(&inv.input).is_ok_and(|m| m.is_file());
    if inv.input.starts_with("/dev") || (inv.input.exists() && !regular) {
        return Err(Error::Usage(format!(
            "INPUT={}: a pipe or a device has no name to put the index beside; give OUTPUT",
            inv.input.display()
        )));
    }
    Ok(default_output(&inv.input))
}

/// Index the BAM at `input` into `output`, written under a temporary name
/// and renamed into place. Returns what the index counts.
pub fn build(input: &Path, output: &Path, threads: usize) -> Result<Stats> {
    registry::check_writable([output], &[input.to_path_buf()])?;
    if crate::mark_duplicates::same_file(input, output) {
        return Err(Error::Usage(format!(
            "OUTPUT={}: it is the input {}; writing it would destroy it",
            output.display(),
            input.display()
        )));
    }
    let source = BamSource::open(input, threads)?;
    let order = &source.header().sort_order;
    if order != "coordinate" {
        return Err(Error::Metric(format!(
            "{}: Input bam file must be sorted by coordinate (its header says the sort order \
             is {order})",
            input.display()
        )));
    }
    let rows: Vec<Row> = source
        .header()
        .references
        .iter()
        .map(|(name, length)| Row {
            name: name.clone(),
            length: *length as i64,
            aligned: 0,
            unaligned: 0,
        })
        .collect();
    let written = |source| Error::Output {
        path: output.to_path_buf(),
        source,
    };
    let file = AtomicFile::new(output);
    let out = BufWriter::with_capacity(1 << 16, file.create().map_err(written)?);
    let (mut out, counts) = source.index(out, output)?;
    out.flush().map_err(written)?;
    drop(out);
    file.commit().map_err(written)?;
    Ok(stats(rows, &counts))
}

/// The rows of the header with the counts the index holds.
fn stats(mut rows: Vec<Row>, counts: &Counts) -> Stats {
    for (row, (aligned, unaligned)) in rows.iter_mut().zip(&counts.references) {
        row.aligned = i64::try_from(*aligned).unwrap_or(i64::MAX);
        row.unaligned = i64::try_from(*unaligned).unwrap_or(i64::MAX);
    }
    Stats {
        rows,
        no_coordinate: i64::try_from(counts.no_coordinate).unwrap_or(i64::MAX),
    }
}

/// The tool writes its own file, so the output holds rows and no text:
/// the registry writer has nothing to print.
fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let output = output_path(inv)?;
    let stats = build(&inv.input, &output, inv.threads)?;
    Ok(MetricsOutput::text(&COLUMNS, stats.values(), String::new()).into())
}

#[cfg(test)]
mod tests {
    use noodles::bgzf;

    use super::*;

    #[test]
    fn the_default_name_is_picards() {
        for (input, index, in_cwd) in [
            ("x.bam", "x.bai", "x.bai"),
            ("dir/a.b.bam", "dir/a.b.bai", "a.b.bai"),
            ("/abs/x.BAM", "/abs/x.BAM.bai", "x.BAM.bai"),
            ("../x.cram", "../x.cram.bai", "x.cram.bai"),
            ("x", "x.bai", "x.bai"),
            ("a b/x.bam.gz", "a b/x.bam.gz.bai", "x.bam.gz.bai"),
            ("/dev/stdin", "/dev/stdin.bai", "stdin.bai"),
        ] {
            let input = Path::new(input);
            assert_eq!(default_output(input), Path::new(index));
            assert_eq!(picard_default_output(input), Path::new(in_cwd));
        }
    }

    #[test]
    fn the_picard_face_names_the_index_in_the_current_directory() {
        let args =
            |list: &[&str]| -> Vec<String> { list.iter().map(|a| (*a).to_owned()).collect() };
        assert_eq!(
            with_default_output(&args(&["I=dir/x.bam", "QUIET=true"])).expect("parses"),
            args(&["I=dir/x.bam", "QUIET=true", "OUTPUT=x.bai"])
        );
        let given = args(&["--INPUT", "dir/x.bam", "-O", "y.bai"]);
        assert_eq!(with_default_output(&given).expect("parses"), given);
        assert!(with_default_output(&args(&["QUIET=true"])).is_err());
    }

    #[test]
    fn the_option_table_is_picards() {
        let parse = |args: &[&str]| {
            let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
            TOOL.parse(&args)
        };
        let inv = parse(&[
            "I=x.bam",
            "O=y.bai",
            "VALIDATION_STRINGENCY=SILENT",
            "R=ref.fa",
        ])
        .expect("parses");
        assert_eq!(inv.input, Path::new("x.bam"));
        assert_eq!(inv.get("OUTPUT"), Some("y.bai"));
        assert_eq!(inv.output, None, "the tool writes OUTPUT itself");
        for refused in [
            "ASSUME_SORTED=true",
            "GA4GH_CLIENT_SECRETS=x",
            "STOP_AFTER=1",
        ] {
            assert!(parse(&["I=x.bam", refused]).is_err(), "{refused}");
        }
        assert!(parse(&["O=y.bai"]).is_err(), "INPUT is required");
        assert!(parse(&["I=x.bam", "O=a", "O=b"]).is_err(), "one OUTPUT");
    }

    /// A scratch directory of this test, readable by this user only.
    fn scratch(name: &str) -> PathBuf {
        use std::os::unix::fs::DirBuilderExt;
        let dir =
            std::env::temp_dir().join(format!("bammetrics-index-{}-{name}", std::process::id()));
        std::fs::DirBuilder::new()
            .mode(0o700)
            .create(&dir)
            .expect("scratch dir");
        dir
    }

    /// A coordinate-sorted BAM at `path`: the references, then the record
    /// stream as given.
    fn write_bam(path: &Path, references: &[(&str, u32)], records: &[u8]) {
        let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
        for (name, length) in references {
            text.push_str(&format!("@SQ\tSN:{name}\tLN:{length}\n"));
        }
        let mut bytes = b"BAM\x01".to_vec();
        bytes.extend_from_slice(&(text.len() as u32).to_le_bytes());
        bytes.extend_from_slice(text.as_bytes());
        bytes.extend_from_slice(&(references.len() as u32).to_le_bytes());
        for (name, length) in references {
            bytes.extend_from_slice(&(name.len() as u32 + 1).to_le_bytes());
            bytes.extend_from_slice(name.as_bytes());
            bytes.push(0);
            bytes.extend_from_slice(&length.to_le_bytes());
        }
        bytes.extend_from_slice(records);
        let file = std::fs::File::create_new(path).expect("a new BAM");
        // Stored blocks: two streams of one length take the same offsets.
        let mut out = bgzf::io::writer::Builder::default()
            .set_compression_level(bgzf::io::writer::CompressionLevel::NONE)
            .build_from_writer(file);
        out.write_all(&bytes).expect("write");
        out.finish().expect("finish");
    }

    /// `count` mapped records of 50 bases on `reference`, 100 bases apart.
    fn reads(reference: i32, count: i32) -> Vec<u8> {
        let mut stream = Vec::new();
        for i in 0..count {
            crate::bam_raw::build("r", 0, reference, i * 100, &[(50, 0)], &[30; 50], b"")
                .write_to(&mut stream)
                .expect("memory");
        }
        stream
    }

    #[test]
    fn a_reference_of_length_0_is_indexed() {
        let dir = scratch("empty-reference");
        let mut records = reads(0, 40);
        records.extend(reads(2, 25));
        // The same bytes but for one digit: the header noodles reads.
        let (empty, one) = (dir.join("empty.bam"), dir.join("one.bam"));
        write_bam(
            &empty,
            &[("a", 70_000), ("zero", 0), ("b", 9_000)],
            &records,
        );
        write_bam(&one, &[("a", 70_000), ("zero", 1), ("b", 9_000)], &records);
        let (empty_bai, one_bai) = (dir.join("empty.bai"), dir.join("one.bai"));
        let stats = build(&empty, &empty_bai, 1).expect("Picard indexes LN:0");
        build(&one, &one_bai, 1).expect("an index");
        let rows: Vec<(&str, i64, i64)> = stats
            .rows
            .iter()
            .map(|row| (row.name.as_str(), row.length, row.aligned))
            .collect();
        assert_eq!(rows, [("a", 70_000, 40), ("zero", 0, 0), ("b", 9_000, 25)]);
        let index = std::fs::read(&empty_bai).expect("the index");
        assert_eq!(index, std::fs::read(&one_bai).expect("the index"));
        assert_eq!(&index[..8], b"BAI\x01\x03\0\0\0");
        std::fs::remove_dir_all(&dir).expect("remove the scratch dir");
    }

    #[test]
    fn a_block_size_past_the_file_is_a_truncated_bam_and_writes_nothing() {
        let dir = scratch("block-size");
        let mut records = reads(0, 3);
        // A fourth record whose block_size claims 0xFFFFFF00 bytes.
        records.extend_from_slice(&0xFFFF_FF00u32.to_le_bytes());
        records.extend_from_slice(&[0u8; 100]);
        let (bam, bai) = (dir.join("x.bam"), dir.join("x.bai"));
        write_bam(&bam, &[("a", 70_000)], &records);
        let error = build(&bam, &bai, 1).expect_err("the record is cut");
        assert!(error.to_string().contains("truncated BAM"), "{error}");
        assert_eq!(std::fs::read_dir(&dir).expect("list").count(), 1);
        std::fs::remove_dir_all(&dir).expect("remove the scratch dir");
    }
}
