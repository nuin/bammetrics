//! `CollectIndexDepthMetrics` (docs/06 section 7): depth estimated from the
//! index of a BAM alone. Not a Picard tool.
//!
//! The `.bai` or `.csi` beside the BAM (`bam_index_stats::find_index`)
//! holds, per reference, a metadata pseudo-bin with the count of records
//! without the unmapped flag (`MAPPED`) and with it (`UNMAPPED_PLACED`),
//! and after the last reference the count of records without a coordinate.
//! [`Summary`] reads those counts and holds no bin. The read length is the
//! mean `l_seq` of the mapped primary records (not secondary, not
//! supplementary, with a sequence) among the first `READ_LENGTH_SAMPLE`
//! records of the file, or `READ_LENGTH` when given; then no record is
//! read. `ESTIMATED_MEAN_DEPTH` is `MAPPED x READ_LENGTH / LENGTH`.
//!
//! The reference names and lengths are the BAM's binary reference list
//! (`doctor::input::bam_header`: the text and the list under one 64 MiB
//! cap), the list the index is built against, so a reference of length 0
//! is kept. Rows: one per reference in header order, one per group
//! ([`Group`], from the names of `genome_build.rs`), then `ALL`. The sex
//! call is `sex_check::row` on `MAPPED` and `LENGTH`: reads per base, the
//! quantity `CollectSexCheckMetrics` compares, with its thresholds.
//! `MT_COPY_NUMBER` is twice the mitochondrial ratio to the autosomes, the
//! estimate of `CollectMtCopyNumberMetrics`.
//!
//! The index counts every record: duplicates, secondary and supplementary
//! records, QC failures and MAPQ 0 included. Nothing filters them, so the
//! estimate runs above a filtered depth by the duplicate rate and the other
//! records (docs/06 section 7.4).

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use crate::bai::{MAX_INDEX_BYTES, Summary};
use crate::bam_index_stats::find_index;
use crate::bam_raw::RawRecord;
use crate::cram::{self, Kind};
use crate::doctor::input::{BamStream, bam_header};
use crate::error::{Error, Result};
use crate::genome_build::{self, Kind as ContigKind};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};
use crate::sex_check;

/// The metrics class string. Not a Picard class: readers key on the last
/// name.
pub const CLASS: &str = "bammetrics.qc.IndexDepthMetrics";

/// Records read for the read length when `READ_LENGTH` is not given.
pub const DEFAULT_READ_LENGTH_SAMPLE: i64 = 10_000;

/// The most records `READ_LENGTH_SAMPLE` may ask for: about 0.3 s of BAM.
pub const MAX_READ_LENGTH_SAMPLE: i64 = 1_000_000;

/// The flags that keep a record out of the read length sample.
const UNMAPPED: u16 = 0x4;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

/// The group a contig is summed into, from its name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Group {
    Autosome,
    X,
    Y,
    Mt,
    Alt,
    Decoy,
    Hla,
    Other,
}

impl Group {
    /// The group rows, in output order.
    pub const ALL: [Group; 8] = [
        Group::Autosome,
        Group::X,
        Group::Y,
        Group::Mt,
        Group::Alt,
        Group::Decoy,
        Group::Hla,
        Group::Other,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Group::Autosome => "AUTOSOME",
            Group::X => "X",
            Group::Y => "Y",
            Group::Mt => "MT",
            Group::Alt => "ALT",
            Group::Decoy => "DECOY",
            Group::Hla => "HLA",
            Group::Other => "OTHER",
        }
    }

    /// `chr1`..`chr22` and `1`..`22`, `chrX` or `X`, `chrY` or `Y`, the
    /// mitochondrial names, then the alt, decoy and HLA kinds of
    /// `genome_build.rs`; every other name is `OTHER`.
    pub fn of(name: &str) -> Self {
        match genome_build::chromosome(name) {
            Some("X") => Group::X,
            Some("Y") => Group::Y,
            Some("M") => Group::Mt,
            Some(_) => Group::Autosome,
            None => match genome_build::kind(name) {
                ContigKind::Alt => Group::Alt,
                ContigKind::Decoy => Group::Decoy,
                ContigKind::Hla => Group::Hla,
                _ => Group::Other,
            },
        }
    }
}

/// The name of the last row.
pub const ALL_ROW: &str = "ALL";

/// One reference: its binary-list name and length and its index counts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub length: u64,
    pub mapped: u64,
    pub unmapped: u64,
}

/// What the index says of a BAM.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Counts {
    pub contigs: Vec<Contig>,
    pub no_coordinate: u64,
}

/// One row of the metrics file.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub contig: Option<String>,
    pub group: &'static str,
    pub length: i64,
    pub mapped: i64,
    pub unmapped_placed: i64,
    pub estimated_mean_depth: f64,
    pub ratio_to_autosomes: f64,
    pub read_length: Option<f64>,
    pub no_coordinate: Option<i64>,
    pub inferred_sex: Option<String>,
    pub mt_copy_number: Option<f64>,
    pub reason: Option<String>,
}

impl MetricRow for Row {
    const CLASS: &'static str = CLASS;
    const COLUMNS: &'static [&'static str] = &[
        "CONTIG",
        "GROUP",
        "LENGTH",
        "MAPPED",
        "UNMAPPED_PLACED",
        "ESTIMATED_MEAN_DEPTH",
        "RATIO_TO_AUTOSOMES",
        "READ_LENGTH",
        "NO_COORDINATE",
        "INFERRED_SEX",
        "MT_COPY_NUMBER",
        "REASON",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.contig.clone().into(),
            self.group.into(),
            self.length.into(),
            self.mapped.into(),
            self.unmapped_placed.into(),
            self.estimated_mean_depth.into(),
            self.ratio_to_autosomes.into(),
            self.read_length.into(),
            self.no_coordinate.into(),
            self.inferred_sex.clone().into(),
            self.mt_copy_number.into(),
            self.reason.clone().into(),
        ]
    }
}

/// A count or length as a metrics integer: the index stores 64-bit
/// unsigned counts, a value past `i64::MAX` saturates.
fn int(value: u64) -> i64 {
    i64::try_from(value).unwrap_or(i64::MAX)
}

/// Mapped records, placed unmapped records and bases of a set of contigs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Sum {
    length: u64,
    mapped: u64,
    unmapped: u64,
}

impl Sum {
    fn add(&mut self, contig: &Contig) {
        self.length = self.length.saturating_add(contig.length);
        self.mapped = self.mapped.saturating_add(contig.mapped);
        self.unmapped = self.unmapped.saturating_add(contig.unmapped);
    }

    /// Mapped records per base, the quantity `sex_check::row` divides.
    fn per_base(self) -> f64 {
        self.mapped as f64 / self.length as f64
    }

    /// `MAPPED x READ_LENGTH / LENGTH`: `?` for a length of 0, 0 for no
    /// mapped record (whatever the read length), `?` for an unknown read
    /// length otherwise.
    fn depth(self, read_length: Option<f64>) -> f64 {
        if self.length == 0 {
            f64::NAN
        } else if self.mapped == 0 {
            0.0
        } else {
            read_length.map_or(f64::NAN, |rl| self.mapped as f64 * rl / self.length as f64)
        }
    }

    /// The row of these counts.
    fn row(
        self,
        contig: Option<String>,
        group: &'static str,
        read_length: Option<f64>,
        autosomes: Sum,
    ) -> Row {
        Row {
            contig,
            group,
            length: int(self.length),
            mapped: int(self.mapped),
            unmapped_placed: int(self.unmapped),
            estimated_mean_depth: self.depth(read_length),
            ratio_to_autosomes: self.per_base() / autosomes.per_base(),
            read_length,
            no_coordinate: None,
            inferred_sex: None,
            mt_copy_number: None,
            reason: None,
        }
    }

    fn pair(self) -> (i64, i64) {
        (int(self.mapped), int(self.length))
    }
}

/// The rows of `counts` at `read_length`: one per contig, one per
/// [`Group`], then `ALL` with the calls. `sample` is `READ_LENGTH_SAMPLE`
/// when the read length was sampled, for the reason of an unknown one.
pub fn rows(counts: &Counts, read_length: Option<f64>, sample: Option<i64>) -> Vec<Row> {
    let mut groups = [Sum::default(); 8];
    let mut all = Sum::default();
    let mut of_contig = Vec::with_capacity(counts.contigs.len());
    for contig in &counts.contigs {
        let group = Group::of(&contig.name);
        groups[group as usize].add(contig);
        all.add(contig);
        of_contig.push(group);
    }
    let autosomes = groups[Group::Autosome as usize];
    let mut rows: Vec<Row> = counts
        .contigs
        .iter()
        .zip(&of_contig)
        .map(|(contig, group)| {
            let mut one = Sum::default();
            one.add(contig);
            one.row(
                Some(contig.name.clone()),
                group.name(),
                read_length,
                autosomes,
            )
        })
        .collect();
    for group in Group::ALL {
        rows.push(groups[group as usize].row(None, group.name(), read_length, autosomes));
    }
    let sex = sex_check::row(
        autosomes.pair(),
        groups[Group::X as usize].pair(),
        groups[Group::Y as usize].pair(),
        None,
        0,
    );
    let mt = groups[Group::Mt as usize];
    let mut reasons: Vec<String> = sex.reason.into_iter().collect();
    if mt.length == 0 {
        reasons.push("no mitochondrial contig".to_owned());
    }
    if read_length.is_none() && all.mapped > 0 {
        reasons.push(match sample {
            Some(n) => format!(
                "no mapped primary record with a sequence among the first {n} records; give \
                 READ_LENGTH"
            ),
            None => "no read length".to_owned(),
        });
    }
    let mut last = all.row(None, ALL_ROW, read_length, autosomes);
    last.no_coordinate = Some(int(counts.no_coordinate));
    last.inferred_sex = Some(sex.inferred_sex);
    last.mt_copy_number = Some(2.0 * mt.per_base() / autosomes.per_base());
    last.reason = (!reasons.is_empty()).then(|| reasons.join("; "));
    rows.push(last);
    rows
}

/// The options of one run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Options {
    /// `READ_LENGTH`: no record is read when given.
    pub read_length: Option<i64>,
    pub read_length_sample: i64,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let read_length: Option<i64> = inv.parse_opt("READ_LENGTH")?;
        if let Some(rl) = read_length
            && !(1..=i64::from(i32::MAX)).contains(&rl)
        {
            return Err(Error::Usage(format!("READ_LENGTH={rl}: 1 to {}", i32::MAX)));
        }
        let read_length_sample = inv.parse("READ_LENGTH_SAMPLE", DEFAULT_READ_LENGTH_SAMPLE)?;
        if !(1..=MAX_READ_LENGTH_SAMPLE).contains(&read_length_sample) {
            return Err(Error::Usage(format!(
                "READ_LENGTH_SAMPLE={read_length_sample}: 1 to {MAX_READ_LENGTH_SAMPLE}"
            )));
        }
        Ok(Self {
            read_length,
            read_length_sample,
        })
    }
}

/// The refusal for a CRAM: its `.crai` holds offsets and spans, no count.
fn cram_refusal(path: &Path) -> Error {
    Error::Metric(format!(
        "{}: a CRAM index (.crai) holds no read counts, so CollectIndexDepthMetrics reads \
         none; samtools idxstats counts a CRAM by reading its containers, and bammetrics \
         sex-check or mt-copy-number read every record",
        path.display()
    ))
}

/// Why `path` cannot be read from its index alone; `None` for a regular
/// file that holds a BAM.
fn refusal(path: &Path) -> Result<Option<Error>> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let name = path.to_string_lossy();
    if [".bai", ".csi", ".crai"].iter().any(|s| name.ends_with(s)) {
        return Ok(Some(Error::Usage(
            "INPUT should be the BAM file name, not its index file".to_owned(),
        )));
    }
    let regular = path != Path::new("-")
        && std::fs::metadata(path).map_err(input)?.is_file()
        && !path.canonicalize().is_ok_and(|p| p.starts_with("/dev"));
    if !regular {
        return Ok(Some(Error::Metric(format!(
            "{name}: not a regular file; a BAM on a pipe has no index, and \
             CollectIndexDepthMetrics never reads a whole file"
        ))));
    }
    let kind = cram::sniff(File::open(path).map_err(input)?)
        .map_err(input)?
        .kind;
    Ok(match kind {
        Kind::Bam => None,
        Kind::Cram => Some(cram_refusal(path)),
        Kind::Sam | Kind::GzipSam => Some(Error::Metric(format!(
            "{name}: SAM text; CollectIndexDepthMetrics reads the index of a BAM"
        ))),
    })
}

/// The index beside `bam`, or the refusal that names how to make one.
fn index_of(bam: &Path) -> Result<PathBuf> {
    find_index(bam).ok_or_else(|| {
        Error::Metric(format!(
            "{}: no .bai or .csi beside the BAM (x.bai, x.csi, x.bam.bai or x.bam.csi); \
             CollectIndexDepthMetrics reads the index alone and never the whole file. \
             samtools index or bammetrics build-bam-index writes one",
            bam.display()
        ))
    })
}

/// The counts of the first `keep` references of the index at `path`, at
/// most [`MAX_INDEX_BYTES`] of it (inflated, for a `.csi`).
fn read_index(path: &Path, keep: usize) -> Result<Summary> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = BufReader::new(File::open(path).map_err(input)?);
    let csi = path.extension().is_some_and(|e| e == "csi");
    let inner: Box<dyn Read> = if csi {
        Box::new(noodles::bgzf::io::Reader::new(file))
    } else {
        Box::new(file)
    };
    let mut capped = inner.take(MAX_INDEX_BYTES + 1);
    let summary = if csi {
        Summary::read_csi(&mut capped, keep)
    } else {
        Summary::read(&mut capped, keep)
    };
    if capped.limit() == 0 {
        return Err(input(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("the index is longer than {MAX_INDEX_BYTES} bytes"),
        )));
    }
    summary.map_err(input)
}

/// A read that ends before the bytes it was promised: the BAM is cut, or a
/// record claims more bytes than the file holds.
fn truncated(error: std::io::Error) -> std::io::Error {
    if error.kind() != std::io::ErrorKind::UnexpectedEof {
        return error;
    }
    std::io::Error::new(
        std::io::ErrorKind::UnexpectedEof,
        "truncated BAM: a record of the read length sample ends past the end of the file; \
         READ_LENGTH skips the sample",
    )
}

/// The smallest BGZF block: an 18-byte header, 2 bytes of deflate and an
/// 8-byte footer. Each block inflates to at most 65,536 bytes.
const MIN_BLOCK_BYTES: u128 = 28;
const MAX_BLOCK_DATA: u128 = 65_536;
/// The smallest BAM record: `block_size` and the 32 fixed bytes.
const MIN_RECORD_BYTES: u128 = 36;

/// The most records a BAM of `bytes` bytes can hold.
fn max_records(bytes: u64) -> u128 {
    u128::from(bytes) / MIN_BLOCK_BYTES * MAX_BLOCK_DATA / MIN_RECORD_BYTES
}

/// The mean `l_seq` of the mapped primary records with a sequence among
/// the first `sample` records of `stream`; `None` when there is none.
fn sample_read_length(stream: &mut BamStream, sample: i64, path: &Path) -> Result<Option<f64>> {
    let mut record = RawRecord::default();
    let (mut bases, mut reads) = (0u64, 0u64);
    for _ in 0..sample {
        let more = record.read_from(stream).map_err(|source| Error::Input {
            path: path.to_path_buf(),
            source: truncated(source),
        })?;
        if !more {
            break;
        }
        let length = record.sequence_length() as u64;
        if record.flags() & (UNMAPPED | SECONDARY | SUPPLEMENTARY) == 0 && length > 0 {
            bases += length;
            reads += 1;
        }
    }
    Ok((reads > 0).then(|| bases as f64 / reads as f64))
}

/// The counts and the read length of the BAM at `bam`, from its index and
/// its first records.
pub fn measure(bam: &Path, opts: &Options) -> Result<(Counts, Option<f64>)> {
    if let Some(refused) = refusal(bam)? {
        return Err(refused);
    }
    let index = index_of(bam)?;
    let (header, mut stream) = bam_header(bam).map_err(|source| Error::Input {
        path: bam.to_path_buf(),
        source,
    })?;
    let read_length = match opts.read_length {
        Some(rl) => Some(rl as f64),
        None => sample_read_length(&mut stream, opts.read_length_sample, bam)?,
    };
    let references = header.binary.len();
    let summary = read_index(&index, references)?;
    if summary.references != references {
        return Err(Error::Metric(format!(
            "{}: the index holds {} reference sequences, the BAM header {references}; it is \
             the index of another file",
            index.display(),
            summary.references
        )));
    }
    let contigs = header
        .binary
        .into_iter()
        .zip(summary.metadata)
        .map(|((name, length), metadata)| {
            let (mapped, unmapped) = metadata.map_or((0, 0), |m| (m.aligned, m.unaligned));
            Contig {
                name,
                length,
                mapped,
                unmapped,
            }
        })
        .collect();
    let counts = Counts {
        contigs,
        no_coordinate: summary.no_coordinate.unwrap_or(0),
    };
    let bytes = std::fs::metadata(bam)
        .map_err(|source| Error::Input {
            path: bam.to_path_buf(),
            source,
        })?
        .len();
    let records: u128 = counts
        .contigs
        .iter()
        .map(|c| u128::from(c.mapped) + u128::from(c.unmapped))
        .sum::<u128>()
        + u128::from(counts.no_coordinate);
    if records > max_records(bytes) {
        return Err(Error::Metric(format!(
            "{}: the index counts {records} records, more than the {} a BAM of {bytes} bytes \
             can hold; it is damaged or the index of another file",
            index.display(),
            max_records(bytes)
        )));
    }
    Ok((counts, read_length))
}

pub const TOOL: Tool = Tool {
    name: "CollectIndexDepthMetrics",
    summary: "one row of the 12 IndexDepthMetrics columns per reference, per contig group \
              (AUTOSOME, X, Y, MT, ALT, DECOY, HLA, OTHER) and for ALL: the mapped and placed \
              unmapped counts of the .bai or .csi, MAPPED x READ_LENGTH / LENGTH as the \
              depth, the ratio to the autosomes, and on ALL the sex call, MT_COPY_NUMBER and \
              NO_COORDINATE. No whole read. Not a Picard tool (docs/06).",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The BAM with its .bai or .csi beside it; a CRAM, SAM text or a pipe is refused.",
        )
        .required(),
        Opt::int(
            &["READ_LENGTH"],
            "The read length in bases; then no record is read (default the mean of the \
             sample).",
        ),
        Opt::int(
            &["READ_LENGTH_SAMPLE"],
            "Records read from the start of the BAM for the read length (default 10000, at \
             most 1000000).",
        ),
        Opt::ignored(
            &["VALIDATION_STRINGENCY"],
            "the index and the sampled records are read as they are",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["REFERENCE_SEQUENCE", "R"]),
        Opt::rejected(&["STOP_AFTER"]),
        Opt::rejected(&["METRIC_ACCUMULATION_LEVEL", "LEVEL"]),
        Opt::rejected(&["MINIMUM_MAPPING_QUALITY", "MQ"]),
    ],
    run: run_tool,
    collector: None,
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    // `OUTPUT` may not replace the index the run reads.
    if let Some(index) = find_index(&inv.input) {
        crate::registry::check_writable(inv.output.as_deref(), &[index])?;
    }
    let (counts, read_length) = measure(&inv.input, &opts)?;
    let sample = opts
        .read_length
        .is_none()
        .then_some(opts.read_length_sample);
    Ok(MetricsOutput::from_rows(&rows(&counts, read_length, sample), None).into())
}

#[cfg(test)]
mod tests;
