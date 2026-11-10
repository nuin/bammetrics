//! The second read of `MarkDuplicates.doWork`: every record copied to
//! `OUTPUT` with its duplicate flag decided, and the metrics counted.
//!
//! Rules copied from Picard 3.5.0 and htsjdk 5.0.0:
//! * Record `i` is a duplicate when `i` is the next index of the sorted
//!   duplicate list. The list advances by one entry whenever `i` passes the
//!   current entry (in query-name order only when the read name changes),
//!   and in query-name order every record after a duplicate with the same
//!   name is a duplicate too, so secondary, supplementary and unmapped
//!   records of a duplicate template are flagged. Optical duplicates are
//!   found the same way in their own list. A flag already set in the input
//!   is cleared when the record is not a duplicate.
//! * `CLEAR_DT` removes `DT`; with `TAGGING_POLICY` a duplicate gets
//!   `DT:Z:SQ` when optical, else `DT:Z:LB` under `All`. With
//!   `TAG_DUPLICATE_SET_MEMBERS` a mapped primary record of a pair set
//!   gets `DI` (the representative's index) and `DS` (the set size). With
//!   `BARCODE_TAG` and `MOLECULAR_IDENTIFIER_TAG`, every record gets
//!   `contig:position/` (the read's start on the reverse strand, else the
//!   mate's start; contig `null` when unmapped), `/A` or `/B` appended for
//!   a duplex UMI of known strand.
//! * `REMOVE_DUPLICATES` drops duplicates and `REMOVE_SEQUENCING_DUPLICATES`
//!   optical duplicates, after they were counted. Each written record's
//!   `PG` becomes the new `@PG` ID chained to its old value, or is removed
//!   when its old value was not seen in the first read.
//! * A record whose tags were set or removed (always, under the default
//!   `CLEAR_DT=true`) is written with its tags in htsjdk order
//!   ([`crate::bam_raw::Tags`]); any other record keeps its bytes.
//! * `@PG`: one record per distinct `PG` value seen in the first read, in
//!   Java `HashSet` order, IDs from `PROGRAM_RECORD_ID` made unique as
//!   htsjdk's `PgIdGenerator` does, with `VN`, `CL`, `PN` and `PP` (the old
//!   value). `VN` defaults to `bammetrics <version>` and `CL` to the
//!   bammetrics command line, where Picard writes its own.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::bgzf;
use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::program::tag;

use super::build::{
    ReadGroups, TagCache, check_order, java_hash, latin1, strand, string_tag, tag_name,
};
use super::inputs::{self, HeaderTag, SortOrder, base36};
use super::mark::Marks;
use super::md_nm::{self, ContigCache};
use super::{Options, TaggingPolicy};
use crate::bam_raw::{RawRecord, TagValue};
use crate::cleanup::AtomicFile;
use crate::duplicate_metrics::Row;
use crate::error::{Error, Result};
use crate::metric::{HistogramSection, HistogramSeries, MetricsOutput};
use crate::metrics_file::KeyClass;
use crate::reference::Reference;
use crate::spill::Sorted;

/// `MarkDuplicates.NO_SUCH_INDEX`.
const NO_SUCH_INDEX: i64 = i64::MAX;

fn next_index(list: &mut Sorted<'_, u64>) -> Result<i64> {
    Ok(list.next_item()?.map_or(NO_SUCH_INDEX, |i| i as i64))
}

/// The new `@PG` ID for each old `PG` value (`None` for no `PG` tag).
type Chained = HashMap<Option<Vec<u8>>, Vec<u8>>;

/// The output header and the new `@PG` ID for each old `PG` value.
fn output_header(
    header: &sam::Header,
    opts: &Options,
    pg_ids: &[Option<Vec<u8>>],
) -> (sam::Header, Option<Chained>) {
    let mut out = header.clone();
    // htsjdk's BAM writer gives `@HD` the version it writes (1.6), whatever
    // the input said.
    if let Some(map) = out.header_mut().as_mut() {
        *map.version_mut() = noodles::sam::header::record::value::map::header::Version::new(1, 6);
    }
    if opts.assume_sort_order == Some(SortOrder::Queryname) {
        inputs::set_header_tag(&mut out, HeaderTag::GroupOrder, "query");
        inputs::set_header_tag(&mut out, HeaderTag::SortOrder, "unknown");
    }
    for comment in &opts.comments {
        let text = comment.strip_prefix("@CO\t").unwrap_or(comment);
        out.comments_mut().push(java_bytes(text).into());
    }
    let Some(record_id) = opts.program_record_id.as_deref() else {
        return (out, None);
    };
    let mut taken: HashSet<Vec<u8>> = out.programs().as_ref().keys().map(|k| k.to_vec()).collect();
    let mut counter = taken.len() as i32;
    let mut next_id = || {
        let wanted = java_bytes(record_id);
        if taken.insert(wanted.clone()) {
            counter += 1;
            return wanted;
        }
        loop {
            let candidate = java_bytes(&format!("{record_id}.{}", base36(counter)));
            counter += 1;
            if taken.insert(candidate.clone()) {
                return candidate;
            }
        }
    };
    let version = opts
        .program_group_version
        .clone()
        .unwrap_or_else(|| format!("bammetrics {}", env!("CARGO_PKG_VERSION")));
    let command_line = opts
        .program_group_command_line
        .clone()
        .unwrap_or_else(|| opts.command.clone());
    let hashes: Vec<i32> = pg_ids
        .iter()
        .map(|id| id.as_deref().map_or(0, java_hash))
        .collect();
    let mut chained = HashMap::new();
    for i in crate::java::hash_map_order_by_hash(&hashes) {
        let existing = &pg_ids[i];
        let id = next_id();
        let mut map = Map::<Program>::default();
        let fields = map.other_fields_mut();
        fields.insert(tag::VERSION, java_bytes(&version).into());
        fields.insert(tag::COMMAND_LINE, java_bytes(&command_line).into());
        if let Some(name) = opts.program_group_name.as_deref() {
            fields.insert(tag::NAME, java_bytes(name).into());
        }
        if let Some(previous) = existing {
            fields.insert(tag::PREVIOUS_PROGRAM_ID, previous.clone().into());
        }
        out.programs_mut().as_mut().insert(id.clone().into(), map);
        chained.insert(existing.clone(), id);
    }
    (out, Some(chained))
}

/// Text from the command line as htsjdk writes a Java `String`: the low
/// byte of each UTF-16 unit (`StringUtil.stringToBytes`), so `café` is
/// written `caf` and the byte 0xE9.
fn java_bytes(text: &str) -> Vec<u8> {
    text.encode_utf16().map(|unit| unit as u8).collect()
}

/// The fields of one header record as `(tag, value)` byte pairs.
fn fields<'a, K, V>(
    other: impl IntoIterator<Item = (&'a K, &'a V)>,
) -> impl Iterator<Item = (&'a [u8], &'a [u8])>
where
    K: AsRef<[u8; 2]> + 'a,
    V: AsRef<[u8]> + 'a,
{
    other
        .into_iter()
        .map(|(t, v)| (t.as_ref().as_slice(), v.as_ref()))
}

/// One header line: `kind`, then a tab and `TAG:value` per field.
fn header_line<'a>(
    text: &mut Vec<u8>,
    kind: &[u8],
    fields: impl IntoIterator<Item = (&'a [u8], &'a [u8])>,
) {
    text.extend_from_slice(kind);
    for (tag, value) in fields {
        text.push(b'\t');
        text.extend_from_slice(tag);
        text.push(b':');
        text.extend_from_slice(value);
    }
    text.push(b'\n');
}

/// The SAM header text: `@HD`, `@SQ`, `@RG`, `@PG` and `@CO` lines in
/// that order, each record's fields in the order noodles keeps them, the
/// order noodles' own writer uses. Values are written byte for byte, as
/// htsjdk writes them; noodles' writer refuses bytes outside printable
/// ASCII, which a valid BAM may hold (`@RG SM:s\xc3\xa9`, a tab in `@CO`).
fn header_text(header: &sam::Header) -> Vec<u8> {
    let mut text = Vec::new();
    if let Some(map) = header.header() {
        let version = map.version();
        let version = format!("{}.{}", version.major(), version.minor());
        let vn = std::iter::once((b"VN".as_slice(), version.as_bytes()));
        header_line(&mut text, b"@HD", vn.chain(fields(map.other_fields())));
    }
    for (name, map) in header.reference_sequences() {
        let length = map.length().to_string();
        let sn_ln = [
            (b"SN".as_slice(), name.as_ref()),
            (b"LN".as_slice(), length.as_bytes()),
        ];
        header_line(
            &mut text,
            b"@SQ",
            sn_ln.into_iter().chain(fields(map.other_fields())),
        );
    }
    for (id, map) in header.read_groups() {
        let id = std::iter::once((b"ID".as_slice(), id.as_ref()));
        header_line(&mut text, b"@RG", id.chain(fields(map.other_fields())));
    }
    for (id, map) in header.programs().as_ref() {
        let id = std::iter::once((b"ID".as_slice(), id.as_ref()));
        header_line(&mut text, b"@PG", id.chain(fields(map.other_fields())));
    }
    for comment in header.comments() {
        text.extend_from_slice(b"@CO\t");
        text.extend_from_slice(comment);
        text.push(b'\n');
    }
    text
}

/// The BAM header: magic, the text, the reference list.
fn write_bam_header(out: &mut impl Write, header: &sam::Header) -> io::Result<()> {
    let invalid = |e| io::Error::new(io::ErrorKind::InvalidInput, e);
    let text = header_text(header);
    out.write_all(b"BAM\x01")?;
    out.write_all(&i32::try_from(text.len()).map_err(invalid)?.to_le_bytes())?;
    out.write_all(&text)?;
    let references = header.reference_sequences();
    out.write_all(
        &i32::try_from(references.len())
            .map_err(invalid)?
            .to_le_bytes(),
    )?;
    for (name, map) in references {
        let length = u32::try_from(name.len() + 1).map_err(invalid)?;
        out.write_all(&length.to_le_bytes())?;
        out.write_all(name)?;
        out.write_all(&[0])?;
        let length = i32::try_from(usize::from(map.length())).map_err(invalid)?;
        out.write_all(&length.to_le_bytes())?;
    }
    Ok(())
}

/// `MarkDuplicates` second read: write `OUTPUT`, return the metrics file.
pub fn run(opts: &Options, inputs: &[PathBuf], mut marks: Marks<'_>) -> Result<MetricsOutput> {
    let rows = inputs::open(
        inputs,
        opts.reference.as_deref(),
        opts.threads,
        opts.assume_sort_order,
        |mut stream| {
            let order = SortOrder::of(&stream.header);
            check_order(order, opts.assume_sort_order)?;
            let (header, chained) = output_header(&stream.header, opts, &marks.pg_ids);
            let groups = ReadGroups::new(&stream.header);
            let mut rows = BTreeMap::new();
            for id in stream.header.read_groups().keys() {
                let library = groups.library(Some(id.as_ref())).to_owned();
                rows.entry(library.clone()).or_insert_with(|| Row {
                    library,
                    ..Row::default()
                });
            }
            let mut writer = Writer::create(opts, &header)?;
            let input_header = stream.header.clone();
            let contigs = match (stream.cram, opts.reference.as_deref()) {
                (true, Some(fasta)) => {
                    Some(ContigCache::new(Reference::open(fasta, &input_header)?))
                }
                _ => None,
            };
            let mut copier = Copier {
                opts,
                header: &input_header,
                queryname: order == SortOrder::Queryname,
                chained: chained.as_ref(),
                tags: TagCache::default(),
                cram: stream.cram,
                contigs,
            };
            copier.copy(&mut stream, &mut marks, &groups, &mut rows, &mut writer)?;
            let output = writer.finish()?;
            let coordinate = SortOrder::of(&header) == SortOrder::Coordinate;
            Ok((rows, coordinate, output))
        },
    )?;
    let (rows, coordinate, output) = rows;
    let output_error = |source| Error::Output {
        path: opts.output.clone(),
        source,
    };
    let index = if opts.create_index && coordinate {
        Some(write_index(
            output.writing(),
            &opts.output.with_extension("bai"),
        )?)
    } else {
        None
    };
    let md5 = if opts.create_md5_file {
        Some(write_md5(output.writing(), &md5_path(&opts.output))?)
    } else {
        None
    };
    output.commit().map_err(output_error)?;
    for file in index.into_iter().chain(md5) {
        let path = file.path().to_path_buf();
        file.commit()
            .map_err(|source| Error::Output { path, source })?;
    }
    Ok(metrics(rows, &marks))
}

/// The BAI of the BAM at `bam`, into a file to be renamed to `target`.
fn write_index(bam: &Path, target: &Path) -> Result<AtomicFile> {
    let index = noodles::bam::fs::index(bam).map_err(|source| Error::Input {
        path: bam.to_path_buf(),
        source,
    })?;
    let noodles::bam::Index::Bai(index) = index else {
        return Err(Error::Metric(
            "CREATE_INDEX: a contig is too long for a BAI index".to_owned(),
        ));
    };
    let file = AtomicFile::new(target);
    file.create()
        .and_then(|out| noodles::bam::bai::io::Writer::new(out).write_index(&index))
        .map_err(|source| Error::Output {
            path: target.to_path_buf(),
            source,
        })?;
    Ok(file)
}

/// `finalizeAndWriteMetrics`.
fn metrics(rows: BTreeMap<String, Row>, marks: &Marks<'_>) -> MetricsOutput {
    let rows: Vec<Row> = rows
        .into_values()
        .map(|mut row| {
            if let Some(count) = marks
                .library_ids
                .get(&row.library)
                .and_then(|id| marks.optical_by_library.get(id))
            {
                row.read_pair_optical_duplicates = *count as i64;
            }
            row.finish()
        })
        .collect();
    let mut series = Vec::new();
    if let [only] = rows.as_slice()
        && let Some(histogram) = only.roi_histogram()
    {
        series.push(HistogramSeries {
            label: "CoverageMult".to_owned(),
            histogram,
        });
    }
    let key_label = if series.is_empty() { "set_size" } else { "BIN" };
    for (label, histogram) in [
        ("all_sets", &marks.all_sets),
        ("optical_sets", &marks.optical_sets),
        ("non_optical_sets", &marks.non_optical_sets),
    ] {
        series.push(HistogramSeries {
            label: label.to_owned(),
            histogram: histogram.clone(),
        });
    }
    let histogram = series
        .iter()
        .any(|s| !s.histogram.is_empty())
        .then_some(HistogramSection {
            key_label,
            key_class: KeyClass::Double,
            series,
        });
    MetricsOutput::from_rows(&rows, histogram)
}

/// Where `CREATE_MD5_FILE` writes: `OUTPUT` with `.md5` appended.
pub(super) fn md5_path(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_owned();
    name.push(".md5");
    PathBuf::from(name)
}

/// The hex MD5 of the BAM at `path`, as htsjdk writes it, into a file to
/// be renamed to `target`.
fn write_md5(path: &Path, target: &Path) -> Result<AtomicFile> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let mut file = File::open(path).map_err(input)?;
    let mut md5 = crate::md5::Md5::new();
    let mut buffer = vec![0u8; 1 << 20];
    loop {
        let n = file.read(&mut buffer).map_err(input)?;
        if n == 0 {
            break;
        }
        md5.update(&buffer[..n]);
    }
    let file = AtomicFile::new(target);
    file.create()
        .and_then(|mut out| out.write_all(crate::md5::hex(&md5.finish()).as_bytes()))
        .map_err(|source| Error::Output {
            path: target.to_path_buf(),
            source,
        })?;
    Ok(file)
}

/// The BGZF output, compressed on worker threads. Records must come in the
/// header's order, as htsjdk's presorted `SAMFileWriterImpl` asserts:
/// coordinate (reference, then start) or query name; any other order is
/// not checked.
struct Writer {
    path: PathBuf,
    /// `OUTPUT` under its temporary name.
    file: AtomicFile,
    inner: noodles::bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<File>>>,
    order: SortOrder,
    /// Reference and start of the last record written, and its name.
    previous: Option<(i32, i32)>,
    previous_name: Vec<u8>,
}

impl Writer {
    /// `assertPresorted`.
    fn check_order(&mut self, record: &RawRecord) -> Result<()> {
        if let Some((reference, start)) = self.previous {
            let out_of_order = match self.order {
                SortOrder::Coordinate => {
                    let next = record.reference_id();
                    let by_reference = match (reference, next) {
                        (-1, -1) => std::cmp::Ordering::Equal,
                        (-1, _) => std::cmp::Ordering::Greater,
                        (_, -1) => std::cmp::Ordering::Less,
                        _ => reference.cmp(&next),
                    };
                    by_reference
                        .then(start.cmp(&record.alignment_start()))
                        .is_gt()
                }
                _ => self.previous_name.as_slice() > record.name(),
            };
            if out_of_order {
                return Err(Error::Metric(format!(
                    "Alignments added out of order in SAMFileWriterImpl.addAlignment for {}. \
                     Sort order is {}. Offending records are at [{}] and [{}]",
                    self.path.display(),
                    self.order.name(),
                    String::from_utf8_lossy(&self.previous_name),
                    String::from_utf8_lossy(record.name())
                )));
            }
        }
        self.previous = Some((record.reference_id(), record.alignment_start()));
        self.previous_name.clear();
        self.previous_name.extend_from_slice(record.name());
        Ok(())
    }

    fn create(opts: &Options, header: &sam::Header) -> Result<Self> {
        let path = opts.output.clone();
        let error = |source| Error::Output {
            path: path.clone(),
            source,
        };
        let atomic = AtomicFile::new(&path);
        let file = atomic.create().map_err(error)?;
        let level = bgzf::io::writer::CompressionLevel::try_from(opts.compression_level)
            .map_err(|e| error(io::Error::new(io::ErrorKind::InvalidInput, e)))?;
        let workers = NonZero::new(opts.threads).unwrap_or(NonZero::<usize>::MIN);
        let bgzf = bgzf::io::multithreaded_writer::Builder::default()
            .set_compression_level(level)
            .set_worker_count(workers)
            .build_from_writer(BufWriter::with_capacity(1 << 20, file));
        let mut inner = noodles::bam::io::Writer::from(bgzf);
        write_bam_header(inner.get_mut(), header).map_err(error)?;
        Ok(Self {
            path,
            file: atomic,
            inner,
            order: SortOrder::of(header),
            previous: None,
            previous_name: Vec::new(),
        })
    }

    fn write(&mut self, record: &RawRecord) -> Result<()> {
        if matches!(self.order, SortOrder::Coordinate | SortOrder::Queryname) {
            self.check_order(record)?;
        }
        record
            .write_to(self.inner.get_mut())
            .map_err(|source| Error::Output {
                path: self.path.clone(),
                source,
            })
    }

    /// Write the last block and the EOF marker. The caller renames the
    /// file into place.
    fn finish(mut self) -> Result<AtomicFile> {
        let error = |source| Error::Output {
            path: self.path.clone(),
            source,
        };
        let mut file = self.inner.get_mut().finish().map_err(error)?;
        file.flush().map_err(error)?;
        Ok(self.file)
    }
}

/// One record's state for the copy.
struct Copier<'a> {
    opts: &'a Options,
    header: &'a sam::Header,
    queryname: bool,
    chained: Option<&'a Chained>,
    tags: TagCache,
    /// Records come from a CRAM: every record is written with sorted tags.
    cram: bool,
    /// Reference bases for the `MD` and `NM` of CRAM records.
    contigs: Option<ContigCache>,
}

/// The index of the next representative entry, the representative and the set size.
#[derive(Clone, Copy)]
struct SetEntry {
    read_index: i64,
    representative: i32,
    size: i32,
}

impl Copier<'_> {
    /// `nextIndexIfNeeded`.
    fn advance(
        &self,
        current: i64,
        index: i64,
        last_name: Option<&[u8]>,
        name: &[u8],
        list: Option<&mut Sorted<'_, u64>>,
    ) -> Result<i64> {
        let need = index > current && (!self.queryname || last_name != Some(name));
        match list {
            Some(list) if need => next_index(list),
            _ => Ok(current),
        }
    }

    fn copy(
        &mut self,
        stream: &mut inputs::Inputs<'_>,
        marks: &mut Marks<'_>,
        groups: &ReadGroups,
        rows: &mut BTreeMap<String, Row>,
        writer: &mut Writer,
    ) -> Result<()> {
        let opts = self.opts;
        let mut record = RawRecord::default();
        let mut index: i64 = 0;
        let mut next_duplicate = next_index(&mut marks.duplicates)?;
        let mut next_optical = match marks.optical.as_mut() {
            Some(list) => next_index(list)?,
            None => NO_SUCH_INDEX,
        };
        let mut set = SetEntry {
            read_index: -1,
            representative: -1,
            size: -1,
        };
        if let Some(list) = marks.representatives.as_mut()
            && let Some(first) = list.next_item()?
        {
            set = SetEntry {
                read_index: i64::from(first.read_index),
                representative: first.representative_index,
                size: first.set_size,
            };
        }
        let mut duplicate_name: Option<Vec<u8>> = None;
        let mut representative_name: Option<Vec<u8>> = None;
        let mi_tag = match (&opts.barcode_tag, &opts.molecular_identifier_tag) {
            (Some(_), Some(tag)) => Some(tag_name(tag)?),
            _ => None,
        };
        while stream.next(&mut record)? {
            self.tags.load(&record)?;
            let name = record.name().to_vec();
            let flags = record.flags();
            let library = groups.library(string_tag(&self.tags, *b"RG")?).to_owned();
            let row = rows.entry(library).or_insert_with_key(|library| Row {
                library: library.clone(),
                ..Row::default()
            });
            row.add_read(flags);
            next_duplicate = self.advance(
                next_duplicate,
                index,
                duplicate_name.as_deref(),
                &name,
                Some(&mut marks.duplicates),
            )?;
            let duplicate = index == next_duplicate
                || (self.queryname
                    && index > next_duplicate
                    && duplicate_name.as_deref() == Some(name.as_slice()));
            let flags = if duplicate {
                flags | 0x400
            } else {
                flags & !0x400
            };
            record.set_flags(flags);
            if duplicate {
                row.add_duplicate(flags);
            }
            next_optical = self.advance(
                next_optical,
                index,
                duplicate_name.as_deref(),
                &name,
                marks.optical.as_mut(),
            )?;
            let optical = (self.queryname
                && index > next_optical
                && duplicate_name.as_deref() == Some(name.as_slice()))
                || index == next_optical;
            let mut edits: Vec<([u8; 2], Option<TagValue>)> = Vec::new();
            if opts.clear_dt {
                edits.push((*b"DT", None));
            }
            if opts.tagging_policy != TaggingPolicy::DontTag && duplicate {
                if optical {
                    edits.push((*b"DT", Some(TagValue::String(b"SQ".to_vec()))));
                } else if opts.tagging_policy == TaggingPolicy::All {
                    edits.push((*b"DT", Some(TagValue::String(b"LB".to_vec()))));
                }
            }
            if let Some(list) = marks.representatives.as_mut() {
                let need = index > set.read_index
                    && (!self.queryname || representative_name.as_deref() != Some(name.as_slice()));
                if need && let Some(next) = list.next_item()? {
                    set = SetEntry {
                        read_index: i64::from(next.read_index),
                        representative: next.representative_index,
                        size: next.set_size,
                    };
                }
                let in_set = index == set.read_index
                    || (self.queryname
                        && index > set.read_index
                        && representative_name.as_deref() == Some(name.as_slice()));
                if in_set && flags & 0x900 == 0 && flags & 0x4 == 0 {
                    edits.push((*b"DI", Some(TagValue::Int(i64::from(set.representative)))));
                    edits.push((*b"DS", Some(TagValue::Int(i64::from(set.size)))));
                    representative_name = Some(name.clone());
                }
            }
            if let Some(tag) = mi_tag {
                edits.push((
                    tag,
                    Some(TagValue::String(self.molecular_identifier(&record)?)),
                ));
            }
            if duplicate {
                duplicate_name = Some(name);
            }
            index += 1;
            if (opts.remove_duplicates && duplicate)
                || (opts.remove_sequencing_duplicates && optical)
            {
                continue;
            }
            if let Some(chained) = self.chained
                && opts.add_pg_tag_to_reads
            {
                let old = string_tag(&self.tags, *b"PG")?.map(<[u8]>::to_vec);
                let new = chained.get(&old).map(|id| TagValue::String(id.clone()));
                edits.push((*b"PG", new));
            }
            if self.cram || !edits.is_empty() {
                let mut tags = record.tags().map_err(Error::Record)?;
                if let Some(contigs) = self.contigs.as_mut() {
                    md_nm::restore(&record, &mut tags, contigs)?;
                }
                for (tag, value) in edits {
                    match value {
                        Some(value) => tags.set(tag, value),
                        None => tags.remove(tag),
                    }
                }
                record.set_tags(&tags).map_err(Error::Record)?;
            }
            writer.write(&record)?;
        }
        Ok(())
    }

    /// `UmiUtil.setMolecularIdentifier` with an empty assigned UMI.
    fn molecular_identifier(&self, record: &RawRecord) -> Result<Vec<u8>> {
        let flags = record.flags();
        let contig = if flags & 0x4 != 0 {
            "null".to_owned()
        } else {
            usize::try_from(record.reference_id())
                .ok()
                .and_then(|i| self.header.reference_sequences().get_index(i))
                .map_or_else(|| "*".to_owned(), |(name, _)| latin1(name))
        };
        let position = if flags & 0x10 != 0 {
            record.alignment_start()
        } else {
            record.mate_alignment_start()
        };
        let mut text = format!("{contig}:{position}/");
        if self.opts.duplex_umi {
            match strand(record, &self.tags)? {
                Some(true) => text.push_str("/A"),
                Some(false) => text.push_str("/B"),
                None => {}
            }
        }
        Ok(text.into_bytes())
    }
}
