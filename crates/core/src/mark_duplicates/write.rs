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
//!
//! [`Copier`] holds the cursor of that read and decides one record per
//! call, handing each kept record to an [`Emit`]: the BAM writer for the
//! tool, the next step's batch for a chain ([`crate::stage`]), where the
//! copier is the [`Transform`] the sealed barrier replays the source
//! through.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::PathBuf;

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
use crate::bam_writer::{self, BamWriter, Companions};
use crate::duplicate_metrics::Row;
use crate::error::{Error, Result};
use crate::metric::{HistogramSection, HistogramSeries, MetricsOutput};
use crate::metrics_file::KeyClass;
use crate::reference::Reference;
use crate::spill::Sorted;
use crate::stage::{Counters, Emit, Products, Transform};

/// `MarkDuplicates.NO_SUCH_INDEX`.
const NO_SUCH_INDEX: i64 = i64::MAX;

fn next_index(list: &mut Sorted<'_, u64>) -> Result<i64> {
    Ok(list.next_item()?.map_or(NO_SUCH_INDEX, |i| i as i64))
}

/// The new `@PG` ID for each old `PG` value (`None` for no `PG` tag).
pub(super) type Chained = HashMap<Option<Vec<u8>>, Vec<u8>>;

/// The output header and the new `@PG` ID for each old `PG` value.
pub(super) fn output_header(
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

/// `MarkDuplicates` second read: write `OUTPUT`, return the metrics file.
pub fn run(opts: &Options, inputs: &[PathBuf], marks: Marks<'_>) -> Result<MetricsOutput> {
    let (copier, header, writer) = inputs::open(
        inputs,
        opts.reference.as_deref(),
        opts.threads,
        opts.assume_sort_order,
        |mut stream| {
            let order = SortOrder::of(&stream.header);
            check_order(order, opts.assume_sort_order)?;
            let (header, chained) = output_header(&stream.header, opts, &marks.pg_ids);
            let mut writer =
                BamWriter::create(&opts.output, opts.compression_level, opts.threads, &header)?;
            let mut copier = Copier::new(opts, &stream.header, chained, stream.cram, marks)?;
            let mut record = RawRecord::default();
            while stream.next(&mut record)? {
                copier.copy(&mut record, &mut writer)?;
            }
            Ok((copier, header, writer))
        },
    )?;
    let companions = Companions {
        index: opts.create_index,
        md5: opts.create_md5_file,
    };
    bam_writer::finish(writer, &header, companions)?;
    Ok(copier.metrics())
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

/// The index of the next representative entry, the representative and the set size.
#[derive(Clone, Copy)]
struct SetEntry {
    read_index: i64,
    representative: i32,
    size: i32,
}

/// `nextIndexIfNeeded`.
fn advance(
    queryname: bool,
    current: i64,
    index: i64,
    names: (Option<&[u8]>, &[u8]),
    list: Option<&mut Sorted<'_, u64>>,
) -> Result<i64> {
    let (last_name, name) = names;
    let need = index > current && (!queryname || last_name != Some(name));
    match list {
        Some(list) if need => next_index(list),
        _ => Ok(current),
    }
}

/// One tag edit: the value to set, or `None` to remove the tag.
type Edit = ([u8; 2], Option<TagValue>);

/// The second read's cursor: the index lists, the record index and the
/// names the query-name rules compare, plus the library rows it counts.
pub(super) struct Copier<'d> {
    opts: Options,
    /// The header the records are read with.
    header: sam::Header,
    queryname: bool,
    chained: Option<Chained>,
    tags: TagCache,
    /// Records come from a CRAM: every record is written with sorted tags.
    cram: bool,
    /// Reference bases for the `MD` and `NM` of CRAM records.
    contigs: Option<ContigCache>,
    groups: ReadGroups,
    rows: BTreeMap<String, Row>,
    marks: Marks<'d>,
    mi_tag: Option<[u8; 2]>,
    index: i64,
    next_duplicate: i64,
    next_optical: i64,
    set: SetEntry,
    duplicate_name: Option<Vec<u8>>,
    representative_name: Option<Vec<u8>>,
    counters: Counters,
}

impl<'d> Copier<'d> {
    /// The cursor before the first record. `header` is the input's (with
    /// `ASSUME_SORT_ORDER` applied), `chained` the second value of
    /// [`output_header`].
    pub(super) fn new(
        opts: &Options,
        header: &sam::Header,
        chained: Option<Chained>,
        cram: bool,
        mut marks: Marks<'d>,
    ) -> Result<Self> {
        let groups = ReadGroups::new(header);
        let mut rows = BTreeMap::new();
        for id in header.read_groups().keys() {
            let library = groups.library(Some(id.as_ref())).to_owned();
            rows.entry(library.clone()).or_insert_with(|| Row {
                library,
                ..Row::default()
            });
        }
        let contigs = match (cram, opts.reference.as_deref()) {
            (true, Some(fasta)) => Some(ContigCache::new(Reference::open(fasta, header)?)),
            _ => None,
        };
        let next_duplicate = next_index(&mut marks.duplicates)?;
        let next_optical = match marks.optical.as_mut() {
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
        let mi_tag = match (&opts.barcode_tag, &opts.molecular_identifier_tag) {
            (Some(_), Some(tag)) => Some(tag_name(tag)?),
            _ => None,
        };
        Ok(Self {
            opts: opts.clone(),
            header: header.clone(),
            queryname: SortOrder::of(header) == SortOrder::Queryname,
            chained,
            tags: TagCache::default(),
            cram,
            contigs,
            groups,
            rows,
            marks,
            mi_tag,
            index: 0,
            next_duplicate,
            next_optical,
            set,
            duplicate_name: None,
            representative_name: None,
            counters: Counters::default(),
        })
    }

    /// Decide the next record of the input and hand it to `emit`, unless
    /// the options remove it. Records must come in input order.
    pub(super) fn copy(&mut self, record: &mut RawRecord, emit: &mut dyn Emit) -> Result<()> {
        let index = self.index;
        self.counters.records_in += 1;
        self.tags.load(record)?;
        let name = record.name().to_vec();
        let input_flags = record.flags();
        let library = self
            .groups
            .library(string_tag(&self.tags, *b"RG")?)
            .to_owned();
        let row = self.rows.entry(library).or_insert_with_key(|library| Row {
            library: library.clone(),
            ..Row::default()
        });
        row.add_read(input_flags);
        self.next_duplicate = advance(
            self.queryname,
            self.next_duplicate,
            index,
            (self.duplicate_name.as_deref(), &name),
            Some(&mut self.marks.duplicates),
        )?;
        let duplicate = index == self.next_duplicate
            || (self.queryname
                && index > self.next_duplicate
                && self.duplicate_name.as_deref() == Some(name.as_slice()));
        let flags = if duplicate {
            input_flags | 0x400
        } else {
            input_flags & !0x400
        };
        record.set_flags(flags);
        if flags != input_flags {
            self.counters.changed += 1;
        }
        if duplicate {
            row.add_duplicate(flags);
        }
        self.next_optical = advance(
            self.queryname,
            self.next_optical,
            index,
            (self.duplicate_name.as_deref(), &name),
            self.marks.optical.as_mut(),
        )?;
        let optical = (self.queryname
            && index > self.next_optical
            && self.duplicate_name.as_deref() == Some(name.as_slice()))
            || index == self.next_optical;
        let mut edits: Vec<Edit> = Vec::new();
        if self.opts.clear_dt {
            edits.push((*b"DT", None));
        }
        if self.opts.tagging_policy != TaggingPolicy::DontTag && duplicate {
            if optical {
                edits.push((*b"DT", Some(TagValue::String(b"SQ".to_vec()))));
            } else if self.opts.tagging_policy == TaggingPolicy::All {
                edits.push((*b"DT", Some(TagValue::String(b"LB".to_vec()))));
            }
        }
        self.set_members(index, &name, flags, &mut edits)?;
        if let Some(tag) = self.mi_tag {
            edits.push((
                tag,
                Some(TagValue::String(self.molecular_identifier(record)?)),
            ));
        }
        if duplicate {
            self.duplicate_name = Some(name);
        }
        self.index += 1;
        if (self.opts.remove_duplicates && duplicate)
            || (self.opts.remove_sequencing_duplicates && optical)
        {
            return Ok(());
        }
        if let Some(chained) = self.chained.as_ref()
            && self.opts.add_pg_tag_to_reads
        {
            let old = string_tag(&self.tags, *b"PG")?.map(<[u8]>::to_vec);
            let new = chained.get(&old).map(|id| TagValue::String(id.clone()));
            edits.push((*b"PG", new));
        }
        self.edit_tags(record, edits)?;
        self.counters.records_out += 1;
        emit.emit(record)
    }

    /// Apply `edits` to the record's tags. A record with no edit keeps its
    /// bytes, unless it came from a CRAM.
    fn edit_tags(&mut self, record: &mut RawRecord, edits: Vec<Edit>) -> Result<()> {
        if !self.cram && edits.is_empty() {
            return Ok(());
        }
        let mut tags = record.tags().map_err(Error::Record)?;
        if let Some(contigs) = self.contigs.as_mut() {
            md_nm::restore(record, &mut tags, contigs)?;
        }
        for (tag, value) in edits {
            match value {
                Some(value) => tags.set(tag, value),
                None => tags.remove(tag),
            }
        }
        record.set_tags(&tags).map_err(Error::Record)
    }

    /// `TAG_DUPLICATE_SET_MEMBERS`: the `DI` and `DS` edits of a mapped
    /// primary record of a pair set.
    fn set_members(
        &mut self,
        index: i64,
        name: &[u8],
        flags: u16,
        edits: &mut Vec<Edit>,
    ) -> Result<()> {
        let Some(list) = self.marks.representatives.as_mut() else {
            return Ok(());
        };
        let need = index > self.set.read_index
            && (!self.queryname || self.representative_name.as_deref() != Some(name));
        if need && let Some(next) = list.next_item()? {
            self.set = SetEntry {
                read_index: i64::from(next.read_index),
                representative: next.representative_index,
                size: next.set_size,
            };
        }
        let in_set = index == self.set.read_index
            || (self.queryname
                && index > self.set.read_index
                && self.representative_name.as_deref() == Some(name));
        if in_set && flags & 0x900 == 0 && flags & 0x4 == 0 {
            edits.push((
                *b"DI",
                Some(TagValue::Int(i64::from(self.set.representative))),
            ));
            edits.push((*b"DS", Some(TagValue::Int(i64::from(self.set.size)))));
            self.representative_name = Some(name.to_vec());
        }
        Ok(())
    }

    /// The metrics of the records copied so far: the whole file after the
    /// last record.
    pub(super) fn metrics(&self) -> MetricsOutput {
        metrics(self.rows.clone(), &self.marks)
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

/// The second read as a chain transform: each batch of the replayed
/// source goes through [`Copier::copy`], and the kept records stay in the
/// batch for the next step.
impl Transform for Copier<'_> {
    fn apply(&mut self, batch: &mut Vec<RawRecord>) -> Result<()> {
        let mut kept = Vec::with_capacity(batch.len());
        for mut record in batch.drain(..) {
            self.copy(&mut record, &mut kept)?;
        }
        *batch = kept;
        Ok(())
    }

    fn counters(&self) -> Counters {
        self.counters
    }

    fn products(&mut self) -> Result<Products> {
        Ok(Products {
            output: self.metrics().into(),
            ..Products::default()
        })
    }
}
