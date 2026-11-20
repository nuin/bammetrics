//! The edited records between the edit and the file, when Picard's writer
//! does not take them as they come (waves W10 and later).
//!
//! * [`WriterOrder::Sort`]: a `SAMFileWriter` made with `presorted=false`
//!   for a header whose `SO` names a sort holds every record in htsjdk's
//!   `SortingCollection` and writes them in that order at close
//!   (`AddOrReplaceReadGroups` with a `SORT_ORDER` that differs from the
//!   input's, `RevertOriginalBaseQualitiesAndAddMateCigar`). The sort is
//!   [`crate::sam_sort`], the one `SortSam` uses: `MAX_RECORDS_IN_RAM`
//!   records per run, spilled to `TMP_DIR`.
//! * [`WriterOrder::Unknown`]: the same writer for `SO:unknown` has no
//!   comparator, so Java's `Arrays.sort` casts the second record to
//!   `Comparable` and Picard fails at close (`ClassCastException`); a run
//!   of 0 or 1 record completes. Nothing is written then.
//! * [`Reorder::regroup`]: a queryname sort first, then each group of
//!   records that share a read name edited as one (`SamPairUtil.
//!   SetMateInfoIterator`), then the writer's order.
//!
//! A [`Gate`] looks at the records as they are read and may decide that
//! the run writes nothing (the early exit of
//! `RevertOriginalBaseQualitiesAndAddMateCigar`).
//!
//! In a chain the same passes make the step a barrier that holds its
//! records ([`sorting_stage`]); a step whose writer is presorted for the
//! stream it reads streams it with the writer's order check
//! ([`presorted_stage`]).

use std::collections::VecDeque;
use std::path::PathBuf;

use noodles::sam;

use super::{Mode, Options, Out, check_stage};
use crate::bam_raw::RawRecord;
use crate::bam_writer::Presorted;
use crate::error::{Error, Result};
use crate::header_merge::SortOrder;
use crate::registry::Invocation;
use crate::sam_order::Order;
use crate::sam_sort::{Limits, RecordSorter, SortedRecords};
use crate::spill::{SortStats, TempDir};
use crate::stage::{
    BATCH_RECORDS, Barrier, Counters, Frame, Next, Products, Records, Sealed, Stage, Transform,
};

/// Picard's `MAX_RECORDS_IN_RAM` default.
pub const DEFAULT_RECORDS_IN_RAM: usize = 500_000;

/// Edit one group of records that share a read name, in queryname order.
pub type Regroup = fn(&mut [RawRecord]) -> Result<()>;

/// What the writer does with the records it is handed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WriterOrder {
    /// The records as they come; a writer whose header names a sort checks
    /// it (htsjdk `assertPresorted`).
    #[default]
    Presorted,
    /// htsjdk's `SortingCollection` in this order.
    Sort(Order),
    /// `SO:unknown` through a writer that sorts: fails on a second record.
    Unknown,
    /// The writer Picard makes for the output order `order` (`None`: the
    /// input's), with `presorted` true when `presorted_if_same` and the
    /// order is the input's. Resolved against the input header
    /// ([`WriterOrder::against`]).
    Header {
        order: Option<SortOrder>,
        presorted_if_same: bool,
    },
}

impl WriterOrder {
    /// The writer for an input whose header has `input`. A writer that
    /// does not sort and is not presorted (`unsorted`) takes the records
    /// as they come.
    pub fn against(self, input: SortOrder) -> Self {
        let Self::Header {
            order,
            presorted_if_same,
        } = self
        else {
            return self;
        };
        let out = order.unwrap_or(input);
        if presorted_if_same && out == input {
            return Self::Presorted;
        }
        match out {
            SortOrder::Unsorted => Self::Presorted,
            SortOrder::Queryname => Self::Sort(Order::Queryname),
            SortOrder::Coordinate => Self::Sort(Order::Coordinate),
            SortOrder::Duplicate => Self::Sort(Order::Duplicate),
            SortOrder::Unknown => Self::Unknown,
        }
    }
}

/// The passes between the edit and the file of a stand-alone run.
#[derive(Clone, Debug)]
pub struct Reorder {
    pub regroup: Option<Regroup>,
    pub writer: WriterOrder,
    /// `MAX_RECORDS_IN_RAM` per sort, and the threads that compress a run.
    pub limits: Limits,
    /// `TMP_DIR`: where sorted runs spill.
    pub tmp_dir: Option<PathBuf>,
}

impl Default for Reorder {
    fn default() -> Self {
        Self {
            regroup: None,
            writer: WriterOrder::Presorted,
            limits: Limits {
                records: DEFAULT_RECORDS_IN_RAM,
                bytes: usize::MAX,
                threads: 1,
            },
            tmp_dir: None,
        }
    }
}

/// `MAX_RECORDS_IN_RAM` (Picard's default 500,000) and the first
/// `TMP_DIR` of a run whose writer may sort. A count below 1 is 0 here,
/// and fails when a sort takes its first record or ends, as htsjdk's
/// `SortingCollection` fails when it is made (`maxRecordsInRam must be >
/// 0`): a run that sorts nothing completes, as in Picard.
pub fn sort_limits(inv: &Invocation) -> Result<(Limits, Option<PathBuf>)> {
    let records = match inv.parse_opt::<i64>("MAX_RECORDS_IN_RAM")? {
        None => DEFAULT_RECORDS_IN_RAM,
        Some(n) => usize::try_from(n).unwrap_or(0),
    };
    let limits = Limits {
        records,
        bytes: usize::MAX,
        threads: inv.threads,
    };
    Ok((limits, inv.all("TMP_DIR").first().map(PathBuf::from)))
}

/// The limits of a chain step's sort: `MAX_RECORDS_IN_RAM` when the step
/// gives it, and `share` of the chain's `--sort-memory`.
pub fn stage_limits(inv: &Invocation, share: usize) -> Result<Limits> {
    let records = match inv.parse_opt::<i64>("MAX_RECORDS_IN_RAM")? {
        None => usize::MAX,
        Some(n) => usize::try_from(n).ok().filter(|n| *n >= 1).ok_or_else(|| {
            Error::Usage(format!(
                "MAX_RECORDS_IN_RAM={n}: maxRecordsInRam must be > 0 for the sort of this step"
            ))
        })?,
    };
    Ok(Limits {
        records,
        bytes: inv.sort_memory.map_or(usize::MAX, |b| b / share.max(1)),
        threads: inv.threads,
    })
}

/// Decides from the first records whether a run writes anything.
pub trait Gate {
    /// The record before its edit: `Some(true)` skips the run, `Some(false)`
    /// runs it, `None` looks at the next record.
    fn look(&mut self, record: &RawRecord) -> Result<Option<bool>>;
    /// The input ended undecided: true skips the run.
    fn at_end(&mut self) -> bool;
}

/// The error of a writer that sorts `SO:unknown` records.
fn unknown_order(records: u64) -> Error {
    Error::Metric(format!(
        "SO:unknown through a writer that sorts: htsjdk holds no comparator for unknown, and \
         Picard fails at close with a ClassCastException on {records} records; nothing is \
         written"
    ))
}

/// The error of a sort whose `MAX_RECORDS_IN_RAM` is below 1.
fn no_records_in_ram() -> Error {
    Error::Usage("MAX_RECORDS_IN_RAM: maxRecordsInRam must be > 0 for this run's sort".to_owned())
}

/// The writer's end of a [`Pipeline`].
enum Last<'d> {
    Direct,
    Sort(Box<RecordSorter<'d>>),
    Unknown(u64),
}

impl<'d> Last<'d> {
    fn new(order: WriterOrder, header: &sam::Header, temp: &'d TempDir, limits: Limits) -> Self {
        match order {
            WriterOrder::Sort(order) => {
                Self::Sort(Box::new(RecordSorter::new(order, header, temp, limits)))
            }
            WriterOrder::Unknown => Self::Unknown(0),
            WriterOrder::Presorted | WriterOrder::Header { .. } => Self::Direct,
        }
    }
}

/// Records from the edit to the file of a stand-alone run.
pub(super) struct Pipeline<'d> {
    group: Option<(RecordSorter<'d>, Regroup)>,
    last: Last<'d>,
    /// The file, made by [`Pipeline::open`] or at the first record written:
    /// a run a [`Gate`] skips makes none.
    out: Option<Out>,
    opts: &'d Options,
    /// The output header.
    header: sam::Header,
    mode: Mode,
    /// `MAX_RECORDS_IN_RAM` below 1: a sort fails.
    no_room: bool,
}

impl<'d> Pipeline<'d> {
    /// `header` is the output header: the writer's sort reads its read
    /// groups (`duplicate` order). `input` is the order of the header
    /// read, which a [`WriterOrder::Header`] is resolved against.
    pub(super) fn new(
        opts: &'d Options,
        input: SortOrder,
        header: &sam::Header,
        temp: &'d TempDir,
        out: Option<Out>,
    ) -> Self {
        let reorder = &opts.reorder;
        let mut limits = reorder.limits;
        let no_room = limits.records == 0;
        limits.records = limits.records.max(1);
        let group = reorder.regroup.map(|regroup| {
            let sorter = RecordSorter::new(Order::Queryname, header, temp, limits);
            (sorter, regroup)
        });
        let writer = reorder.writer.against(input);
        Self {
            group,
            last: Last::new(writer, header, temp, limits),
            out,
            opts,
            header: header.clone(),
            mode: opts.mode,
            no_room,
        }
    }

    /// Make the file now, when it is not made yet.
    pub(super) fn open(&mut self) -> Result<&mut Out> {
        let out = match self.out.take() {
            Some(out) => out,
            None => Out::create(self.opts, &self.header)?,
        };
        Ok(self.out.insert(out))
    }

    fn write(&mut self, record: &mut RawRecord) -> Result<()> {
        let mode = self.mode;
        self.open()?.write(record, mode)
    }

    fn sorts(&self) -> bool {
        self.group.is_some() || matches!(self.last, Last::Sort(_))
    }

    /// Hand on `record`; its content afterwards is unspecified.
    pub(super) fn push(&mut self, record: &mut RawRecord) -> Result<()> {
        if self.no_room && self.sorts() {
            return Err(no_records_in_ram());
        }
        match &mut self.group {
            Some((sorter, _)) => sorter.push(std::mem::take(record)),
            None => self.feed_last(record),
        }
    }

    fn feed_last(&mut self, record: &mut RawRecord) -> Result<()> {
        match &mut self.last {
            Last::Direct => {}
            Last::Sort(sorter) => return sorter.push(std::mem::take(record)),
            Last::Unknown(n) => *n += 1,
        }
        self.write(record)
    }

    /// After the last record: every record in the file, the file not yet
    /// closed.
    pub(super) fn finish(mut self) -> Result<Out> {
        if self.no_room && self.sorts() {
            return Err(no_records_in_ram());
        }
        if let Some((sorter, regroup)) = self.group.take() {
            let (sorted, _) = sorter.finish()?;
            let mut grouped = Grouped::new(sorted, regroup);
            while let Some(mut record) = grouped.next_record()? {
                self.feed_last(&mut record)?;
            }
        }
        match std::mem::replace(&mut self.last, Last::Direct) {
            Last::Direct => {}
            Last::Unknown(n) if n >= 2 => return Err(unknown_order(n)),
            Last::Unknown(_) => {}
            Last::Sort(sorter) => {
                let (mut sorted, _) = (*sorter).finish()?;
                while let Some(mut record) = sorted.next_record()? {
                    self.write(&mut record)?;
                }
            }
        }
        match self.out.take() {
            Some(out) => Ok(out),
            None => Out::create(self.opts, &self.header),
        }
    }
}

fn add(a: SortStats, b: SortStats) -> SortStats {
    SortStats {
        spilled_bytes: a.spilled_bytes + b.spilled_bytes,
        merge_levels: a.merge_levels.max(b.merge_levels),
        ..b
    }
}

/// Queryname-sorted records, one read name at a time through a
/// [`Regroup`].
struct Grouped<'d> {
    sorted: SortedRecords<'d>,
    regroup: Regroup,
    /// The first record of the next group.
    next: Option<RawRecord>,
    ready: VecDeque<RawRecord>,
}

impl<'d> Grouped<'d> {
    fn new(sorted: SortedRecords<'d>, regroup: Regroup) -> Self {
        Self {
            sorted,
            regroup,
            next: None,
            ready: VecDeque::new(),
        }
    }

    fn next_record(&mut self) -> Result<Option<RawRecord>> {
        if let Some(record) = self.ready.pop_front() {
            return Ok(Some(record));
        }
        let first = match self.next.take() {
            Some(record) => record,
            None => match self.sorted.next_record()? {
                Some(record) => record,
                None => return Ok(None),
            },
        };
        let mut group = vec![first];
        while let Some(record) = self.sorted.next_record()? {
            if record.name() == group[0].name() {
                group.push(record);
            } else {
                self.next = Some(record);
                break;
            }
        }
        (self.regroup)(&mut group)?;
        self.ready.extend(group);
        Ok(self.ready.pop_front())
    }
}

impl Records for Grouped<'_> {
    fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool> {
        out.clear();
        while out.len() < BATCH_RECORDS {
            match self.next_record()? {
                Some(record) => out.push(record),
                None => break,
            }
        }
        Ok(!out.is_empty())
    }
}

/// A record-editing step whose stand-alone writer sorts: the edit on each
/// batch, every record held, then handed on in order under the edited
/// header ([`Next::Records`]). `transform` is the step's edit
/// ([`super::stage`]). A [`WriterOrder::Unknown`] is a usage error: the
/// tool on its own fails on any stream of two records.
pub fn sorting_stage<'t>(
    inv: &Invocation,
    frame: &Frame,
    temp: &'t TempDir,
    reorder: &Reorder,
    transform: Box<dyn Transform + 't>,
) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let header = transform.header(&frame.header)?;
    let group = reorder.regroup.map(|regroup| {
        let sorter = RecordSorter::new(Order::Queryname, &header, temp, reorder.limits);
        (sorter, regroup)
    });
    let last = match reorder.writer.against(SortOrder::of(&frame.header)) {
        WriterOrder::Sort(order) => Some(RecordSorter::new(order, &header, temp, reorder.limits)),
        WriterOrder::Presorted | WriterOrder::Header { .. } => None,
        WriterOrder::Unknown => {
            return Err(Error::Usage(
                "SO:unknown in a chain: the tool on its own sorts the records with no \
                 comparator and fails; give SORT_ORDER"
                    .to_owned(),
            ));
        }
    };
    Ok(Stage::Barrier(Box::new(Sorting {
        transform,
        header,
        group,
        last,
        batch: Vec::new(),
    })))
}

/// A record-editing step whose stand-alone writer is presorted for this
/// stream (`AddOrReplaceReadGroups` with a `SORT_ORDER` that is the
/// input's): the edit on each batch, then htsjdk's `assertPresorted` on
/// every record it hands on, in the order its output header names. The
/// records go on as they come, so the step stays in its segment.
pub fn presorted_stage<'t>(
    frame: &Frame,
    transform: Box<dyn Transform + 't>,
    tool: &str,
) -> Result<Stage<'t>> {
    let header = transform.header(&frame.header)?;
    Ok(Stage::Transform(Box::new(Checked {
        transform,
        order: Presorted::of(&header),
        stream: PathBuf::from(format!("the stream of {tool}")),
    })))
}

/// The transform of [`presorted_stage`].
struct Checked<'t> {
    transform: Box<dyn Transform + 't>,
    order: Presorted,
    /// What the order error names in place of a file.
    stream: PathBuf,
}

impl Checked<'_> {
    fn check(&mut self, records: &[RawRecord]) -> Result<()> {
        records
            .iter()
            .try_for_each(|record| self.order.check(record, &self.stream))
    }
}

impl Transform for Checked<'_> {
    fn header(&self, up: &sam::Header) -> Result<sam::Header> {
        self.transform.header(up)
    }

    fn apply(&mut self, batch: &mut Vec<RawRecord>) -> Result<()> {
        self.transform.apply(batch)?;
        self.check(batch)
    }

    fn finish(&mut self, out: &mut Vec<RawRecord>) -> Result<()> {
        self.transform.finish(out)?;
        self.check(out)
    }

    fn counters(&self) -> Counters {
        self.transform.counters()
    }

    fn products(&mut self) -> Result<Products> {
        self.transform.products()
    }
}

/// The barrier of [`sorting_stage`].
struct Sorting<'t> {
    transform: Box<dyn Transform + 't>,
    header: sam::Header,
    group: Option<(RecordSorter<'t>, Regroup)>,
    last: Option<RecordSorter<'t>>,
    batch: Vec<RawRecord>,
}

impl Sorting<'_> {
    fn take(&mut self, record: RawRecord) -> Result<()> {
        match (&mut self.group, &mut self.last) {
            (Some((sorter, _)), _) | (None, Some(sorter)) => sorter.push(record),
            (None, None) => Err(Error::Metric(
                "a sorting step without a sort: an internal error".to_owned(),
            )),
        }
    }
}

impl<'t> Barrier<'t> for Sorting<'t> {
    fn push(&mut self, batch: &[RawRecord]) -> Result<()> {
        let mut edited = std::mem::take(&mut self.batch);
        edited.clear();
        edited.extend_from_slice(batch);
        self.transform.apply(&mut edited)?;
        for record in edited.drain(..) {
            self.take(record)?;
        }
        self.batch = edited;
        Ok(())
    }

    fn seal(mut self: Box<Self>) -> Result<Sealed<'t>> {
        let mut tail = Vec::new();
        self.transform.finish(&mut tail)?;
        for record in tail {
            self.take(record)?;
        }
        let counters = self.transform.counters();
        let Sorting {
            header,
            group,
            last,
            ..
        } = *self;
        let mut stats = SortStats::default();
        let next: Box<dyn Records + 't> = match (group, last) {
            (Some((sorter, regroup)), Some(mut last)) => {
                let (sorted, first) = sorter.finish()?;
                let mut grouped = Grouped::new(sorted, regroup);
                while let Some(record) = grouped.next_record()? {
                    last.push(record)?;
                }
                let (sorted, second) = last.finish()?;
                stats = add(add(stats, first), second);
                Box::new(sorted)
            }
            (Some((sorter, regroup)), None) => {
                let (sorted, first) = sorter.finish()?;
                stats = add(stats, first);
                Box::new(Grouped::new(sorted, regroup))
            }
            (None, Some(last)) => {
                let (sorted, second) = last.finish()?;
                stats = add(stats, second);
                Box::new(sorted)
            }
            (None, None) => {
                return Err(Error::Metric(
                    "a sorting step without a sort: an internal error".to_owned(),
                ));
            }
        };
        let mut notes = Vec::new();
        if stats.spilled_bytes > 0 {
            notes.push(format!(
                "spilled {} bytes in sorted runs, {} merge level{}",
                stats.spilled_bytes,
                stats.merge_levels,
                if stats.merge_levels == 1 { "" } else { "s" }
            ));
        }
        Ok(Sealed {
            header,
            next: Next::Records(next),
            products: Products {
                counters: Some(counters),
                notes,
                ..Products::default()
            },
        })
    }
}
