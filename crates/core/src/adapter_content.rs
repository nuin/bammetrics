//! `CollectAdapterContentMetrics` (docs/06 section 7): adapter content as
//! FastQC's Adapter Content module reports it, read from a BAM. Not a
//! Picard tool.
//!
//! One pass. A record counts when it is primary (not secondary, not
//! supplementary); unmapped records count, since adapter dimers rarely
//! map. `PF_READS_ONLY` drops QC-failed records and `INCLUDE_DUPLICATES`
//! false drops duplicates. The bases are taken in sequencing order through
//! `read_bases.rs`, as `CollectTelomereMetrics` takes them: a record with
//! flag 0x10 is reverse complemented, as `samtools fastq` does.
//! Hard-clipped bases are not in the record and cannot be restored.
//!
//! For each adapter the search key is its first `ADAPTER_KMER` bases
//! (default 12, FastQC's length). A read carries the adapter when the key
//! occurs in it exactly, as FastQC's `String.indexOf`; the start is the
//! first occurrence. The by-cycle table gives, per cycle c, the percent of
//! the read end's reads whose start is at or before c: FastQC's cumulative
//! curve. Reads are split by read end: `FIRST` (flags 0x1 and 0x40 without
//! 0x80), `SECOND` (0x1 and 0x80 without 0x40), else `UNPAIRED`.
//!
//! `PAIR_READ_THROUGH` needs no sequence: the fraction of proper pairs
//! (flags 0x1 and 0x2, both mates mapped, TLEN not 0) whose |TLEN| is
//! below the record's read length, so the read runs into the adapter.
//!
//! Search: every key is a 2-bit code of at most 32 bases. Each read is
//! rolled base by base into a code, a 65,536-bit filter rejects most
//! positions, and a binary search over the sorted keys finds the adapters
//! a code belongs to. Cost is linear in the read length.

use std::collections::{BTreeMap, HashSet};

use noodles::{bam, sam};

use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::levels::{Labels, Level, Levels, ReadGroupInfo, parse_levels, read_groups};
use crate::metric::{
    HistogramSection, HistogramSeries, MetricRow, MetricsOutput, ToolOutput, Value,
};
use crate::metrics_file::KeyClass;
use crate::pass::{self, Accumulator};
use crate::read_bases;
use crate::registry::{Invocation, Opt, Tool};

/// The metrics class string. Not a Picard class: readers key on the last
/// name.
pub const CLASS: &str = "bammetrics.qc.AdapterContentMetrics";

/// FastQC's k-mer length (`Configuration/adapter_list.txt`: "we've gone
/// with 12bp fragments").
pub const DEFAULT_KMER: usize = 12;
/// The longest key a 64-bit code holds.
pub const MAX_KMER: usize = 32;
/// Adapters one run searches: one bit each in a 64-bit mask.
pub const MAX_ADAPTERS: usize = 64;
/// The by-cycle table stops at this cycle. A start after it still counts
/// in `READS_WITH_ADAPTER` and `MEAN_ADAPTER_START`.
pub const MAX_TABLE_CYCLES: usize = 1_000;
/// The cells one run may write: 12 per row plus one per table column and
/// cycle. One unit holds at most 3 x 64 x (12 + 1,000) = 194,304; a level
/// with many units can ask for more, which is refused.
pub const MAX_OUTPUT_CELLS: usize = 5_000_000;
/// The histogram key label. `CYCLE` would take the glossary rows of
/// `MeanQualityByCycle` under `CollectMultipleMetrics`.
pub const KEY_LABEL: &str = "READ_CYCLE";

/// FastQC 0.13.0 `Configuration/adapter_list.txt`, in file order.
pub const DEFAULT_ADAPTERS: [(&str, &str); 6] = [
    ("Illumina Universal Adapter", "AGATCGGAAGAG"),
    ("Illumina Small RNA 3' Adapter", "TGGAATTCTCGG"),
    ("Illumina Small RNA 5' Adapter", "GATCGTCGGACT"),
    ("Nextera Transposase Sequence", "CTGTCTCTTATA"),
    ("PolyA", "AAAAAAAAAAAA"),
    ("PolyG", "GGGGGGGGGGGG"),
];

/// The read end of a record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReadEnd {
    First,
    Second,
    Unpaired,
}

impl ReadEnd {
    pub const ALL: [ReadEnd; 3] = [ReadEnd::First, ReadEnd::Second, ReadEnd::Unpaired];

    pub fn name(self) -> &'static str {
        match self {
            ReadEnd::First => "FIRST",
            ReadEnd::Second => "SECOND",
            ReadEnd::Unpaired => "UNPAIRED",
        }
    }

    fn of(flags: sam::alignment::record::Flags) -> Self {
        let (first, last) = (flags.is_first_segment(), flags.is_last_segment());
        match (flags.is_segmented(), first, last) {
            (true, true, false) => ReadEnd::First,
            (true, false, true) => ReadEnd::Second,
            _ => ReadEnd::Unpaired,
        }
    }

    fn index(self) -> usize {
        self as usize
    }
}

/// One adapter: its name and the key searched.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Adapter {
    pub name: String,
    pub kmer: String,
}

/// An `ADAPTER_SEQUENCE` value: `NAME=SEQ`, or `SEQ` named by itself. The
/// key is the first `k` bases, upper-cased, which must be A, C, G or T.
pub fn parse_adapter(text: &str, k: usize) -> Result<Adapter> {
    let (name, sequence) = text.rsplit_once('=').unwrap_or((text, text));
    let bad = |why: String| Error::Usage(format!("ADAPTER_SEQUENCE={text}: {why}"));
    if name.is_empty() {
        return Err(bad("the name is empty".to_owned()));
    }
    if name.contains(['\t', '\n', '\r']) {
        return Err(bad("the name holds a tab or a line end".to_owned()));
    }
    let kmer = sequence
        .chars()
        .take(k)
        .collect::<String>()
        .to_ascii_uppercase();
    if kmer.chars().count() < k {
        return Err(bad(format!(
            "{} bases, ADAPTER_KMER is {k}",
            sequence.chars().count()
        )));
    }
    if !kmer.bytes().all(|b| base_code(b).is_some()) {
        return Err(bad(format!(
            "the first {k} bases must be A, C, G or T, got {kmer}"
        )));
    }
    Ok(Adapter {
        name: name.to_owned(),
        kmer,
    })
}

/// The 2-bit code of an upper-case base, `None` for any other byte.
fn base_code(base: u8) -> Option<u64> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// The options of one run.
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub adapters: Vec<Adapter>,
    pub kmer: usize,
    pub pf_reads_only: bool,
    pub include_duplicates: bool,
    pub levels: Vec<Level>,
}

impl Options {
    pub(crate) fn from_invocation(inv: &Invocation) -> Result<Self> {
        let value: i64 = inv.parse("ADAPTER_KMER", DEFAULT_KMER as i64)?;
        let kmer = usize::try_from(value)
            .ok()
            .filter(|k| (1..=MAX_KMER).contains(k))
            .ok_or_else(|| Error::Usage(format!("ADAPTER_KMER={value}: 1 to {MAX_KMER}")))?;
        let given = inv.all("ADAPTER_SEQUENCE");
        if given.len() > MAX_ADAPTERS {
            return Err(Error::Usage(format!(
                "ADAPTER_SEQUENCE: {} values, at most {MAX_ADAPTERS}",
                given.len()
            )));
        }
        let adapters = if given.is_empty() {
            if kmer > DEFAULT_KMER {
                return Err(Error::Usage(format!(
                    "ADAPTER_KMER={kmer}: FastQC's adapters have {DEFAULT_KMER} bases; give \
                     ADAPTER_SEQUENCE values of at least {kmer} bases"
                )));
            }
            DEFAULT_ADAPTERS
                .iter()
                .map(|(name, seq)| parse_adapter(&format!("{name}={seq}"), kmer))
                .collect::<Result<Vec<_>>>()?
        } else {
            given
                .iter()
                .map(|text| parse_adapter(text, kmer))
                .collect::<Result<Vec<_>>>()?
        };
        let mut names = HashSet::new();
        if let Some(twice) = adapters.iter().find(|a| !names.insert(a.name.as_str())) {
            return Err(Error::Usage(format!(
                "ADAPTER_SEQUENCE: the name {} is given twice",
                twice.name
            )));
        }
        Ok(Self {
            adapters,
            kmer,
            pf_reads_only: inv.flag("PF_READS_ONLY", false)?,
            include_duplicates: inv.flag("INCLUDE_DUPLICATES", true)?,
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
        })
    }
}

/// Bits of the code filter: one per value of the low 16 bits of a code.
const FILTER_WORDS: usize = 1 << 10;

/// Exact search of every key in a read, one pass over its bases.
#[derive(Clone, Debug)]
pub struct Matcher {
    k: usize,
    mask: u64,
    filter: Vec<u64>,
    /// Sorted codes, each with the bit of every adapter whose key it is.
    codes: Vec<(u64, u64)>,
}

impl Matcher {
    pub fn new(adapters: &[Adapter], k: usize) -> Self {
        let mask = if k >= MAX_KMER {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let mut by_code: BTreeMap<u64, u64> = BTreeMap::new();
        for (index, adapter) in adapters.iter().enumerate() {
            let code = adapter
                .kmer
                .bytes()
                .filter_map(base_code)
                .fold(0, |code, b| (code << 2) | b);
            *by_code.entry(code).or_default() |= 1 << index;
        }
        let mut filter = vec![0u64; FILTER_WORDS];
        for &code in by_code.keys() {
            let low = (code & 0xFFFF) as usize;
            filter[low / 64] |= 1 << (low % 64);
        }
        Self {
            k,
            mask,
            filter,
            codes: by_code.into_iter().collect(),
        }
    }

    fn adapters_of(&self, code: u64) -> u64 {
        self.codes
            .binary_search_by_key(&code, |&(c, _)| c)
            .map_or(0, |i| self.codes[i].1)
    }

    /// The 0-based start of the first occurrence of each key in `bases`,
    /// as `(adapter index, start)` pairs in `out`.
    pub fn first_starts(&self, bases: &[u8], out: &mut Vec<(usize, usize)>) {
        let codes = bases
            .iter()
            .map(|&b| base_code(b).map_or(NO_BASE, |c| c as u8));
        self.scan(codes, out);
    }

    /// [`Matcher::first_starts`] over base codes in sequencing order,
    /// [`NO_BASE`] for a base other than A, C, G or T.
    fn scan(&self, codes: impl Iterator<Item = u8>, out: &mut Vec<(usize, usize)>) {
        out.clear();
        let mut scan = Scan::new(self);
        for base in codes {
            if scan.step(base, out) {
                return;
            }
        }
    }

    /// [`Matcher::first_starts`] over a BAM record's packed bases, taken
    /// in sequencing order by [`read_bases::Bases::each`]. No decoded copy is
    /// made.
    fn scan_record(
        &self,
        packed: &[u8],
        length: usize,
        reverse: bool,
        out: &mut Vec<(usize, usize)>,
    ) {
        out.clear();
        let mut search = Search {
            scan: Scan::new(self),
            out,
        };
        SEQUENCED.each(packed, length, reverse, &mut search);
    }
}

/// The rolling state of one read's search.
struct Scan<'m> {
    matcher: &'m Matcher,
    /// Every adapter's bit.
    all: u64,
    code: u64,
    /// Bases since the last base other than A, C, G or T.
    run: usize,
    found: u64,
    /// The 0-based position of the next base.
    position: usize,
}

impl<'m> Scan<'m> {
    fn new(matcher: &'m Matcher) -> Self {
        Self {
            matcher,
            all: matcher.codes.iter().fold(0, |m, &(_, bits)| m | bits),
            code: 0,
            run: 0,
            found: 0,
            position: 0,
        }
    }

    /// Take the next base; true once every key is found. A base costs one
    /// branch, taken only when the code passes the filter.
    #[inline(always)]
    fn step(&mut self, base: u8, out: &mut Vec<(usize, usize)>) -> bool {
        let m = self.matcher;
        self.code = ((self.code << 2) | u64::from(base & 3)) & m.mask;
        self.run = if base == NO_BASE { 0 } else { self.run + 1 };
        self.position += 1;
        let low = (self.code & 0xFFFF) as usize;
        if self.run < m.k || m.filter[low / 64] & (1 << (low % 64)) == 0 {
            return false;
        }
        let mut new = m.adapters_of(self.code) & !self.found;
        self.found |= new;
        while new != 0 {
            out.push((new.trailing_zeros() as usize, self.position - m.k));
            new &= new - 1;
        }
        self.found == self.all
    }
}

/// A [`Scan`] and the hits it writes, as [`read_bases::Bases::each`] takes
/// them.
struct Search<'m, 'o> {
    scan: Scan<'m>,
    out: &'o mut Vec<(usize, usize)>,
}

impl read_bases::Take for Search<'_, '_> {
    #[inline(always)]
    fn take(&mut self, base: u8) -> bool {
        self.scan.step(base, self.out)
    }
}

/// A base other than A, C, G or T in [`Matcher::scan`].
const NO_BASE: u8 = 4;

/// The 2-bit code of each 4-bit BAM base (`=ACMGRSVTWYHKDBN`): A 1, C 2,
/// G 4 and T 8 have one, every other code is [`NO_BASE`].
const NIBBLE_CODE: [u8; 16] = {
    let mut table = [NO_BASE; 16];
    table[1] = 0;
    table[2] = 1;
    table[4] = 2;
    table[8] = 3;
    table
};

/// [`NIBBLE_CODE`] of each base in sequencing order.
static SEQUENCED: read_bases::Bases = read_bases::Bases::new(NIBBLE_CODE);

/// The counts of one read end in one unit.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct End {
    pub reads: i64,
    /// The longest read, in bases.
    pub longest: usize,
    /// Per adapter: reads by 0-based start of the first occurrence. Empty
    /// until the end's first read, so a unit without reads holds no map.
    pub starts: Vec<BTreeMap<usize, i64>>,
    pub proper_pairs: i64,
    pub read_through: i64,
}

/// The counts of one accumulation unit.
#[derive(Clone, Debug, PartialEq)]
pub struct Unit {
    pub ends: [End; 3],
}

impl Unit {
    fn new() -> Self {
        Self {
            ends: [End::default(), End::default(), End::default()],
        }
    }
}

/// One row: an adapter in one read end of one unit.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub adapter: String,
    pub kmer: String,
    pub read_end: ReadEnd,
    pub reads: i64,
    pub reads_with_adapter: i64,
    pub pct_reads_with_adapter: f64,
    pub mean_adapter_start: f64,
    pub pct_at_last_cycle: f64,
    pub pair_read_through: f64,
    pub labels: Labels,
}

impl MetricRow for Row {
    const CLASS: &'static str = CLASS;
    const COLUMNS: &'static [&'static str] = &[
        "ADAPTER",
        "KMER",
        "READ_END",
        "READS",
        "READS_WITH_ADAPTER",
        "PCT_READS_WITH_ADAPTER",
        "MEAN_ADAPTER_START",
        "PCT_AT_LAST_CYCLE",
        "PAIR_READ_THROUGH",
        "SAMPLE",
        "LIBRARY",
        "READ_GROUP",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.adapter.clone().into(),
            self.kmer.clone().into(),
            self.read_end.name().into(),
            self.reads.into(),
            self.reads_with_adapter.into(),
            self.pct_reads_with_adapter.into(),
            self.mean_adapter_start.into(),
            self.pct_at_last_cycle.into(),
            self.pair_read_through.into(),
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]
    }
}

/// `100 * count / reads`, NaN without a read.
fn percent(count: i64, reads: i64) -> f64 {
    if reads == 0 {
        f64::NAN
    } else {
        (count as f64 * 100.0) / reads as f64
    }
}

/// `part / whole`, NaN when `whole` is 0.
fn fraction(part: i64, whole: i64) -> f64 {
    if whole == 0 {
        f64::NAN
    } else {
        part as f64 / whole as f64
    }
}

/// The row of one adapter in one read end; `last_cycle` is the table's
/// last cycle.
fn row(
    adapter: &Adapter,
    read_end: ReadEnd,
    end: &End,
    starts: &BTreeMap<usize, i64>,
    last_cycle: usize,
) -> Row {
    let with: i64 = starts.values().sum();
    let start_sum: f64 = starts
        .iter()
        .map(|(&s, &n)| (s + 1) as f64 * n as f64)
        .sum();
    let by_last: i64 = starts.range(..last_cycle).map(|(_, &n)| n).sum();
    Row {
        adapter: adapter.name.clone(),
        kmer: adapter.kmer.clone(),
        read_end,
        reads: end.reads,
        reads_with_adapter: with,
        pct_reads_with_adapter: percent(with, end.reads),
        mean_adapter_start: if with == 0 {
            f64::NAN
        } else {
            start_sum / with as f64
        },
        pct_at_last_cycle: percent(by_last, end.reads),
        pair_read_through: fraction(end.read_through, end.proper_pairs),
        labels: Labels::default(),
    }
}

/// The cumulative percent at every cycle 1 to `last_cycle`.
fn cumulative(starts: &BTreeMap<usize, i64>, reads: i64, last_cycle: usize) -> Histogram {
    let mut histogram = Histogram::new();
    let mut counted = 0;
    let mut next = starts.iter().peekable();
    for cycle in 1..=last_cycle {
        while let Some((_, &n)) = next.next_if(|&(&s, _)| s < cycle) {
            counted += n;
        }
        histogram.increment_by(cycle as i64, percent(counted, reads));
    }
    histogram
}

/// The usage error of a run whose output would exceed [`MAX_OUTPUT_CELLS`].
fn check_cells(cells: usize, units: usize, adapters: usize) -> Result<()> {
    if cells <= MAX_OUTPUT_CELLS {
        return Ok(());
    }
    Err(Error::Usage(format!(
        "{units} units and {adapters} adapters make {cells} cells, at most {MAX_OUTPUT_CELLS}; \
         choose a coarser METRIC_ACCUMULATION_LEVEL or fewer ADAPTER_SEQUENCE values"
    )))
}

/// The units the header alone gives: one for `ALL_READS`, one per distinct
/// sample, library or platform unit (a missing field is one more) at each
/// other level. The `unknown` unit of a record is not counted.
fn header_units(header: &sam::Header, levels: &[Level]) -> usize {
    let groups = read_groups(header);
    let distinct = |field: fn(&ReadGroupInfo) -> Option<&str>| {
        groups.iter().map(field).collect::<HashSet<_>>().len()
    };
    levels
        .iter()
        .map(|level| match level {
            Level::AllReads => 1,
            Level::Sample => distinct(|g| g.sample.as_deref()),
            Level::Library => distinct(|g| g.library.as_deref()),
            Level::ReadGroup => distinct(|g| g.platform_unit.as_deref()),
        })
        .sum()
}

/// What a run returns: the rows and the by-cycle table.
#[derive(Clone, Debug)]
pub struct Report {
    pub rows: Vec<Row>,
    pub table: Option<HistogramSection>,
}

/// `CollectAdapterContentMetrics` as a pass [`Accumulator`].
pub struct AdapterContent {
    adapters: Vec<Adapter>,
    matcher: Matcher,
    pf_reads_only: bool,
    include_duplicates: bool,
    levels: Levels<Unit>,
    /// Read ends with at least one counted record.
    seen: [bool; 3],
    hits: Vec<(usize, usize)>,
}

impl AdapterContent {
    pub fn new(opts: &Options, header: &sam::Header) -> Self {
        Self {
            matcher: Matcher::new(&opts.adapters, opts.kmer),
            adapters: opts.adapters.clone(),
            pf_reads_only: opts.pf_reads_only,
            include_duplicates: opts.include_duplicates,
            levels: Levels::new(header, &opts.levels, Unit::new()),
            seen: [false; 3],
            hits: Vec::new(),
        }
    }

    fn counted(&self, flags: sam::alignment::record::Flags) -> bool {
        !(flags.is_secondary()
            || flags.is_supplementary()
            || (self.pf_reads_only && flags.is_qc_fail())
            || (!self.include_duplicates && flags.is_duplicate()))
    }
}

/// True for a proper pair whose insert is shorter than the read: `Some`
/// for a record that counts in `PAIR_READ_THROUGH`.
fn read_through(record: &bam::Record, length: usize) -> Option<bool> {
    let flags = record.flags();
    let proper = flags.is_segmented()
        && flags.is_properly_segmented()
        && !flags.is_unmapped()
        && !flags.is_mate_unmapped();
    let tlen = i64::from(record.template_length()).abs();
    (proper && tlen != 0).then_some(tlen < length as i64)
}

impl Accumulator for AdapterContent {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if !self.counted(flags) {
            return Ok(());
        }
        let read_end = ReadEnd::of(flags);
        self.seen[read_end.index()] = true;
        let sequence = record.sequence();
        let length = sequence.len();
        let reverse = flags.is_reverse_complemented();
        self.matcher
            .scan_record(sequence.as_bytes(), length, reverse, &mut self.hits);
        let through = read_through(record, length);
        let hits = &self.hits;
        let adapters = self.adapters.len();
        self.levels.route(record, |unit| {
            let end = &mut unit.ends[read_end.index()];
            if end.starts.is_empty() {
                end.starts.resize(adapters, BTreeMap::new());
            }
            end.reads += 1;
            end.longest = end.longest.max(length);
            for &(adapter, start) in hits {
                *end.starts[adapter].entry(start).or_default() += 1;
            }
            if let Some(through) = through {
                end.proper_pairs += 1;
                end.read_through += i64::from(through);
            }
        });
        Ok(())
    }

    fn finish(self) -> Result<Report> {
        let k = self.matcher.k;
        let buckets = self.levels.into_buckets();
        let longest = buckets
            .iter()
            .flat_map(|b| b.state.ends.iter().map(|e| e.longest))
            .max()
            .unwrap_or(0);
        let last_cycle = (longest + 1).saturating_sub(k).min(MAX_TABLE_CYCLES);
        let mut ends: Vec<ReadEnd> = ReadEnd::ALL
            .into_iter()
            .filter(|e| self.seen[e.index()])
            .collect();
        if ends.is_empty() {
            ends.push(ReadEnd::Unpaired);
        }
        let row_count = buckets.len() * ends.len() * self.adapters.len();
        let read_ends = buckets
            .iter()
            .flat_map(|b| &b.state.ends)
            .filter(|e| e.reads > 0)
            .count();
        let cells = row_count * Row::COLUMNS.len() + read_ends * self.adapters.len() * last_cycle;
        check_cells(cells, buckets.len(), self.adapters.len())?;
        let no_start = BTreeMap::new();
        let mut rows = Vec::with_capacity(row_count);
        let mut series = Vec::new();
        for bucket in &buckets {
            for &read_end in &ends {
                let end = &bucket.state.ends[read_end.index()];
                for (index, adapter) in self.adapters.iter().enumerate() {
                    let starts = end.starts.get(index).unwrap_or(&no_start);
                    rows.push(Row {
                        labels: bucket.labels.clone(),
                        ..row(adapter, read_end, end, starts, last_cycle)
                    });
                    // A read end without a read in its unit gets its row and
                    // no column: every cell would be `?`.
                    if end.reads > 0 {
                        series.push(HistogramSeries {
                            label: format!("{}.{}.{}", bucket.key, read_end.name(), adapter.name),
                            histogram: cumulative(starts, end.reads, last_cycle),
                        });
                    }
                }
            }
        }
        let table = (last_cycle > 0).then_some(HistogramSection {
            key_label: KEY_LABEL,
            key_class: KeyClass::Integer,
            series,
        });
        Ok(Report { rows, table })
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectAdapterContentMetrics",
    summary: "one row of the 12 AdapterContentMetrics columns per adapter, read end and \
              accumulation unit (reads, reads with the adapter's first ADAPTER_KMER bases, \
              their mean start, PAIR_READ_THROUGH from TLEN), and a by-cycle histogram of \
              the cumulative percent of reads with the adapter, as FastQC's Adapter Content \
              module. Not a Picard tool (docs/06).",
    options: &[
        Opt::text(
            &["ADAPTER_SEQUENCE"],
            "An adapter as NAME=SEQ or SEQ; repeat for each. Any value replaces FastQC's six \
             (Illumina Universal, Small RNA 3' and 5', Nextera, PolyA, PolyG).",
        )
        .repeated(),
        Opt::int(
            &["ADAPTER_KMER"],
            "Bases of each adapter searched, exact match, 1 to 32 (default 12, as FastQC).",
        ),
        Opt::bool(
            &["PF_READS_ONLY"],
            "Count only reads that pass the vendor filter (default false).",
        ),
        Opt::bool(
            &["INCLUDE_DUPLICATES"],
            "Count reads flagged as duplicates (default true).",
        ),
        crate::registry::LEVEL_OPTION,
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn pass::Collector>> {
    let opts = Options::from_invocation(inv)?;
    // Every unit gets at least one row per adapter: refuse before the first
    // record when the header alone is over the bound.
    let units = header_units(header, &opts.levels);
    let adapters = opts.adapters.len();
    check_cells(units * adapters * Row::COLUMNS.len(), units, adapters)?;
    Ok(pass::boxed(AdapterContent::new(&opts, header), |report| {
        Ok(MetricsOutput::from_rows(&report.rows, report.table).into())
    }))
}

#[cfg(test)]
mod tests;
