//! `CollectTelomereMetrics` (docs/telomere.md): telomeric read content and
//! a telomere length estimate as TelSeq defines them (Ding et al. 2014,
//! Nucleic Acids Research 42:e75). Not a Picard tool.
//!
//! One pass. A record counts when it is primary (not secondary, not
//! supplementary) and not QC-failed; unmapped records count, since most
//! telomeric reads are unmapped. Duplicates are left out unless
//! `COUNT_DUPLICATES=true`. A record with no bases (`SEQ` `*`) is not
//! counted. Each counted read is restored to the strand it
//! was sequenced on (reverse complemented when flag 0x10 is set), then
//! [`count_repeats`] counts the non-overlapping occurrences of `REPEAT` and
//! of its reverse complement, scanning left to right, and keeps the larger
//! count, as TelSeq's `countMotif` does. A read is telomeric at
//! `TELOMERE_REPEAT_THRESHOLD` (TelSeq's k) or more. One row per
//! accumulation unit (`levels.rs`), empty units included.

use noodles::{bam, sam};

use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::levels::{Labels, Level, Levels, parse_levels};
use crate::metric::{
    HistogramSection, HistogramSeries, MetricRow, MetricsOutput, ToolOutput, Value,
};
use crate::metrics_file::KeyClass;
use crate::pass::{self, Accumulator};
use crate::registry::{Invocation, Opt, Tool};

/// The metrics class string. Not a Picard class: readers key on the last
/// name.
pub const CLASS: &str = "bammetrics.qc.TelomereMetrics";

/// TelSeq's motif and its k.
pub const DEFAULT_REPEAT: &str = "TTAGGG";
pub const DEFAULT_THRESHOLD: i64 = 7;
/// TelSeq's `TELOMERE_ENDS`: 23 chromosome pairs, two ends each.
pub const DEFAULT_TELOMERE_ENDS: i64 = 46;
/// TelSeq's `GENOME_LENGTH_AT_TEL_GC`: the bases of hg19 in regions of 48
/// to 52 percent GC.
pub const DEFAULT_GENOME_LENGTH_AT_GC: i64 = 332_720_800;
/// TelSeq's `LENGTH_UNIT`: the estimate is in kilobases.
const LENGTH_UNIT: f64 = 1000.0;
/// The longest `REPEAT` taken.
pub const MAX_REPEAT_LENGTH: usize = 100;
/// A read is off length when its length differs from `READ_LENGTH` by more
/// than this percentage of `READ_LENGTH`.
pub const LENGTH_TOLERANCE_PERCENT: i64 = 10;
/// `ESTIMATED_LENGTH_KB` is empty when more than this percentage of the
/// counted reads are off length.
pub const MAX_OFF_LENGTH_PERCENT: i64 = 5;
/// The histogram prints every repeat count from 0 to `READ_LENGTH` over the
/// repeat length, at most to this one; a count a read has is printed above.
/// Long reads with a short repeat would otherwise print a row per base.
pub const MAX_ZERO_FILLED_KEY: i64 = 1000;
/// The zero-filled bins of all units together: with more units, each unit
/// is filled to `MAX_ZERO_FILLED_CELLS / units` at most, so a header with
/// many read groups cannot multiply the section by 1,000.
pub const MAX_ZERO_FILLED_CELLS: i64 = 1_000_000;

/// One row: a unit's counts.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub total_reads: i64,
    pub telomeric_reads: i64,
    pub pct_telomeric: f64,
    pub reads_gc_48_52: i64,
    pub read_length: Option<i64>,
    pub pct_off_length: f64,
    pub estimated_length_kb: Option<f64>,
    pub labels: Labels,
}

impl MetricRow for Row {
    const CLASS: &'static str = CLASS;
    const COLUMNS: &'static [&'static str] = &[
        "TOTAL_READS",
        "TELOMERIC_READS",
        "PCT_TELOMERIC",
        "READS_GC_48_52",
        "READ_LENGTH",
        "PCT_OFF_LENGTH",
        "ESTIMATED_LENGTH_KB",
        "SAMPLE",
        "LIBRARY",
        "READ_GROUP",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.total_reads.into(),
            self.telomeric_reads.into(),
            self.pct_telomeric.into(),
            self.reads_gc_48_52.into(),
            self.read_length.into(),
            self.pct_off_length.into(),
            self.estimated_length_kb.into(),
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]
    }
}

/// The counts of one accumulation unit.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Counts {
    pub total_reads: i64,
    pub telomeric_reads: i64,
    pub reads_gc_48_52: i64,
    /// Reads by repeat count: the histogram section.
    pub repeats: Histogram,
    /// Reads by length, for `READ_LENGTH` and `PCT_OFF_LENGTH`.
    pub lengths: Histogram,
}

/// The constants of TelSeq's length formula.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Formula {
    pub telomere_ends: i64,
    pub genome_length_at_gc: i64,
}

impl Default for Formula {
    fn default() -> Self {
        Self {
            telomere_ends: DEFAULT_TELOMERE_ENDS,
            genome_length_at_gc: DEFAULT_GENOME_LENGTH_AT_GC,
        }
    }
}

/// TelSeq's length estimate in kilobases, `l = t_k x s / c`:
/// `telomeric / gc x genome_length_at_gc / 1000 / telomere_ends`, in this
/// order. `None` when no read is in the GC bin.
pub fn estimate_kb(telomeric: i64, gc: i64, formula: Formula) -> Option<f64> {
    (gc > 0).then(|| {
        telomeric as f64 / gc as f64 * formula.genome_length_at_gc as f64
            / LENGTH_UNIT
            / formula.telomere_ends as f64
    })
}

/// True when a read of `length` bases with `gc` G or C bases is in
/// TelSeq's 48 to 52 percent bin: `48 x length < 100 x gc <= 52 x length`.
/// TelSeq bins `gc / length` with `floor((gc - 0.4f) / 0.02f)` and keeps
/// bins 4 and 5; in that float arithmetic 0.48 falls in bin 3 and 0.52 in
/// bin 5, and for reads under 2,000,000 bases the bins hold exactly the
/// fractions this integer test keeps (docs/telomere.md section 2).
pub fn in_gc_bin(gc: usize, length: usize) -> bool {
    let (gc, length) = (gc as u128 * 100, length as u128);
    length > 0 && 48 * length < gc && gc <= 52 * length
}

/// Non-overlapping occurrences of `pattern` in `read`, scanning left to
/// right and resuming after each match (TelSeq's `countMotif`).
pub fn count_occurrences(read: &[u8], pattern: &[u8]) -> usize {
    let m = pattern.len();
    let Some(&first) = pattern.first() else {
        return 0;
    };
    let Some(last_start) = read.len().checked_sub(m) else {
        return 0;
    };
    let (mut count, mut i) = (0, 0);
    while i <= last_start {
        // The next position holding the pattern's first base.
        let Some(skip) = read[i..=last_start].iter().position(|&b| b == first) else {
            break;
        };
        i += skip;
        if read[i..i + m].iter().zip(pattern).all(|(a, b)| a == b) {
            count += 1;
            i += m;
        } else {
            i += 1;
        }
    }
    count
}

/// The repeat count of a read as sequenced: the larger of the counts of
/// `repeat` and of `reverse`, its reverse complement, both as
/// [`count_occurrences`] counts them.
///
/// A repeat of at most 8 bases is matched in one pass: the last bases read
/// are kept in a `u64`, and a match ending at `i` is taken when it starts at
/// or after the end of the last match taken for the same pattern, which is
/// the left-to-right greedy count.
pub fn count_repeats(read: &[u8], repeat: &[u8], reverse: &[u8]) -> usize {
    let m = repeat.len();
    if m == 0 || m > 8 || reverse.len() != m {
        return count_occurrences(read, repeat).max(count_occurrences(read, reverse));
    }
    let mask = u64::MAX >> (64 - 8 * m);
    let word = |p: &[u8]| p.iter().fold(0u64, |w, &b| (w << 8) | u64::from(b));
    let (forward_word, reverse_word) = (word(repeat), word(reverse));
    let (mut window, mut forward, mut backward) = (0u64, 0, 0);
    // The first start a next match of each pattern may have.
    let (mut forward_free, mut backward_free) = (0, 0);
    for (i, &base) in read.iter().enumerate() {
        window = ((window << 8) | u64::from(base)) & mask;
        let Some(start) = (i + 1).checked_sub(m) else {
            continue;
        };
        if window == forward_word && start >= forward_free {
            forward += 1;
            forward_free = i + 1;
        }
        if window == reverse_word && start >= backward_free {
            backward += 1;
            backward_free = i + 1;
        }
    }
    forward.max(backward)
}

/// The bases of the 4-bit codes of a BAM record, and their complements
/// (`N` for every code that is not A, C, G or T).
const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
const COMPLEMENTS: [u8; 16] = *b"NTGNCNNNANNNNNNN";

/// The reverse complement of A, C, G, T text; any other byte becomes `N`.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|&b| complement(b)).collect()
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// The modal length (the shorter on a tie) and the reads more than
/// [`LENGTH_TOLERANCE_PERCENT`] away from it.
fn length_summary(lengths: &Histogram) -> Option<(i64, i64)> {
    let mode = lengths.mode()?;
    let off = lengths
        .keys()
        .filter(|&len| (len - mode).abs() * 100 > LENGTH_TOLERANCE_PERCENT * mode)
        .map(|len| lengths.get(len).unwrap_or(0.0) as i64)
        .sum();
    Some((mode, off))
}

/// The row of one unit. `ESTIMATED_LENGTH_KB` is empty with no read in the
/// GC bin, or when more than [`MAX_OFF_LENGTH_PERCENT`] of the reads are off
/// length.
pub fn row(counts: &Counts, formula: Formula) -> Row {
    let ratio = |n: i64| {
        if counts.total_reads == 0 {
            f64::NAN
        } else {
            n as f64 / counts.total_reads as f64
        }
    };
    let summary = length_summary(&counts.lengths);
    let off = summary.map_or(0, |(_, off)| off);
    let estimated_length_kb = if off * 100 > MAX_OFF_LENGTH_PERCENT * counts.total_reads {
        None
    } else {
        estimate_kb(counts.telomeric_reads, counts.reads_gc_48_52, formula)
    };
    Row {
        total_reads: counts.total_reads,
        telomeric_reads: counts.telomeric_reads,
        pct_telomeric: ratio(counts.telomeric_reads),
        reads_gc_48_52: counts.reads_gc_48_52,
        read_length: summary.map(|(mode, _)| mode),
        pct_off_length: ratio(off),
        estimated_length_kb,
        labels: Labels::default(),
    }
}

/// The options of one run.
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub repeat: Vec<u8>,
    pub threshold: usize,
    pub count_duplicates: bool,
    pub formula: Formula,
    pub levels: Vec<Level>,
}

/// A positive integer option.
fn positive(inv: &Invocation, name: &str, default: i64) -> Result<i64> {
    let value: i64 = inv.parse(name, default)?;
    if value < 1 {
        return Err(Error::Usage(format!("{name}={value}: must be at least 1")));
    }
    Ok(value)
}

impl Options {
    pub(crate) fn from_invocation(inv: &Invocation) -> Result<Self> {
        let text = inv.get("REPEAT").unwrap_or(DEFAULT_REPEAT);
        let repeat = text.to_ascii_uppercase().into_bytes();
        if repeat.is_empty()
            || repeat.len() > MAX_REPEAT_LENGTH
            || !repeat.iter().all(|b| b"ACGT".contains(b))
        {
            return Err(Error::Usage(format!(
                "REPEAT={text}: 1 to {MAX_REPEAT_LENGTH} bases of A, C, G and T"
            )));
        }
        let threshold = positive(inv, "TELOMERE_REPEAT_THRESHOLD", DEFAULT_THRESHOLD)?;
        Ok(Self {
            repeat,
            threshold: usize::try_from(threshold).unwrap_or(usize::MAX),
            count_duplicates: inv.flag("COUNT_DUPLICATES", false)?,
            formula: Formula {
                telomere_ends: positive(inv, "TELOMERE_ENDS", DEFAULT_TELOMERE_ENDS)?,
                genome_length_at_gc: positive(
                    inv,
                    "GENOME_LENGTH_AT_GC",
                    DEFAULT_GENOME_LENGTH_AT_GC,
                )?,
            },
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
        })
    }
}

/// What one read adds to its units.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scan {
    pub length: usize,
    pub repeats: usize,
    pub gc_bin: bool,
}

/// `CollectTelomereMetrics` as a pass [`Accumulator`].
pub struct Telomere {
    repeat: Vec<u8>,
    reverse: Vec<u8>,
    threshold: usize,
    count_duplicates: bool,
    formula: Formula,
    levels: Levels<Counts>,
    /// The read's bases, reused from record to record.
    bases: Vec<u8>,
}

impl Telomere {
    pub fn new(opts: &Options, header: &sam::Header) -> Self {
        Self {
            reverse: reverse_complement(&opts.repeat),
            repeat: opts.repeat.clone(),
            threshold: opts.threshold,
            count_duplicates: opts.count_duplicates,
            formula: opts.formula,
            levels: Levels::new(header, &opts.levels, Counts::default()),
            bases: Vec::new(),
        }
    }

    fn counted(&self, record: &bam::Record) -> bool {
        let flags = record.flags();
        !(flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_qc_fail()
            || (flags.is_duplicate() && !self.count_duplicates))
    }

    /// The read as sequenced, scanned.
    fn scan(&mut self, record: &bam::Record) -> Scan {
        let sequence = record.sequence();
        let reverse = record.flags().is_reverse_complemented();
        let table = if reverse { &COMPLEMENTS } else { &BASES };
        self.bases.clear();
        for &pair in sequence.as_bytes() {
            self.bases.push(table[usize::from(pair >> 4)]);
            self.bases.push(table[usize::from(pair & 0x0F)]);
        }
        self.bases.truncate(sequence.len());
        if reverse {
            self.bases.reverse();
        }
        let gc = self
            .bases
            .iter()
            .filter(|&&b| b == b'G' || b == b'C')
            .count();
        Scan {
            length: self.bases.len(),
            repeats: count_repeats(&self.bases, &self.repeat, &self.reverse),
            gc_bin: in_gc_bin(gc, self.bases.len()),
        }
    }
}

impl Accumulator for Telomere {
    type Output = Vec<Bucket>;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if !self.counted(record) {
            return Ok(());
        }
        let scan = self.scan(record);
        if scan.length == 0 {
            return Ok(());
        }
        let telomeric = scan.repeats >= self.threshold;
        self.levels.route(record, |unit| {
            unit.total_reads += 1;
            unit.telomeric_reads += i64::from(telomeric);
            unit.reads_gc_48_52 += i64::from(scan.gc_bin);
            unit.repeats.increment(scan.repeats as i64);
            unit.lengths.increment(scan.length as i64);
        });
        Ok(())
    }

    fn finish(self) -> Result<Vec<Bucket>> {
        let buckets = self.levels.into_buckets();
        let filled_units = buckets.iter().filter(|b| b.state.total_reads > 0).count();
        let fill_cap = zero_fill_cap(filled_units);
        let motif = self.repeat.len() as i64;
        Ok(buckets
            .into_iter()
            .map(|bucket| {
                let mut row = row(&bucket.state, self.formula);
                row.labels = bucket.labels;
                let mut repeats = bucket.state.repeats;
                if let Some(length) = row.read_length {
                    for key in 0..=(length / motif).min(fill_cap) {
                        repeats.increment_by(key, 0.0);
                    }
                }
                Bucket {
                    row,
                    label: format!("{}.read_count", bucket.key),
                    repeats,
                }
            })
            .collect())
    }
}

/// The last key zero-filled in each of `units` units with reads: at most
/// [`MAX_ZERO_FILLED_KEY`], and at most [`MAX_ZERO_FILLED_CELLS`] bins over
/// all units.
pub fn zero_fill_cap(units: usize) -> i64 {
    let units = i64::try_from(units.max(1)).unwrap_or(i64::MAX);
    (MAX_ZERO_FILLED_CELLS / units - 1).clamp(0, MAX_ZERO_FILLED_KEY)
}

/// One unit's row and its histogram column.
#[derive(Clone, Debug)]
pub struct Bucket {
    pub row: Row,
    pub label: String,
    pub repeats: Histogram,
}

/// The metrics file: the rows, then the reads by repeat count per unit.
pub fn output(buckets: Vec<Bucket>) -> MetricsOutput {
    let rows: Vec<Row> = buckets.iter().map(|b| b.row.clone()).collect();
    let series: Vec<HistogramSeries> = buckets
        .into_iter()
        .filter(|b| !b.repeats.is_empty())
        .map(|b| HistogramSeries {
            label: b.label,
            histogram: b.repeats,
        })
        .collect();
    let histogram = (!series.is_empty()).then_some(HistogramSection {
        key_label: "repeat_count",
        key_class: KeyClass::Integer,
        series,
    });
    MetricsOutput::from_rows(&rows, histogram)
}

pub const TOOL: Tool = Tool {
    name: "CollectTelomereMetrics",
    summary: "one row of the 10 TelomereMetrics columns per accumulation unit: reads holding at \
              least TELOMERE_REPEAT_THRESHOLD copies of TTAGGG or CCCTAA, reads of 48 to 52 \
              percent GC, and TelSeq's telomere length estimate in kb; the histogram counts \
              reads by repeat count. Not a Picard tool (docs/telomere.md).",
    options: &[
        Opt::int(
            &["TELOMERE_REPEAT_THRESHOLD", "K"],
            "A read with at least this many repeats is telomeric, TelSeq's k (default 7).",
        ),
        Opt::text(
            &["REPEAT"],
            "The telomere repeat, A, C, G and T; its reverse complement is counted too \
             (default TTAGGG).",
        ),
        Opt::bool(
            &["COUNT_DUPLICATES"],
            "Count reads flagged as duplicates (default false).",
        ),
        Opt::int(
            &["TELOMERE_ENDS"],
            "Telomere ends per genome in the length formula (default 46).",
        ),
        Opt::int(
            &["GENOME_LENGTH_AT_GC"],
            "Reference bases in regions of 48 to 52 percent GC in the length formula (default \
             332720800, TelSeq's hg19 value).",
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
    Ok(pass::boxed(Telomere::new(&opts, header), |buckets| {
        Ok(output(buckets).into())
    }))
}

#[cfg(test)]
mod tests;
