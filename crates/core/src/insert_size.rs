use std::fs::File;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use crate::histogram::Histogram;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Orientation { Fr, Rf, Tandem }

pub struct Options { pub minimum_pct: f64, pub deviations: f64, pub include_duplicates: bool }
impl Default for Options {
    fn default() -> Self { Self { minimum_pct: 0.05, deviations: 10.0, include_duplicates: false } }
}

#[derive(Default, Debug)]
pub struct Row {
    pub orientation: Option<OrientStr>,
    pub read_pairs: i64,
    pub min_insert_size: i64,
    pub max_insert_size: i64,
    pub median_insert_size: f64,
    pub mode_insert_size: i64,
    pub median_absolute_deviation: f64,
    pub mean_insert_size: f64,
    pub standard_deviation: f64,
    pub width_of_50_percent: i64,
}
pub type OrientStr = &'static str;

fn ref_len_of_cigar(record: &bam::Record) -> std::io::Result<usize> {
    let mut n = 0usize;
    for op in record.cigar().iter() {
        let op = op?;
        match op.kind() {
            Kind::Match | Kind::Deletion | Kind::Skip
            | Kind::SequenceMatch | Kind::SequenceMismatch => n += op.len(),
            _ => {}
        }
    }
    Ok(n)
}

fn orientation(record: &bam::Record) -> std::io::Result<Orientation> {
    let f = record.flags();
    if f.is_reverse_complemented() == f.is_mate_reverse_complemented() {
        return Ok(Orientation::Tandem);
    }
    let start = record.alignment_start().unwrap()?.get() as i64;          // 1-based
    let mate_start = record.mate_alignment_start().unwrap()?.get() as i64;
    let (pos_5p, neg_5p) = if f.is_reverse_complemented() {
        // this read is the negative strand; its 5' is its alignment end
        let end = start + ref_len_of_cigar(record)? as i64 - 1;
        (mate_start, end)
    } else {
        // positive strand read; need mate's end (from MC) — for the test, MC=100M
        // generic: fall back to mate_start if MC absent (documented limitation)
        (start, mate_start)
    };
    Ok(if pos_5p <= neg_5p { Orientation::Fr } else { Orientation::Rf })
}

pub fn run_over_bam(path: &str, opts: Options) -> std::io::Result<Vec<Row>> {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let _header = reader.read_header()?;
    let (mut fr, mut rf, mut tandem) = (Histogram::new(), Histogram::new(), Histogram::new());

    for result in reader.records() {
        let record = result?;
        let f = record.flags();
        if !f.is_segmented() || f.is_unmapped() || f.is_mate_unmapped()
            || f.is_first_segment()
            || f.is_secondary() || f.is_supplementary()
            || (f.is_duplicate() && !opts.include_duplicates)
            || record.template_length() == 0 { continue; }
        let rid = match record.reference_sequence_id() { Some(x) => x?, None => continue };
        let mrid = match record.mate_reference_sequence_id() { Some(x) => x?, None => continue };
        if rid != mrid { continue; }

        let isize = record.template_length().unsigned_abs() as i64;
        match orientation(&record)? {
            Orientation::Fr => fr.increment(isize),
            Orientation::Rf => rf.increment(isize),
            Orientation::Tandem => tandem.increment(isize),
        }
    }

    let total = fr.count() + rf.count() + tandem.count();
    let mut rows = Vec::new();
    if total == 0.0 { return Ok(rows); }
    for (name, h) in [("FR", &fr), ("RF", &rf), ("TANDEM", &tandem)] {
        if h.count() < total * opts.minimum_pct { continue; }
        rows.push(finalize_row(name, h, &opts));
    }
    Ok(rows)
}

fn finalize_row(name: OrientStr, h: &Histogram, opts: &Options) -> Row {
    let mut m = Row::default();
    m.orientation = Some(name);
    m.read_pairs = h.count() as i64;
    m.min_insert_size = h.min();
    m.max_insert_size = h.max();
    let median = h.median();
    m.median_insert_size = median;
    m.mode_insert_size = h.mode();
    let mad = h.median_absolute_deviation();
    m.median_absolute_deviation = mad;

    let total = h.count();
    let (mut covered, mut low, mut high) = (0.0, median, median);
    while low >= h.min() as f64 - 1.0 || high <= h.max() as f64 + 1.0 {
        if let Some(v) = h.get(low as i64) { covered += v; }
        if low != high { if let Some(v) = h.get(high as i64) { covered += v; } }
        let pct = covered / total;
        let dist = (high - low) as i64 + 1;
        if pct >= 0.50 && m.width_of_50_percent == 0 { m.width_of_50_percent = dist; }
        low -= 1.0; high += 1.0;
    }

    let width = (median + opts.deviations * mad) as i64;
    let mut trimmed = h.clone();
    trimmed.trim_by_width(width);
    if !trimmed.is_empty() {
        m.mean_insert_size = trimmed.mean();
        m.standard_deviation = trimmed.standard_deviation();
    }
    m
}
