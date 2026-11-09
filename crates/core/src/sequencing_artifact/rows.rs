//! The rows of `CollectSequencingArtifactMetrics` and how one library's
//! counts become them (`ArtifactCounter.finish`,
//! `ContextAccumulator.calculateMetrics`, the error summary of `finish`).

use std::collections::BTreeMap;

use crate::metric::{MetricRow, Value};

/// The smallest error rate reported (`MIN_ERROR`).
const MIN_ERROR: f64 = 1e-10;
pub(super) const BASES: [u8; 4] = *b"ACGT";

/// `Transition.ALT_VALUES`: every substitution, reference then call base
/// in `ACGT` order.
pub(super) const ALT_TRANSITIONS: [(usize, usize); 12] = [
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 0),
    (1, 2),
    (1, 3),
    (2, 0),
    (2, 1),
    (2, 3),
    (3, 0),
    (3, 1),
    (3, 2),
];

/// `PreAdapterDetailMetrics`.
#[derive(Clone, Debug, PartialEq)]
pub struct PreAdapterDetail {
    pub sample_alias: String,
    pub library: String,
    pub ref_base: u8,
    pub alt_base: u8,
    pub context: String,
    pub pro_ref_bases: i64,
    pub pro_alt_bases: i64,
    pub con_ref_bases: i64,
    pub con_alt_bases: i64,
    pub error_rate: f64,
    pub qscore: f64,
}

impl MetricRow for PreAdapterDetail {
    const CLASS: &'static str =
        "picard.analysis.artifacts.SequencingArtifactMetrics$PreAdapterDetailMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "SAMPLE_ALIAS",
        "LIBRARY",
        "REF_BASE",
        "ALT_BASE",
        "CONTEXT",
        "PRO_REF_BASES",
        "PRO_ALT_BASES",
        "CON_REF_BASES",
        "CON_ALT_BASES",
        "ERROR_RATE",
        "QSCORE",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.sample_alias.as_str().into(),
            self.library.as_str().into(),
            char_value(self.ref_base),
            char_value(self.alt_base),
            self.context.as_str().into(),
            self.pro_ref_bases.into(),
            self.pro_alt_bases.into(),
            self.con_ref_bases.into(),
            self.con_alt_bases.into(),
            self.error_rate.into(),
            self.qscore.into(),
        ]
    }
}

impl PreAdapterDetail {
    /// `calculateDerivedStatistics`.
    fn derive(&mut self) {
        self.error_rate = MIN_ERROR;
        let total =
            self.pro_ref_bases + self.pro_alt_bases + self.con_ref_bases + self.con_alt_bases;
        if total > 0 {
            let raw = (self.pro_alt_bases - self.con_alt_bases) as f64 / total as f64;
            self.error_rate = java_max(MIN_ERROR, raw);
        }
        self.qscore = phred(self.error_rate);
    }
}

/// `BaitBiasDetailMetrics`.
#[derive(Clone, Debug, PartialEq)]
pub struct BaitBiasDetail {
    pub sample_alias: String,
    pub library: String,
    pub ref_base: u8,
    pub alt_base: u8,
    pub context: String,
    pub fwd_cxt_ref_bases: i64,
    pub fwd_cxt_alt_bases: i64,
    pub rev_cxt_ref_bases: i64,
    pub rev_cxt_alt_bases: i64,
    pub fwd_error_rate: f64,
    pub rev_error_rate: f64,
    pub error_rate: f64,
    pub qscore: f64,
}

impl MetricRow for BaitBiasDetail {
    const CLASS: &'static str =
        "picard.analysis.artifacts.SequencingArtifactMetrics$BaitBiasDetailMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "SAMPLE_ALIAS",
        "LIBRARY",
        "REF_BASE",
        "ALT_BASE",
        "CONTEXT",
        "FWD_CXT_REF_BASES",
        "FWD_CXT_ALT_BASES",
        "REV_CXT_REF_BASES",
        "REV_CXT_ALT_BASES",
        "FWD_ERROR_RATE",
        "REV_ERROR_RATE",
        "ERROR_RATE",
        "QSCORE",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.sample_alias.as_str().into(),
            self.library.as_str().into(),
            char_value(self.ref_base),
            char_value(self.alt_base),
            self.context.as_str().into(),
            self.fwd_cxt_ref_bases.into(),
            self.fwd_cxt_alt_bases.into(),
            self.rev_cxt_ref_bases.into(),
            self.rev_cxt_alt_bases.into(),
            self.fwd_error_rate.into(),
            self.rev_error_rate.into(),
            self.error_rate.into(),
            self.qscore.into(),
        ]
    }
}

impl BaitBiasDetail {
    /// `calculateDerivedStatistics`.
    fn derive(&mut self) {
        let rate = |alt: i64, reference: i64| {
            let bases = reference + alt;
            if bases > 0 {
                java_max(MIN_ERROR, alt as f64 / bases as f64)
            } else {
                MIN_ERROR
            }
        };
        self.fwd_error_rate = rate(self.fwd_cxt_alt_bases, self.fwd_cxt_ref_bases);
        self.rev_error_rate = rate(self.rev_cxt_alt_bases, self.rev_cxt_ref_bases);
        self.error_rate = java_max(MIN_ERROR, self.fwd_error_rate - self.rev_error_rate);
        self.qscore = phred(self.error_rate);
    }
}

/// `PreAdapterSummaryMetrics` or `BaitBiasSummaryMetrics`: the same
/// columns, told apart by the artifact names.
#[derive(Clone, Debug, PartialEq)]
pub struct Summary {
    pub sample_alias: String,
    pub library: String,
    pub ref_base: u8,
    pub alt_base: u8,
    pub total_qscore: f64,
    pub worst_cxt: String,
    pub worst_cxt_qscore: f64,
    pub worst_pre_cxt: String,
    pub worst_pre_cxt_qscore: f64,
    pub worst_post_cxt: String,
    pub worst_post_cxt_qscore: f64,
    pub artifact_name: &'static str,
}

const SUMMARY_COLUMNS: &[&str] = &[
    "SAMPLE_ALIAS",
    "LIBRARY",
    "REF_BASE",
    "ALT_BASE",
    "TOTAL_QSCORE",
    "WORST_CXT",
    "WORST_CXT_QSCORE",
    "WORST_PRE_CXT",
    "WORST_PRE_CXT_QSCORE",
    "WORST_POST_CXT",
    "WORST_POST_CXT_QSCORE",
    "ARTIFACT_NAME",
];

impl Summary {
    fn values(&self) -> Vec<Value> {
        vec![
            self.sample_alias.as_str().into(),
            self.library.as_str().into(),
            char_value(self.ref_base),
            char_value(self.alt_base),
            self.total_qscore.into(),
            self.worst_cxt.as_str().into(),
            self.worst_cxt_qscore.into(),
            self.worst_pre_cxt.as_str().into(),
            self.worst_pre_cxt_qscore.into(),
            self.worst_post_cxt.as_str().into(),
            self.worst_post_cxt_qscore.into(),
            self.artifact_name.into(),
        ]
    }
}

/// A pre-adapter summary row.
#[derive(Clone, Debug, PartialEq)]
pub struct PreAdapterSummary(pub Summary);

impl MetricRow for PreAdapterSummary {
    const CLASS: &'static str =
        "picard.analysis.artifacts.SequencingArtifactMetrics$PreAdapterSummaryMetrics";
    const COLUMNS: &'static [&'static str] = SUMMARY_COLUMNS;

    fn values(&self) -> Vec<Value> {
        self.0.values()
    }
}

/// A bait-bias summary row.
#[derive(Clone, Debug, PartialEq)]
pub struct BaitBiasSummary(pub Summary);

impl MetricRow for BaitBiasSummary {
    const CLASS: &'static str =
        "picard.analysis.artifacts.SequencingArtifactMetrics$BaitBiasSummaryMetrics";
    const COLUMNS: &'static [&'static str] = SUMMARY_COLUMNS;

    fn values(&self) -> Vec<Value> {
        self.0.values()
    }
}

/// `ErrorSummaryMetrics`.
#[derive(Clone, Debug, PartialEq)]
pub struct ErrorSummary {
    pub ref_base: u8,
    pub alt_base: u8,
    pub substitution: String,
    pub ref_count: i64,
    pub alt_count: i64,
    pub substitution_rate: f64,
}

impl MetricRow for ErrorSummary {
    const CLASS: &'static str = "picard.analysis.artifacts.ErrorSummaryMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "REF_BASE",
        "ALT_BASE",
        "SUBSTITUTION",
        "REF_COUNT",
        "ALT_COUNT",
        "SUBSTITUTION_RATE",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            char_value(self.ref_base),
            char_value(self.alt_base),
            self.substitution.as_str().into(),
            self.ref_count.into(),
            self.alt_count.into(),
            self.substitution_rate.into(),
        ]
    }
}

/// The rows of the five files, detail rows already cut to
/// `CONTEXTS_TO_PRINT`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Report {
    pub pre_adapter_details: Vec<PreAdapterDetail>,
    pub pre_adapter_summaries: Vec<PreAdapterSummary>,
    pub bait_bias_details: Vec<BaitBiasDetail>,
    pub bait_bias_summaries: Vec<BaitBiasSummary>,
    pub error_summaries: Vec<ErrorSummary>,
}

/// A one-character field.
fn char_value(base: u8) -> Value {
    Value::Str(char::from(base).to_string())
}

/// `Math.max`: NaN when either is NaN.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// `QualityUtil.getPhredScoreFromErrorProbability`: rounded to an integer.
fn phred(error_rate: f64) -> f64 {
    (-10.0 * error_rate.log10()).round()
}

/// Index of an upper-case `ACGT` base.
pub(super) fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// `SequenceUtil.complement` of a base or of `N`.
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// `SequenceUtil.reverseComplement` of a context.
fn reverse_complement(context: &str) -> String {
    context
        .bytes()
        .rev()
        .map(|b| char::from(complement(b)))
        .collect()
}

/// `AlignmentAccumulator`: bases seen on read one forward, read one
/// reverse, read two forward and read two reverse.
pub(super) type Strands = [i64; 4];
/// One context's strands by called base `A`, `C`, `G`, `T`.
pub(super) type Calls = [Strands; 4];

pub(super) const R1_POS: usize = 0;
pub(super) const R1_NEG: usize = 1;
pub(super) const R2_POS: usize = 2;
pub(super) const R2_NEG: usize = 3;

fn merge(into: &mut Calls, from: &Calls) {
    for (a, b) in into.iter_mut().zip(from) {
        for (x, y) in a.iter_mut().zip(b) {
            *x += y;
        }
    }
}

/// Counts by context string, sorted as Picard's `TreeSet` sorts them.
type Table = BTreeMap<String, Calls>;

/// Detail rows of one table, one list per transition (`ref * 4 + alt`).
type DetailLists = Vec<Vec<(PreAdapterDetail, BaitBiasDetail)>>;

/// One library's detail rows and what its summaries need.
pub(super) struct LibraryRows {
    sample_alias: String,
    library: String,
    context_size: usize,
    /// Full-context rows per `ALT_TRANSITIONS` entry.
    pub(super) details: Vec<Vec<(PreAdapterDetail, BaitBiasDetail)>>,
    full: DetailLists,
    half: DetailLists,
    zero: DetailLists,
}

impl LibraryRows {
    pub(super) fn new(
        sample_alias: &str,
        library: &str,
        k: usize,
        tandem: bool,
        counts: &[Calls],
    ) -> Self {
        let length = 2 * k + 1;
        let mut full = Table::new();
        for (code, calls) in counts.iter().enumerate() {
            let mut context = vec![b'A'; length];
            let mut rest = code;
            for slot in context.iter_mut().rev() {
                *slot = BASES[rest % 4];
                rest /= 4;
            }
            full.insert(String::from_utf8(context).expect("ACGT"), *calls);
        }
        let padding = "N".repeat(k);
        let leading = |c: &str| format!("{}{padding}", &c[..=k]);
        let trailing = |c: &str| format!("{padding}{}", &c[k..]);
        let zero_key = |c: &str| format!("{padding}{}{padding}", &c[k..=k]);
        // `fillHalfRecords` and `fillZeroRecords`: every key exists first,
        // then each full context is merged in. With CONTEXT_SIZE=0 the
        // leading and trailing keys are one, merged twice.
        let mut half = Table::new();
        let mut zero = Table::new();
        for context in full.keys() {
            half.insert(leading(context), [[0; 4]; 4]);
            half.insert(trailing(context), [[0; 4]; 4]);
            zero.insert(zero_key(context), [[0; 4]; 4]);
        }
        for (context, calls) in &full {
            merge(half.get_mut(&trailing(context)).expect("key"), calls);
            merge(half.get_mut(&leading(context)).expect("key"), calls);
            merge(zero.get_mut(&zero_key(context)).expect("key"), calls);
        }
        let full = calculate(&full, tandem, sample_alias, library);
        let half = calculate(&half, tandem, sample_alias, library);
        let zero = calculate(&zero, tandem, sample_alias, library);
        let details = ALT_TRANSITIONS
            .iter()
            .map(|&(r, a)| full[r * 4 + a].clone())
            .collect();
        Self {
            sample_alias: sample_alias.to_owned(),
            library: library.to_owned(),
            context_size: k,
            details,
            full,
            half,
            zero,
        }
    }

    /// `getSummaryMetrics` for one transition: the pre-adapter and the
    /// bait-bias summary rows.
    pub(super) fn summary(&self, (r, a): (usize, usize)) -> (Summary, Summary) {
        let t = r * 4 + a;
        let k = self.context_size;
        let padding = "N".repeat(k);
        let is_leading = |c: &str| k == 0 || c.ends_with(&padding);
        let is_trailing = |c: &str| k == 0 || c.starts_with(&padding);
        let total = &self.zero[t][0];
        let full: Vec<&_> = self.full[t].iter().collect();
        let leading: Vec<&_> = self.half[t]
            .iter()
            .filter(|(p, _)| is_leading(&p.context))
            .collect();
        let trailing: Vec<&_> = self.half[t]
            .iter()
            .filter(|(p, _)| is_trailing(&p.context))
            .collect();
        let worst_pre = |rows: &[&(PreAdapterDetail, BaitBiasDetail)]| {
            rows.iter()
                .map(|(p, _)| p)
                .min_by(|x, y| order(x.qscore, &x.context, y.qscore, &y.context))
                .map(|p| (p.context.clone(), p.qscore))
                .expect("one context at least")
        };
        let worst_bait = |rows: &[&(PreAdapterDetail, BaitBiasDetail)]| {
            rows.iter()
                .map(|(_, b)| b)
                .min_by(|x, y| order(x.qscore, &x.context, y.qscore, &y.context))
                .map(|b| (b.context.clone(), b.qscore))
                .expect("one context at least")
        };
        let summary = |total_qscore, worst: [(String, f64); 3], artifact_name| {
            let [cxt, pre, post] = worst;
            Summary {
                sample_alias: self.sample_alias.clone(),
                library: self.library.clone(),
                ref_base: BASES[r],
                alt_base: BASES[a],
                total_qscore,
                worst_cxt: cxt.0,
                worst_cxt_qscore: cxt.1,
                worst_pre_cxt: pre.0,
                worst_pre_cxt_qscore: pre.1,
                worst_post_cxt: post.0,
                worst_post_cxt_qscore: post.1,
                artifact_name,
            }
        };
        let pre_name = match (BASES[r], BASES[a]) {
            (b'G', b'T') => "OxoG",
            (b'C', b'T') => "Deamination",
            _ => "NA",
        };
        let bait_name = match (BASES[r], BASES[a]) {
            (b'G', b'T') => "Gref",
            (b'C', b'A') => "Cref",
            _ => "NA",
        };
        let pre = summary(
            total.0.qscore,
            [worst_pre(&full), worst_pre(&leading), worst_pre(&trailing)],
            pre_name,
        );
        let bait = summary(
            total.1.qscore,
            [
                worst_bait(&full),
                worst_bait(&leading),
                worst_bait(&trailing),
            ],
            bait_name,
        );
        (pre, bait)
    }
}

/// `compareTo` of two detail rows of one transition: `QSCORE`, then the
/// context.
fn order(q1: f64, c1: &str, q2: f64, c2: &str) -> std::cmp::Ordering {
    q1.total_cmp(&q2).then_with(|| c1.cmp(c2))
}

/// `ContextAccumulator.calculateMetrics`: for each context in sorted order
/// and each called base, the pre-adapter and bait-bias row of the
/// transition from the centre base to it, ref-to-ref included.
fn calculate(table: &Table, tandem: bool, sample_alias: &str, library: &str) -> DetailLists {
    let mut lists: DetailLists = vec![Vec::new(); 16];
    for (context, calls) in table {
        let centre = context.as_bytes()[context.len() / 2];
        let r = base_index(centre).expect("contexts are centred on ACGT");
        let reverse = &table[&reverse_complement(context)];
        for a in 0..4 {
            let fwd_ref = &calls[r];
            let fwd_alt = &calls[a];
            let rev_ref = &reverse[3 - r];
            let rev_alt = &reverse[3 - a];
            let (pro_ref, pro_alt, con_ref, con_alt) = if tandem {
                (
                    fwd_ref[R1_POS] + fwd_ref[R2_POS] + rev_ref[R1_NEG] + rev_ref[R2_NEG],
                    fwd_alt[R1_POS] + fwd_alt[R2_POS] + rev_alt[R1_NEG] + rev_alt[R2_NEG],
                    fwd_ref[R1_NEG] + fwd_ref[R2_NEG] + rev_ref[R1_POS] + rev_ref[R2_POS],
                    fwd_alt[R1_NEG] + fwd_alt[R2_NEG] + rev_alt[R1_POS] + rev_alt[R2_POS],
                )
            } else {
                (
                    fwd_ref[R1_POS] + fwd_ref[R2_NEG] + rev_ref[R1_NEG] + rev_ref[R2_POS],
                    fwd_alt[R1_POS] + fwd_alt[R2_NEG] + rev_alt[R1_NEG] + rev_alt[R2_POS],
                    fwd_ref[R1_NEG] + fwd_ref[R2_POS] + rev_ref[R1_POS] + rev_ref[R2_NEG],
                    fwd_alt[R1_NEG] + fwd_alt[R2_POS] + rev_alt[R1_POS] + rev_alt[R2_NEG],
                )
            };
            let mut pre = PreAdapterDetail {
                sample_alias: sample_alias.to_owned(),
                library: library.to_owned(),
                ref_base: BASES[r],
                alt_base: BASES[a],
                context: context.clone(),
                pro_ref_bases: pro_ref,
                pro_alt_bases: pro_alt,
                con_ref_bases: con_ref,
                con_alt_bases: con_alt,
                error_rate: 0.0,
                qscore: 0.0,
            };
            pre.derive();
            let sum = |s: &Strands| s.iter().sum::<i64>();
            let mut bait = BaitBiasDetail {
                sample_alias: sample_alias.to_owned(),
                library: library.to_owned(),
                ref_base: BASES[r],
                alt_base: BASES[a],
                context: context.clone(),
                fwd_cxt_ref_bases: sum(fwd_ref),
                fwd_cxt_alt_bases: sum(fwd_alt),
                rev_cxt_ref_bases: sum(rev_ref),
                rev_cxt_alt_bases: sum(rev_alt),
                fwd_error_rate: 0.0,
                rev_error_rate: 0.0,
                error_rate: 0.0,
                qscore: 0.0,
            };
            bait.derive();
            lists[r * 4 + a].push((pre, bait));
        }
    }
    lists
}

/// The error summary of the printed pre-adapter rows: a `G` or `T`
/// reference folded to its complement, then one row per substitution in
/// sorted order.
pub(super) fn error_summaries(details: &[PreAdapterDetail]) -> Vec<ErrorSummary> {
    let mut groups: BTreeMap<String, (u8, u8, i64, i64)> = BTreeMap::new();
    for m in details {
        let (r, a) = if matches!(m.ref_base, b'G' | b'T') {
            (complement(m.ref_base), complement(m.alt_base))
        } else {
            (m.ref_base, m.alt_base)
        };
        let key = format!("{}>{}", char::from(r), char::from(a));
        let entry = groups.entry(key).or_insert((r, a, 0, 0));
        entry.2 += m.pro_ref_bases + m.con_ref_bases;
        entry.3 += m.pro_alt_bases + m.con_alt_bases;
    }
    groups
        .into_iter()
        .map(|(substitution, (r, a, ref_count, alt_count))| {
            let total = (ref_count + alt_count) as f64;
            ErrorSummary {
                ref_base: r,
                alt_base: a,
                substitution,
                ref_count,
                alt_count,
                substitution_rate: if total == 0.0 {
                    0.0
                } else {
                    alt_count as f64 / total
                },
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred_scores_round_and_floor_at_100() {
        assert_eq!(phred(MIN_ERROR), 100.0);
        assert_eq!(phred(0.001588), 28.0);
        assert_eq!(phred(0.5), 3.0);
        assert_eq!(phred(1.0), 0.0);
    }

    #[test]
    fn java_max_keeps_nan() {
        assert!(java_max(f64::NAN, 1.0).is_nan());
        assert!(java_max(1e-10, f64::NAN).is_nan());
        assert_eq!(java_max(1e-10, -0.5), 1e-10);
    }

    #[test]
    fn reverse_complement_keeps_padding() {
        assert_eq!(reverse_complement("AAN"), "NTT");
        assert_eq!(reverse_complement("NCN"), "NGN");
        assert_eq!(reverse_complement("ACGT"), "ACGT");
    }

    /// One pre-adapter row checked by hand: A>C at `CAG`, read one forward
    /// counts pro, read one reverse con; the reverse complement `CTG`
    /// counts the other way round.
    #[test]
    fn pro_and_con_follow_read_end_and_strand() {
        let mut counts = vec![[[0i64; 4]; 4]; 64];
        let code = |c: &str| c.bytes().fold(0, |n, b| n * 4 + base_index(b).unwrap());
        counts[code("CAG")][1][R1_POS] = 3; // C called, read one forward
        counts[code("CAG")][1][R1_NEG] = 1;
        counts[code("CAG")][0][R2_NEG] = 10; // A (reference) on read two reverse
        counts[code("CTG")][3][R1_NEG] = 7; // T, the reference complement
        counts[code("CTG")][2][R1_NEG] = 2; // G, the call complement
        let rows = LibraryRows::new("S", "L", 1, false, &counts);
        let (pre, bait) = rows.details[0]
            .iter()
            .find(|(p, _)| p.context == "CAG")
            .expect("CAG row");
        assert_eq!(
            (
                pre.pro_ref_bases,
                pre.pro_alt_bases,
                pre.con_ref_bases,
                pre.con_alt_bases
            ),
            (17, 5, 0, 1)
        );
        assert_eq!(
            (
                bait.fwd_cxt_ref_bases,
                bait.fwd_cxt_alt_bases,
                bait.rev_cxt_ref_bases,
                bait.rev_cxt_alt_bases
            ),
            (10, 4, 7, 2)
        );
        let tandem = LibraryRows::new("S", "L", 1, true, &counts);
        let (pre, _) = tandem.details[0]
            .iter()
            .find(|(p, _)| p.context == "CAG")
            .expect("CAG row");
        assert_eq!(
            (
                pre.pro_ref_bases,
                pre.pro_alt_bases,
                pre.con_ref_bases,
                pre.con_alt_bases
            ),
            (7, 5, 10, 1)
        );
    }

    #[test]
    fn empty_counts_give_quality_100_everywhere() {
        let counts = vec![[[0i64; 4]; 4]; 4];
        let rows = LibraryRows::new("S", "L", 0, false, &counts);
        let (pre, bait) = rows.summary((2, 3));
        assert_eq!(pre.artifact_name, "OxoG");
        assert_eq!(bait.artifact_name, "Gref");
        assert_eq!((pre.total_qscore, pre.worst_cxt.as_str()), (100.0, "G"));
        assert_eq!(
            (pre.worst_pre_cxt.as_str(), pre.worst_post_cxt.as_str()),
            ("G", "G")
        );
    }

    #[test]
    fn error_summary_folds_g_and_t() {
        let row = |r: u8, a: u8, pro_ref, pro_alt| PreAdapterDetail {
            sample_alias: String::new(),
            library: String::new(),
            ref_base: r,
            alt_base: a,
            context: String::new(),
            pro_ref_bases: pro_ref,
            pro_alt_bases: pro_alt,
            con_ref_bases: 0,
            con_alt_bases: 0,
            error_rate: 0.0,
            qscore: 0.0,
        };
        let rows = error_summaries(&[
            row(b'G', b'T', 10, 2),
            row(b'C', b'A', 5, 1),
            row(b'A', b'C', 0, 0),
        ]);
        let got: Vec<(&str, i64, i64, f64)> = rows
            .iter()
            .map(|r| {
                (
                    r.substitution.as_str(),
                    r.ref_count,
                    r.alt_count,
                    r.substitution_rate,
                )
            })
            .collect();
        assert_eq!(got, [("A>C", 0, 0, 0.0), ("C>A", 15, 3, 3.0 / 18.0)]);
    }
}
