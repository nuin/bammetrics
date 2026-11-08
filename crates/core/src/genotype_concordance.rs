//! `GenotypeConcordance` (plan row P7, machinery M21): the scheme table,
//! the allele normalization and the three metrics files.
//!
//! The truth and call VCFs are walked together in position order of the
//! truth dictionary (`PairedVariantSubContextIterator`): a record present
//! on one side only pairs with nothing. Each record is cut down to the
//! chosen sample (`VariantContext.subContextFromSample`), which re-derives
//! its alleles from that genotype, so a hom-ref or no-call genotype makes a
//! `NO_VARIATION` record. The pair's types pick the counter
//! (`classifyVariants`): a SNP truth with a SNP, mixed or non-variant call
//! goes to the SNP counter, an indel truth likewise to the indel counter, a
//! mixed or non-variant truth follows the call's type, and every other
//! combination is counted nowhere.
//!
//! Each side gets a state (`determineState`). A missing record is
//! `MISSING`; a mixed record `IS_MIXED`; a filtered record `VC_FILTERED`
//! (the call side only unless `IGNORE_FILTER_STATUS`); then the genotype:
//! all no-call `NO_CALL`, `FT` set `GT_FILTERED`, `GQ` under `MIN_GQ`
//! `LOW_GQ`, `DP` under `MIN_DP` `LOW_DP`, half called `NO_CALL`. The
//! remaining genotypes are named by their alleles' positions in one ordered
//! list built by `normalizeAlleles`: the reference, the truth alleles, then
//! the call alleles, after the shorter reference of the two records has the
//! longer one's suffix spliced into its alleles. When a call allele already
//! sits past position 1, the truth alleles are re-added in swapped order
//! first. `TruthState` reads positions 0 to 2 (`HOM_REF`, `HET_REF_VAR1`,
//! `HET_VAR1_VAR2`, `HOM_VAR1`); `CallState` reads 0 to 4 (adding
//! `HET_REF_VAR2`, `HET_REF_VAR3`, `HET_VAR1_VAR3`, `HET_VAR3_VAR4`,
//! `HOM_VAR2`, `HOM_VAR3`).
//!
//! The scheme (`GA4GHScheme`, or `GA4GHSchemeWithMissingAsHomRef` under
//! `MISSING_SITES_HOM_REF`) maps each (truth, call) cell to its contingency
//! values; a cell can hold two or three (`TP,FN`), `EMPTY`, or `NA` for a
//! cell the states cannot reach. Sensitivity is TP over TP plus FN over the
//! het, hom-var or all variant truth states; PPV is TP over TP plus FP over
//! the call states; specificity TN over TN plus FP; the two concordances
//! divide the cells whose states share a name by all counted cells (the
//! non-reference one over cells with a variant state on either side), with
//! `MISSING` cells skipped unless `MISSING_SITES_HOM_REF`, which also adds
//! the interval bases no record covered to the `MISSING`/`MISSING` cell.

use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};
use crate::variant_calling_metrics::output_path;
use crate::vcf::{self, Dictionary, Regions, Variant, VariantType, VcfReader};

pub const SUMMARY_EXTENSION: &str = "genotype_concordance_summary_metrics";
pub const DETAIL_EXTENSION: &str = "genotype_concordance_detail_metrics";
pub const CONTINGENCY_EXTENSION: &str = "genotype_concordance_contingency_metrics";

/// `GenotypeConcordanceStates.TruthState`, in declaration order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TruthState {
    Missing,
    HomRef,
    HetRefVar1,
    HetVar1Var2,
    HomVar1,
    NoCall,
    LowGq,
    LowDp,
    VcFiltered,
    GtFiltered,
    IsMixed,
}

pub const TRUTH_STATES: [TruthState; 11] = [
    TruthState::Missing,
    TruthState::HomRef,
    TruthState::HetRefVar1,
    TruthState::HetVar1Var2,
    TruthState::HomVar1,
    TruthState::NoCall,
    TruthState::LowGq,
    TruthState::LowDp,
    TruthState::VcFiltered,
    TruthState::GtFiltered,
    TruthState::IsMixed,
];

impl TruthState {
    pub fn as_str(self) -> &'static str {
        match self {
            TruthState::Missing => "MISSING",
            TruthState::HomRef => "HOM_REF",
            TruthState::HetRefVar1 => "HET_REF_VAR1",
            TruthState::HetVar1Var2 => "HET_VAR1_VAR2",
            TruthState::HomVar1 => "HOM_VAR1",
            TruthState::NoCall => "NO_CALL",
            TruthState::LowGq => "LOW_GQ",
            TruthState::LowDp => "LOW_DP",
            TruthState::VcFiltered => "VC_FILTERED",
            TruthState::GtFiltered => "GT_FILTERED",
            TruthState::IsMixed => "IS_MIXED",
        }
    }

    fn index(self) -> usize {
        self as usize
    }

    /// `TruthState.getHom`.
    fn hom(allele: usize) -> Result<Self> {
        match allele {
            0 => Ok(TruthState::HomRef),
            1 => Ok(TruthState::HomVar1),
            _ => Err(state_error()),
        }
    }

    /// `TruthState.getVar`.
    fn var(a: usize, b: usize) -> Result<Self> {
        match (a, b) {
            (0, 1) | (1, 0) => Ok(TruthState::HetRefVar1),
            (1, 2) | (2, 1) => Ok(TruthState::HetVar1Var2),
            _ => Err(state_error()),
        }
    }

    /// `getCode`: the code the call state must share for a concordant cell.
    fn code(self) -> &'static str {
        self.as_str()
    }
}

/// `GenotypeConcordanceStates.CallState`, in declaration order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CallState {
    Missing,
    HomRef,
    HetRefVar1,
    HetVar1Var2,
    HomVar1,
    HetRefVar2,
    HetRefVar3,
    HetVar1Var3,
    HetVar3Var4,
    HomVar2,
    HomVar3,
    NoCall,
    LowGq,
    LowDp,
    VcFiltered,
    GtFiltered,
    IsMixed,
}

pub const CALL_STATES: [CallState; 17] = [
    CallState::Missing,
    CallState::HomRef,
    CallState::HetRefVar1,
    CallState::HetVar1Var2,
    CallState::HomVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
    CallState::HomVar2,
    CallState::HomVar3,
    CallState::NoCall,
    CallState::LowGq,
    CallState::LowDp,
    CallState::VcFiltered,
    CallState::GtFiltered,
    CallState::IsMixed,
];

impl CallState {
    pub fn as_str(self) -> &'static str {
        match self {
            CallState::Missing => "MISSING",
            CallState::HomRef => "HOM_REF",
            CallState::HetRefVar1 => "HET_REF_VAR1",
            CallState::HetVar1Var2 => "HET_VAR1_VAR2",
            CallState::HomVar1 => "HOM_VAR1",
            CallState::HetRefVar2 => "HET_REF_VAR2",
            CallState::HetRefVar3 => "HET_REF_VAR3",
            CallState::HetVar1Var3 => "HET_VAR1_VAR3",
            CallState::HetVar3Var4 => "HET_VAR3_VAR4",
            CallState::HomVar2 => "HOM_VAR2",
            CallState::HomVar3 => "HOM_VAR3",
            CallState::NoCall => "NO_CALL",
            CallState::LowGq => "LOW_GQ",
            CallState::LowDp => "LOW_DP",
            CallState::VcFiltered => "VC_FILTERED",
            CallState::GtFiltered => "GT_FILTERED",
            CallState::IsMixed => "IS_MIXED",
        }
    }

    fn index(self) -> usize {
        self as usize
    }

    /// `CallState.getHom`.
    fn hom(allele: usize) -> Result<Self> {
        match allele {
            0 => Ok(CallState::HomRef),
            1 => Ok(CallState::HomVar1),
            2 => Ok(CallState::HomVar2),
            3 => Ok(CallState::HomVar3),
            _ => Err(state_error()),
        }
    }

    /// `CallState.getHet`: the lower index first; (2, 3) is
    /// `HET_VAR3_VAR4`, as Picard's special case has it.
    fn het(a: usize, b: usize) -> Result<Self> {
        let (a, b) = if a > b { (b, a) } else { (a, b) };
        match (a, b) {
            (0, 1) => Ok(CallState::HetRefVar1),
            (0, 2) => Ok(CallState::HetRefVar2),
            (0, 3) => Ok(CallState::HetRefVar3),
            (1, 2) => Ok(CallState::HetVar1Var2),
            (1, 3) => Ok(CallState::HetVar1Var3),
            (2, 3) | (3, 4) => Ok(CallState::HetVar3Var4),
            _ => Err(state_error()),
        }
    }

    /// `getCode`: the six allele states beyond `HOM_VAR1` share
    /// `INCOMPARABLE_CODE`, which no truth state has.
    fn code(self) -> &'static str {
        match self {
            CallState::HetRefVar2
            | CallState::HetRefVar3
            | CallState::HetVar1Var3
            | CallState::HetVar3Var4
            | CallState::HomVar2
            | CallState::HomVar3 => "INCOMPARABLE",
            other => other.as_str(),
        }
    }

    fn is_var(self) -> bool {
        matches!(
            self,
            CallState::HetRefVar1
                | CallState::HetRefVar2
                | CallState::HetRefVar3
                | CallState::HetVar1Var2
                | CallState::HetVar1Var3
                | CallState::HetVar3Var4
                | CallState::HomVar1
                | CallState::HomVar2
                | CallState::HomVar3
        )
    }
}

fn state_error() -> Error {
    Error::Metric("Shouldn't be here.".to_owned())
}

/// `ContingencyState`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Contingency {
    Tp,
    Fp,
    Tn,
    Fn,
    Na,
    Empty,
}

impl Contingency {
    fn as_str(self) -> &'static str {
        match self {
            Contingency::Tp => "TP",
            Contingency::Fp => "FP",
            Contingency::Tn => "TN",
            Contingency::Fn => "FN",
            Contingency::Na => "NA",
            Contingency::Empty => "EMPTY",
        }
    }
}

use Contingency::{Empty as E, Fn as FN, Fp as FP, Na as N, Tn as TN, Tp as TP};

const NA: &[Contingency] = &[N];
const EMPTY: &[Contingency] = &[E];
const TP_ONLY: &[Contingency] = &[TP];
const FP_ONLY: &[Contingency] = &[FP];
const TN_ONLY: &[Contingency] = &[TN];
const FN_ONLY: &[Contingency] = &[FN];
const TP_FN: &[Contingency] = &[TP, FN];
const TP_FP: &[Contingency] = &[TP, FP];
const TP_TN: &[Contingency] = &[TP, TN];
const FP_FN: &[Contingency] = &[FP, FN];
const FP_TN: &[Contingency] = &[FP, TN];
const FP_TN_FN: &[Contingency] = &[FP, TN, FN];
const TP_FP_FN: &[Contingency] = &[TP, FP, FN];
const TN_FN: &[Contingency] = &[TN, FN];

/// One scheme row: the cell for each truth state, in `TRUTH_STATES` order.
type Row = [&'static [Contingency]; 11];

/// `GA4GHScheme.initiateScheme`, one row per call state in `CALL_STATES`
/// order. Columns: MISSING, HOM_REF, HET_REF_VAR1, HET_VAR1_VAR2, HOM_VAR1,
/// NO_CALL, LOW_GQ, LOW_DP, VC_FILTERED, GT_FILTERED, IS_MIXED.
const GA4GH: [Row; 17] = [
    [
        NA, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [
        TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [
        FP_TN, FP_TN, TP_TN, TP_FN, TP_FN, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [
        FP_ONLY, FP_ONLY, TP_FP, TP_ONLY, TP_FP_FN, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [
        FP_ONLY, FP_ONLY, TP_FP, TP_FN, TP_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [NA, NA, FP_TN_FN, NA, FP_FN, NA, NA, NA, NA, NA, NA],
    [NA, NA, NA, FP_FN, NA, NA, NA, NA, NA, NA, NA],
    [NA, NA, NA, TP_FP_FN, NA, NA, NA, NA, NA, NA, NA],
    [
        FP_ONLY, FP_ONLY, FP_FN, FP_FN, FP_FN, NA, NA, NA, NA, NA, NA,
    ],
    [NA, NA, FP_FN, TP_FN, FP_FN, NA, NA, NA, NA, NA, NA],
    [NA, NA, NA, FP_FN, NA, NA, NA, NA, NA, NA, NA],
    [EMPTY; 11],
    [
        EMPTY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [
        EMPTY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [
        EMPTY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [
        EMPTY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
    ],
    [EMPTY; 11],
];

/// `GA4GHSchemeWithMissingAsHomRef`: the GA4GH table with a missing truth
/// read as hom-ref (`TN`) for a missing, filtered, low-GQ or low-DP call.
const GA4GH_MISSING_HOM_REF: [Row; 17] = {
    let mut rows = GA4GH;
    rows[0][0] = TN_ONLY;
    rows[12][0] = TN_ONLY;
    rows[13][0] = TN_ONLY;
    rows[14][0] = TN_ONLY;
    rows[15][0] = TN_ONLY;
    rows
};

/// The scheme in use.
#[derive(Clone, Copy, Debug)]
pub struct Scheme {
    rows: &'static [Row; 17],
}

impl Scheme {
    /// `GenotypeConcordanceSchemeFactory.getScheme`.
    pub fn new(missing_sites_hom_ref: bool) -> Self {
        Self {
            rows: if missing_sites_hom_ref {
                &GA4GH_MISSING_HOM_REF
            } else {
                &GA4GH
            },
        }
    }

    pub fn cell(&self, truth: TruthState, call: CallState) -> &'static [Contingency] {
        self.rows[call.index()][truth.index()]
    }

    /// `getContingencyStateString`.
    pub fn cell_text(&self, truth: TruthState, call: CallState) -> String {
        self.cell(truth, call)
            .iter()
            .map(|c| c.as_str())
            .collect::<Vec<_>>()
            .join(",")
    }
}

/// `GenotypeConcordanceCounts`: a count per (truth, call) cell, kept as
/// doubles like Picard's `Histogram`.
#[derive(Clone, Debug, Default)]
pub struct Counts {
    cells: [[f64; 17]; 11],
}

const HET_TRUTH: [TruthState; 2] = [TruthState::HetRefVar1, TruthState::HetVar1Var2];
const HOM_VAR_TRUTH: [TruthState; 1] = [TruthState::HomVar1];
const VAR_TRUTH: [TruthState; 5] = [
    TruthState::HetRefVar1,
    TruthState::HetVar1Var2,
    TruthState::HomVar1,
    TruthState::HomRef,
    TruthState::Missing,
];
const HET_CALL: [CallState; 6] = [
    CallState::HetRefVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var2,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
];
const HOM_VAR_CALL: [CallState; 3] = [CallState::HomVar1, CallState::HomVar2, CallState::HomVar3];
const VAR_CALL: [CallState; 9] = [
    CallState::HetRefVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var2,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
    CallState::HomVar1,
    CallState::HomVar2,
    CallState::HomVar3,
];

impl Counts {
    pub fn increment(&mut self, truth: TruthState, call: CallState, by: f64) {
        self.cells[truth.index()][call.index()] += by;
    }

    /// `getCount`: the bin value as a long.
    pub fn count(&self, truth: TruthState, call: CallState) -> i64 {
        self.cells[truth.index()][call.index()] as i64
    }

    /// `getCounterSize`: every count summed.
    pub fn size(&self) -> f64 {
        self.cells.iter().flatten().sum()
    }

    /// `getSensitivity` (`TP` over `TP` plus `FN`), `Ppv` (`TP` over `TP`
    /// plus `FP`) and `getSpecificity` (`TN` over `TN` plus `FP`), each a
    /// double division that is NaN over no cells.
    fn rate(
        &self,
        scheme: &Scheme,
        truths: &[TruthState],
        calls: &[CallState],
        hit: Contingency,
        miss: Contingency,
    ) -> f64 {
        let (mut numerator, mut denominator) = (0.0, 0.0);
        for &truth in truths {
            for &call in calls {
                let count = self.count(truth, call) as f64;
                for &c in scheme.cell(truth, call) {
                    if c == hit {
                        numerator += count;
                        denominator += count;
                    } else if c == miss {
                        denominator += count;
                    }
                }
            }
        }
        numerator / denominator
    }

    fn sensitivity(&self, scheme: &Scheme, truths: &[TruthState]) -> f64 {
        self.rate(scheme, truths, &CALL_STATES, TP, FN)
    }

    fn ppv(&self, scheme: &Scheme, calls: &[CallState]) -> f64 {
        self.rate(scheme, &TRUTH_STATES, calls, TP, FP)
    }

    fn specificity(&self, scheme: &Scheme, truths: &[TruthState]) -> f64 {
        self.rate(scheme, truths, &CALL_STATES, TN, FP)
    }

    /// `calculateGenotypeConcordanceUtil`.
    fn concordance(&self, missing_sites: bool, include_hom_ref: bool) -> f64 {
        let (mut numerator, mut denominator) = (0.0, 0.0);
        for truth in TRUTH_STATES {
            for call in CALL_STATES {
                let missing = truth == TruthState::Missing || call == CallState::Missing;
                if !missing_sites && missing {
                    continue;
                }
                let is_var = matches!(
                    truth,
                    TruthState::HomVar1 | TruthState::HetRefVar1 | TruthState::HetVar1Var2
                ) || call.is_var();
                if include_hom_ref || is_var {
                    let count = self.count(truth, call) as f64;
                    if truth.code() == call.code() {
                        numerator += count;
                    }
                    denominator += count;
                }
            }
        }
        if denominator > 0.0 {
            numerator / denominator
        } else {
            f64::NAN
        }
    }

    /// `getContingencyStateCounts`: TP, TN, FP, FN and EMPTY totals.
    fn contingency_counts(&self, scheme: &Scheme) -> [i64; 6] {
        let mut totals = [0i64; 6];
        for truth in TRUTH_STATES {
            for call in CALL_STATES {
                let count = self.count(truth, call);
                for &c in scheme.cell(truth, call) {
                    totals[c as usize] += count;
                }
            }
        }
        totals
    }
}

/// One `GenotypeConcordanceSummaryMetrics` row.
#[derive(Clone, Debug, PartialEq)]
pub struct SummaryRow {
    pub variant_type: &'static str,
    pub truth_sample: String,
    pub call_sample: String,
    pub het_sensitivity: f64,
    pub het_ppv: f64,
    pub homvar_sensitivity: f64,
    pub homvar_ppv: f64,
    pub var_sensitivity: f64,
    pub var_ppv: f64,
    pub var_specificity: f64,
    pub genotype_concordance: f64,
    pub non_ref_genotype_concordance: f64,
}

impl SummaryRow {
    fn new(
        variant_type: &'static str,
        counts: &Counts,
        scheme: &Scheme,
        truth_sample: &str,
        call_sample: &str,
        missing_sites: bool,
    ) -> Self {
        Self {
            variant_type,
            truth_sample: truth_sample.to_owned(),
            call_sample: call_sample.to_owned(),
            het_sensitivity: counts.sensitivity(scheme, &HET_TRUTH),
            het_ppv: counts.ppv(scheme, &HET_CALL),
            homvar_sensitivity: counts.sensitivity(scheme, &HOM_VAR_TRUTH),
            homvar_ppv: counts.ppv(scheme, &HOM_VAR_CALL),
            var_sensitivity: counts.sensitivity(scheme, &VAR_TRUTH),
            var_ppv: counts.ppv(scheme, &VAR_CALL),
            var_specificity: counts.specificity(scheme, &VAR_TRUTH),
            genotype_concordance: counts.concordance(missing_sites, true),
            non_ref_genotype_concordance: counts.concordance(missing_sites, false),
        }
    }
}

impl MetricRow for SummaryRow {
    const CLASS: &'static str = "picard.vcf.GenotypeConcordanceSummaryMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "VARIANT_TYPE",
        "TRUTH_SAMPLE",
        "CALL_SAMPLE",
        "HET_SENSITIVITY",
        "HET_PPV",
        "HET_SPECIFICITY",
        "HOMVAR_SENSITIVITY",
        "HOMVAR_PPV",
        "HOMVAR_SPECIFICITY",
        "VAR_SENSITIVITY",
        "VAR_PPV",
        "VAR_SPECIFICITY",
        "GENOTYPE_CONCORDANCE",
        "NON_REF_GENOTYPE_CONCORDANCE",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.variant_type.into(),
            self.truth_sample.clone().into(),
            self.call_sample.clone().into(),
            self.het_sensitivity.into(),
            self.het_ppv.into(),
            // `HET_SPECIFICITY` and `HOMVAR_SPECIFICITY` are NaN by definition.
            f64::NAN.into(),
            self.homvar_sensitivity.into(),
            self.homvar_ppv.into(),
            f64::NAN.into(),
            self.var_sensitivity.into(),
            self.var_ppv.into(),
            self.var_specificity.into(),
            self.genotype_concordance.into(),
            self.non_ref_genotype_concordance.into(),
        ]
    }
}

/// One `GenotypeConcordanceDetailMetrics` row: a (truth, call) cell.
#[derive(Clone, Debug, PartialEq)]
pub struct DetailRow {
    pub variant_type: &'static str,
    pub truth_sample: String,
    pub call_sample: String,
    pub truth_state: TruthState,
    pub call_state: CallState,
    pub count: i64,
    pub contingency_values: String,
}

impl MetricRow for DetailRow {
    const CLASS: &'static str = "picard.vcf.GenotypeConcordanceDetailMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "VARIANT_TYPE",
        "TRUTH_SAMPLE",
        "CALL_SAMPLE",
        "TRUTH_STATE",
        "CALL_STATE",
        "COUNT",
        "CONTINGENCY_VALUES",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.variant_type.into(),
            self.truth_sample.clone().into(),
            self.call_sample.clone().into(),
            self.truth_state.as_str().into(),
            self.call_state.as_str().into(),
            self.count.into(),
            self.contingency_values.clone().into(),
        ]
    }
}

/// One `GenotypeConcordanceContingencyMetrics` row.
#[derive(Clone, Debug, PartialEq)]
pub struct ContingencyRow {
    pub variant_type: &'static str,
    pub truth_sample: String,
    pub call_sample: String,
    pub tp_count: i64,
    pub tn_count: i64,
    pub fp_count: i64,
    pub fn_count: i64,
    pub empty_count: i64,
}

impl MetricRow for ContingencyRow {
    const CLASS: &'static str = "picard.vcf.GenotypeConcordanceContingencyMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "VARIANT_TYPE",
        "TRUTH_SAMPLE",
        "CALL_SAMPLE",
        "TP_COUNT",
        "TN_COUNT",
        "FP_COUNT",
        "FN_COUNT",
        "EMPTY_COUNT",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.variant_type.into(),
            self.truth_sample.clone().into(),
            self.call_sample.clone().into(),
            self.tp_count.into(),
            self.tn_count.into(),
            self.fp_count.into(),
            self.fn_count.into(),
            self.empty_count.into(),
        ]
    }
}

/// The normalized alleles of one pair: the ordered allele list and each
/// side's two alleles as positions in it (`GenotypeConcordance.Alleles`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alleles {
    pub all: Vec<String>,
    pub truth: Option<(usize, usize)>,
    pub call: Option<(usize, usize)>,
}

/// `spliceOrAppendString`: `insert` goes in at `at`, or at the end of a
/// shorter string; the spanning-deletion allele stays `*`.
fn splice(destination: &str, insert: &str, at: usize) -> String {
    if destination == "*" {
        return destination.to_owned();
    }
    if at <= destination.len() {
        format!("{}{insert}{}", &destination[..at], &destination[at..])
    } else {
        format!("{destination}{insert}")
    }
}

/// The two allele strings of a context's only genotype (`getBaseString`:
/// `.` for a no-call).
fn genotype_alleles(ctx: &Variant) -> Result<(String, String)> {
    let g = &ctx.genotypes[0];
    if g.alleles.len() != 2 {
        return Err(Error::Metric(format!(
            "Genotype for Variant Context: {}:{} does not have exactly 2 alleles",
            ctx.contig, ctx.start
        )));
    }
    let text = |a: Option<usize>| a.map_or(".".to_owned(), |i| ctx.alleles[i].bases.clone());
    Ok((text(g.alleles[0]), text(g.alleles[1])))
}

/// `getStringSuffix`.
fn suffix(longer: &str, shorter: &str, truth: &Variant, call: &Variant) -> Result<String> {
    longer
        .strip_prefix(shorter)
        .map(str::to_owned)
        .ok_or_else(|| ref_mismatch(truth, call))
}

fn ref_mismatch(truth: &Variant, call: &Variant) -> Error {
    Error::Metric(format!(
        "Ref alleles mismatch between: {}:{} {} and {}:{} {}",
        truth.contig,
        truth.start,
        truth.reference().bases,
        call.contig,
        call.start,
        call.reference().bases
    ))
}

/// `normalizeAlleles`. A truth context that is mixed or filtered, or a
/// call context that is mixed or (unless `ignore_filtered`) filtered,
/// contributes no alleles.
pub fn normalize_alleles(
    truth: Option<&Variant>,
    call: Option<&Variant>,
    ignore_filtered: bool,
) -> Result<Alleles> {
    let truth = truth.filter(|t| t.kind() != VariantType::Mixed && !t.filtered);
    let call = call.filter(|c| c.kind() != VariantType::Mixed && (ignore_filtered || !c.filtered));
    let mut truth_ref = truth.map(|t| t.reference().bases.clone());
    let mut call_ref = call.map(|c| c.reference().bases.clone());
    let mut truth_alleles = truth.map(genotype_alleles).transpose()?;
    let mut call_alleles = call.map(genotype_alleles).transpose()?;
    if let (Some(tr), Some(cr), Some(t), Some(c)) = (&truth_ref, &call_ref, truth, call)
        && tr != cr
    {
        let splice_pair = |pair: &mut Option<(String, String)>, suffix: &str, at: usize| {
            if let Some((a, b)) = pair {
                if a != "." {
                    *a = splice(a, suffix, at);
                }
                if b != "." {
                    *b = splice(b, suffix, at);
                }
            }
        };
        if tr.len() < cr.len() {
            let s = suffix(cr, tr, t, c)?;
            splice_pair(&mut truth_alleles, &s, tr.len());
            truth_ref = Some(format!("{tr}{s}"));
        } else if tr.len() > cr.len() {
            let s = suffix(tr, cr, t, c)?;
            splice_pair(&mut call_alleles, &s, cr.len());
            call_ref = Some(format!("{cr}{s}"));
        } else {
            return Err(ref_mismatch(t, c));
        }
    }
    let mut all: Vec<String> = Vec::new();
    let smart_add = |all: &mut Vec<String>, allele: &str| {
        if !all.iter().any(|a| a == allele) {
            all.push(allele.to_owned());
        }
    };
    if truth_alleles.is_some() || call_alleles.is_some() {
        let zeroth = if truth_alleles.is_some() {
            truth_ref.as_deref()
        } else {
            call_ref.as_deref()
        };
        smart_add(&mut all, zeroth.unwrap_or_default());
    }
    if let Some((a, b)) = &truth_alleles {
        smart_add(&mut all, a);
        smart_add(&mut all, b);
    }
    if let Some((a, b)) = &call_alleles {
        let index = |all: &[String], x: &str| all.iter().position(|v| v == x);
        let past_one = |i: Option<usize>| i.is_some_and(|i| i > 1);
        if past_one(index(&all, a)) || past_one(index(&all, b)) {
            // The call carries a truth allele from position 2: re-add the
            // truth alleles in swapped order so it comes first.
            all.remove(2);
            all.remove(1);
            if let Some((ta, tb)) = &truth_alleles {
                smart_add(&mut all, tb);
                smart_add(&mut all, ta);
            }
        }
        smart_add(&mut all, a);
        smart_add(&mut all, b);
    }
    let position = |pair: &Option<(String, String)>| {
        pair.as_ref().map(|(a, b)| {
            let find = |x: &str| all.iter().position(|v| v == x).expect("added above");
            (find(a), find(b))
        })
    };
    Ok(Alleles {
        truth: position(&truth_alleles),
        call: position(&call_alleles),
        all,
    })
}

/// `getStateCode` as the states it maps to: `None` when the genotype's
/// alleles decide.
enum Code {
    IsMixed,
    VcFiltered,
    NoCall,
    GtFiltered,
    LowGq,
    LowDp,
}

fn state_code(ctx: Option<&Variant>, min_gq: i64, min_dp: i64) -> Option<Code> {
    let ctx = ctx?;
    if ctx.kind() == VariantType::Mixed {
        return Some(Code::IsMixed);
    }
    if ctx.filtered {
        return Some(Code::VcFiltered);
    }
    let g = &ctx.genotypes[0];
    if g.is_no_call() {
        Some(Code::NoCall)
    } else if g.filtered {
        Some(Code::GtFiltered)
    } else if g.gq != -1 && g.gq < min_gq {
        Some(Code::LowGq)
    } else if g.dp != -1 && g.dp < min_dp {
        Some(Code::LowDp)
    } else if g.is_mixed() {
        Some(Code::NoCall)
    } else {
        None
    }
}

/// `determineState`: the truth and call states of one pair of
/// sub-contexts.
pub fn determine_state(
    truth: Option<&Variant>,
    call: Option<&Variant>,
    min_gq: i64,
    min_dp: i64,
    ignore_filtered: bool,
) -> Result<(TruthState, CallState)> {
    let truth_state = match truth {
        None => Some(TruthState::Missing),
        Some(_) => state_code(truth, min_gq, min_dp).map(|c| match c {
            Code::IsMixed => TruthState::IsMixed,
            Code::VcFiltered => TruthState::VcFiltered,
            Code::NoCall => TruthState::NoCall,
            Code::GtFiltered => TruthState::GtFiltered,
            Code::LowGq => TruthState::LowGq,
            Code::LowDp => TruthState::LowDp,
        }),
    };
    let call_state = match call {
        None => Some(CallState::Missing),
        Some(_) => state_code(call, min_gq, min_dp)
            .filter(|c| !(ignore_filtered && matches!(c, Code::VcFiltered)))
            .map(|c| match c {
                Code::IsMixed => CallState::IsMixed,
                Code::VcFiltered => CallState::VcFiltered,
                Code::NoCall => CallState::NoCall,
                Code::GtFiltered => CallState::GtFiltered,
                Code::LowGq => CallState::LowGq,
                Code::LowDp => CallState::LowDp,
            }),
    };
    let alleles = normalize_alleles(
        if truth_state.is_none() { truth } else { None },
        if call_state.is_none() { call } else { None },
        ignore_filtered,
    )?;
    let truth_state = match truth_state {
        Some(s) => s,
        None => {
            let (a, b) = alleles.truth.ok_or_else(state_error)?;
            if a == b {
                TruthState::hom(a)?
            } else {
                TruthState::var(a, b)?
            }
        }
    };
    let call_state = match call_state {
        Some(s) => s,
        None => {
            let (a, b) = alleles.call.ok_or_else(state_error)?;
            if a == b {
                CallState::hom(a)?
            } else {
                CallState::het(a, b)?
            }
        }
    };
    Ok((truth_state, call_state))
}

/// The SNP and indel counters.
#[derive(Clone, Debug, Default)]
pub struct Counters {
    pub snp: Counts,
    pub indel: Counts,
}

impl Counters {
    /// `classifyVariants`: true when the pair was counted.
    pub fn classify(
        &mut self,
        truth: Option<&Variant>,
        call: Option<&Variant>,
        min_gq: i64,
        min_dp: i64,
        ignore_filtered: bool,
    ) -> Result<bool> {
        use VariantType::{Indel, Mixed, NoVariation, Snp};
        let truth_type = truth.map_or(NoVariation, Variant::kind);
        let call_type = call.map_or(NoVariation, Variant::kind);
        let (t, c) = determine_state(truth, call, min_gq, min_dp, ignore_filtered)?;
        let counter = match (truth_type, call_type) {
            (Snp, Snp | Mixed | NoVariation) => &mut self.snp,
            (Indel, Indel | Mixed | NoVariation) => &mut self.indel,
            (Mixed | NoVariation, Snp) => &mut self.snp,
            (Mixed | NoVariation, Indel) => &mut self.indel,
            _ => return Ok(false),
        };
        counter.increment(t, c, 1.0);
        Ok(true)
    }
}

/// `GenotypeConcordance` in the tool registry. `OUTPUT_VCF` is not
/// implemented and so stays unknown.
pub const TOOL: Tool = Tool {
    name: "GenotypeConcordance",
    summary: "the truth and call VCFs compared by position. OUTPUT is a prefix: \
              PREFIX.genotype_concordance_summary_metrics, _detail_metrics and \
              _contingency_metrics; Python returns one row per file: file and rows.",
    options: &[
        Opt::path(&["TRUTH_VCF", "TV"], "The truth VCF.").required(),
        crate::registry::REFERENCE_NOT_DECODED,
        Opt::path(&["CALL_VCF", "CV"], "The call VCF.").required(),
        Opt::output(&["OUTPUT", "O"])
            .with_help("The output prefix of the three files.")
            .required(),
        Opt::text(
            &["TRUTH_SAMPLE", "TS"],
            "Sample of the truth VCF; required with more than one.",
        ),
        Opt::text(
            &["CALL_SAMPLE", "CS"],
            "Sample of the call VCF; required with more than one.",
        ),
        Opt::path(
            &["INTERVALS"],
            "An interval list; repeat for several. Both VCFs must be indexed.",
        )
        .repeated(),
        Opt::bool(
            &["INTERSECT_INTERVALS"],
            "Use the intersection of INTERVALS, false their union (default true).",
        ),
        Opt::int(
            &["MIN_GQ"],
            "Genotypes below this GQ count as no call (default 0).",
        ),
        Opt::int(
            &["MIN_DP"],
            "Genotypes below this depth count as no call (default 0).",
        ),
        Opt::bool(
            &["OUTPUT_ALL_ROWS"],
            "Write rows whose counts are all zero (default false).",
        ),
        Opt::bool(
            &["USE_VCF_INDEX"],
            "Read the VCFs through their index (default false).",
        ),
        Opt::bool(
            &["MISSING_SITES_HOM_REF", "MISSING_HOM"],
            "Sites absent from a VCF are hom-ref; needs INTERVALS (default false).",
        ),
        Opt::bool(
            &["IGNORE_FILTER_STATUS"],
            "Count filtered sites (default false).",
        ),
        Opt::rejected(&["INPUT", "I"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

pub struct Options {
    pub truth_vcf: PathBuf,
    pub call_vcf: PathBuf,
    pub prefix: PathBuf,
    pub truth_sample: Option<String>,
    pub call_sample: Option<String>,
    pub intervals: Vec<PathBuf>,
    pub intersect_intervals: bool,
    pub min_gq: i64,
    pub min_dp: i64,
    pub output_all_rows: bool,
    pub use_vcf_index: bool,
    pub missing_sites_hom_ref: bool,
    pub ignore_filter_status: bool,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let required = |name: &str, value: Option<&str>| {
            value
                .map(PathBuf::from)
                .ok_or_else(|| Error::Usage(format!("Option '{name}' is required.")))
        };
        let opts = Self {
            truth_vcf: required("TRUTH_VCF", inv.get("TRUTH_VCF"))?,
            call_vcf: required("CALL_VCF", inv.get("CALL_VCF"))?,
            prefix: required("OUTPUT", inv.output.as_deref().and_then(Path::to_str))?,
            truth_sample: inv.get("TRUTH_SAMPLE").map(str::to_owned),
            call_sample: inv.get("CALL_SAMPLE").map(str::to_owned),
            intervals: inv.all("INTERVALS").iter().map(PathBuf::from).collect(),
            intersect_intervals: inv.flag("INTERSECT_INTERVALS", true)?,
            min_gq: inv.parse("MIN_GQ", 0)?,
            min_dp: inv.parse("MIN_DP", 0)?,
            output_all_rows: inv.flag("OUTPUT_ALL_ROWS", false)?,
            use_vcf_index: inv.flag("USE_VCF_INDEX", false)?,
            missing_sites_hom_ref: inv.flag("MISSING_SITES_HOM_REF", false)?,
            ignore_filter_status: inv.flag("IGNORE_FILTER_STATUS", false)?,
        };
        // `customCommandLineValidation`.
        let using_intervals = !opts.intervals.is_empty();
        let mut errors = Vec::new();
        if opts.use_vcf_index || using_intervals {
            for (path, side) in [(&opts.truth_vcf, "TRUTH"), (&opts.call_vcf, "CALL")] {
                if !vcf::index_exists(path) {
                    errors.push(format!(
                        "The index file was not found for the {side} VCF.  Note that if intervals are specified, the VCF files must be indexed."
                    ));
                }
            }
        }
        if opts.missing_sites_hom_ref && !using_intervals {
            errors.push(
                "You cannot use the MISSING_HOM option without also supplying an interval list over which missing sites are considered confident homozygous reference calls."
                    .to_owned(),
            );
        }
        if !errors.is_empty() {
            return Err(Error::Usage(errors.join("\n")));
        }
        Ok(opts)
    }
}

/// The sample to compare: the named one, or the only one.
fn choose_sample(
    reader: &VcfReader,
    wanted: Option<&str>,
    option: &str,
    side: &str,
) -> Result<String> {
    let samples = reader.samples();
    let name = match wanted {
        Some(name) => name.to_owned(),
        None if samples.len() > 1 => {
            return Err(Error::Usage(format!(
                "{option} is required when the {side} has more than one sample"
            )));
        }
        None => samples.first().cloned().ok_or_else(|| {
            Error::Usage(format!(
                "{side} {} has no sample (Picard fails with an IndexOutOfBoundsException)",
                reader.path().display()
            ))
        })?,
    };
    if !samples.contains(&name) {
        return Err(Error::Usage(format!(
            "File {} does not contain genotypes for sample {name}",
            reader.path().display()
        )));
    }
    Ok(name)
}

/// A VCF walked record by record with one record of look-ahead, limited to
/// the regions when given (`ByIntervalListVariantContextIterator`).
struct Stream {
    reader: VcfReader,
    regions: Option<Regions>,
    peeked: Option<Variant>,
}

impl Stream {
    fn new(reader: VcfReader, regions: Option<Regions>) -> Result<Self> {
        let mut stream = Self {
            reader,
            regions,
            peeked: None,
        };
        stream.advance()?;
        Ok(stream)
    }

    fn advance(&mut self) -> Result<()> {
        loop {
            let next = self.reader.next_variant()?;
            match (&self.regions, &next) {
                (Some(regions), Some(v)) if !regions.overlaps(&v.contig, v.start, v.end) => {
                    continue;
                }
                _ => {
                    self.peeked = next;
                    return Ok(());
                }
            }
        }
    }

    fn take(&mut self) -> Result<Option<Variant>> {
        let taken = self.peeked.take();
        if taken.is_some() {
            self.advance()?;
        }
        Ok(taken)
    }
}

/// `VariantContextComparator`: the contig's index in the truth dictionary,
/// then the start.
fn position_key(v: &Variant, dictionary: &Dictionary) -> Result<(usize, i64)> {
    let index = dictionary
        .iter()
        .position(|c| c.name == v.contig)
        .ok_or_else(|| {
            Error::Metric(format!(
                "contig {} is not in the truth sequence dictionary (Picard fails with a NullPointerException)",
                v.contig
            ))
        })?;
    Ok((index, v.start))
}

/// The intervals of `INTERVALS`: the first list's dictionary, the lists
/// intersected or united, then `uniqued`.
fn load_intervals(opts: &Options) -> Result<Option<Regions>> {
    let mut combined: Option<Regions> = None;
    for path in &opts.intervals {
        let list = Regions::read(path)?;
        combined = Some(match combined {
            None => list,
            Some(current) if opts.intersect_intervals => current.intersection(&list),
            Some(current) => current.union(&list),
        });
    }
    Ok(combined.map(|r| r.uniqued()))
}

/// `doWork`: the three metrics files' rows.
pub fn run(opts: &Options) -> Result<(Vec<SummaryRow>, Vec<DetailRow>, Vec<ContingencyRow>)> {
    let scheme = Scheme::new(opts.missing_sites_hom_ref);
    let intervals = load_intervals(opts)?;
    let truth_reader = VcfReader::open(&opts.truth_vcf)?;
    let call_reader = VcfReader::open(&opts.call_vcf)?;
    let truth_sample = choose_sample(
        &truth_reader,
        opts.truth_sample.as_deref(),
        "TRUTH_SAMPLE",
        "TRUTH_VCF",
    )?;
    let call_sample = choose_sample(
        &call_reader,
        opts.call_sample.as_deref(),
        "CALL_SAMPLE",
        "CALL_VCF",
    )?;
    let truth_dictionary = truth_reader.dictionary();
    vcf::assert_dictionaries_equal(
        truth_dictionary.as_ref(),
        &opts.truth_vcf.display().to_string(),
        call_reader.dictionary().as_ref(),
        &opts.call_vcf.display().to_string(),
    )?;
    if let Some(intervals) = &intervals {
        vcf::assert_dictionaries_equal(
            Some(&intervals.dictionary),
            "provided intervals",
            truth_dictionary.as_ref(),
            &opts.truth_vcf.display().to_string(),
        )?;
    }
    let dictionary = truth_dictionary.ok_or_else(|| Error::Intervals {
        path: opts.truth_vcf.clone(),
        message:
            "no ##contig line (Picard's VariantContextComparator fails with a NullPointerException)"
                .to_owned(),
    })?;
    let truth_index = truth_reader
        .samples()
        .iter()
        .position(|s| *s == truth_sample)
        .expect("checked");
    let call_index = call_reader
        .samples()
        .iter()
        .position(|s| *s == call_sample)
        .expect("checked");
    let mut truth = Stream::new(truth_reader, intervals.clone())?;
    let mut call = Stream::new(call_reader, intervals.clone())?;
    let mut counters = Counters::default();
    loop {
        let (take_truth, take_call) = match (&truth.peeked, &call.peeked) {
            (None, None) => break,
            (Some(_), None) => (true, false),
            (None, Some(_)) => (false, true),
            (Some(t), Some(c)) => {
                let ordering = position_key(t, &dictionary)?.cmp(&position_key(c, &dictionary)?);
                (!ordering.is_gt(), !ordering.is_lt())
            }
        };
        let left = if take_truth {
            truth.take()?.map(|v| v.sub_context(truth_index))
        } else {
            None
        };
        let right = if take_call {
            call.take()?.map(|v| v.sub_context(call_index))
        } else {
            None
        };
        counters.classify(
            left.as_ref(),
            right.as_ref(),
            opts.min_gq,
            opts.min_dp,
            opts.ignore_filter_status,
        )?;
    }
    if opts.missing_sites_hom_ref {
        let base_count = intervals.as_ref().map_or_else(
            || {
                dictionary
                    .iter()
                    .map(|c| c.length.unwrap_or(0))
                    .sum::<i64>()
            },
            Regions::base_count,
        );
        for counter in [&mut counters.snp, &mut counters.indel] {
            let missing = base_count as f64 - counter.size();
            counter.increment(TruthState::Missing, CallState::Missing, missing);
        }
    }
    let types = [("SNP", &counters.snp), ("INDEL", &counters.indel)];
    let summary = types
        .iter()
        .map(|(name, counts)| {
            SummaryRow::new(
                name,
                counts,
                &scheme,
                &truth_sample,
                &call_sample,
                opts.missing_sites_hom_ref,
            )
        })
        .collect();
    let mut detail = Vec::new();
    for (name, counts) in &types {
        for truth_state in TRUTH_STATES {
            for call_state in CALL_STATES {
                let count = counts.count(truth_state, call_state);
                if count > 0 || opts.output_all_rows {
                    detail.push(DetailRow {
                        variant_type: name,
                        truth_sample: truth_sample.clone(),
                        call_sample: call_sample.clone(),
                        truth_state,
                        call_state,
                        count,
                        contingency_values: scheme.cell_text(truth_state, call_state),
                    });
                }
            }
        }
    }
    let contingency = types
        .iter()
        .map(|(name, counts)| {
            let totals = counts.contingency_counts(&scheme);
            ContingencyRow {
                variant_type: name,
                truth_sample: truth_sample.clone(),
                call_sample: call_sample.clone(),
                tp_count: totals[Contingency::Tp as usize],
                tn_count: totals[Contingency::Tn as usize],
                fp_count: totals[Contingency::Fp as usize],
                fn_count: totals[Contingency::Fn as usize],
                empty_count: totals[Contingency::Empty as usize],
            }
        })
        .collect();
    Ok((summary, detail, contingency))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let (summary, detail, contingency) = run(&opts)?;
    Ok(ToolOutput {
        main: None,
        extra: vec![
            (
                output_path(&opts.prefix, SUMMARY_EXTENSION),
                MetricsOutput::from_rows(&summary, None),
            ),
            (
                output_path(&opts.prefix, DETAIL_EXTENSION),
                MetricsOutput::from_rows(&detail, None),
            ),
            (
                output_path(&opts.prefix, CONTINGENCY_EXTENSION),
                MetricsOutput::from_rows(&contingency, None),
            ),
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{Allele, Genotype};

    /// A one-sample record already in sub-context form: `gt` indexes
    /// `alleles`, the reference first.
    fn ctx(alleles: &[&str], gt: &[Option<usize>], filtered: bool, gq: i64, dp: i64) -> Variant {
        Variant {
            contig: "c1".to_owned(),
            start: 5,
            end: 5 + alleles[0].len() as i64 - 1,
            alleles: alleles
                .iter()
                .map(|a| Allele {
                    bases: (*a).to_owned(),
                    symbolic: false,
                })
                .collect(),
            filtered,
            genotypes: vec![Genotype {
                alleles: gt.to_vec(),
                gq,
                dp,
                ad: None,
                filtered: false,
            }],
        }
    }

    fn state(truth: Option<&Variant>, call: Option<&Variant>) -> (TruthState, CallState) {
        determine_state(truth, call, 0, 0, false).expect("state")
    }

    #[test]
    fn states_follow_allele_positions() {
        let het = ctx(&["A", "G"], &[Some(0), Some(1)], false, 30, 10);
        let hom = ctx(&["A", "G"], &[Some(1), Some(1)], false, 30, 10);
        let het2 = ctx(&["A", "G", "T"], &[Some(1), Some(2)], false, 30, 10);
        let other = ctx(&["A", "T"], &[Some(0), Some(1)], false, 30, 10);
        let third = ctx(&["A", "C"], &[Some(1), Some(1)], false, 30, 10);
        let two_new = ctx(&["A", "T", "C"], &[Some(1), Some(2)], false, 30, 10);
        assert_eq!(
            state(Some(&het), Some(&het)),
            (TruthState::HetRefVar1, CallState::HetRefVar1)
        );
        assert_eq!(
            state(Some(&het), None),
            (TruthState::HetRefVar1, CallState::Missing)
        );
        assert_eq!(
            state(None, Some(&hom)),
            (TruthState::Missing, CallState::HomVar1)
        );
        assert_eq!(
            state(Some(&het), Some(&other)),
            (TruthState::HetRefVar1, CallState::HetRefVar2)
        );
        assert_eq!(
            state(Some(&hom), Some(&two_new)),
            (TruthState::HomVar1, CallState::HetVar3Var4)
        );
        assert_eq!(
            state(Some(&het2), Some(&third)),
            (TruthState::HetVar1Var2, CallState::HomVar3)
        );
        // The call's T sits at position 2 of [A, G, T]: the truth alleles
        // are re-added as T, G and the call T/T becomes HOM_VAR1.
        let tt = ctx(&["A", "T"], &[Some(1), Some(1)], false, 30, 10);
        assert_eq!(
            state(Some(&het2), Some(&tt)),
            (TruthState::HetVar1Var2, CallState::HomVar1)
        );
        let ref_t = ctx(&["A", "T"], &[Some(0), Some(1)], false, 30, 10);
        assert_eq!(
            state(Some(&het2), Some(&ref_t)),
            (TruthState::HetVar1Var2, CallState::HetRefVar1)
        );
        let tc = ctx(&["A", "T", "C"], &[Some(1), Some(2)], false, 30, 10);
        assert_eq!(
            state(Some(&het2), Some(&tc)),
            (TruthState::HetVar1Var2, CallState::HetVar1Var3)
        );
        let gc = ctx(&["A", "G", "C"], &[Some(1), Some(2)], false, 30, 10);
        assert_eq!(
            state(Some(&het2), Some(&gc)),
            (TruthState::HetVar1Var2, CallState::HetVar1Var3)
        );
        let ref_c = ctx(&["A", "C"], &[Some(0), Some(1)], false, 30, 10);
        assert_eq!(
            state(Some(&het2), Some(&ref_c)),
            (TruthState::HetVar1Var2, CallState::HetRefVar3)
        );
    }

    #[test]
    fn quality_filter_and_no_call_states() {
        let het = ctx(&["A", "G"], &[Some(0), Some(1)], false, 30, 10);
        let low = ctx(&["A", "G"], &[Some(0), Some(1)], false, 5, 3);
        let no_gq = ctx(&["A", "G"], &[Some(0), Some(1)], false, -1, -1);
        let filtered = ctx(&["A", "G"], &[Some(0), Some(1)], true, 30, 10);
        let mixed_gt = ctx(&["A", "G"], &[None, Some(1)], false, 30, 10);
        let no_call = ctx(&["A"], &[None, None], false, 30, 10);
        let mixed = ctx(&["A", "G", "AT"], &[Some(1), Some(2)], false, 30, 10);
        let by = |t: &Variant, c: &Variant, gq, dp, ignore| {
            determine_state(Some(t), Some(c), gq, dp, ignore).expect("state")
        };
        assert_eq!(by(&low, &het, 20, 10, false).0, TruthState::LowGq);
        assert_eq!(by(&het, &low, 1, 10, false).1, CallState::LowDp);
        assert_eq!(
            by(&no_gq, &no_gq, 20, 10, false),
            (TruthState::HetRefVar1, CallState::HetRefVar1)
        );
        assert_eq!(
            by(&filtered, &filtered, 0, 0, false),
            (TruthState::VcFiltered, CallState::VcFiltered)
        );
        assert_eq!(
            by(&filtered, &filtered, 0, 0, true),
            (TruthState::VcFiltered, CallState::HetRefVar1)
        );
        assert_eq!(
            by(&mixed_gt, &no_call, 0, 0, false),
            (TruthState::NoCall, CallState::NoCall)
        );
        assert_eq!(
            by(&mixed, &het, 0, 0, false),
            (TruthState::IsMixed, CallState::HetRefVar1)
        );
        let mut gt_filtered = het.clone();
        gt_filtered.genotypes[0].filtered = true;
        assert_eq!(by(&het, &gt_filtered, 0, 0, false).1, CallState::GtFiltered);
    }

    #[test]
    fn indel_references_are_normalized() {
        let truth = ctx(&["AT", "A"], &[Some(0), Some(1)], false, 30, 10);
        let call = ctx(&["ATC", "AC"], &[Some(0), Some(1)], false, 30, 10);
        let alleles = normalize_alleles(Some(&truth), Some(&call), false).expect("normalized");
        assert_eq!(alleles.all, ["ATC", "AC"]);
        assert_eq!(
            state(Some(&truth), Some(&call)),
            (TruthState::HetRefVar1, CallState::HetRefVar1)
        );
        let ins = ctx(&["A", "AT"], &[Some(0), Some(1)], false, 30, 10);
        let ins_long = ctx(&["AC", "ATC"], &[Some(0), Some(1)], false, 30, 10);
        let alleles = normalize_alleles(Some(&ins), Some(&ins_long), false).expect("normalized");
        assert_eq!(alleles.all, ["AC", "ACT", "ATC"]);
        let bad = ctx(&["AG", "A"], &[Some(0), Some(1)], false, 30, 10);
        assert!(normalize_alleles(Some(&truth), Some(&bad), false).is_err());
        let haploid = ctx(&["A", "G"], &[Some(1)], false, 30, 10);
        assert!(determine_state(Some(&haploid), None, 0, 0, false).is_err());
    }

    #[test]
    fn classification_by_type_and_summary_rates() {
        let het = ctx(&["A", "G"], &[Some(0), Some(1)], false, 30, 10);
        let hom = ctx(&["A", "G"], &[Some(1), Some(1)], false, 30, 10);
        let indel = ctx(&["A", "AT"], &[Some(0), Some(1)], false, 30, 10);
        let mut counters = Counters::default();
        let classify = |counters: &mut Counters, t: Option<&Variant>, c: Option<&Variant>| {
            counters.classify(t, c, 0, 0, false).expect("classified")
        };
        assert!(classify(&mut counters, Some(&het), Some(&het)));
        assert!(classify(&mut counters, Some(&het), Some(&hom)));
        assert!(classify(&mut counters, Some(&hom), None));
        assert!(classify(&mut counters, None, Some(&indel)));
        assert!(!classify(&mut counters, Some(&het), Some(&indel)));
        let scheme = Scheme::new(false);
        let row = SummaryRow::new("SNP", &counters.snp, &scheme, "t", "c", false);
        // Het truth: TP,TN and TP,FP (both TP); the missing hom call is FN.
        assert_eq!(row.het_sensitivity, 1.0);
        assert_eq!(row.homvar_sensitivity, 0.0);
        assert!((row.var_sensitivity - 2.0 / 3.0).abs() < 1e-12);
        assert!((row.genotype_concordance - 0.5).abs() < 1e-12);
        assert!((row.var_specificity - 0.5).abs() < 1e-12);
        let totals = counters.indel.contingency_counts(&scheme);
        assert_eq!(totals[Contingency::Fp as usize], 1);
        assert_eq!(totals[Contingency::Tn as usize], 1);
        assert_eq!(
            scheme.cell_text(TruthState::Missing, CallState::Missing),
            "NA"
        );
        assert_eq!(
            Scheme::new(true).cell_text(TruthState::Missing, CallState::Missing),
            "TN"
        );
        assert_eq!(
            scheme.cell_text(TruthState::NoCall, CallState::HomVar1),
            "EMPTY"
        );
    }
}
