use std::path::PathBuf;

use bammetrics_core::alignment_summary::{self, DEFAULT_MAX_INSERT_SIZE};
use bammetrics_core::insert_size::{self, Options};
use bammetrics_core::quality_yield;
use bammetrics_core::record::Orientation;
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;

#[pyclass(get_all)]
#[derive(Clone)]
struct InsertSizeMetrics {
    pair_orientation: String,
    read_pairs: i64,
    min_insert_size: i64,
    max_insert_size: i64,
    median_insert_size: f64,
    mode_insert_size: i64,
    median_absolute_deviation: f64,
    mean_insert_size: f64,
    standard_deviation: f64,
    width_of_10_percent: i64,
    width_of_20_percent: i64,
    width_of_30_percent: i64,
    width_of_40_percent: i64,
    width_of_50_percent: i64,
    width_of_60_percent: i64,
    width_of_70_percent: i64,
    width_of_80_percent: i64,
    width_of_90_percent: i64,
    width_of_95_percent: i64,
    width_of_99_percent: i64,
}

#[pymethods]
impl InsertSizeMetrics {
    fn __repr__(&self) -> String {
        format!(
            "InsertSizeMetrics(orientation={}, read_pairs={}, median={}, mean={:.4}, sd={:.4})",
            self.pair_orientation,
            self.read_pairs,
            self.median_insert_size,
            self.mean_insert_size,
            self.standard_deviation
        )
    }
}

impl From<insert_size::Row> for InsertSizeMetrics {
    fn from(r: insert_size::Row) -> Self {
        Self {
            pair_orientation: r.orientation.as_str().to_owned(),
            read_pairs: r.read_pairs,
            min_insert_size: r.min_insert_size,
            max_insert_size: r.max_insert_size,
            median_insert_size: r.median_insert_size,
            mode_insert_size: r.mode_insert_size,
            median_absolute_deviation: r.median_absolute_deviation,
            mean_insert_size: r.mean_insert_size,
            standard_deviation: r.standard_deviation,
            width_of_10_percent: r.widths[0],
            width_of_20_percent: r.widths[1],
            width_of_30_percent: r.widths[2],
            width_of_40_percent: r.widths[3],
            width_of_50_percent: r.widths[4],
            width_of_60_percent: r.widths[5],
            width_of_70_percent: r.widths[6],
            width_of_80_percent: r.widths[7],
            width_of_90_percent: r.widths[8],
            width_of_95_percent: r.widths[9],
            width_of_99_percent: r.widths[10],
        }
    }
}

#[pyfunction]
#[pyo3(name = "insert_size", signature = (
    path,
    minimum_pct=0.05,
    deviations=10.0,
    histogram_width=None,
    min_histogram_width=None,
    include_duplicates=false,
))]
fn py_insert_size(
    path: &str,
    minimum_pct: f32,
    deviations: f64,
    histogram_width: Option<i64>,
    min_histogram_width: Option<i64>,
    include_duplicates: bool,
) -> PyResult<Vec<InsertSizeMetrics>> {
    let opts = Options {
        minimum_pct,
        deviations,
        histogram_width,
        min_histogram_width,
        include_duplicates,
    };
    let rows =
        insert_size::run_over_bam(path, &opts).map_err(|e| PyIOError::new_err(e.to_string()))?;
    Ok(rows.into_iter().map(Into::into).collect())
}

/// One `AlignmentSummaryMetrics` row; field names are Picard's columns in
/// lower case.
#[pyclass(get_all)]
#[derive(Clone)]
struct AlignmentSummaryMetrics {
    category: String,
    total_reads: i64,
    pf_reads: i64,
    pct_pf_reads: f64,
    pf_noise_reads: i64,
    pf_reads_aligned: i64,
    pct_pf_reads_aligned: f64,
    pf_aligned_bases: i64,
    pf_hq_aligned_reads: i64,
    pf_hq_aligned_bases: i64,
    pf_hq_aligned_q20_bases: i64,
    pf_hq_median_mismatches: f64,
    pf_mismatch_rate: f64,
    pf_hq_error_rate: f64,
    pf_indel_rate: f64,
    mean_read_length: f64,
    sd_read_length: f64,
    median_read_length: f64,
    mad_read_length: f64,
    min_read_length: f64,
    max_read_length: f64,
    mean_aligned_read_length: f64,
    reads_aligned_in_pairs: i64,
    pct_reads_aligned_in_pairs: f64,
    pf_reads_improper_pairs: i64,
    pct_pf_reads_improper_pairs: f64,
    bad_cycles: i64,
    strand_balance: f64,
    pct_chimeras: f64,
    pct_adapter: f64,
    pct_softclip: f64,
    pct_hardclip: f64,
    avg_pos_3prime_softclip_length: f64,
}

#[pymethods]
impl AlignmentSummaryMetrics {
    fn __repr__(&self) -> String {
        format!(
            "AlignmentSummaryMetrics(category={}, total_reads={}, pf_reads_aligned={}, pf_mismatch_rate={:.6})",
            self.category, self.total_reads, self.pf_reads_aligned, self.pf_mismatch_rate
        )
    }
}

impl From<alignment_summary::Row> for AlignmentSummaryMetrics {
    fn from(r: alignment_summary::Row) -> Self {
        Self {
            category: r.category.as_str().to_owned(),
            total_reads: r.total_reads,
            pf_reads: r.pf_reads,
            pct_pf_reads: r.pct_pf_reads,
            pf_noise_reads: r.pf_noise_reads,
            pf_reads_aligned: r.pf_reads_aligned,
            pct_pf_reads_aligned: r.pct_pf_reads_aligned,
            pf_aligned_bases: r.pf_aligned_bases,
            pf_hq_aligned_reads: r.pf_hq_aligned_reads,
            pf_hq_aligned_bases: r.pf_hq_aligned_bases,
            pf_hq_aligned_q20_bases: r.pf_hq_aligned_q20_bases,
            pf_hq_median_mismatches: r.pf_hq_median_mismatches,
            pf_mismatch_rate: r.pf_mismatch_rate,
            pf_hq_error_rate: r.pf_hq_error_rate,
            pf_indel_rate: r.pf_indel_rate,
            mean_read_length: r.mean_read_length,
            sd_read_length: r.sd_read_length,
            median_read_length: r.median_read_length,
            mad_read_length: r.mad_read_length,
            min_read_length: r.min_read_length,
            max_read_length: r.max_read_length,
            mean_aligned_read_length: r.mean_aligned_read_length,
            reads_aligned_in_pairs: r.reads_aligned_in_pairs,
            pct_reads_aligned_in_pairs: r.pct_reads_aligned_in_pairs,
            pf_reads_improper_pairs: r.pf_reads_improper_pairs,
            pct_pf_reads_improper_pairs: r.pct_pf_reads_improper_pairs,
            bad_cycles: r.bad_cycles,
            strand_balance: r.strand_balance,
            pct_chimeras: r.pct_chimeras,
            pct_adapter: r.pct_adapter,
            pct_softclip: r.pct_softclip,
            pct_hardclip: r.pct_hardclip,
            avg_pos_3prime_softclip_length: r.avg_pos_3prime_softclip_length,
        }
    }
}

/// `CollectAlignmentSummaryMetrics` at `ALL_READS`. A list given for
/// `expected_pair_orientations` or `adapter_sequences` replaces Picard's
/// default (`["FR"]` and its six Illumina adapters); `None` keeps it.
#[pyfunction]
#[pyo3(name = "alignment_summary", signature = (
    path,
    reference=None,
    max_insert_size=DEFAULT_MAX_INSERT_SIZE,
    expected_pair_orientations=None,
    adapter_sequences=None,
    collect_alignment_information=true,
))]
fn py_alignment_summary(
    path: &str,
    reference: Option<PathBuf>,
    max_insert_size: i32,
    expected_pair_orientations: Option<Vec<String>>,
    adapter_sequences: Option<Vec<String>>,
    collect_alignment_information: bool,
) -> PyResult<Vec<AlignmentSummaryMetrics>> {
    let defaults = alignment_summary::Options::default();
    let expected_orientations = match expected_pair_orientations {
        None => defaults.expected_orientations,
        Some(names) => names
            .iter()
            .map(|n| n.parse::<Orientation>())
            .collect::<Result<_, _>>()
            .map_err(|e| PyValueError::new_err(e.to_string()))?,
    };
    let opts = alignment_summary::Options {
        reference,
        max_insert_size,
        expected_orientations,
        adapter_sequences: adapter_sequences.unwrap_or(defaults.adapter_sequences),
        collect_alignment_information,
    };
    let report = alignment_summary::run_over_bam(path, &opts)
        .map_err(|e| PyIOError::new_err(e.to_string()))?;
    Ok(report.rows.into_iter().map(Into::into).collect())
}

/// The `QualityYieldMetrics` row; field names are Picard's columns in lower case.
#[pyclass(get_all)]
#[derive(Clone)]
struct QualityYieldMetrics {
    total_reads: i64,
    pf_reads: i64,
    read_length: i32,
    total_bases: i64,
    pf_bases: i64,
    q20_bases: i64,
    pf_q20_bases: i64,
    q30_bases: i64,
    pf_q30_bases: i64,
    q20_equivalent_yield: i64,
    pf_q20_equivalent_yield: i64,
}

#[pymethods]
impl QualityYieldMetrics {
    fn __repr__(&self) -> String {
        format!(
            "QualityYieldMetrics(total_reads={}, total_bases={}, q30_bases={})",
            self.total_reads, self.total_bases, self.q30_bases
        )
    }
}

impl From<quality_yield::Row> for QualityYieldMetrics {
    fn from(r: quality_yield::Row) -> Self {
        Self {
            total_reads: r.total_reads,
            pf_reads: r.pf_reads,
            read_length: r.read_length,
            total_bases: r.total_bases,
            pf_bases: r.pf_bases,
            q20_bases: r.q20_bases,
            pf_q20_bases: r.pf_q20_bases,
            q30_bases: r.q30_bases,
            pf_q30_bases: r.pf_q30_bases,
            q20_equivalent_yield: r.q20_equivalent_yield,
            pf_q20_equivalent_yield: r.pf_q20_equivalent_yield,
        }
    }
}

/// `CollectQualityYieldMetrics`: one `QualityYieldMetrics` for the BAM.
#[pyfunction]
#[pyo3(name = "quality_yield", signature = (
    path,
    use_original_qualities=true,
    include_secondary_alignments=false,
    include_supplemental_alignments=false,
))]
fn py_quality_yield(
    path: &str,
    use_original_qualities: bool,
    include_secondary_alignments: bool,
    include_supplemental_alignments: bool,
) -> PyResult<QualityYieldMetrics> {
    let opts = quality_yield::Options {
        use_original_qualities,
        include_secondary_alignments,
        include_supplemental_alignments,
    };
    quality_yield::run_over_bam(path, &opts)
        .map(Into::into)
        .map_err(|e| PyIOError::new_err(e.to_string()))
}

#[pymodule]
fn bammetrics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<InsertSizeMetrics>()?;
    m.add_function(wrap_pyfunction!(py_insert_size, m)?)?;
    m.add_class::<AlignmentSummaryMetrics>()?;
    m.add_function(wrap_pyfunction!(py_alignment_summary, m)?)?;
    m.add_class::<QualityYieldMetrics>()?;
    m.add_function(wrap_pyfunction!(py_quality_yield, m)?)?;
    Ok(())
}
