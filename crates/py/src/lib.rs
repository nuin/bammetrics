use std::collections::BTreeMap;
use std::path::PathBuf;

use bammetrics_core::alignment_summary::{self, DEFAULT_MAX_INSERT_SIZE};
use bammetrics_core::collect::{self, Outputs, Request};
use bammetrics_core::histogram::Histogram;
use bammetrics_core::hs_metrics::{self, COLUMNS, Field};
use bammetrics_core::insert_size::{self, Options};
use bammetrics_core::levels::{Level, parse_levels};
use bammetrics_core::metric::{MetricsOutput, Value};
use bammetrics_core::pass;
use bammetrics_core::quality_yield;
use bammetrics_core::record::Orientation;
use bammetrics_core::registry::{self, Tool};
use pyo3::exceptions::{PyAttributeError, PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyCFunction, PyDict, PyList, PyTuple};

/// Declares a row class: every field readable as an attribute, `to_dict()`
/// returning them in declaration order, and a `__repr__` that shows the
/// fields named after `repr`.
macro_rules! row_class {
    (
        $(#[$meta:meta])*
        $name:ident { $($field:ident: $ty:ty),+ $(,)? }
        repr [$($shown:ident),+ $(,)?]
    ) => {
        $(#[$meta])*
        #[pyclass(get_all)]
        #[derive(Clone)]
        struct $name {
            $($field: $ty),+
        }

        #[pymethods]
        impl $name {
            /// Every field as {name: value}, in declaration order.
            fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                let dict = PyDict::new(py);
                $(dict.set_item(stringify!($field), &self.$field)?;)+
                Ok(dict)
            }

            fn __repr__(&self) -> String {
                let shown: Vec<String> = vec![$(format!("{}={:?}", stringify!($shown), self.$shown)),+];
                format!("{}({})", stringify!($name), shown.join(", "))
            }
        }
    };
}

row_class! {
    InsertSizeMetrics {
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
    sample: Option<String>,
    library: Option<String>,
    read_group: Option<String>,
    }
    repr [pair_orientation, read_pairs, median_insert_size, mean_insert_size, standard_deviation]
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
            sample: r.labels.sample,
            library: r.labels.library,
            read_group: r.labels.read_group,
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
    metric_accumulation_level=None,
))]
fn py_insert_size(
    path: &str,
    minimum_pct: f32,
    deviations: f64,
    histogram_width: Option<i64>,
    min_histogram_width: Option<i64>,
    include_duplicates: bool,
    metric_accumulation_level: Option<Vec<String>>,
) -> PyResult<Vec<InsertSizeMetrics>> {
    let levels = metric_accumulation_level
        .map_or_else(|| Ok(vec![Level::AllReads]), |values| parse_levels(&values));
    let opts = Options {
        minimum_pct,
        deviations,
        histogram_width,
        min_histogram_width,
        include_duplicates,
        levels: levels.map_err(|e| py_error(&e))?,
    };
    let rows = insert_size::run_over_bam(path, &opts).map_err(|e| py_error(&e))?;
    Ok(rows.into_iter().map(Into::into).collect())
}

row_class! {
    /// One `AlignmentSummaryMetrics` row; field names are Picard's columns in
    /// lower case.
    AlignmentSummaryMetrics {
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
    repr [category, total_reads, pf_reads_aligned, pf_mismatch_rate]
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
    let report = alignment_summary::run_over_bam(path, &opts).map_err(|e| py_error(&e))?;
    Ok(report.rows.into_iter().map(Into::into).collect())
}

row_class! {
    /// The `QualityYieldMetrics` row; field names are Picard's columns in lower case.
    QualityYieldMetrics {
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
    repr [total_reads, total_bases, q30_bases]
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
        .map_err(|e| py_error(&e))
}

/// The `HsMetrics` row. Every Picard column is an attribute in lower case
/// (`m.mean_target_coverage`); `HS_LIBRARY_SIZE` and the read-group columns
/// may be None. The two histograms are dicts keyed by depth or quality.
#[pyclass]
struct HsMetrics {
    fields: Vec<(String, Field)>,
    /// `high_quality_coverage_count`: target bases per high-quality depth.
    #[pyo3(get)]
    high_quality_coverage: BTreeMap<i64, f64>,
    /// `unfiltered_baseq_count`: on-target bases above Q2 per base quality.
    #[pyo3(get)]
    unfiltered_base_quality: BTreeMap<i64, f64>,
}

fn histogram_dict(histogram: &Histogram) -> BTreeMap<i64, f64> {
    histogram
        .keys()
        .map(|k| (k, histogram.get(k).unwrap_or(0.0)))
        .collect()
}

fn field_object(py: Python<'_>, field: &Field) -> PyResult<Py<PyAny>> {
    Ok(match field {
        Field::Text(s) => s.into_pyobject(py)?.into_any().unbind(),
        Field::Long(v) => v.into_pyobject(py)?.into_any().unbind(),
        Field::Double(v) => v.into_pyobject(py)?.into_any().unbind(),
        Field::Null => py.None(),
    })
}

#[pymethods]
impl HsMetrics {
    fn __getattr__(&self, py: Python<'_>, name: &str) -> PyResult<Py<PyAny>> {
        let field = self
            .fields
            .iter()
            .find(|(column, _)| column == name)
            .map(|(_, field)| field)
            .ok_or_else(|| PyAttributeError::new_err(format!("HsMetrics has no field {name:?}")))?;
        field_object(py, field)
    }

    /// Every column, in Picard's order, as {lower-case name: value}. The two
    /// histograms are attributes, not entries.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (name, field) in &self.fields {
            dict.set_item(name, field_object(py, field)?)?;
        }
        Ok(dict)
    }

    fn __repr__(&self) -> String {
        let text = |name: &str| {
            self.fields
                .iter()
                .find(|(column, _)| column == name)
                .map_or_else(String::new, |(_, f)| f.picard_text())
        };
        format!(
            "HsMetrics(bait_set={}, mean_target_coverage={}, pct_target_bases_20x={})",
            text("bait_set"),
            text("mean_target_coverage"),
            text("pct_target_bases_20x")
        )
    }
}

impl From<hs_metrics::Report> for HsMetrics {
    fn from(report: hs_metrics::Report) -> Self {
        Self {
            fields: COLUMNS
                .iter()
                .map(|c| c.to_ascii_lowercase())
                .zip(report.row.fields())
                .collect(),
            high_quality_coverage: histogram_dict(&report.depth),
            unfiltered_base_quality: histogram_dict(&report.base_quality),
        }
    }
}

/// One path or a list of paths, as Picard's repeatable interval options.
fn path_list(value: &Bound<'_, PyAny>) -> PyResult<Vec<PathBuf>> {
    match value.extract::<PathBuf>() {
        Ok(path) => Ok(vec![path]),
        Err(_) => value.extract::<Vec<PathBuf>>(),
    }
}

/// `CollectHsMetrics`: one `HsMetrics` for the BAM. `bait_intervals` and
/// `target_intervals` take a path or a list of paths (interval lists or BED).
/// `per_base_coverage` names a file for Picard's PER_BASE_COVERAGE output.
#[pyfunction]
#[pyo3(name = "hs_metrics", signature = (
    path,
    bait_intervals,
    target_intervals,
    reference=None,
    bait_set_name=None,
    minimum_mapping_quality=hs_metrics::DEFAULT_MINIMUM_MAPPING_QUALITY,
    minimum_base_quality=hs_metrics::DEFAULT_MINIMUM_BASE_QUALITY,
    clip_overlapping_reads=true,
    include_indels=false,
    coverage_cap=hs_metrics::DEFAULT_COVERAGE_CAP,
    near_distance=hs_metrics::DEFAULT_NEAR_DISTANCE,
    sample_size=hs_metrics::DEFAULT_SAMPLE_SIZE,
    per_base_coverage=None,
))]
#[allow(clippy::too_many_arguments)]
fn py_hs_metrics(
    path: &str,
    bait_intervals: &Bound<'_, PyAny>,
    target_intervals: &Bound<'_, PyAny>,
    reference: Option<PathBuf>,
    bait_set_name: Option<String>,
    minimum_mapping_quality: i32,
    minimum_base_quality: i32,
    clip_overlapping_reads: bool,
    include_indels: bool,
    coverage_cap: i32,
    near_distance: i32,
    sample_size: usize,
    per_base_coverage: Option<PathBuf>,
) -> PyResult<HsMetrics> {
    let opts = hs_metrics::Options {
        bait_intervals: path_list(bait_intervals)?,
        target_intervals: path_list(target_intervals)?,
        bait_set_name,
        reference,
        minimum_mapping_quality,
        minimum_base_quality,
        clip_overlapping_reads,
        include_indels,
        coverage_cap,
        near_distance,
        sample_size,
        per_base_coverage,
    };
    hs_metrics::run_over_bam(path, &opts)
        .map(Into::into)
        .map_err(|e| py_error(&e))
}

/// An error from the core as a Python exception: bad arguments as
/// ValueError, everything else as IOError.
fn py_error(e: &bammetrics_core::Error) -> PyErr {
    match e {
        bammetrics_core::Error::Usage(message) => PyValueError::new_err(message.clone()),
        other => PyIOError::new_err(other.to_string()),
    }
}

/// Several metrics from one pass over the BAM. `metrics` lists any of
/// "quality_yield", "alignment_summary", "insert_size" and "hs_metrics"; by
/// default the first three, plus "hs_metrics" when `hs` is given. `hs` is a
/// list of (bait_intervals, target_intervals) pairs, one `CollectHsMetrics`
/// run each; each side takes a path or a list of paths. `reference` goes to
/// alignment_summary and hs_metrics. Every other option is Picard's default.
/// `threads` is the BGZF worker count (default: cores, at most 8).
///
/// Returns {name: result}, each result as the single-metric function
/// returns it, and a list of HsMetrics for "hs_metrics". With
/// `output_prefix`, every metric is also written as a Picard file, named as
/// by `bammetrics collect -O PREFIX`. `hs_per_base` (needs `output_prefix`)
/// also writes each HS run's per-base coverage, as `--hs-per-base` does. A
/// failed metric raises, after the files of the others are written.
#[pyfunction]
#[pyo3(name = "collect", signature = (
    path,
    metrics=None,
    reference=None,
    hs=None,
    output_prefix=None,
    threads=None,
    hs_per_base=false,
))]
#[allow(clippy::too_many_arguments)]
fn py_collect<'py>(
    py: Python<'py>,
    path: PathBuf,
    metrics: Option<Vec<String>>,
    reference: Option<PathBuf>,
    hs: Option<Vec<(Bound<'py, PyAny>, Bound<'py, PyAny>)>>,
    output_prefix: Option<PathBuf>,
    threads: Option<usize>,
    hs_per_base: bool,
) -> PyResult<Bound<'py, PyDict>> {
    const NAMES: [&str; 4] = [
        "quality_yield",
        "alignment_summary",
        "insert_size",
        "hs_metrics",
    ];
    let hs = hs.unwrap_or_default();
    let wanted = metrics.unwrap_or_else(|| {
        let n = if hs.is_empty() { 3 } else { 4 };
        NAMES[..n].iter().map(|s| (*s).to_owned()).collect()
    });
    if let Some(unknown) = wanted.iter().find(|m| !NAMES.contains(&m.as_str())) {
        return Err(PyValueError::new_err(format!(
            "unknown metric {unknown:?}; expected one of {NAMES:?}"
        )));
    }
    let has = |name: &str| wanted.iter().any(|m| m == name);
    if has("hs_metrics") == hs.is_empty() {
        return Err(PyValueError::new_err(
            "\"hs_metrics\" in metrics and a non-empty hs list go together",
        ));
    }
    if hs_per_base && (hs.is_empty() || output_prefix.is_none()) {
        return Err(PyValueError::new_err(
            "hs_per_base needs a non-empty hs list and output_prefix",
        ));
    }
    let mut request = Request {
        quality_yield: has("quality_yield").then(quality_yield::Options::default),
        alignment_summary: has("alignment_summary").then(|| alignment_summary::Options {
            reference: reference.clone(),
            ..alignment_summary::Options::default()
        }),
        insert_size: has("insert_size").then(Options::default),
        hs_metrics: hs
            .iter()
            .map(|(baits, targets)| {
                Ok(hs_metrics::Options {
                    bait_intervals: path_list(baits)?,
                    target_intervals: path_list(targets)?,
                    reference: reference.clone(),
                    ..hs_metrics::Options::default()
                })
            })
            .collect::<PyResult<_>>()?,
    };
    let outputs = output_prefix
        .as_deref()
        .map(|prefix| Outputs::new(prefix, &request))
        .transpose()
        .map_err(|e| py_error(&e))?;
    if let (true, Some(outputs)) = (hs_per_base, &outputs) {
        for (opts, path) in request.hs_metrics.iter_mut().zip(&outputs.hs_metrics) {
            opts.per_base_coverage = Some(collect::per_base_coverage_path(path));
        }
    }
    let threads = threads.unwrap_or_else(pass::default_threads);
    let results = py
        .detach(|| {
            let results = collect::run_over_bam(&path, &request, threads)?;
            if let Some(outputs) = &outputs {
                collect::write_files(&request, &results, &path, outputs)?;
            }
            Ok::<_, bammetrics_core::Error>(results)
        })
        .map_err(|e| py_error(&e))?;
    if let Some((tool, e)) = results.errors().first() {
        return Err(PyIOError::new_err(format!("{tool}: {e}")));
    }
    let dict = PyDict::new(py);
    if let Some(Ok(row)) = results.quality_yield {
        dict.set_item("quality_yield", QualityYieldMetrics::from(row))?;
    }
    if let Some(Ok(report)) = results.alignment_summary {
        let rows: Vec<AlignmentSummaryMetrics> = report.rows.into_iter().map(Into::into).collect();
        dict.set_item("alignment_summary", rows)?;
    }
    if let Some(Ok(rows)) = results.insert_size {
        let rows: Vec<InsertSizeMetrics> = rows.into_iter().map(Into::into).collect();
        dict.set_item("insert_size", rows)?;
    }
    if !request.hs_metrics.is_empty() {
        let reports: Vec<HsMetrics> = results
            .hs_metrics
            .into_iter()
            .filter_map(Result::ok)
            .map(Into::into)
            .collect();
        dict.set_item("hs_metrics", reports)?;
    }
    Ok(dict)
}

/// A metrics value as a Python object: int, float, str, bool or None.
fn value_object(py: Python<'_>, value: &Value) -> PyResult<Py<PyAny>> {
    Ok(match value {
        Value::Int(v) => v.into_pyobject(py)?.into_any().unbind(),
        Value::Float(v) => v.into_pyobject(py)?.into_any().unbind(),
        Value::Str(v) => v.into_pyobject(py)?.into_any().unbind(),
        Value::Bool(v) => v.into_pyobject(py)?.to_owned().into_any().unbind(),
        Value::Null => py.None(),
    })
}

/// The generic converter: every registered tool's rows as a list of dicts
/// keyed by Picard's column names, in column order.
fn rows_as_dicts<'py>(py: Python<'py>, output: &MetricsOutput) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty(py);
    for record in output.records() {
        let dict = PyDict::new(py);
        for (column, value) in record {
            dict.set_item(column, value_object(py, value)?)?;
        }
        list.append(dict)?;
    }
    Ok(list)
}

/// Python keyword arguments as Picard `KEY=VALUE` arguments: upper-case
/// names, bools as true/false, lists as repeated options, `None` skipped.
fn keyword_arguments(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<String>> {
    let mut args = Vec::new();
    let Some(kwargs) = kwargs else {
        return Ok(args);
    };
    for (key, value) in kwargs.iter() {
        let name = key.extract::<String>()?.to_uppercase();
        if value.is_none() {
            continue;
        }
        let texts: Vec<String> = if let Ok(items) = value.extract::<Vec<String>>() {
            items
        } else if let Ok(flag) = value.extract::<bool>() {
            vec![flag.to_string()]
        } else {
            vec![value.str()?.to_string()]
        };
        args.extend(texts.into_iter().map(|t| format!("{name}={t}")));
    }
    Ok(args)
}

/// Run a registered tool: the input path plus Picard options as keyword
/// arguments (`minimum_pct=0.1`, `metric_accumulation_level=["SAMPLE"]`).
/// Returns one dict per metrics row, empty when Picard would write no file.
/// With `output=` the Picard metrics file is written as well.
fn run_registered<'py>(
    py: Python<'py>,
    tool: &'static Tool,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyList>> {
    let input: String = args.get_item(0)?.str()?.to_string();
    let mut options = vec![format!("INPUT={input}")];
    options.extend(keyword_arguments(kwargs)?);
    let (invocation, output) = py
        .detach(|| tool.execute(&options))
        .map_err(|e| py_error(&e))?;
    let Some(output) = output else {
        return Ok(PyList::empty(py));
    };
    if let Some(path) = invocation.output.as_deref() {
        let headers = bammetrics_core::metrics_file::headers(invocation.command.clone());
        registry::write_to(Some(path), |out| output.write(headers, out))
            .map_err(|e| py_error(&e))?;
    }
    rows_as_dicts(py, &output)
}

/// `CollectInsertSizeMetrics` -> `collect_insert_size_metrics`.
fn snake_case(name: &str) -> String {
    let mut out = String::new();
    for (i, c) in name.chars().enumerate() {
        if c.is_ascii_uppercase() {
            if i > 0 {
                out.push('_');
            }
            out.push(c.to_ascii_lowercase());
        } else {
            out.push(c);
        }
    }
    out
}

/// Add one function per registered tool, named after its Picard class in
/// snake case.
fn add_registered_tools(m: &Bound<'_, PyModule>) -> PyResult<()> {
    for &tool in registry::TOOLS {
        let function = PyCFunction::new_closure(
            m.py(),
            None,
            None,
            move |args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>| {
                run_registered(args.py(), tool, args, kwargs).map(Bound::unbind)
            },
        )?;
        m.add(snake_case(tool.name).as_str(), function)?;
    }
    Ok(())
}

/// Picard names of the registered tools.
#[pyfunction]
fn tools() -> Vec<&'static str> {
    registry::TOOLS.iter().map(|t| t.name).collect()
}

#[pymodule]
fn bammetrics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_registered_tools(m)?;
    m.add_function(wrap_pyfunction!(tools, m)?)?;
    m.add_class::<InsertSizeMetrics>()?;
    m.add_function(wrap_pyfunction!(py_insert_size, m)?)?;
    m.add_class::<AlignmentSummaryMetrics>()?;
    m.add_function(wrap_pyfunction!(py_alignment_summary, m)?)?;
    m.add_class::<QualityYieldMetrics>()?;
    m.add_function(wrap_pyfunction!(py_quality_yield, m)?)?;
    m.add_class::<HsMetrics>()?;
    m.add_function(wrap_pyfunction!(py_hs_metrics, m)?)?;
    m.add_function(wrap_pyfunction!(py_collect, m)?)?;
    Ok(())
}
