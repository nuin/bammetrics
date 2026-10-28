use bammetrics_core::insert_size::{self, Options};
use pyo3::exceptions::PyIOError;
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

#[pymodule]
fn bammetrics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<InsertSizeMetrics>()?;
    m.add_function(wrap_pyfunction!(py_insert_size, m)?)?;
    Ok(())
}
