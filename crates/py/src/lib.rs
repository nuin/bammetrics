mod histogram;
mod insert_size;

use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
struct InsertSizeMetrics {
    #[pyo3(get)] pair_orientation: String,
    #[pyo3(get)] read_pairs: i64,
    #[pyo3(get)] min_insert_size: i64,
    #[pyo3(get)] max_insert_size: i64,
    #[pyo3(get)] median_insert_size: f64,
    #[pyo3(get)] mode_insert_size: i64,
    #[pyo3(get)] median_absolute_deviation: f64,
    #[pyo3(get)] mean_insert_size: f64,
    #[pyo3(get)] standard_deviation: f64,
    #[pyo3(get)] width_of_50_percent: i64,
}

#[pymethods]
impl InsertSizeMetrics {
    fn __repr__(&self) -> String {
        format!("InsertSizeMetrics(orientation={}, read_pairs={}, median={}, mean={:.4}, sd={:.4})",
            self.pair_orientation, self.read_pairs, self.median_insert_size,
            self.mean_insert_size, self.standard_deviation)
    }
}

#[pyfunction]
#[pyo3(name = "insert_size", signature = (path, minimum_pct=0.05, deviations=10.0, include_duplicates=false))]
fn py_insert_size(path: &str, minimum_pct: f64, deviations: f64, include_duplicates: bool)
    -> PyResult<Vec<InsertSizeMetrics>>
{
    let opts = insert_size::Options { minimum_pct, deviations, include_duplicates };
    let rows = insert_size::run_over_bam(path, opts)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    Ok(rows.into_iter().map(|r| InsertSizeMetrics {
        pair_orientation: r.orientation.unwrap_or("NA").to_string(),
        read_pairs: r.read_pairs,
        min_insert_size: r.min_insert_size,
        max_insert_size: r.max_insert_size,
        median_insert_size: r.median_insert_size,
        mode_insert_size: r.mode_insert_size,
        median_absolute_deviation: r.median_absolute_deviation,
        mean_insert_size: r.mean_insert_size,
        standard_deviation: r.standard_deviation,
        width_of_50_percent: r.width_of_50_percent,
    }).collect())
}

#[pymodule]
fn bammetrics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<InsertSizeMetrics>()?;
    m.add_function(wrap_pyfunction!(py_insert_size, m)?)?;
    Ok(())
}
