mod distance;
mod needle;

use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

pub use crate::distance::AlignmentStats;
pub use crate::needle::Aligner;

/// Returns two strings that represent aligned `target` and aligned `query` respectively.
#[pyfunction]
fn align_to_str(target: &str, query: &str) -> PyResult<(String, String)> {
    let aligner = Aligner::default();
    aligner
        .align(target.as_bytes(), query.as_bytes())
        .as_strings()
        .map_err(exceptions::PyUnicodeEncodeError::new_err)
}

/// Returns 4 distances between `target` and `query`.
///
/// Performs alignment.
#[allow(clippy::unnecessary_wraps)]
#[pyfunction]
fn seq_distances(target: &str, query: &str) -> PyResult<[f64; 4]> {
    let aligner = Aligner::default();
    let alignment = aligner.align(target.as_bytes(), query.as_bytes());
    let mut alignment_stats = AlignmentStats::new();
    alignment
        .common_path_iter()
        .for_each(|pair| alignment_stats.update(pair));
    Ok([
        alignment_stats.pdistance(),
        alignment_stats.jukes_cantor_distance(),
        alignment_stats.kimura2p_distance(),
        alignment_stats.pdistance_counting_gaps(),
    ])
}

/// A Python module implemented in Rust.
#[pymodule]
fn calculate_distances(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(align_to_str, m)?)?;
    m.add_function(wrap_pyfunction!(seq_distances, m)?)?;

    Ok(())
}
