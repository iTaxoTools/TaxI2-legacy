mod distance;
mod needle;

use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

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

/// A Python module implemented in Rust.
#[pymodule]
fn calculate_distances(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(align_to_str, m)?)?;

    Ok(())
}
