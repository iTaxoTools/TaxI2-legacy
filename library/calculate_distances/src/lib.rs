mod needle;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

pub use crate::needle::Aligner;

/// Returns two strings that represent aligned `target` and aligned `query` respectively.
#[pyfunction]
fn align_to_str(target: &str, query: &str) -> PyResult<(String, String)> {
    let aligner = Aligner::default();
    Ok(aligner.align_to_str(target, query))
}

/// A Python module implemented in Rust.
#[pymodule]
fn calculate_distances(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(align_to_str, m)?)?;

    Ok(())
}
