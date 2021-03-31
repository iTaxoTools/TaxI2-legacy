mod distance;
mod needle;

use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

pub use crate::distance::AlignmentStats;
pub use crate::needle::Aligner;

/// Makes an Aligner with given scores
#[pyfunction]
#[text_signature = "(match_score, mismatch_score, end_open_gap_score, end_extend_gap_score, internal_open_gap_score, internal_extend_gap_score, /)"]
pub fn make_aligner(
    match_score: i16,
    mismatch_score: i16,
    end_open_gap_score: i16,
    end_extend_gap_score: i16,
    internal_open_gap_score: i16,
    internal_extend_gap_score: i16,
) -> Aligner {
    Aligner {
        match_score,
        mismatch_score,
        end_gap_penalty: end_open_gap_score,
        end_gap_extend_penalty: end_extend_gap_score,
        gap_penalty: internal_open_gap_score,
        gap_extend_penalty: internal_extend_gap_score,
    }
}

/// Returns two strings that represent aligned `target` and aligned `query` respectively.
#[pyfunction]
#[text_signature = "(target, query, /)"]
fn align_to_str(target: &str, query: &str) -> PyResult<(String, String)> {
    let aligner = Aligner::default();
    aligner
        .align(target.as_bytes(), query.as_bytes())
        .as_strings()
        .map_err(exceptions::PyUnicodeEncodeError::new_err)
}

/// Returns two strings that represent aligned `target` and aligned `query` respectively.
#[pyfunction]
#[text_signature = "(aligner, target, query, /)"]
fn show_alignment(aligner: &Aligner, target: &str, query: &str) -> PyResult<String> {
    aligner
        .align(target.as_bytes(), query.as_bytes())
        .show_alignment()
        .map_err(exceptions::PyUnicodeEncodeError::new_err)
}

/// Returns 4 distances between `target` and `query`.
///
/// Performs alignment.
#[pyfunction]
#[text_signature = "(aligner, target, query, /)"]
fn seq_distances(aligner: &Aligner, target: &str, query: &str) -> [f64; 4] {
    let alignment = aligner.align(target.as_bytes(), query.as_bytes());
    let mut alignment_stats = AlignmentStats::new();
    alignment
        .common_path_iter()
        .for_each(|pair| alignment_stats.update(pair));
    [
        alignment_stats.pdistance(),
        alignment_stats.jukes_cantor_distance(),
        alignment_stats.kimura2p_distance(),
        alignment_stats.pdistance_counting_gaps(),
    ]
}

// Returns true if the character is part of a meaningful part of a sequences
fn is_nucleotide(c: char) -> bool {
    !matches!(c, '-' | 'n' | 'N' | '?')
}

// Returns the inclusive boundaries of the common non-gap part of given sequences
fn common_content(target: &str, query: &str) -> Option<(usize, usize)> {
    let target_start = target.find(is_nucleotide)?;
    let query_start = query.find(is_nucleotide)?;
    let target_end = target.rfind(is_nucleotide)?;
    let query_end = query.rfind(is_nucleotide)?;
    let start = usize::max(target_start, query_start);
    let end = usize::min(target_end, query_end);
    Some((start, end))
}

/// Returns 4 distances between `target` and `query`.
///
/// Expects aligned sequences.
#[pyfunction]
#[text_signature = "(target, query, /)"]
fn seq_distances_aligned(target: &str, query: &str) -> [f64; 4] {
    let (start, end) = match common_content(target, query) {
        None => return [f64::NAN; 4],
        Some(x) => x,
    };
    let target = &target[start..=end];
    let query = &query[start..=end];
    let mut alignment_stats = AlignmentStats::new();
    target
        .bytes()
        .zip(query.bytes())
        .for_each(|pair| alignment_stats.update(pair));
    [
        alignment_stats.pdistance(),
        alignment_stats.jukes_cantor_distance(),
        alignment_stats.kimura2p_distance(),
        alignment_stats.pdistance_counting_gaps(),
    ]
}

/// A Python module implemented in Rust.
#[pymodule]
fn calculate_distances(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(align_to_str, m)?)?;
    m.add_function(wrap_pyfunction!(make_aligner, m)?)?;
    m.add_function(wrap_pyfunction!(seq_distances, m)?)?;
    m.add_function(wrap_pyfunction!(seq_distances_aligned, m)?)?;
    m.add_function(wrap_pyfunction!(show_alignment, m)?)?;

    Ok(())
}
