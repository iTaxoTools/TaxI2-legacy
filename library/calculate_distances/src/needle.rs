//! Aligmnent using Needleman-Wunsch algorithm.
mod score;
mod table;

use score::{Dir, Score};
use table::Table;

/// Contains parameters for the Needleman-Wunsch algorithm.
pub struct Aligner {
    match_score: i16,
    mismatch_score: i16,
    gap_penalty: i16,
    gap_extend_penalty: i16,
    end_gap_penalty: i16,
    end_gap_extend_penalty: i16,
}

impl Aligner {
    /// Returns two strings representing aligned `target` and `query` respectively.
    ///
    /// `target` and `query` should be ASCII strings.
    ///
    /// # Panics
    /// May panic if non-ASCII characters are encountered.
    pub fn align_to_str(&self, target: &str, query: &str) -> (String, String) {
        // Convert arguments into byte arrays
        let target = target.as_bytes();
        let query = query.as_bytes();
        // Allocate space for the result
        let (mut target_align, mut query_align) = {
            let len = target.len().max(query.len());
            (Vec::with_capacity(len), Vec::with_capacity(len))
        };
        // Perform alignment
        let table = self.align(target, query);
        // Start from the bottom-right corner
        let (mut i, mut j) = (target.len(), query.len());
        while i > 0 || j > 0 {
            // iterate until the top-left corner
            let current = table[[i, j]];
            // Push corresponding characters to the result
            let (target_c, query_c) = match current.dir {
                Dir::Diagonal => (target[i - 1], query[j - 1]),
                Dir::Up => (b'-', query[j - 1]),
                Dir::Left => (target[i - 1], b'-'),
            };
            target_align.push(target_c);
            query_align.push(query_c);
            // Move to the next cell
            match current.dir {
                Dir::Diagonal => {
                    i -= 1;
                    j -= 1;
                }
                Dir::Up => {
                    j -= 1;
                }
                Dir::Left => {
                    i -= 1;
                }
            };
        }
        // Reverse the result, since the iteration was backwards
        target_align.reverse();
        query_align.reverse();

        (
            String::from_utf8(target_align).unwrap(),
            String::from_utf8(query_align).unwrap(),
        )
    }

    /// Construct the Needleman-Wunsch table for `target` and `query`.
    ///
    /// `target` is horizontal, `query` is vertical.
    fn align(&self, target: &[u8], query: &[u8]) -> Table<Score> {
        // Allocate and initialize table
        let mut table = Table::<Score>::new(target.len() + 1, query.len() + 1);
        table.fill_default();
        // Fill the left column and the top row
        table[[0, 0]] = Score::new(Dir::Diagonal, 0);
        for i in 0..target.len() {
            table[[i + 1, 0]] = Score::new(
                Dir::Left,
                self.end_gap_penalty + i as i16 * self.end_gap_extend_penalty,
            );
        }
        for j in 0..query.len() {
            table[[0, j + 1]] = Score::new(
                Dir::Up,
                self.end_gap_penalty + j as i16 * self.end_gap_extend_penalty,
            );
        }

        // Fill the rest of the table
        for j in 0..query.len() {
            // row iteration for efficiency
            for i in 0..target.len() {
                // Calculate score from the left
                let left_cell = table[[i, j + 1]];
                let left_score = left_cell.score
                    + match (left_cell.dir, j == query.len() - 1) {
                        // `j == query.len() - 1` means that we are currently in the bottom row,
                        // therefore the gap is an end gap.
                        // if the left cell has the left direction, the gap is being extended
                        // otherwise it's being opened.
                        (Dir::Left, false) => self.gap_extend_penalty,
                        (Dir::Left, true) => self.end_gap_extend_penalty,
                        (_, false) => self.gap_penalty,
                        (_, true) => self.end_gap_penalty,
                    };
                // Calculate score from up
                let up_cell = table[[i + 1, j]];
                // `i == target.len() - 1` means that we are currently in the right column,
                // therefore the gap is an end gap.
                // if the left cell has the up direction, the gap is being extended
                // otherwise it's being opened.
                let up_score = up_cell.score
                    + match (up_cell.dir, i == target.len() - 1) {
                        (Dir::Up, false) => self.gap_extend_penalty,
                        (Dir::Up, true) => self.end_gap_extend_penalty,
                        (_, false) => self.gap_penalty,
                        (_, true) => self.end_gap_penalty,
                    };
                // Calculate the diagonal score
                let diag_score = table[[i, j]].score
                    // It's either match or mismatch
                    + if target[i] == query[j] {
                        self.match_score
                    } else {
                        self.mismatch_score
                    };
                // Find the maximal score and record the direction from which it came
                let new_score = if left_score <= diag_score && up_score <= diag_score {
                    Score::new(Dir::Diagonal, diag_score)
                } else if left_score <= up_score {
                    Score::new(Dir::Up, up_score)
                } else {
                    Score::new(Dir::Left, left_score)
                };
                table[[i + 1, j + 1]] = new_score;
            }
        }

        table
    }
}

impl Default for Aligner {
    fn default() -> Self {
        Aligner {
            match_score: 1,
            mismatch_score: -1,
            gap_penalty: -100,
            gap_extend_penalty: -10,
            end_gap_penalty: -2,
            end_gap_extend_penalty: -1,
        }
    }
}

#[cfg(test)]
mod test_super {
    use super::*;
}
