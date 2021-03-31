//! Aligmnent using Needleman-Wunsch algorithm.
mod score;
mod table;

use pyo3::prelude::pyclass;

use score::{Dir, Score};
use std::string::FromUtf8Error;
use table::Table;

/// Contains parameters for the Needleman-Wunsch algorithm.
#[pyclass]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Aligner {
    pub(super) match_score: i16,
    pub(super) mismatch_score: i16,
    pub(super) gap_penalty: i16,
    pub(super) gap_extend_penalty: i16,
    pub(super) end_gap_penalty: i16,
    pub(super) end_gap_extend_penalty: i16,
}

impl Aligner {
    /// Construct the Needleman-Wunsch table for `target` and `query`.
    pub fn align<'target, 'query>(
        &self,
        target: &'target [u8],
        query: &'query [u8],
    ) -> Alignment<'target, 'query> {
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

        Alignment {
            target,
            query,
            table,
        }
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

/// Contains two aligned sequences
pub struct Alignment<'target, 'query> {
    target: &'target [u8],
    query: &'query [u8],
    table: Table<Score>,
}

impl<'target, 'query> Alignment<'target, 'query> {
    /// Iterate backwards over the path in the Needleman-Wunsch table.
    pub fn iter<'alignment>(&'alignment self) -> AlignmentIter<'target, 'query, 'alignment> {
        AlignmentIter {
            current_row: self.table.height - 1,
            current_column: self.table.width - 1,
            target: self.target,
            query: self.query,
            table: &self.table,
        }
    }

    /// Iterate backwards over the path in the Needleman-Wunsch table, skips the end gaps
    pub fn common_path_iter<'alignment>(
        &'alignment self,
    ) -> AlignmentCommonIter<'target, 'query, 'alignment> {
        AlignmentCommonIter {
            back_gap: true,
            inner: self.iter(),
        }
    }

    /// Returns two strings representing aligned target and query respectively.
    /// # Errors
    /// Returns [Err] if invalid UTF-8 has been constructed
    pub fn as_strings(&self) -> Result<(String, String), FromUtf8Error> {
        // Allocate space for the result
        let (mut target_align, mut query_align) = {
            let len = self.target.len().max(self.query.len());
            (Vec::with_capacity(len), Vec::with_capacity(len))
        };
        // Collect the symbols
        for (target_c, query_c) in self.iter() {
            target_align.push(target_c);
            query_align.push(query_c);
        }

        // Reverse the result, since the iteration was backwards
        target_align.reverse();
        query_align.reverse();

        Ok((
            String::from_utf8(target_align)?,
            String::from_utf8(query_align)?,
        ))
    }

    /// Returns a string showing alignment in the same way as Biopython
    /// # Errors
    /// Returns [Err] if invalid UTF-8 has been constructed
    pub fn show_alignment(&self) -> Result<String, FromUtf8Error> {
        // Allocate space for the result
        let (mut target_align, mut correspondence, mut query_align) = {
            let len = self.target.len().max(self.query.len());
            (
                Vec::with_capacity(len),
                Vec::with_capacity(len),
                Vec::with_capacity(len),
            )
        };
        // Collect the symbols
        for (target_c, query_c) in self.iter() {
            target_align.push(target_c);
            query_align.push(query_c);
            let correspondence_c = match (target_c, query_c) {
                (b'-', _) => b'-',
                (_, b'-') => b'-',
                _ if target_c == query_c => b'|',
                _ => b'.',
            };
            correspondence.push(correspondence_c);
        }

        // Reverse the result, since the iteration was backwards
        target_align.reverse();
        correspondence.reverse();
        query_align.reverse();

        Ok(String::from_utf8(target_align)?
            + "\n"
            + &String::from_utf8(correspondence)?
            + "\n"
            + &String::from_utf8(query_align)?)
    }
}

/// Iterator over the alignment path
pub struct AlignmentIter<'target, 'query, 'alignment> {
    current_row: usize,
    current_column: usize,
    target: &'target [u8],
    query: &'query [u8],
    table: &'alignment Table<Score>,
}

impl<'target, 'query, 'alignment> AlignmentIter<'target, 'query, 'alignment> {
    // Returns the current cell of the alignment table
    fn current_cell(&self) -> &Score {
        &self.table[[self.current_column, self.current_row]]
    }

    // Returns the current char of target sequence
    //
    // # Panics
    // Panics if self.current_column == 0
    fn current_target(&self) -> u8 {
        self.target[self.current_column - 1]
    }

    // Returns the current char of query sequence
    //
    // # Panics
    // Panics if self.current_row == 0
    fn current_query(&self) -> u8 {
        self.query[self.current_row - 1]
    }

    // Moves to the previous cell of the alignment sequence
    // according to the current direction
    fn shift_back(&mut self) {
        match self.current_cell().dir {
            Dir::Diagonal => {
                self.current_row -= 1;
                self.current_column -= 1;
            }
            Dir::Up => {
                self.current_row -= 1;
            }
            Dir::Left => {
                self.current_column -= 1;
            }
        }
    }
}

impl<'target, 'query, 'alignment> Iterator for AlignmentIter<'target, 'query, 'alignment> {
    type Item = (u8, u8);
    // Yields the previous pair of aligned characters
    fn next(&mut self) -> Option<(u8, u8)> {
        // If current cell is the top-left cell,
        // the iteration is finished
        if self.current_row == 0 && self.current_column == 0 {
            return None;
        }
        let (target_c, query_c) = match self.current_cell().dir {
            Dir::Diagonal => (self.current_target(), self.current_query()),
            Dir::Up => (b'-', self.current_query()),
            Dir::Left => (self.current_target(), b'-'),
        };
        self.shift_back();
        Some((target_c, query_c))
    }
}

/// Iterates over the alignment path, skipping
pub struct AlignmentCommonIter<'target, 'query, 'alignment> {
    back_gap: bool,
    inner: AlignmentIter<'target, 'query, 'alignment>,
}

impl<'target, 'query, 'alignment> AlignmentCommonIter<'target, 'query, 'alignment> {
    /// Returns the next pair of symbols after the back gap
    fn next_after_back_gap(&mut self) -> Option<<AlignmentIter as Iterator>::Item> {
        if !self.back_gap {
            return self.inner.next();
        }
        loop {
            if let Dir::Diagonal = self.inner.current_cell().dir {
                self.back_gap = false;
                return self.inner.next();
            }
            self.inner.next();
        }
    }
}

impl<'target, 'query, 'alignment> Iterator for AlignmentCommonIter<'target, 'query, 'alignment> {
    type Item = <AlignmentIter<'target, 'query, 'alignment> as Iterator>::Item;
    fn next(&mut self) -> Option<Self::Item> {
        if self.inner.current_row == 0 || self.inner.current_column == 0 {
            None
        } else {
            self.next_after_back_gap()
        }
    }
}

#[cfg(test)]
mod test_super {
    use super::*;

    fn test_aligner() -> Aligner {
        Aligner {
            match_score: 1,
            mismatch_score: -1,
            gap_penalty: -10,
            gap_extend_penalty: -5,
            end_gap_penalty: -2,
            end_gap_extend_penalty: -1,
        }
    }

    #[test]
    fn test_align() -> Result<(), FromUtf8Error> {
        let aligned = vec![
            ("ttcctcgt", "cattctcgt", "-ttcctcgt", "cattctcgt"),
            ("cccgtgcg", "acgtccg", "cccgtgcg", "-acgtccg"),
            ("gtcattag", "gtcatttag", "gtcattag-", "gtcatttag"),
            ("cccaaggt", "cacaaggt", "cccaaggt", "cacaaggt"),
            ("cctagtag", "cctatt", "cctagtag", "cctatt--"),
            ("tgatcagg", "tgactgag", "tgatcagg", "tgactgag"),
            ("tcagttct", "tcagttcgct", "tcagttct--", "tcagttcgct"),
            ("catccaac", "cttcca", "catccaac", "cttcca--"),
        ];
        let aligner = test_aligner();

        for (target, query, target_align_test, query_align_test) in aligned {
            let (target_align, query_align) = aligner
                .align(target.as_bytes(), query.as_bytes())
                .as_strings()?;
            assert_eq!(target_align, target_align_test);
            assert_eq!(query_align, query_align_test);
        }
        Ok(())
    }

    #[test]
    fn test_next_after_back_gap_no_gap() {
        let target = [1, 2, 3, 4, 5, 6, 7];
        let query = [1, 2, 3, 11, 5, 6, 8];
        let alignment = test_aligner().align(&target, &query);
        assert_eq!(
            alignment.common_path_iter().next_after_back_gap(),
            Some((7, 8))
        )
    }

    #[test]
    fn test_next_after_back_gap() {
        let target = [1, 2, 3, 4, 5, 6, 7];
        let query = [2, 3, 4, 11];
        let alignment = test_aligner().align(&target, &query);
        assert_eq!(
            alignment.common_path_iter().next_after_back_gap(),
            Some((5, 11))
        )
    }

    #[test]
    fn test_common_path_iter() {
        let target = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17];
        let query = [4, 5, 6, 7, 8, 10, 11, 12, 13, 14];
        let alignment = test_aligner().align(&target, &query);
        assert_eq!(
            alignment.common_path_iter().collect::<Vec<_>>(),
            vec![
                (14, 14),
                (13, 13),
                (12, 12),
                (11, 11),
                (10, 10),
                (9, b'-'),
                (8, 8),
                (7, 7),
                (6, 6),
                (5, 5),
                (4, 4),
            ]
        )
    }
}
