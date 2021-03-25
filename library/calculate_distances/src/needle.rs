use core::fmt::{self, Debug};
use core::ops::Index;
use core::ops::IndexMut;

pub struct Aligner {
    match_score: i16,
    mismatch_score: i16,
    gap_penalty: i16,
    gap_extend_penalty: i16,
    end_gap_penalty: i16,
    end_gap_extend_penalty: i16,
}

impl Aligner {
    pub fn align_to_str(&self, target: &str, query: &str) -> (String, String) {
        let target = target.as_bytes();
        let query = query.as_bytes();
        let (mut target_align, mut query_align) = {
            let len = target.len().max(query.len());
            (Vec::with_capacity(len), Vec::with_capacity(len))
        };
        let table = self.align(target, query);
        let (mut i, mut j) = (target.len(), query.len());
        while i > 0 || j > 0 {
            let current = table[[i, j]];
            let (target_c, query_c) = match current.dir {
                Dir::Diagonal => (target[i - 1], query[j - 1]),
                Dir::Up => (b'-', query[j - 1]),
                Dir::Left => (target[i - 1], b'-'),
            };
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
            target_align.push(target_c);
            query_align.push(query_c);
        }
        target_align.reverse();
        query_align.reverse();
        (
            String::from_utf8(target_align).unwrap(),
            String::from_utf8(query_align).unwrap(),
        )
    }

    fn align(&self, target: &[u8], query: &[u8]) -> Table<Score> {
        let mut table = Table::<Score>::new(target.len() + 1, query.len() + 1);
        table.fill_default();
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

        for j in 0..query.len() {
            for i in 0..target.len() {
                let left_cell = table[[i, j + 1]];
                let left_score = left_cell.score
                    + match (left_cell.dir, j == query.len() - 1) {
                        (Dir::Left, false) => self.gap_extend_penalty,
                        (Dir::Left, true) => self.end_gap_extend_penalty,
                        (_, false) => self.gap_penalty,
                        (_, true) => self.end_gap_penalty,
                    };
                let up_cell = table[[i + 1, j]];
                let up_score = up_cell.score
                    + match (up_cell.dir, i == target.len() - 1) {
                        (Dir::Up, false) => self.gap_extend_penalty,
                        (Dir::Up, true) => self.end_gap_extend_penalty,
                        (_, false) => self.gap_penalty,
                        (_, true) => self.end_gap_penalty,
                    };
                let diag_score = table[[i, j]].score
                    + if target[i] == query[j] {
                        self.match_score
                    } else {
                        self.mismatch_score
                    };
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

#[derive(Copy, Clone, PartialEq, Eq)]
enum Dir {
    Left,
    Up,
    Diagonal,
}

impl Debug for Dir {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let arrow = match self {
            Dir::Up => "↑",
            Dir::Left => "←",
            Dir::Diagonal => "↖",
        };
        f.write_str(arrow)
    }
}

#[derive(Copy, Clone, PartialEq, Eq)]
struct Score {
    dir: Dir,
    score: i16,
}

impl Debug for Score {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}{:>4}", self.dir, self.score)
    }
}

impl Default for Score {
    fn default() -> Self {
        Score {
            dir: Dir::Diagonal,
            score: 0,
        }
    }
}

impl Score {
    fn new(dir: Dir, score: i16) -> Self {
        Score { dir, score }
    }
}

struct Table<T> {
    width: usize,
    height: usize,
    data: Vec<T>,
}

impl<T> Index<[usize; 2]> for Table<T> {
    type Output = T;

    fn index(&self, [x, y]: [usize; 2]) -> &T {
        if x < self.width && y < self.height {
            &self.data[self.width * y + x]
        } else {
            panic!(format!("{:?} is not a valid index", [x, y]))
        }
    }
}

impl<T> IndexMut<[usize; 2]> for Table<T> {
    fn index_mut(&mut self, [x, y]: [usize; 2]) -> &mut T {
        if x < self.width && y < self.height {
            &mut self.data[self.width * y + x]
        } else {
            panic!(format!("{:?} is not a valid index", [x, y]))
        }
    }
}

impl<T> Table<T> {
    fn new(width: usize, height: usize) -> Self {
        Table {
            width,
            height,
            data: Vec::with_capacity(width * height),
        }
    }
}

impl<T: Default> Table<T> {
    fn fill_default(&mut self) {
        self.data.resize_with(self.width * self.height, T::default);
    }
}

impl<T: Debug> Debug for Table<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for j in 0..self.height {
            for i in 0..self.width {
                write!(f, "{:?} ", self[[i, j]])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test_super {
    use super::*;

    #[test]
    fn test_align() {
        let aligner = Aligner::default();
        let table = aligner.align(b"gctagc", b"gctgc");
        println!("{:?}", table);
    }

    #[test]
    fn test_align_to_str() {
        let aligner = Aligner::default();
        let (target_align, query_align) = aligner.align_to_str("gctagc", "gctgc");
        dbg!(target_align);
        dbg!(query_align);
    }
}
