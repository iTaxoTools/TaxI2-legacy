//! Scores for the Needleman-Wunsch algorithm
use core::fmt::{self, Debug};

/// For marking direction in the Needleman-Wunsch algorithm.
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum Dir {
    /// The optimal score comes from the left.
    Left,
    /// The optimal score comes from the top.
    Up,
    /// The optimal score comes from the topleft.
    Diagonal,
}

impl Debug for Dir {
    /// Print the direction as an arrow.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let arrow = match self {
            Dir::Up => "↑",
            Dir::Left => "←",
            Dir::Diagonal => "↖",
        };
        f.write_str(arrow)
    }
}

/// Contents of a cell of the Needleman-Wunsch table.
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Score {
    /// Direction of the previous optimal score.
    pub dir: Dir,
    /// Current optimal score.
    pub score: i16,
}

impl Debug for Score {
    /// Formats the [Score] with the width 5.
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
    /// Contructs a [Score]
    pub fn new(dir: Dir, score: i16) -> Self {
        Score { dir, score }
    }
}
