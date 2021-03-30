//! Calculating distances between sequences

/// State for the distance calculation
pub struct AlignmentStats {
    total_length: usize,
    common_length: usize,
    total_gap_length: usize,
    transitions: usize,
    transversions: usize,
}

impl AlignmentStats {
    /// Zeroed state.
    pub fn new() -> Self {
        AlignmentStats {
            total_length: 0,
            common_length: 0,
            total_gap_length: 0,
            transitions: 0,
            transversions: 0,
        }
    }

    // Number of substitions calculated so far
    fn substitions(&self) -> usize {
        self.transversions + self.transitions
    }

    /// pairwise uncorrelated distance
    pub fn pdistance(&self) -> f64 {
        (self.substitions() as f64 / self.common_length as f64).abs()
    }

    /// pairwise uncorrelated distance with gaps
    pub fn pdistance_counting_gaps(&self) -> f64 {
        f64::abs(
            (self.substitions() as f64 + self.total_gap_length as f64) / self.total_length as f64,
        )
    }

    /// Jukes-Cantor distance
    pub fn jukes_cantor_distance(&self) -> f64 {
        let p = self.substitions() as f64 / self.common_length as f64;
        if p > 3.0 / 4.0 {
            f64::INFINITY
        } else {
            f64::abs(-(3.0 / 4.0) * f64::ln_1p(-(4.0 / 3.0) * p))
        }
    }

    /// Kimura's two parameter distance
    pub fn kimura2p_distance(&self) -> f64 {
        let p = self.transitions as f64 / self.common_length as f64;
        let q = self.transversions as f64 / self.common_length as f64;
        let distance =
            f64::abs(-(1.0 / 2.0) * f64::ln((1.0 - 2.0 * p - q) * f64::sqrt(1.0 - 2.0 * q)));
        if distance.is_nan() {
            f64::INFINITY
        } else {
            distance
        }
    }

    fn count_gap(&mut self) {
        self.total_length += 1;
        self.total_gap_length += 1;
    }

    fn count_match(&mut self) {
        self.total_length += 1;
        self.common_length += 1;
    }

    fn count_transition(&mut self) {
        self.total_length += 1;
        self.common_length += 1;
        self.transitions += 1;
    }
    fn count_transversion(&mut self) {
        self.total_length += 1;
        self.common_length += 1;
        self.transversions += 1;
    }

    /// Count `(x, y)` pair.
    pub fn update(&mut self, (x, y): (u8, u8)) {
        use NucleotideType::*;
        use SymbolType::*;
        match (classify(x), classify(y)) {
            (Gap, Nucleotide(_)) => self.count_gap(),
            (Nucleotide(_), Gap) => self.count_gap(),
            (Nucleotide(_), Nucleotide(_)) if x == y => self.count_match(),
            (Nucleotide(Purine), Nucleotide(Purine))
            | (Nucleotide(Pyrimidine), Nucleotide(Pyrimidine)) => self.count_transition(),
            (Nucleotide(_), Nucleotide(_)) => self.count_transversion(),
            _ => {}
        }
    }
}

enum SymbolType {
    Gap,
    Missing,
    Nucleotide(NucleotideType),
}

enum NucleotideType {
    Purine,
    Pyrimidine,
    Unknown,
}

fn classify(x: u8) -> SymbolType {
    use NucleotideType::*;
    use SymbolType::*;
    match x {
        b'-' => Gap,
        b'n' | b'N' | b'?' => Missing,
        b'a' | b'A' | b'g' | b'G' => Nucleotide(Purine),
        b'c' | b'C' | b't' | b'T' => Nucleotide(Pyrimidine),
        _ => Nucleotide(Unknown),
    }
}

#[cfg(test)]
mod test_super {
    use super::*;

    #[test]
    fn test_distance() {
        let target = "gg-ccnccta";
        let query = "ggaccaccaa";
        let mut alignment_stats = AlignmentStats::new();
        target
            .bytes()
            .zip(query.bytes())
            .for_each(|pair| alignment_stats.update(pair));
        assert_eq!(alignment_stats.pdistance(), 1.0 / 8.0);
        assert_eq!(alignment_stats.pdistance_counting_gaps(), 2.0 / 9.0);
    }
}
