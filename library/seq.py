from typing import Tuple, Optional
from Bio.Align import PairwiseAligner
import os
import math
import numpy as np


with open(os.path.join('data', 'scores.tab')) as scores_file:
    scores_dict = {}
    for line in scores_file:
        score_name, _, val = line.partition('\t')
        try:
            scores_dict[score_name] = int(val)
        except ValueError as ex:
            raise ValueError(
                f"The value for '{score_name}' in data/scores.tab is not a number") from ex

try:
    GAP_PENALTY = scores_dict['gap penalty']
    GAP_EXTEND_PENALTY = scores_dict['gap extend penalty']
    END_GAP_PENALTY = scores_dict['end gap penalty']
    END_GAP_EXTEND_PENALTY = scores_dict['end gap extend penalty']
    MATCH_SCORE = scores_dict['match score']
    MISMATCH_SCORE = scores_dict['mismatch score']
except KeyError as ex:
    raise ValueError(f"'{ex.args[0]}' is missing in data/scores.tab") from ex


aligner = PairwiseAligner(match_score=MATCH_SCORE, mismatch_score=MISMATCH_SCORE, end_open_gap_score=END_GAP_PENALTY,
                          end_extend_gap_score=END_GAP_EXTEND_PENALTY, internal_open_gap_score=GAP_PENALTY, internal_extend_gap_score=GAP_EXTEND_PENALTY)


class Seq():
    """
    Wrapper for sequences
    """

    def __init__(self, sequence: str) -> None:
        self.seq = sequence

    def __len__(self) -> int:
        return len(self.seq)

    def align(self, other: 'Seq') -> 'Alignment':
        return Alignment(aligner.align(self.seq, other.seq)[0].aligned)


class Alignment():
    """
    Represents the alignment
    """

    def __init__(self, alignment: Tuple[Tuple[Tuple[int, int], ...], Tuple[Tuple[int, int], ...]]) -> None:
        self.alignment = alignment

    @classmethod
    def already_aligned(cls, seq: Seq) -> 'Alignment':
        return cls((((0, len(seq)),), ((0, len(seq)),)))


Fragment = Tuple[Tuple[int, int], Tuple[int, int]]


class AlignmentStats():
    """
    Summary of alignment information
    """

    def __init__(self) -> None:
        self.total_length = 0
        self.common_length = 0
        self.total_gap_length = 0
        self.transitions = 0
        self.transversions = 0
        self._previous_target_end: Optional[int] = None
        self._previous_query_end: Optional[int] = None

    @property
    def substitions(self) -> int:
        return self.transitions + self.transversions

    def calculate(self, alignment: Alignment, target: Seq, query: Seq) -> None:
        """
        Calculates the stats, when alignment is the alignment of target and query
        """
        (target_fragments, query_fragments) = alignment.alignment
        for fragment in zip(target_fragments, query_fragments):
            self.update(fragment, target.seq, query.seq)

    def pdistance(self) -> float:
        return self.substitions / self.common_length

    def pdistance_counting_gaps(self) -> float:
        return (self.substitions + self.total_gap_length) / self.common_length

    def jukes_cantor_distance(self) -> float:
        p = self.substitions / self.common_length
        return - (3 / 4) * math.log(1 - (4 / 3) * p)

    def kimura2p_distance(self) -> float:
        p = self.transitions / self.common_length
        q = self.transversions / self.common_length
        return - (1 / 2) * math.log((1 - 2 * p - q) * math.sqrt(1 - 2 * q))

    def update(self, frag: Fragment, target: str, query: str) -> None:
        self.update_total_length(frag)
        self.update_common_length(frag)
        self.update_total_gap_length(frag)
        self.update_transitions_and_transversions(frag, target, query)
        self._update_previous_target_end(frag)
        self._update_previous_query_end(frag)

    def update_total_length(self, frag: Fragment) -> None:
        ((target_start, target_end), (query_start, query_end)) = frag
        self.total_length += target_end - target_start + \
            self.target_gap(target_start) + self.query_gap(query_start)

    def update_common_length(self, frag: Fragment) -> None:
        ((target_start, target_end), _) = frag
        self.common_length += target_end - target_start

    def update_total_gap_length(self, frag: Fragment) -> None:
        ((target_start, _), (query_start, _)) = frag
        self.total_gap_length += self.target_gap(
            target_start) + self.query_gap(query_start)

    def update_transitions_and_transversions(self, frag: Fragment, target: str, query: str) -> None:
        (target_frag, query_frag) = frag
        for pair in zip(target[slice(target_frag)], query[slice(query_frag)]):
            if pair[0] == pair[0]:
                continue
            else:
                self.transitions += AlignmentStats.is_transition(pair)
                self.transversions += 1 - AlignmentStats.is_transition(pair)

    def _update_previous_target_end(self, frag: Fragment) -> None:
        self._previous_target_end = frag[0][1]

    def _update_previous_query_end(self, frag: Fragment) -> None:
        self._previous_query_end = frag[1][1]

    @staticmethod
    def is_transition(pair: Tuple[str, str]) -> bool:
        return pair in {('a', 'g'), ('g', 'a'), ('c', 't'), ('t', 'c')}

    def target_gap(self, target_start: int) -> int:
        if self._previous_target_end:
            return target_start - self._previous_target_end
        else:
            return 0

    def query_gap(self, query_start: int) -> int:
        if self._previous_query_end:
            return query_start - self._previous_query_end
        else:
            return 0


# the constants representing indices to extract corresponding distance from the result of seq_distances
PDISTANCE = 0
JUKES_CANTOR = 1
KIMURA_2P = 2
PDISTANCE_GAPS = 3
# number of distance options
NDISTANCES = 4


def seq_distances(target: str, query: str) -> np.array:
    """
    Returns array of 4 floats representing various distance between sequences.

    Index with constants to extract the distances:
    PDISTANCE - pairwise uncorrected distance
    JUKES_CANTOR - pairwise Jukes-Cantor distance
    KIMURA_2P - pairwise Kimura-2-Parameter distance
    PDISTANCE_GAPS - pairwise uncorrected distance including gaps
    """
    seq_target = Seq(target)
    seq_query = Seq(query)
    alignment = seq_target.align(seq_query)
    stats = AlignmentStats()
    stats.calculate(alignment, seq_target, seq_query)
    return np.array([stats.pdistance(), stats.jukes_cantor_distance(), stats.kimura2p_distance(), stats.pdistance_counting_gaps()])


def seq_distance_aligned(target: str, query: str) -> np.array:
    """
    Returns array of 4 floats representing various distance between sequences.

    Index with constants to extract the distances:
    PDISTANCE - pairwise uncorrected distance
    JUKES_CANTOR - pairwise Jukes-Cantor distance
    KIMURA_2P - pairwise Kimura-2-Parameter distance
    PDISTANCE_GAPS - pairwise uncorrected distance including gaps

    Expects aligned sequences
    """
    seq_target = Seq(target)
    seq_query = Seq(query)
    alignment = Alignment.already_aligned(seq_target)
    stats = AlignmentStats()
    stats.calculate(alignment, seq_target, seq_query)
    return np.array([stats.pdistance(), stats.jukes_cantor_distance(), stats.kimura2p_distance(), stats.pdistance_counting_gaps()])


seq_distances_ufunc: np.ufunc = np.frompyfunc(seq_distances, 2, 1)
seq_distances_aligned_ufunc: np.ufunc = np.frompyfunc(seq_distances, 2, 1)
