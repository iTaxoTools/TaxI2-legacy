from Bio.Align import PairwiseAligner
import os


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
