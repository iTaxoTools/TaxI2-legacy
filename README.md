# TaxI2
Calculates genetic differences between DNA sequences

## Scores for alignment
The file `data/scores.tab` contains the scores used in the sequence alignment.
Each line has the format:
```
score_identifier<Tab>value
```

The possible scores are:
* `gap penalty`: Score to open a gap in the middle of a sequence
* `gap extend penalty`: Score to extend an existing gap in the middle of a sequence
* `end gap penalty`: Score to create a gap at an end of a sequence.
* `end gap extend penalty`: Score to extend a gap at an end of a sequence.
* `match score`: Score for matching nucleotides
* `mismatch score`: Score for non-matching nucleotides
