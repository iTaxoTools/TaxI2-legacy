from typing import Union, TextIO
from library.seq import PDISTANCE, JUKES_CANTOR, KIMURA_2P, PDISTANCE_GAPS, NDISTANCES, seq_distances_ufunc, seq_distances_aligned_ufunc
import tkinter as tk
import pandas as pd
import numpy as np


class FileFormat():
    """
    Interface for file formats supported by the program
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        raise NotImplementedError


class TabFormat(FileFormat):
    """
    Format for delimiter-separated table
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        return pd.read_csv(filepath_or_buffer, sep='\t').rename(
            columns=str.casefold).rename(columns={'organism': 'species'})['seqid', 'specimen_voucher', 'species']


class ProgramState():
    """
    Encapsulates the state of TaxI2
    """

    formats = dict(
        Tabfile=TabFormat
    )

    def __init__(self, root: tk.Tk) -> None:
        self.input_format_name = tk.StringVar(root, value="Tabfile")
        self.already_aligned = tk.BooleanVar(root, value=False)
        self.distance_options = tuple(tk.BooleanVar(root, value=False)
                                      for _ in range(NDISTANCES))
        self.distance_options[PDISTANCE].set(True)

    @property
    def input_format(self) -> FileFormat:
        return ProgramState.formats[self.input_format_name.get()]()

    def process(self, input_file: str, output_file: str) -> None:
        table = self.input_format.load_table(input_file).pivot(
            index=["seqid", "specimen_voucher", "species"]).squeeze()
        distance_table = distance_table(table, self.already_aligned.get())
        distance_table.pipe(select_distance, PDISTANCE).pipe(
            seqid_distance_table).to_csv(output_file, sep='\t', line_terminator='\n')


def distance_table(sequences: pd.Series, already_aligned: bool) -> pd.DataFrame:
    """
    Takes a series of sequences with a multi-index and returns a square dataframe

    Index and columns of the dataframe are the same as the index of the series

    The entries are arrays of pairwise distances
    """
    if already_aligned:
        distance_array = seq_distances_ufunc.outer(
            np.asarray(sequences), np.asarray(sequences))
    else:
        distance_array = seq_distances_aligned_ufunc.outer(
            np.asarray(sequences), np.asarray(sequences))
    return pd.DataFrame(distance_array, index=sequences.index, columns=sequences.index)


def select_distance(distance_table: pd.DataFrame, kind: int) -> pd.DataFrame:
    """
    Transforms table of arrays of distances into the table of distances

    kind should be one of the distance selection constants
    """
    if kind >= NDISTANCES:
        raise ValueError(f"{kind} doesn't corresponds to a valid distance")
    return distance_table.applymap(lambda arr: arr[kind])


def seqid_distance_table(distance_table: pd.DataFrame) -> pd.DataFrame:
    """
    Changes the index of the table to only seqid
    """
    return distance_table.reindex(index=distance_table.index.get_level_values('seqid'), columns=distance_table.columns.get_level_values('seqid'))
