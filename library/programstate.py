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
