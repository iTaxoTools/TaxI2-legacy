from typing import Union, TextIO
from library.seq import PDISTANCE, JUKES_CANTOR, KIMURA_2P, PDISTANCE_GAPS, NDISTANCES, seq_distances_ufunc, seq_distances_aligned_ufunc
import tkinter as tk
import pandas as pd
import numpy as np

distances_names = ["pairwise uncorrelated distance", "Jukes-Cantor distance",
                   "Kimura-2-Parameter distance", "pairwise uncorrelated distance counting gaps"]


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
        try:
            return pd.read_csv(filepath_or_buffer, sep='\t').rename(
                columns=str.casefold).rename(columns={'organism': 'species'})[['seqid', 'specimen_voucher', 'species', 'sequence']]
        except KeyError as ex:
            raise ValueError(
                "'seqid', 'specimen_voucher', 'species' or 'organism', or 'sequence' column is missing") from ex


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
        with open(output_file, "w") as outfile:
            table = self.input_format.load_table(input_file)
            sequences = table.set_index(
                [column for column in table.columns if column != 'sequence']).squeeze()
            if not isinstance(sequences, pd.Series):
                raise ValueError(
                    "Extracting sequence from the table failed for some reason")
            distance_table = make_distance_table(
                sequences, self.already_aligned.get())
            for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
                print(
                    f"{distances_names[kind]} between sequences", file=outfile)
                distance_table.pipe(select_distance, kind).pipe(
                    seqid_distance_table).to_csv(outfile, sep='\t', line_terminator='\n', float_format="%.4g")
                outfile.write('\n')


def make_distance_table(sequences: pd.Series, already_aligned: bool) -> pd.DataFrame:
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
    result = distance_table.copy()
    result.index = distance_table.index.get_level_values('seqid')
    result.columns = distance_table.columns.get_level_values('seqid')
    return result
