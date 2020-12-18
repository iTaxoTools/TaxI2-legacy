from typing import Union, TextIO, Iterator, Tuple, Any
from library.fasta import Fastafile
from library.genbank import GenbankFile
from library.record import Record
from library.seq import PDISTANCE, JUKES_CANTOR, KIMURA_2P, PDISTANCE_GAPS, NDISTANCES, seq_distances_ufunc, seq_distances_aligned_ufunc, aligner
import tkinter as tk
import pandas as pd
import numpy as np
import os

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


class FastaFormat(FileFormat):
    """
    Format for fasta files
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer) as infile:
                return self._load_table(infile)
        else:
            return self._load_table(infile)

    def _load_table(self, file: TextIO) -> pd.DataFrame:
        _, records = Fastafile.read(file)
        return pd.DataFrame(([record['seqid'], record['sequence']] for record in records()), columns=['seqid', 'sequence'])


class GenbankFormat(FileFormat):
    """
    Format for fasta files
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer) as infile:
                return self._load_table(infile)
        else:
            return self._load_table(infile)

    def _load_table(self, file: TextIO) -> pd.DataFrame:
        _, records = GenbankFile.read(file)
        try:
            return pd.DataFrame((record._fields for record in records())).rename(columns=str.casefold).rename(columns={'organism': 'species'})[['seqid', 'specimen_voucher', 'species', 'sequence']]
        except KeyError as ex:
            raise ValueError(f"{str(ex)} is missing") from ex


class ProgramState():
    """
    Encapsulates the state of TaxI2
    """

    formats = dict(
        Tabfile=TabFormat,
        Fasta=FastaFormat,
        Genbank=GenbankFormat
    )

    def __init__(self, root: tk.Tk) -> None:
        self.input_format_name = tk.StringVar(root, value="Tabfile")
        self.already_aligned = tk.BooleanVar(root, value=False)
        self.distance_options = tuple(tk.BooleanVar(root, value=False)
                                      for _ in range(NDISTANCES))
        self.distance_options[PDISTANCE].set(True)
        self.print_alignments = tk.BooleanVar(root, value=False)

    @property
    def input_format(self) -> FileFormat:
        return ProgramState.formats[self.input_format_name.get()]()

    def process(self, input_file: str, output_file: str) -> None:
        if self.input_format_name.get() == "Genbank" and self.already_aligned.get():
            raise ValueError(
                "'Already aligned' option is not allowed for the Genbank format.")
        table = self.input_format.load_table(input_file)
        sequences = table.set_index(
            [column for column in table.columns if column != 'sequence']).squeeze()
        sequences = normalize_sequences(sequences)
        if not isinstance(sequences, pd.Series):
            raise ValueError(
                "Extracting sequence from the table failed for some reason")
        distance_table = make_distance_table(
            sequences, self.already_aligned.get())
        with open(output_file, "w") as outfile:
            # The table of most similar sequences
            print(f"Most similar sequences", file=outfile)
            table_closest(distance_table).to_csv(
                outfile, sep='\t', line_terminator='\n', float_format="%.4g")
            outfile.write('\n')

            # The matrix of distance between seqids (input order)
            for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
                print(
                    f"{distances_names[kind]} between sequences", file=outfile)
                distance_table.pipe(select_distance, kind).pipe(
                    seqid_distance_table).to_csv(outfile, sep='\t', line_terminator='\n', float_format="%.4g")
                outfile.write('\n')

            # The matrix of distance between seqids (alphabetical order)
            for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
                print(
                    f"{distances_names[kind]} between sequences (Alphabetical order)", file=outfile)
                distance_table.pipe(select_distance, kind).pipe(
                    seqid_distance_table).sort_index().sort_index(axis='columns').to_csv(outfile, sep='\t', line_terminator='\n', float_format="%.4g")
                outfile.write('\n')

        if self.print_alignments.get():
            with open(alignment_file_name(output_file), "w") as alignment_file:
                print_alignments(sequences, alignment_file)


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


def normalize_sequences(sequences: pd.Series) -> pd.Series:
    return sequences.str.upper().str.replace("?", "N").str.replace("-", "")


def alignment_file_name(output_file: str) -> str:
    output_file_base, output_file_ext = os.path.splitext(output_file)
    return output_file_base + "_alignments" + output_file_ext


def print_alignments(sequences: pd.Series, alignment_file: TextIO) -> None:
    sequences = sequences.copy()
    sequences.index = sequences.index.get_level_values('seqid')
    for (seqid_target, target) in sequences.items():
        for (seqid_query, query) in sequences.items():
            print(f"{seqid_target} <-> {seqid_query}", file=alignment_file)
            alignment = aligner.align(target, query)[0]
            print(alignment, file=alignment_file)


def table_closest(distance_table: pd.DataFrame) -> pd.DataFrame:
    columns = tuple(map(lambda s: s + " (most similar sequence in the dataset)", distance_table.index.names)
                    ) + ("p-distance", "JC-distance", "K2P-distance", "p-distance with gaps")
    index_rename = {name: (name + " (query)")
                    for name in distance_table.index.names}
    return pd.DataFrame(iterate_closest(distance_table), index=distance_table.index, columns=columns).reset_index().rename(columns=index_rename)


def iterate_closest(distance_table: pd.DataFrame) -> Iterator[Tuple[Any]]:
    for column in distance_table.columns:
        idx = distance_table[column].drop(labels=column).map(
            lambda arr: arr[PDISTANCE]).idxmin()
        yield idx + tuple(distance_table.loc[(column, idx)])
