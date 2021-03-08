from typing import Union, TextIO, Iterator, Tuple, Any, Dict
from library.fasta import Fastafile
from library.genbank import GenbankFile
from library.seq import PDISTANCE, NDISTANCES, seq_distances_ufunc, seq_distances_aligned_ufunc, aligner
import tkinter as tk
import pandas as pd
import numpy as np
import networkx as nx
import os
import math
import re
import warnings

distances_names = ["pairwise uncorrected distance", "Jukes-Cantor distance",
                   "Kimura-2-Parameter distance", "pairwise uncorrected distance counting gaps"]
distances_short_names = ['p-distance', 'JC distance',
                         'K2P distance', 'p-distance with gaps']


class FileFormat():
    """
    Interface for file formats supported by the program
    """

    rename_dict = {
            'organism': 'species',
            'specimen_identifier': 'specimen_voucher',
            'specimen identifier': 'specimen_voucher',
            'specimen voucher': 'specimen_voucher'
            }

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        raise NotImplementedError


    @staticmethod
    def rename_columns(name: str) -> str:
        if "sequence" in name:
            return "sequence"
        else:
            try:
                return FileFormat.rename_dict[name]
            except KeyError:
                return name

class TabFormat(FileFormat):
    """
    Format for delimiter-separated table
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        try:
            with open(filepath_or_buffer, errors='replace') as infile:
                return pd.read_csv(infile, sep='\t', dtype=str).rename(
                    columns=str.casefold).rename(columns=FileFormat.rename_columns)[['seqid', 'specimen_voucher', 'species', 'sequence']].drop_duplicates()
        except KeyError as ex:
            raise ValueError(
                "'seqid', 'specimen_voucher', 'species' or 'organism', or 'sequence' column is missing") from ex


class FastaFormat(FileFormat):
    """
    Format for fasta files
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer, errors='replace') as infile:
                return self._load_table(infile)
        else:
            return self._load_table(filepath_or_buffer)

    def _load_table(self, file: TextIO) -> pd.DataFrame:
        _, records = Fastafile.read(file)
        return pd.DataFrame(([record['seqid'], record['sequence']] for record in records()), columns=['seqid', 'sequence']).drop_duplicates()


class GenbankFormat(FileFormat):
    """
    Format for fasta files
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer, errors='replace') as infile:
                return self._load_table(infile)
        else:
            return self._load_table(filepath_or_buffer)

    def _load_table(self, file: TextIO) -> pd.DataFrame:
        _, records = GenbankFile.read(file)
        try:
            return pd.DataFrame((record._fields for record in records())).rename(columns=str.casefold).rename(columns=FileFormat.rename_columns)[['seqid', 'specimen_voucher', 'species', 'sequence']].drop_duplicates()
        except KeyError as ex:
            raise ValueError(f"{str(ex)} is missing") from ex


class ProgramState():
    """
    Encapsulates the state of TaxI2
    """

    SUMMARY_STATISTICS_NAME = "Summary_statistics.txt"

    formats = dict(
        Tabfile=TabFormat,
        Fasta=FastaFormat,
        Genbank=GenbankFormat
    )

    def __init__(self, root: tk.Misc, output_dir: str) -> None:
        self.input_format_name = tk.StringVar(root, value="Tabfile")
        self.already_aligned = tk.BooleanVar(root, value=False)
        self.distance_options = tuple(tk.BooleanVar(root, value=False)
                                      for _ in range(NDISTANCES))
        self.distance_options[PDISTANCE].set(True)
        self.print_alignments = tk.BooleanVar(root, value=False)
        self.cluster_distance = tk.StringVar(root, value=distances_names[PDISTANCE])
        self.cluster_size = tk.StringVar(root, value='0.3')
        self.output_dir = output_dir

    @property
    def input_format(self) -> FileFormat:
        return ProgramState.formats[self.input_format_name.get()]()

    def output_name(self, description: str) -> str:
        return os.path.join(self.output_dir, description.replace(" ", "_") + ".txt")

    def output(self, description: str, table: pd.DataFrame, index: bool = True) -> None:
        out_name = self.output_name(description)
        with open(out_name, mode='w') as outfile:
            print(description, file=outfile)
            table.to_csv(
                outfile, sep='\t', line_terminator='\n', float_format="%.4g", index=index)

    def process(self, input_file: str) -> None:
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

        # The table of most similar sequences
        self.output(f"Most similar sequences", table_closest(distance_table))

        # The matrix of distances between seqids (input order)
        for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
            self.output(f"{distances_names[kind]} between sequences", distance_table.pipe(select_distance, kind).pipe(
                seqid_distance_table))

        # The matrix of distances between seqids (alphabetical order)
        for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
            self.output(f"{distances_names[kind]} between sequences (Alphabetical order)", distance_table.pipe(select_distance, kind).pipe(
                seqid_distance_table).sort_index().sort_index(axis='columns'))

        # clustering
        self.cluster_analysis(distance_table)

        if 'species' in distance_table.index.names:
            # The matrix of distances between seqids (order by species)
            for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
                self.output(f"{distances_names[kind]} between sequences (Ordered by species)", distance_table.pipe(select_distance, kind).pipe(
                    species_distance_table).sort_index(level=['species', 'seqid']).sort_index(axis='columns', level=['species', 'seqid']))

            # The matrices of distances between species
            for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
                this_distance_table = distance_table.pipe(
                    select_distance, kind).pipe(species_distance_table)
                mean_distances = this_distance_table.groupby(level='species').mean(
                ).groupby(axis=1, level='species').mean().sort_index().sort_index(axis='columns')
                min_distances = this_distance_table.groupby(level='species').min(
                ).groupby(axis=1, level='species').min().sort_index().sort_index(axis='columns')
                max_distances = this_distance_table.groupby(level='species').max(
                ).groupby(axis=1, level='species').max().sort_index().sort_index(axis='columns')

                species_distances = mean_distances.applymap(lambda mean: (mean,)).combine(
                    min_distances, series_append).combine(max_distances, series_append)

                self.output(f"Mean {distances_names[kind]} between species", species_distances.applymap(
                    lambda mean_min_max: mean_min_max[0]))

                self.output(
                    f"Minimum and maximum {distances_names[kind]} between species", species_distances.applymap(show_min_max))

                self.output(f"Mean, minimum and maximum {distances_names[kind]} between species", species_distances.applymap(
                    show_mean_min_max))

                with open(self.output_name(f"Mean, minimum and maximum intra-species {distances_names[kind]}"), mode='w') as outfile:
                    print(
                        f"Mean, minimum and maximum intra-species {distances_names[kind]}", file=outfile)
                    for species in mean_distances.index:
                        mean_min_max = (
                            mean_distances.at[species, species], min_distances.at[species, species], max_distances.at[species, species])
                        print(
                            f"{species}\t{show_mean_min_max(mean_min_max)}", file=outfile)
                    outfile.write('\n')

                with open(self.output_name(f"Closest sequence from different species with {distances_names[kind]}"), mode='w') as outfile:
                    print(
                        f"Closest sequence from different species with {distances_names[kind]}", file=outfile)
                    print(
                        "species\tdistance (closest sequence of different species)\tseqid (closest sequence of different species)", file=outfile)
                    for species, idxmin in find_closest_from_another(this_distance_table).items():
                        other_seqid, other_species, seqid_self = idxmin
                        distance = this_distance_table.at[(
                            seqid_self, species), (other_seqid, other_species)]
                        print(species, f"{distance:.4g}", other_seqid,
                              sep='\t', file=outfile)
                    outfile.write('\n')

                mean_distances_genera = mean_distances.groupby(select_genus).mean().groupby(
                    select_genus, axis=1).mean().sort_index().sort_index(axis=1)
                min_distances_genera = min_distances.groupby(select_genus).min().groupby(
                    select_genus, axis=1).min().sort_index().sort_index(axis=1)
                max_distances_genera = max_distances.groupby(select_genus).max().groupby(
                    select_genus, axis=1).max().sort_index().sort_index(axis=1)

                genera_distances = mean_distances_genera.applymap(lambda mean: (mean,)).combine(
                    min_distances_genera, series_append).combine(max_distances_genera, series_append)

                self.output(f"Mean, minimum and maximum {distances_names[kind]} between genera", genera_distances.applymap(
                    show_mean_min_max))

            distance_table.index.names = map(
                lambda s: s + ' (query 1)', distance_table.index.names)
            distance_table.columns.names = map(
                lambda s: s + ' (query 2)', distance_table.columns.names)
            distance_table = distance_table.stack(
                distance_table.columns.names)
            distance_table = distance_table.reset_index(name='distances')

            distance_table[['genus (query 1)', 'species (query 1)']] = distance_table['species (query 1)'].str.split(
                r' |_', expand=True, n=1)
            genus1_pos = distance_table.columns.get_loc(
                'species (query 1)')
            distance_table.insert(
                genus1_pos, 'genus (query 1)', distance_table.pop('genus (query 1)'))
            distance_table[['genus (query 2)', 'species (query 2)']] = distance_table['species (query 2)'].str.split(
                r' |_', expand=True, n=1)
            genus2_pos = distance_table.columns.get_loc(
                'species (query 2)')
            distance_table.insert(
                genus2_pos, 'genus (query 2)', distance_table.pop('genus (query 2)'))
            for kind in range(NDISTANCES):
                distance_table[distances_short_names[kind]] = distance_table['distances'].map(
                    lambda arr: arr[kind])
            distance_table.pop('distances')
            distance_table = distance_table.loc[distance_table['seqid (query 1)']
                                                != distance_table['seqid (query 2)']]
            same_species = distance_table['species (query 1)'] == distance_table['species (query 2)']
            same_genus = distance_table['genus (query 1)'] == distance_table['genus (query 2)']

            def comparison_type(same_species: bool, same_genus: bool) -> str:
                if same_genus:
                    if same_species:
                        return 'intra-species'
                    else:
                        return 'inter-species'
                else:
                    return 'inter-genus'
            comparison_type_pos = int((
                len(distance_table.columns) - NDISTANCES) / 2)
            distance_table.insert(comparison_type_pos, 'comparison_type', same_species.combine(
                same_genus, comparison_type))

            self.output("Summary statistics", distance_table, index=False)


        if self.print_alignments.get():
            with open(os.path.join(self.output_dir, "taxi2_alignments.txt"), "w") as alignment_file:
                print_alignments(sequences, alignment_file)

    def cluster_analysis(self, distance_table: pd.DataFrame) -> None:
        with open(self.output_name("Cluster analysis"), mode='w') as output_file:
            # extracting options
            distance_kind = distances_names.index(self.cluster_distance.get())
            try:
                cluster_threshold = float(self.cluster_size.get())
            except Exception:
                warnings.warn(f"Invalid cluster threshold {self.cluster_size.get()}.\nUsing default: 0.3")
                cluster_threshold = 0.3

            # preparing the table
            distance_table = seqid_distance_table(distance_table).stack()
            distance_table.index.names = ['seqid1', 'seqid2']
            distance_table = distance_table.map(lambda distances: distances[distance_kind]).reset_index(name="distance")

            # calculating components
            connected_table = distance_table.loc[(distance_table['distance'] < cluster_threshold) | (distance_table["seqid1"].eq(distance_table["seqid2"]))]
            components = nx.connected_components(nx.from_pandas_edgelist(connected_table, source="seqid1", target="seqid2"))

            # add cluster classification to the table
            cluster_of: Dict[str, int] = {}
            for i, component in enumerate(components):
                print(f'Cluster{i+1}: {", ".join(component)}', file=output_file)
                for seqid in component:
                    cluster_of[seqid] = i
            distance_table["cluster1"] = distance_table["seqid1"].map(cluster_of)
            distance_table["cluster2"] = distance_table["seqid2"].map(cluster_of)
            print("\n", file=output_file)

            max_in_cluster_distances = distance_table.loc[distance_table["cluster1"] == distance_table["cluster2"]][["cluster1", "distance"]].groupby("cluster1").max()["distance"]
            print(max_in_cluster_distances)

            print("Maximum intra-sample distance within clusters (marked with # if above specified threshold):", file=output_file)
            for cluster_i, distance in max_in_cluster_distances.items():
                if math.isnan(distance):
                    distance = 0
                print(f'Cluster{cluster_i+1}: {distance:.4g}{" #" if distance > cluster_threshold else ""}', file=output_file)


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
    for i in range(len(sequences)):
        distance_array[(i, i)] = np.full(NDISTANCES, np.nan)
    return pd.DataFrame(distance_array, index=sequences.index.copy(), columns=sequences.index.copy())


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


def species_distance_table(distance_table: pd.DataFrame) -> pd.DataFrame:
    """
    Changes the index of the table to seqid and species
    """
    result = distance_table.copy()
    to_drop = [level for level in result.index.names if level not in {
        'seqid', 'species'}]
    result.index = result.index.droplevel(to_drop)
    result.columns = result.columns.droplevel(to_drop)
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


def series_append(series_tuple: pd.Series, series_elem: pd.Series) -> pd.Series:
    """
    Combines a Series of tuples with a Series of scalars by appending componentwise
    """
    return series_tuple.combine(series_elem, lambda tuple, elem: tuple + (elem,))


def show_min_max(mean_min_max: Tuple[float, float, float]) -> str:
    if np.isnan(mean_min_max[0]):
        return ""
    else:
        return f"{mean_min_max[1]:.4g}-{mean_min_max[2]:.4g}"


def show_mean_min_max(mean_min_max: Tuple[float, float, float]) -> str:
    if np.isnan(mean_min_max[0]):
        return ""
    else:
        return f"{mean_min_max[0]:.4g} ({mean_min_max[1]:.4g}-{mean_min_max[2]:.4g})"


def find_closest_from_another(table: pd.DataFrame) -> pd.Series:
    """
    Returns as Series of tuples:
    species | (seqid_of_closest, species_of_closest, seqid_of_self)
    """
    table = table.copy()
    for lbl in table.index.levels[1]:
        table.loc[(slice(None), lbl), (slice(None), lbl)] = np.nan
    return table.stack(level=0).idxmin()


def select_genus(species: str) -> str:
    genus, _ = re.split(r'[ _]', species, maxsplit=1)
    return genus
