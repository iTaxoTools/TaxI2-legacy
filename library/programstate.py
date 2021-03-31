from typing import Union, TextIO, Iterator, Tuple, Any, Dict, Optional
from library.fasta import Fastafile
from library.genbank import GenbankFile
from library.seq import PDISTANCE, NDISTANCES, seq_distances_ufunc, seq_distances_aligned_ufunc, make_aligner, show_alignment
from multiprocessing import Pool, cpu_count
import tkinter as tk
import pandas as pd
import numpy as np
import networkx as nx
import os
import math
import re
import warnings
import datetime
import time
import gc
import sys
import itertools

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

    chunk_size = 100

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        raise NotImplementedError


    @staticmethod
    def rename_columns(name: str) -> str:
        regex = "([a-zA-Z0-9_-]+)"
        name = ''.join(re.findall(regex, name))
        if "sequence" in name:
            return "sequence"
        else:
            try:
                return FileFormat.rename_dict[name]
            except KeyError:
                return name

    def load_chunks(self, filepath_or_buffer: Union[str, TextIO]) -> Iterator[pd.DataFrame]:
        raise NotImplementedError

class TabFormat(FileFormat):
    """
    Format for delimiter-separated table
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        with open(filepath_or_buffer, errors='replace') as infile:
            return self.process_table(pd.read_csv(infile, sep='\t', dtype=str))

    @staticmethod
    def process_table(table: pd.DataFrame) -> pd.DataFrame:
        try:
            return table.rename(
                    columns=str.casefold).rename(columns=FileFormat.rename_columns)[['seqid', 'specimen_voucher', 'species', 'sequence']].drop_duplicates()
        except KeyError as ex:
            raise ValueError(
                "'seqid', 'specimen_voucher', 'species' or 'organism', or 'sequence' column is missing") from ex

    def load_chunks(self, filepath_or_buffer: Union[str, TextIO]) -> Iterator[pd.DataFrame]:
        with open(filepath_or_buffer, errors='replace') as infile:
            tables = pd.read_csv(infile, sep='\t', dtype=str, chunksize=FileFormat.chunk_size)
            for table in tables:
                yield self.process_table(table)


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
        return pd.DataFrame(([record['seqid'], record['sequence']] for record in records()), columns=['seqid', 'sequence']).drop_duplicates(subset='seqid')

    def load_chunks(self, filepath_or_buffer: Union[str, TextIO]) -> Iterator[pd.DataFrame]:
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer, errors='replace') as infile:
                for chunk in self._load_chunks(infile):
                    yield chunk
        else:
            for chunk in self._load_chunks(filepath_or_buffer):
                yield chunk

    def _load_chunks(self, file: TextIO) -> Iterator[pd.DataFrame]:
        _, records_gen = Fastafile.read(file)
        records = records_gen()
        while True:
            chunk = itertools.islice(records, FileFormat.chunk_size)
            table = pd.DataFrame(([record['seqid'], record['sequence']] for record in chunk), columns=['seqid', 'sequence'])
            if len(table) == 0:
                break
            yield table.drop_duplicates(subset='seqid').copy()




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
        self.root = root
        self.input_format_name = tk.StringVar(root, value="Tabfile")
        self.already_aligned = tk.BooleanVar(root, value=False)
        self.distance_options = tuple(tk.BooleanVar(root, value=False)
                                      for _ in range(NDISTANCES))
        self.distance_options[PDISTANCE].set(True)
        self.reference_comparison = tk.BooleanVar(root, value=True)
        self.print_alignments = tk.BooleanVar(root, value=False)
        self.perform_clustering = tk.BooleanVar(root, value=False)
        self.cluster_distance = tk.StringVar(root, value=distances_names[PDISTANCE])
        self.cluster_size = tk.StringVar(root, value='0.05')
        self.output_dir = output_dir

    def show_progress(self, message: str) -> None:
        full_message = f"{time.monotonic() - self.start_time:.1f}s: {message}\n"
        print(full_message)
        self.root.show_progress(full_message)

    @property
    def input_format(self) -> FileFormat:
        return ProgramState.formats[self.input_format_name.get()]()

    def output_name(self, description: str) -> str:
        return os.path.join(self.output_dir, description.replace(" ", "_") + ".txt")

    def output(self, description: str, table: pd.DataFrame, **kwargs) -> None:
        out_name = self.output_name(description)
        with open(out_name, mode='w') as outfile:
            print(description, file=outfile)
            table.to_csv(
                outfile, sep='\t', line_terminator='\n', float_format="%.4g", **kwargs)

    def process(self, input_file: str) -> None:
        self.start_time = time.monotonic()
        if self.input_format_name.get() == "Genbank" and self.already_aligned.get():
            raise ValueError(
                "'Already aligned' option is not allowed for the Genbank format.")
        table = self.input_format.load_table(input_file)
        self.species_analysis = "species" in table.columns
        table.set_index("seqid", inplace=True)
        if not self.already_aligned.get():
            table["sequence"] = normalize_sequences(table["sequence"])

        if self.reference_comparison.get():
            self.reference_comparison_process(table, reference_file)
            return

        if self.print_alignments.get():
            with open(os.path.join(self.output_dir, "taxi2_alignments.txt"), "w") as alignment_file:
                print_alignments(table["sequence"], alignment_file)
            self.show_progress("Alignment is printed")

        distance_table = make_distance_table(
            table, self.already_aligned.get())
        if self.species_analysis:
            species_table = pd.DataFrame(table["species"])
        else:
            species_table = None
        del table
        self.show_progress("Distance calcution")

        # The table of most similar sequences
        self.output(f"Most similar sequences", table_closest(distance_table))
        self.show_progress("Table of closest sequence")

        # The matrix of distances between seqids (input order)
        for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
            kind_name = distances_short_names[kind]
            out_table = distance_table[["seqid (query 1)", "seqid (query 2)", kind_name]].set_index(["seqid (query 1)", "seqid (query 2)"]).squeeze().unstack()
            out_table.columns.names = [''] * len(out_table.columns.names)
            self.output(f"{distances_names[kind]} between sequences", out_table, index_label=False)
            del out_table

        self.show_progress("Seqid distance table 1")

        # The matrix of distances between seqids (alphabetical order)
        for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
            kind_name = distances_short_names[kind]
            out_table = distance_table[["seqid (query 1)", "seqid (query 2)", kind_name]].set_index(["seqid (query 1)", "seqid (query 2)"]).squeeze().unstack().sort_index().sort_index(axis=1)
            out_table.columns.names = [''] * len(out_table.columns.names)
            self.output(f"{distances_names[kind]} between sequences (Alphabetical order)", out_table, index_label=False)
            del out_table

        self.show_progress("Seqid distance table 2")

        # clustering
        if self.perform_clustering.get():
            self.cluster_analysis(distance_table, species_table, os.path.basename(input_file))
            self.show_progress("Cluster analysis")

        if self.species_analysis:
            # The matrix of distances between seqids (order by species)
            for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
                kind_name = distances_short_names[kind]
                out_table = pd.Series(distance_table[kind_name].array, index=pd.MultiIndex.from_frame(distance_table[["species (query 1)", "species (query 2)", "seqid (query 1)", "seqid (query 2)"]])).squeeze().unstack(level=["species (query 2)", "seqid (query 2)"]).sort_index().sort_index(axis=1)
                out_table.columns.names = [''] * len(out_table.columns.names)
                self.output(f"{distances_names[kind]} between sequences (Ordered by species)", out_table, index_label=False)
                del out_table
            self.show_progress("Seqid distance table 3")

            genus1 = distance_table["species (query 1)"].str.split(pat=r' |_', n=1, expand=True).iloc[:,0]
            species1_index = distance_table.columns.get_loc("species (query 1)")
            distance_table.insert(loc=species1_index, column="genus (query 1)", value=genus1)

            genus2 = distance_table["species (query 2)"].str.split(pat=r' |_', n=1, expand=True).iloc[:,0]
            species2_index = distance_table.columns.get_loc("species (query 2)")
            distance_table.insert(loc=species2_index, column="genus (query 2)", value=genus2)

            # The matrices of distances between species
            for kind in (kind for kind in range(NDISTANCES) if self.distance_options[kind].get()):
                self.show_progress(distances_names[kind])
                kind_name = distances_short_names[kind]
                species_statistics = distance_table[
                        ["species (query 1)", "species (query 2)", "genus (query 1)", "genus (query 2)", kind_name]
                        ].groupby(
                                ["species (query 1)", "species (query 2)"]
                                ).agg(
                                        genus_1 = pd.NamedAgg(column="genus (query 1)", aggfunc="first"),
                                        genus_2 = pd.NamedAgg(column="genus (query 2)", aggfunc="first"),
                                        min=pd.NamedAgg(column=kind_name, aggfunc='min'),
                                        mean=pd.NamedAgg(column=kind_name, aggfunc='mean'),
                                        max=pd.NamedAgg(column=kind_name, aggfunc='max'))
                species_statistics['minmax'] = species_statistics['min'].apply(format_float).str.cat(species_statistics['max'].apply(format_float), sep='-')
                species_statistics['mean_minmax'] = species_statistics['mean'].apply(format_float).str.cat(species_statistics['minmax'] + ")", sep=' (')

                square_species_statistics = species_statistics.unstack("species (query 2)")

                self.output(f"Mean {distances_names[kind]} between species", square_species_statistics['mean'], index_label=False)
                self.show_progress("Mean distance between species")

                self.output(
                    f"Minimum and maximum {distances_names[kind]} between species", square_species_statistics['minmax'], index_label=False)
                self.show_progress("Minimum and maximum distance between species")

                self.output(f"Mean, minimum and maximum {distances_names[kind]} between species", square_species_statistics['mean_minmax'], index_label=False)
                self.show_progress("Mean, minimum and maximum distance between species")

                del square_species_statistics

                with open(self.output_name(f"Mean, minimum and maximum intra-species {distances_names[kind]}"), mode='w') as outfile:
                    print(
                        f"Mean, minimum and maximum intra-species {distances_names[kind]}", file=outfile)
                    species_statistics.reset_index(inplace=True)
                    species_statistics.loc[species_statistics["species (query 1)"] == species_statistics["species (query 2)"]].to_csv(outfile, sep="\t", columns=["species (query 1)", "mean_minmax"], header=False, index=False, line_terminator="\n")
                    outfile.write('\n')
                self.show_progress("Mean, minimum and maximum intra-species distance")

                closest_sequence_indices = distance_table.loc[distance_table["species (query 1)"] != distance_table["species (query 2)"], ["species (query 1)", kind_name]].groupby("species (query 1)").idxmin()[kind_name].dropna()
                closest_sequences = distance_table.loc[closest_sequence_indices, ["species (query 1)", kind_name, "seqid (query 2)"]]

                with open(self.output_name(f"Closest sequence from different species with {distances_names[kind]}"), mode='w') as outfile:
                    print(
                        f"Closest sequence from different species with {distances_names[kind]}", file=outfile)
                    print(
                        "species\tdistance (closest sequence of different species)\tseqid (closest sequence of different species)", file=outfile)
                    closest_sequences.to_csv(outfile, sep='\t', float_format='%.4g', header=False, index=False, line_terminator='\n')
                    outfile.write('\n')

                del closest_sequences
                del closest_sequence_indices

                self.show_progress("Closest sequences")

                genera_statistics = species_statistics.groupby(["genus_1", "genus_2"]).agg({'min': 'min', 'mean': 'mean', 'max': 'max'})
                genera_mean_minmax = genera_statistics["mean"].map(format_float) + " (" + genera_statistics["min"].map(format_float) + "-" + genera_statistics["max"].map(format_float) + ")"

                self.output(f"Mean, minimum and maximum {distances_names[kind]} between genera", genera_mean_minmax.unstack(), index_label=False)
                self.show_progress("Mean, minimum and maximum distances between genera")

                genera_mean_minmax = genera_mean_minmax.reset_index(name="mean_minmax")
                intra_genus = genera_mean_minmax.loc[genera_mean_minmax["genus_1"] == genera_mean_minmax["genus_2"], ["genus_1", "mean_minmax"]]

                with open(self.output_name(f"Mean, minimum and maximum intra-genus {distances_names[kind]}"), mode='w') as outfile:
                    print(
                        f"Mean, minimum and maximum intra-genus {distances_names[kind]}", file=outfile)
                    intra_genus.to_csv(outfile, sep="\t", header=False, index=False, line_terminator="\n")
                    outfile.write('\n')
                self.show_progress("Mean, minimum and maximum intra-genus distance")

                del species_statistics
                del genera_statistics
                del genera_mean_minmax
                del intra_genus


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

            distance_table["species (query 1)"] = distance_table["species (query 1)"].str.split(pat=r' |_', n=1, expand=True).iloc[:,1]
            distance_table["species (query 2)"] = distance_table["species (query 2)"].str.split(pat=r' |_', n=1, expand=True).iloc[:,1]

        self.output("Summary statistics", distance_table, index=False)
        self.show_progress("Final table")


    def cluster_analysis(self, distance_table: pd.DataFrame, species_table: Optional[pd.DataFrame], input_file: str) -> None:
        with open(self.output_name("Cluster analysis"), mode='w') as output_file:
            # extracting options
            distance_kind = distances_short_names[distances_names.index(self.cluster_distance.get())]
            try:
                cluster_threshold = float(self.cluster_size.get())
            except Exception:
                warnings.warn(f"Invalid cluster threshold {self.cluster_size.get()}.\nUsing default: 0.3")
                cluster_threshold = 0.3

            # preparing the table
            nodes = distance_table["seqid (query 1)"].unique()
            distance_table = distance_table[["seqid (query 1)", "seqid (query 2)", distance_kind]].copy()
            distance_table.columns = ["seqid1", "seqid2", "distance"]

            # calculating components
            connected_table = distance_table.loc[(distance_table['distance'] < cluster_threshold)]
            graph = nx.from_pandas_edgelist(connected_table, source="seqid1", target="seqid2")
            graph.add_nodes_from(nodes)
            components = nx.connected_components(graph)

            print(f"Samples were clustered at a threshold of {cluster_threshold:.3g} ({cluster_threshold*100:.2g}%) uncorrected p-distance", file=output_file)

            # add cluster classification to the table
            cluster_of: Dict[str, int] = {}
            max_samples = 0
            min_samples = len(distance_table)
            for i, component in enumerate(components):
                print(f'Cluster{i+1}: {", ".join(component)}', file=output_file)
                min_samples = min(min_samples, len(component))
                max_samples = max(max_samples, len(component))
                for seqid in component:
                    cluster_of[seqid] = i
            num_clusters = i + 1
            distance_table["cluster1"] = distance_table["seqid1"].map(cluster_of)
            distance_table["cluster2"] = distance_table["seqid2"].map(cluster_of)
            print("\n", file=output_file)

            max_in_cluster_distances = distance_table.loc[distance_table["cluster1"] == distance_table["cluster2"]][["cluster1", "distance"]].groupby("cluster1").max()["distance"]

            print("Maximum intra-sample distance within clusters (marked with # if above specified threshold):", file=output_file)
            big_clusters = 0
            for cluster_i, distance in max_in_cluster_distances.items():
                if math.isnan(distance):
                    distance = 0
                if distance > cluster_threshold:
                    big_clusters += 1
                    print(f'Cluster{cluster_i+1}: {distance:.4g} #', file=output_file)
                else:
                    print(f'Cluster{cluster_i+1}: {distance:.4g}', file=output_file)

            output_file.write("\n")

            min_between_cluster_distance = distance_table.loc[distance_table["cluster1"] > distance_table["cluster2"]][["cluster1", "cluster2", "distance"]].groupby(["cluster1", "cluster2"]).min()["distance"].unstack()
            for i in range(num_clusters):
                min_between_cluster_distance.at[(i, i)] = 0
            min_between_cluster_distance.sort_index(axis=0, inplace=True)
            min_between_cluster_distance.index.name = ''
            min_between_cluster_distance.rename(index=(lambda i: f"Cluster{i+1}"), columns=(lambda i: f"Cluster{i+1}"), inplace=True)

            print("Minimum distance between clusters:\n", file=output_file)

            min_between_cluster_distance.to_csv(
                output_file, sep='\t', line_terminator='\n', float_format="%.4g", index=True)

            output_file.write('\n')
            print("Total number of clusters:", num_clusters, file=output_file)
            print("Number of clusters violating threshold for intra-cluster distances:", big_clusters, file=output_file)

            output_file.write('\n')
            print(f"A total of {num_clusters} clusters were found, containing between {min_samples} and {max_samples} samples.", file=output_file)

            output_file.write('\n')

            if species_table is not None:
                species_table["cluster"] = species_table.index.to_series().map(cluster_of) + 1
                species_table.drop_duplicates(inplace=True)

                print("Comparison of cluster assignment with species assignment in the input file:", file=output_file)
                exists_multicluster_species = False
                for (species, clusters) in species_table.groupby("species")["cluster"]:
                    if len(clusters) > 1:
                        print(f"Sequences of {species} are included in {len(clusters)} clusters:", ", ".join(clusters.astype(str)), file=output_file)
                        exists_multicluster_species = True
                if exists_multicluster_species:
                    print("Sequences of all other species are included in only one cluster, respectively", file=output_file)
                else:
                    print("Sequences of all species are included in only one cluster, respectively.", file=output_file)

                output_file.write('\n')

                print("List of clusters containing sequences of more than one species (according to species assignment in the input file):", file=output_file)
                exists_multispecies_cluster = False
                for (cluster, species) in species_table.groupby("cluster")["species"]:
                    if len(species) > 1:
                        print(f"Cluster {cluster} contains sequences of {len(species)} species:", ", ".join(species), file=output_file)
                        exists_multispecies_cluster = True
                if exists_multispecies_cluster:
                    print("All other clusters contain sequences of only a single species, respectively.", file=output_file)
                else:
                    print("All clusters contain sequences of only one species, respectively.", file=output_file)

                output_file.write('\n')


        time = datetime.datetime.now()
        with open(os.path.join(self.output_dir, "taxi2_cluster_" + time.strftime("%Y-%m-%dT%H%M%S") + ".spart"), mode='w') as spart_file:
            print("begin spart;", file=spart_file)
            print("project_name = taxi2_clustering;", file=spart_file)
            print("Date = " + time.astimezone().isoformat(), file=spart_file)
            print("N_spartitions = "+str(num_clusters)+" : "+spart_form(input_file) + ";", file=spart_file)
            print("N_individuals = "+str(len(cluster_of))+";", file=spart_file)
            print("N_subsets = "+str(num_clusters) + ";", file=spart_file)
            print(f"[Generated by a simple {cluster_threshold*100:.2g}% threshold clustering in taxi2]", file=spart_file)
            print("[WARNING: The sample names below may have been changed to fit SPART specification (only alphanumeric characters and _ )]", file=spart_file)
            print(f"[The following clusters included sequences differing by a distance above the threshold: {cluster_threshold:.3g}]", file=spart_file)
            print("Individual_assignment = ", file=spart_file, end='')
            for specimen, cluster in cluster_of.items():
                print(f"\n{spart_form(specimen)}: {cluster + 1}", file=spart_file, end='')
            print(";\n", file=spart_file)
            print("end;", file=spart_file)

    def reference_comparison_process(self, input_file: str, reference_file: str) -> None:
        self.start_time = time.monotonic()
        if self.input_format_name.get() == "Genbank":
            raise ValueError(f"Comparison with reference database is not implemented for the Genbank format")
        reference_table = self.input_format.load_table(reference_file)
        reference_table.set_index("seqid", inplace=True)
        if not self.already_aligned.get():
            reference_table["sequence"] = normalize_sequences(reference_table["sequence"])

        with open(self.output_name("Closest reference sequences"), mode="w") as outfile:
            header = True
            sequences_num = 0
            for table in self.input_format.load_chunks(input_file):
                table.set_index("seqid", inplace=True)
                if not self.already_aligned.get():
                    table["sequence"] = normalize_sequences(table["sequence"])
                distance_table = make_distance_table2(table, reference_table, self.already_aligned.get())
                pdistance_name = distances_short_names[PDISTANCE]
                indices_closest = distance_table[["seqid (query 1)", pdistance_name]].groupby("seqid (query 1)").idxmin()[pdistance_name].squeeze().dropna()
                closest_table = distance_table.loc[indices_closest].rename(columns=(lambda col: col.replace("query 2", "closest reference sequence")))
                closest_table = distance_table.loc[indices_closest].rename(columns=(lambda col: col.replace("query 1", "query")))
                closest_table.to_csv(outfile, sep='\t', line_terminator='\n', float_format="%.4g", header=header)
                header=False
                sequences_num += self.input_format.chunk_size
                outfile.flush()
                del closest_table
                del distance_table
                del table
                gc.collect()
                self.show_progress(f"{sequences_num} sequences processed")




def format_float(x: float) -> str:
    return f"{x:.4g}"


def spart_form(s: str) -> str:
    return "_".join(re.compile("[A-Za-z0-9_]+").findall(s))

def _seq_par(x):
    return seq_distances_ufunc.outer(x[0], x[1])

def _seq_par_aligned(x):
    return seq_distances_aligned_ufunc.outer(x[0], x[1])

def _lines_for_data_len(tar, ref):
    _max_comparisons = 10000
    _cores = cpu_count()
    _lines = math.ceil(_max_comparisons / ref)
    _opt_split = math.ceil(tar / _cores)
    _lines = min(_lines, _opt_split)
    return _lines

def make_distance_table(table: pd.DataFrame, already_aligned: bool) -> pd.DataFrame:
    """
    Takes a series of sequences with a multi-index and returns a square dataframe

    Index and columns of the dataframe are the same as the index of the series

    The entries are arrays of pairwise distances
    """
    print('> Calculating distance table:')
    _tasks = 1
    _ufunc = _seq_par_aligned if already_aligned else _seq_par
    _len = len(table["sequence"])
    _lines = _lines_for_data_len(_len, _len)
    _sections = math.ceil(_len / _lines)
    print('> Number of Chunks:', _sections)
    chunks = np.array_split(np.asarray(table["sequence"]), _sections)
    ref_table = np.asarray(table["sequence"])
    x = [(i, ref_table) for i in chunks]
    res_list = []
    with Pool(maxtasksperchild=_tasks) as p:
        print('> {} / {}'.format(0, _sections))
        res = p.imap(_ufunc, x)
        for i, res in enumerate(res):
            print('> {} / {}'.format(i+1, _sections))
            res_list.append(res)
            # print(res) # Write this to a temp file
    print('> Merging Chunks:', _sections)
    distance_array = np.concatenate(res_list)

    # prepare indices
    seqid1 = table.index.copy()
    seqid1.name = "seqid (query 1)"
    seqid2 = table.index.copy()
    seqid2.name = "seqid (query 2)"
    # create distance table
    distance_table = pd.DataFrame(distance_array, index=seqid1, columns=seqid2).stack().reset_index(name="distances")
    distance_table = distance_table[distance_table["seqid (query 1)"] != distance_table["seqid (query 2)"]]

    # add other columns
    table.drop(columns="sequence", inplace=True)
    distance_table = distance_table.join(table, on="seqid (query 1)")
    distance_table.rename(columns={col:(col + " (query 1)") for col in table.columns}, inplace=True)
    distance_table = distance_table.join(table, on="seqid (query 2)")
    distance_table.rename(columns={col:(col + " (query 2)") for col in table.columns}, inplace=True)

    # reorder columns
    col1 = [col for col in distance_table.columns if "query 1" in col]
    col2 = [col for col in distance_table.columns if "query 2" in col]
    distance_table = distance_table[col1 + col2 + ["distances"]]
    distance_table[distances_short_names] = pd.DataFrame(distance_table["distances"].to_list(), index=distance_table.index)
    distance_table.drop(columns="distances", inplace=True)
    gc.collect()
    return distance_table

def make_distance_table2(table: pd.DataFrame, reference_table: pd.DataFrame, already_aligned: bool) -> pd.DataFrame:
    """
    Takes a series of sequences with a multi-index and returns a square dataframe

    Index and columns of the dataframe are the same as the index of the series

    The entries are arrays of pairwise distances
    """
    print('> Calculating distance table:')
    _tasks = 1
    _map_chunksize = 1
    _ufunc = _seq_par_aligned if already_aligned else _seq_par
    _len_tar = len(table["sequence"])
    _len_ref = len(reference_table["sequence"])
    _lines = _lines_for_data_len(_len_tar, _len_ref)
    _sections = math.ceil(_len_tar / _lines)
    print('> Number of Chunks:', _sections)
    chunks = np.array_split(np.asarray(table["sequence"]), _sections)
    ref_table = np.asarray(reference_table["sequence"])
    x = [(i, ref_table) for i in chunks]
    res_list = []
    with Pool(maxtasksperchild=_tasks) as p:
        print('> {} / {}'.format(0, _sections))
        res = p.imap(_ufunc, x, chunksize=_map_chunksize)
        for i, res in enumerate(res):
            print('> {} / {}'.format(i+1, _sections))
            res_list.append(res)
            # print(res) # Write this to a temp file
    print('> Merging Chunks:', _sections)
    distance_array = np.concatenate(res_list)

    # prepare indices
    seqid1 = table.index.copy()
    seqid1.name = "seqid (query 1)"
    seqid2 = reference_table.index.copy()
    seqid2.name = "seqid (query 2)"
    # create distance table
    distance_table = pd.DataFrame(distance_array, index=seqid1, columns=seqid2).stack().reset_index(name="distances")
    distance_table = distance_table[distance_table["seqid (query 1)"] != distance_table["seqid (query 2)"]]

    # add other columns
    table.drop(columns="sequence", inplace=True)
    reference_table = reference_table.drop(columns="sequence")
    distance_table = distance_table.join(table, on="seqid (query 1)")
    distance_table.rename(columns={col:(col + " (query 1)") for col in table.columns}, inplace=True)
    distance_table = distance_table.join(reference_table, on="seqid (query 2)")
    distance_table.rename(columns={col:(col + " (query 2)") for col in reference_table.columns}, inplace=True)

    # reorder columns
    col1 = [col for col in distance_table.columns if "query 1" in col]
    col2 = [col for col in distance_table.columns if "query 2" in col]
    distance_table = distance_table[col1 + col2 + ["distances"]]
    distance_table[distances_short_names] = pd.DataFrame(distance_table["distances"].to_list(), index=distance_table.index)
    distance_table.drop(columns="distances", inplace=True)
    gc.collect()
    return distance_table



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
    aligner = make_aligner()
    for (seqid_target, target) in sequences.items():
        for (seqid_query, query) in sequences.items():
            print(f"{seqid_target} <-> {seqid_query}", file=alignment_file)
            show_alignment(aligner, target, query, alignment_file)


def table_closest(distance_table: pd.DataFrame) -> pd.DataFrame:
    pdistance_name = distances_short_names[PDISTANCE]
    indices_closest = distance_table[["seqid (query 1)", pdistance_name]].groupby("seqid (query 1)").idxmin()[pdistance_name].squeeze().dropna()
    closest_table = distance_table.loc[indices_closest].rename(columns=(lambda col: col.replace("query 2", "most similar sequence in the dataset")))
    return closest_table

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
