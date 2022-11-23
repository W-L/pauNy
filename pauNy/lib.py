import os
import pathlib
import gzip
from itertools import groupby
from typing import List, Dict, TextIO, Tuple, Union

import numpy as np
import pandas as pd

# type aliases
NX_DICT = Dict[str, List[int]]
AUN_DICT = Dict[str, float]



class AssemblyCollection:

    def __init__(
        self,
        paths: List[Union[pathlib.Path, str]],
        reference_path: Union[pathlib.Path, str] = None,
        genome_size: int = None
    ):
        # make sure that input paths are pathlib.Paths
        self._assembly_paths = convert_paths(paths)
        self.ref_name = ''

        if reference_path:
            # if reference file is given, save its name and load the genome size
            self.ref_path = pathlib.Path(reference_path)
            self.ref_name = self.ref_path.name
            self._assembly_paths.append(self.ref_path)
            self.genome_size = Reference(path=self.ref_path).genome_size()
        elif genome_size:
            # if genome size estimate is given
            self.genome_size = genome_size
        else:
            self.genome_size = 0
        # initialise an assembly object for all input files
        self.assemblies = self._initialise_assemblies(assembly_paths=self._assembly_paths)




    def calculate_metrics(self) -> Tuple[NX_DICT, AUN_DICT]:
        """
        Calculate Nx and auN values for all assemblies in the collection

        :return: Dictionaries of Nx and auN values per input assembly
        """
        nx_values = self._calculate_Nx()
        aun_values = self._calculate_auN()
        return nx_values, aun_values


    def generate_dataframes(self, nx_values: NX_DICT = None, aun_values: AUN_DICT = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Generate pandas.DataFrames from the dictionaries of calculated values

        :param nx_values: Dict of Nx values per assembly
        :param aun_values: Dict of auN values per assembly
        :return: Dataframes for both metrics
        """
        if not nx_values:
            nx_values, aun_values = self.calculate_metrics()
        if not aun_values:
            nx_values, aun_values = self.calculate_metrics()
        nx_frame = self._generate_nx_frame(nx_values)
        aun_frame = self._generate_aun_frame(aun_values)
        return nx_frame, aun_frame



    def metric_dataframes(self, nx_frame: pd.DataFrame = None, aun_frame: pd.DataFrame = None, out_name: str = None) -> None:
        """
        Write the pandas.DataFrames of Nx and auN values to csv files

        :param nx_frame: Dataframe of Nx values for all assemblies
        :param aun_frame: Dataframe of auN values for all assemblies
        :param out_name: Base name for output files
        :return: None
        """
        if nx_frame is None:
            nx_frame, aun_frame = self.generate_dataframes()
        if aun_frame is None:
            nx_frame, aun_frame = self.generate_dataframes()

        print_frame(df=nx_frame, out_file=f'{out_name}.nx.csv')
        print_frame(df=aun_frame, out_file=f'{out_name}.aun.csv')


    def _calculate_Nx(self) -> NX_DICT:
        """
        Calculate Nx values for each assembly

        :return: Dict of Nx values per assembly
        """
        nx_values = dict()
        for asm in self.assemblies:
            nx_vector = asm.calculate_Nx()
            nx_values[asm.name] = nx_vector
        return nx_values


    def _calculate_auN(self) -> AUN_DICT:
        """
        Calculate auN values for each assembly

        :return: Dict of auN values per assembly
        """
        aun_values = dict()
        for asm in self.assemblies:
            aun = asm.calculate_auN()
            aun_values[asm.name] = aun
        return aun_values


    def _generate_nx_frame(self, nx_values: NX_DICT) -> pd.DataFrame:
        """
        Generate pandas frame from values of all assemblies

        :param nx_values: Dict of Nx values per assembly
        :return: Dataframe of Nx values for all assemblies
        """
        ndf = {'Nx': [], 'val': [], 'name': [], 'reference': []}
        for name, vector in nx_values.items():
            ndf['Nx'].extend(np.arange(1, 101))
            ndf['val'].extend(vector)
            ndf['name'].extend([name] * 100)
            # mark the reference assembly for downstream analyses
            if name == self.ref_name:
                refv = np.ones(100, dtype="bool")
            else:
                refv = np.zeros(100, dtype="bool")
            ndf['reference'].extend(refv)
        # transform to dataframe
        nx_df = pd.DataFrame(ndf)
        return nx_df


    def _generate_aun_frame(self, aun_values: AUN_DICT) -> pd.DataFrame:
        """
        Generate pandas frame from values of all assemblies

        :param aun_values: Dict of auN values per assembly
        :return: Dataframe of auN values for all assemblies
        """
        ndf = {'auN': [], 'name': [], 'reference': []}
        for name, value in aun_values.items():
            ndf['auN'].append(value)
            ndf['name'].append(name)
            # mark the reference assembly for downstream analyses
            refi = True if name == self.ref_name else False
            ndf['reference'].append(refi)
        # transform to dataframe
        aun_df = pd.DataFrame(ndf)
        return aun_df


    def _initialise_assemblies(self, assembly_paths: List[pathlib.Path]) -> List['Assembly']:
        """
        Load an Assembly object for each input file

        :param assembly_paths: List of input Paths
        :return: List of initialised Assembly objects
        """
        assemblies = []
        for ap in assembly_paths:
            assemblies.append(Assembly(path=ap, gsize=self.genome_size))
        return assemblies




class Assembly:

    def __init__(self, path: Union[pathlib.Path, str], gsize: int = None):
        self.path = pathlib.Path(path)
        self.name = self.path.name
        self.gsize = 0 if not gsize else gsize
        self.gzipped = is_gzipped(self.path)
        self.lengths = self._get_sequence_lengths()


    def _get_sequence_lengths(self) -> np.ndarray:
        """
        Load the sequence lengths from a file

        :return: array of sequence lengths
        """
        assert os.path.getsize(self.path) != 0
        seq_lengths = []
        open_func = gzip.open if self.gzipped else open
        read_func = read_fa_gz if self.gzipped else read_fa
        with open_func(self.path, 'r') as fa:
            for header, seq in read_func(fa):
                seq_lengths.append(len(seq))
        # get lengths of all sequences
        lengths = np.array(seq_lengths)
        return lengths


    def calculate_Nx(self) -> List[int]:
        """
        Calculate Nx values for this assembly

        :return: Nx values of the assembly
        """
        if not self.gsize:
            gsize = np.sum(self.lengths)
        else:
            gsize = self.gsize
        nx = []
        # sort sequence length and calc cumulative sum
        seq_lengths_sorted = np.sort(self.lengths)[::-1]
        seq_lengths_sorted_cuml = np.cumsum(seq_lengths_sorted)
        asm_perc = np.arange(0.01, 1.01, 0.01)
        # multiply either by total contig length
        # or by reference length/genome estimate
        asm_p = asm_perc * gsize
        for i in range(len(asm_p)):
            j = 0
            try:
                while seq_lengths_sorted_cuml[j] < asm_p[i]:
                    j += 1
                nx.append(seq_lengths_sorted[j])
            except IndexError:
                nx.append(0)
        return nx


    def calculate_auN(self) -> float:
        """
        Calculate area under Nx for the assembly

        :return: auN value
        """
        if not self.gsize:
            aun = np.sum(np.power(self.lengths, 2)) / np.sum(self.lengths)
        else:
            aun = np.sum(self.lengths * (self.lengths / self.gsize))
        return aun


class Reference(Assembly):

    def genome_size(self) -> np.ndarray:
        """
        Sum of reference sequence lengths used for NGx and auNG values

        :return: Total sequence length of reference file
        """
        return np.sum(self.lengths)




def print_frame(df: pd.DataFrame, out_file: str = None) -> None:
    """
    Print a dataframe in csv format either to stdout or to file

    :param df: Pandas data frame
    :param out_file: name of an output file
    :return: None
    """
    if not out_file:
        print(df.to_csv())
    else:
        df.to_csv(out_file)



def read_fa(fh: TextIO) -> Tuple[str, str]:
    """
    Yield headers and sequences of a fasta file

    :param fh: File handle of an open file connection
    :return: Tuple of fasta header and sequence
    """
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        headerStr = header.__next__().strip().split(' ')[0]
        seq = "".join(s.strip() for s in faiter.__next__())
        yield headerStr, seq


def read_fa_gz(fh: TextIO) -> Tuple[str, str]:
    """
    Yield headers and sequences of a gzipped fasta file

    :param fh: File handle of an open file connection
    :return: Tuple of fasta header and sequence
    """
    faiter = (x[1] for x in groupby(fh, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        headerStr = str(header.__next__(), 'utf-8').strip().replace('>', '').split()[0]
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield headerStr, seq


def convert_paths(paths: List[Union[pathlib.Path, str]]) -> List[pathlib.Path]:
    """
    Make sure that all paths in a list are pathlib.Path objects

    :param paths: List of paths
    :return: List of pathlib.Path objects
    """
    conv_paths = []
    for p in paths:
        conv_paths.append(pathlib.Path(p))
    return conv_paths



def is_gzipped(f: Union[pathlib.Path, str]) -> bool:
    """
    Check if a file is gzipped

    :param f: File path
    :return: bool
    """
    isgz = True
    with gzip.open(f, 'r') as fh:
        try:
            fh.read(1)
        except gzip.BadGzipFile:
            isgz = False
    return isgz










