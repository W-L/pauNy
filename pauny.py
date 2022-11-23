#!/usr/bin/env python3
import argparse
import os
import glob
from pathlib import Path

import pauNy



def setup_parser():
    parser = argparse.ArgumentParser(prog='pauNy',
                                     description='Nx curves and area under Nx in python')
    parser.add_argument('-i', '--input', dest='input', type=str, required=True, nargs='+',
                        help='input fasta file(s) or director(ies) of files. Can be multiple (space-separated).')
    parser.add_argument('-o', '--out', dest='out', type=str, default="pauNy", help='base name for output files')
    reference_group = parser.add_mutually_exclusive_group(required=False)
    reference_group.add_argument('-r', '--ref', dest='ref', type=str, default=None, help='path to reference fasta file')
    reference_group.add_argument('-g', '--genomesize', dest='genomesize', type=int, default=0, help='genome size or estimate')
    return parser


def parse_input_arg(input_arg):
    """
    Parse the input argument to grab file paths. Globs fasta files if dorectories are given.

    :param input_arg: List of paths to files or directories
    :return: List of pathlib.Path paths
    """
    # use first element to check if input are files or dirs
    first_element = input_arg[0]
    if os.path.isfile(first_element):
        paths = [Path(p) for p in input_arg]
    elif os.path.isdir(first_element):
        path_set = set()
        # glob for all possible fasta file extensions in given dir
        patterns = ["*.fa.gz", "*.fasta.gz", "*.fasta.gzip", "*.fa.gzip", "*.fasta", "*.fa"]
        for i in input_arg:
            for p in patterns:
                path_set.update(glob.glob(f'{i}/{p}'))
        paths = [Path(p) for p in path_set]
    else:
        paths = []
    return sorted(paths)


def main():
    # setup parser and grab arguments
    parser = setup_parser()
    args = parser.parse_args()
    # parse input file paths
    paths = parse_input_arg(args.input)
    # initialise a collection of assemblies from all input files
    asm = pauNy.AssemblyCollection(
        paths=paths,
        reference_path=args.ref,
        genome_size=args.genomesize
    )
    # calculate Nx and auN values for all input files
    nx_values, aun_values = asm.calculate_metrics()
    # generate pandas frames from the values - passing is optional
    nx_frame, aun_frame = asm.generate_dataframes(nx_values, aun_values)
    # write dataframes to files - passing is optional
    asm.metric_dataframes(
        nx_frame=nx_frame,
        aun_frame=aun_frame,
        out_name=args.out
    )
    # generate and save plots - passing data frames is optional
    asm.metric_plots(
        nx_frame=nx_frame,
        aun_frame=aun_frame,
        out_name=args.out
    )
    # produce pngs for the example in repository
    # asm.metric_plots(
    #     nx_frame=nx_frame,
    #     aun_frame=aun_frame,
    #     out_name=args.out,
    #     format="png"
    # )



if __name__ == "__main__":
    main()


