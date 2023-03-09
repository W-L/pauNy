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
                        help='input fasta/fastq file(s) or director(ies) of files. Can be multiple (space-separated).')
    parser.add_argument('-o', '--out', dest='out', type=str, default="pauNy", help='base name for output files')
    parser.add_argument('-f', '--format', dest='format', type=str, default='pdf', help='output format for plots')
    reference_group = parser.add_mutually_exclusive_group(required=False)
    reference_group.add_argument('-r', '--ref', dest='ref', type=str, default=None, help='path to reference sequence file')
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
        patterns_fa = ["*.fa.gz", "*.fasta.gz", "*.fasta.gzip", "*.fa.gzip", "*.fasta", "*.fa"]
        patterns_fq = ["*.fq.gz", "*.fastq.gz", "*.fastq.gzip", "*.fq.gzip", "*.fastq", "*.fq"]
        patterns = patterns_fa + patterns_fq
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
        genome_size=args.genomesize,
        out_name=args.out,
    )
    # calculate Nx and auN values for all input files
    # this also generates pandas frames
    asm.calculate_metrics()
    # visualise results
    asm.plot(format=args.format)

    # produce pngs for the example in repository
    # asm.plot(format="png")


if __name__ == "__main__":
    main()


