#! /usr/bin/env python3

import argparse

from BSB_Index.RRBSGenomeBuild import RRBSGenomeIndexBuild
from BSB_Index.WholeGenomeBuild import WholeGenomeIndexBuild
from BSB_Utils.BSB_UtilityFunctions import check_bowtie2_path

parser = argparse.ArgumentParser(description='BSBot Index Generation Module')

parser.add_argument('-G', type=str, required=True,
                    help='Path for reference genome fasta file, fasta file should contain all contigs')
parser.add_argument('-DB', type=str, required=True,
                    help='Path to index directory, will create directory if folder does not exist')
parser.add_argument('-BT2', type=str, default='bowtie2', help='Path to bowtie2 executable')
parser.add_argument('-BT2-p', type=int, default=2, help='Number of threads for Bowtie2 to use')
parser.add_argument('-rrbs', action="store_true", default=False, help='Generate Reduced Representative'
                                                                      ' Bisulfite Sequencing Index')
parser.add_argument('-rrbs-cut-format', default='C-CGG', help='Cut format to use for generation of RRBS database, '
                                                              'default= C-CGG (MSPI), input multiple enzymes as a '
                                                              'comma seperate string, C-CGG,C-CGG,...')
parser.add_argument('-rrbs-lower', type=int, default=40, help='Lower bound fragment size to consider RRBS index'
                                                              'generation, default = 40')
parser.add_argument('-rrbs-upper', type=int, default=500, help='Upper bound fragment size to consider RRBS index'
                                                               'generation, default = 500')

if __name__ == "__main__":
    arguments = parser.parse_args()
    check_bowtie2_path(bowtie2_path=arguments.BT2)
    if arguments.rrbs:
        print(f'Generating RRBS Database at {arguments.DB}: '
              f'lower bound {arguments.rrbs_lower}, upper bound {arguments.rrbs_upper}: '
              f'Cut Format {arguments.cut_format}')
        index = RRBSGenomeIndexBuild(reference_file=arguments.G,
                                     genome_database=arguments.DB,
                                     bowtie2_path=arguments.BT2,
                                     bowtie2_threads=arguments.BT2_p,
                                     cut_format=arguments.rrbs_cut_format,
                                     lower_bound=arguments.rrbs_lower,
                                     upper_bound=arguments.rrbs_upper)
        index.generate_rrbs_database()
    else:
        print(f'Generating WGBS Database at {arguments.DB}')
        index = WholeGenomeIndexBuild(reference_file=arguments.G,
                                      genome_database=arguments.DB,
                                      bowtie2_path=arguments.BT2,
                                      bowtie2_threads=arguments.BT2_p)
        index.generate_bsseeker_database()
