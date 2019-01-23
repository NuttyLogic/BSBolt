#! /usr/bin/env python3

import argparse
import os
from BSB_Simulate.SimulateMethylatedReads import SimulateMethylatedReads
from BSB_Utils.BSB_UtilityFunctions import check_python_version


parser = argparse.ArgumentParser(description='BSBolt Modules to Simulate Bisulfite Treated Illumina Reads')

parser.add_argument('-G', type=str, required=True,
                    help='Path for reference genome fasta file, fasta file should contain all contigs')
parser.add_argument('-A', type=str, required=True, help='Path to ART executable')
parser.add_argument('-O', type=str, required=True, help='Output prefix')
parser.add_argument('-PE', default=False, action='store_true', help='Simulate Paired End Reads, default Single End')
parser.add_argument('-RL', type=int, default=125, help='Simulated Read Lenghth')
parser.add_argument('-RD', type=int, default=20, help='Simulated Read Depth')
parser.add_argument('-U', default=False, action='store_true', help='Simulate Undirectional Reads, default=Directional')

if __name__ == "__main__":
    check_python_version()
    arguments = parser.parse_args()
    if os.path.isfile(arguments.A):
        read_simulation = SimulateMethylatedReads(reference_file=arguments.G, art_path=arguments.A,
                                                  output_path=arguments.O, paired_end=arguments.PE,
                                                  read_length=arguments.RL, read_depth=arguments.RD,
                                                  undirectional=arguments.U)
        read_simulation.run_simulation()
    else:
        print('ART Executable Path not Valid')
