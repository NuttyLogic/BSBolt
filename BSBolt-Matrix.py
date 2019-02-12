#! /usr/bin/env python3

import argparse
from BSB_Matrix.MatrixAggregator import AggregateMatrix
from BSB_Utils.BSB_UtilityFunctions import check_python_version


def get_sample_info(file_path):
    file_list = []
    with open(file_path, 'r') as file:
        for line in file:
            file_list.append(line.replace('\n', ''))
    return file_list


parser = argparse.ArgumentParser(description='BSBot Matrix Aggregation Module')

parser.add_argument('-F', type=lambda file: [file_path for file_path in file.split(',')], required=True,
                    help='Comma separated list of CGmap file paths, or path to text file with list of line separated '
                         'CGmap file paths')
parser.add_argument('-S', type=lambda sample_labels: [sample for sample in sample_labels.split(',')], default=None,
                    help='Comma separated list of samples labels. If sample labels are not provided sample labels '
                         'are extracted from CGmap file paths. Can also pass path to txt for line separated sample '
                         'labels.')
parser.add_argument('-min-coverage', type=int, default=10, help='Minimum site read depth coverage for a '
                                                                'site to be included in the aggregate matrix')
parser.add_argument('-min-sample', type=float, default=0.80, help='Proportion of samples that must have a valid site '
                                                                  '(above minimum coverage threshold), for a site to be'
                                                                  'included in the aggregate matrix.')
parser.add_argument('-O', type=str, default=None, required=True, help='Aggregate matrix output path')
parser.add_argument('-CG', action="store_true", default=False, help='Only output CG sites')
parser.add_argument('-verbose', action="store_true", default=False, help='Verbose aggregation')


if __name__ == "__main__":
    check_python_version()
    arguments = parser.parse_args()
    if len(arguments.F) == 1:
        arguments.F = get_sample_info(arguments.F[0])
    if arguments.S:
        if len(arguments.S) == 1:
            arguments.S = get_sample_info(arguments.S[0])
    matrix = AggregateMatrix(file_list=arguments.F,
                             sample_list=arguments.S,
                             min_site_coverage=arguments.min_coverage,
                             site_proportion_threshold=arguments.min_sample,
                             output_path=arguments.O,
                             cg_only=arguments.CG,
                             verbose=arguments.verbose)
