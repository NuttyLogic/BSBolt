import random
import numpy as np
import pandas as pd
from tqdm import tqdm
from BSB.BSB_Impute.Imputation.GenomeImputation import GenomeImputation
from BSB.BSB_Utils.MatrixIterator import OpenMatrix


class ImputeMissingValues:
    """

            :param methylation_data:
            :param batch_size:
            :param imputation_window_size:
            :param k:
            :param threads:
            :param verbose:
            :param sep:
            :param output_path:
    """

    def __init__(self, input_matrix_file=None, batch_size=None, imputation_window_size=3000000, k=5,
                 threads=4, verbose=False, sep='\t', output_path=None, randomize_batch=False):
        self.input_matrix = input_matrix_file
        self.batch_size = batch_size
        self.output_path = output_path
        if self.output_path:
            self.output_matrix = open(self.output_path, 'w')
        self.tqdm_disable = False if verbose else True
        self.sep = sep
        self.imputation_kwargs = dict(imputation_window_size=imputation_window_size,
                                      k=k,
                                      threads=threads,
                                      verbose=verbose)

    def impute_values(self):
        if self.complete_matrix and not self.batch_size:
            knn_imputation = self.launch_genome_imputation(self.methylation_dataframe)
            self.methylation_dataframe = pd.DataFrame(data=knn_imputation.genomic_array,
                                                      index=list(self.methylation_dataframe.index),
                                                      columns=list(self.methylation_dataframe))
        else:
            self.batch_imputation()

    def import_matrix(self, bsb_matrix=None):
        methylation_matrix = {}
        sample_number = None
        matrix_samples = None
        with open(bsb_matrix, 'r') as methylation_matrix:
            for line in methylation_matrix:
                if not matrix_samples:
                    matrix_samples = line.replace('\n', '').split('\t')

    def randomly_sample_list(self, batch_size):
        try:
            batch_label = [0] + random.sample(sample_list, batch_size)
        except ValueError:
            batch_label = [0] + sample_list
        return batch_label

    def update_output_matrix(self, batch_knn_imputation):
        if not self.output_header:
            self.output_matrix.write('\t'.join(['Samples'] + batch_knn_imputation.row_labels) + '\n')
            self.output_header = True
        for sample, values in zip(batch_knn_imputation.sample_labels, batch_knn_imputation.genomic_array.T):
            self.output_matrix.write('\t'.join([sample] + [str(value) for value in values]) + '\n')


    def launch_genome_imputation(self, methylation_dataframe):
        imputation_kwargs = dict(self.imputation_kwargs)
        assert isinstance(methylation_dataframe.values, (np.ndarray, np.generic))
        imputation_kwargs.update(dict(genomic_array=methylation_dataframe.values,
                                      sample_labels=list(methylation_dataframe),
                                      row_labels=list(methylation_dataframe.index)))
        knn_imputation: GenomeImputation = GenomeImputation(**imputation_kwargs)
        knn_imputation.impute_windows()
        return knn_imputation


