import random
import numpy as np
import pandas as pd
from tqdm import tqdm
from BSB.BSB_Impute.Imputation.GenomeImputation import GenomeImputation


class ImputeMissingValues:

    def __init__(self, methylation_data=None, batch_size=None, imputation_window_size=3000000, k=5,
                 threads=4, verbose=False, sep='\t', output_path=None):
        self.complete_matrix = False
        self.batch_size = batch_size
        self.output_path = output_path
        if self.output_path:
            self.output_matrix = open(self.output_path, 'w')
            self.output_header = False
        if isinstance(methylation_data, pd.DataFrame):
            self.methylation_dataframe = methylation_data
            self.complete_matrix = True
            self.sample_labels = list(self.methylation_dataframe)
        else:
            if not self.batch_size:
                self.methylation_dataframe = pd.read_csv(methylation_data, index_col=0, sep=sep)
                self.sample_labels = list(self.methylation_dataframe)
                self.complete_matrix = True
            else:
                self.methylation_dataframe = methylation_data
                self.sample_labels = list(pd.read_csv(methylation_data, index_col=0, nrows=1, sep=sep))
        self.verbose = verbose
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

    def batch_imputation(self):
        batches = self.get_batch_split
        sample_labels = [sample_index + 1 for sample_index in range(len(self.sample_labels))]
        batch_order = []
        for batch in tqdm(batches, desc='Processing Batches', total=len(batches),
                          disable=True if not self.verbose else False):
            batch_labels = self.randomly_sample_list(sample_labels, batch)
            batch_order.extend(batch_labels)
            self.remove_samples(sample_labels, batch_labels)
            if self.complete_matrix:
                batch_df: pd.DataFrame = self.methylation_dataframe[batch_labels]
                batch_knn_imputation = self.launch_genome_imputation(batch_df)
                self.update_methylation_dataframe(batch_knn_imputation)
            else:
                batch_df = pd.read_csv(self.methylation_dataframe, index_col=0, sep=self.sep, usecols=batch_labels)
                batch_knn_imputation = self.launch_genome_imputation(batch_df)
                self.update_output_matrix(batch_knn_imputation)
        if self.output_path:
            self.output_matrix.close()

    @staticmethod
    def randomly_sample_list(sample_list, batch_size):
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

    def update_methylation_dataframe(self, batch_imputation: GenomeImputation):
        batch_dataframe = pd.DataFrame(data=batch_imputation.genomic_array,
                                       index=batch_imputation.row_labels,
                                       columns=batch_imputation.sample_labels)
        for sample in list(batch_dataframe):
            self.methylation_dataframe[sample] = batch_dataframe[sample].values

    def launch_genome_imputation(self, methylation_dataframe):
        imputation_kwargs = dict(self.imputation_kwargs)
        assert isinstance(methylation_dataframe.values, (np.ndarray, np.generic))
        imputation_kwargs.update(dict(genomic_array=methylation_dataframe.values,
                                      sample_labels=list(methylation_dataframe),
                                      row_labels=list(methylation_dataframe.index)))
        knn_imputation: GenomeImputation = GenomeImputation(**imputation_kwargs)
        knn_imputation.impute_windows()
        return knn_imputation

    @staticmethod
    def remove_samples(sample_list, sample_removal_list):
        for sample in sample_removal_list[1:]:
            sample_list.remove(sample)

    @property
    def get_batch_split(self):
        # set number of batches according to batch size
        batches = [self.batch_size for _ in range(int(len(self.sample_labels) / self.batch_size))]
        # add remainder of last batch
        batch_remainder = len(self.sample_labels) % self.batch_size
        if batch_remainder < self.batch_size / 2:
            batch_addition = int(batch_remainder / len(batches))
            batch_addition_remainder = batch_remainder % len(batches) + batch_addition
            batches[:-1] = [batch + batch_addition for batch in batches]
            batches[-1] += batch_addition_remainder
        else:
            batches.append(batch_remainder)
        return batches
