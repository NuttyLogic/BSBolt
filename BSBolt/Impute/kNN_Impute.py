import gzip
import io
import random
from typing import List, Tuple
import numpy as np
from BSBolt.Impute.Imputation.GenomeImputation import GenomeImputation
from BSBolt.Impute.Impute_Utils.ImputationFunctions import get_bsb_matrix


class ImputeMissingValues:
    """
    Launch and knn imputation task. This wrapper imports data for imputation, split data for batch imputation,
    and combines data after imputation. Data is held in memory for access during imputation. If closest neighbors null,
    null imputed value returned.

    Params:

    * *input_matrix_file (str)*: Path to BSBolt matrix
    * *batch_size (int)*: Batch size for batch imputation
    * *imputation_window_size (int)*: Size (bp) for imputation window, [3,000,000]
    * *k (int)*: Nearest neighbors used for imputation, [5]
    * *threads (int)*: Number of threads available for imputation, [1]
    * *verbose (bool)*: Verbose imputation, [False]
    * *sep (str)*: separator character used in methylation matrix, [\t]
    * *output_path (str)*: output path
    * *randomize_batch (bool)*: randomize batch, [False]
    * *meth_matrix (np.ndarray)*: imputed methylation matrix
    * *meth_site_order (list)*: ordered methylation sites, sorted by contig then position
    * *sample_ids (list)*: sample names
    """

    def __init__(self, input_matrix_file: str = None, batch_size: int = None, imputation_window_size: int = 3000000,
                 k: int = 5, threads: int = 4, verbose: bool = False, sep: str = '\t', output_path: str = None,
                 randomize_batch: bool = False, meth_matrix: np.ndarray = None, meth_site_order: list = None,
                 sample_ids: list = None):
        self.input_matrix_file = input_matrix_file
        self.meth_matrix = meth_matrix
        self.meth_site_order = meth_site_order
        self.sample_ids = sample_ids
        self.batch_commands = [0, False]
        if batch_size:
            self.batch_commands = [batch_size, randomize_batch]
        self.output_matrix = None
        if output_path:
            self.output_matrix = self.get_output_matrix(output_path)
        self.tqdm_disable = False if verbose else True
        self.sep = sep
        self.imputation_kwargs = dict(imputation_window_size=imputation_window_size,
                                      k=k,
                                      threads=threads,
                                      verbose=verbose)

    def impute_values(self):
        """
        Launch kNN imputation for each batch and set values in original matrix.
        """
        imputation_order = [count for count in range(len(self.sample_ids[1]))]
        if self.batch_commands[1]:
            random.shuffle(imputation_order)
        if not self.batch_commands[0]:
            self.batch_commands[0] = len(self.sample_ids[1])
        imputation_batches = self.process_batch(imputation_order)
        for batch in imputation_batches:
            batch_array, sample_labels = self.get_batch_data(batch)
            batch_impute = self.launch_genome_imputation(batch_array, sample_labels)
            for count, sample in enumerate(batch):
                self.meth_matrix[:, sample] = batch_impute.genomic_array[:, count]

    def process_batch(self, imputation_order: List[int]) -> List[List[int]]:
        """Generate sample batches"""
        batches = []
        if not self.batch_commands[0]:
            self.batch_commands[0] = len(self.sample_ids)
        batch_number = int(len(imputation_order) / self.batch_commands[0])
        batch_remainder = len(imputation_order) % self.batch_commands[0]
        if float(batch_remainder) / float(self.batch_commands[0]) <= .6 and batch_remainder != 0:
            batch_addition = int(batch_remainder / batch_number)
            self.batch_commands[0] += batch_addition + 1
            print(f'Adjusting batch size, new batch size = {self.batch_commands[0]}')
        start, end = 0, self.batch_commands[0]
        while True:
            batch_order = imputation_order[start: end]
            if not batch_order:
                break
            batches.append(batch_order)
            start += self.batch_commands[0]
            end += self.batch_commands[0]
        return batches

    def get_batch_data(self, batch: List[int]) -> Tuple[np.ndarray, List[str]]:
        """Return methylation value for batch imputation

        Returns:

        * *batch_array (np.ndarry)*: array of methylation values
        * *sample_labels (list)*: list of samples in batch
            """
        batch_array = self.meth_matrix[:, batch]
        sample_labels = [self.sample_ids[1][sample] for sample in batch]
        return batch_array, sample_labels

    def import_matrix(self):
        self.meth_matrix, self.meth_site_order, self.sample_ids = get_bsb_matrix(self.input_matrix_file)

    @staticmethod
    def get_output_matrix(output_path: str):
        """Get output object"""
        if output_path.endswith('.gz'):
            out = io.BufferedWriter(gzip.open(output_path, 'wb'))
        else:
            out = open(output_path, 'w')
        return out

    def output_imputed_matrix(self):
        """Write imputed values"""
        self.output_matrix.write('\t'.join([self.sample_ids[0]] + self.sample_ids[1]) + '\n')
        for site_label, values in zip(self.meth_site_order, self.meth_matrix):
            str_values = '\t'.join([str(value) for value in values])
            self.output_matrix.write(f'{site_label}\t{str_values}\n')
        self.output_matrix.close()

    def launch_genome_imputation(self, meth_array: np.ndarray, sample_labels: List) -> object:
        imputation_kwargs = dict(self.imputation_kwargs)
        imputation_kwargs.update(dict(genomic_array=meth_array,
                                      sample_labels=sample_labels,
                                      row_labels=self.meth_site_order))
        knn_imputation: GenomeImputation = GenomeImputation(**imputation_kwargs)
        knn_imputation.impute_windows()
        return knn_imputation
