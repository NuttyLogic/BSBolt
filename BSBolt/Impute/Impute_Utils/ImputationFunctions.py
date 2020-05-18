import os
from typing import List, Optional, Tuple
import numpy as np
from BSBolt.Utils.MatrixIterator import OpenMatrix


def get_bsb_matrix(bsb_matrix: str) -> Tuple[np.ndarray, List[str], List[str]]:
    """Import bsb matrix file
    Arguments:
        bsb_matrix (str): path to matrix file generated using BSBolt
    Returns:
        site_values (np.array): array of methylation matrix values
        site_order (list): list of methylation site identifiers
        sample_ids ([site_labels, [sample_labels]]: site label, list of sample labels
        """
    assert os.path.exists(bsb_matrix), f'{bsb_matrix} does not exist, please change path'
    site_values, site_order, sample_ids = [], [], None
    for line_label, line_values in OpenMatrix(bsb_matrix):
        if not sample_ids:
            sample_ids = line_label, line_values
        else:
            site_order.append(line_label)
            site_values.append(line_values)
    return np.asarray(site_values), site_order, sample_ids
