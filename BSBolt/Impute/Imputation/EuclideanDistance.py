#! /usr/bin/env python3

from typing import List
import numpy as np


def get_euclidean(array: np.ndarray = None, global_neighbors: np.ndarray = None) -> List[List[float]]:
    """ Pairwise euclidean distance
    ----------------------------------------
    input: numpy array with samples as rows
    output: pairwise euclidean distance array"""
    # split numpy array by columns
    assert isinstance(array, (np.ndarray, np.generic))
    samples = np.split(array, 1)
    # output to ensure correct matrix orientation
    # initialize empty list to hold pairwise distance values
    pairwise_distance = []
    # iterate through sample combinations
    for count, vector in enumerate(samples[0]):
        # add list for all sample pairwise comparison
        pairwise_distance.append([])
        for comparison_count, comparison_vector in enumerate(samples[0]):
            # only perform calculation for upper half of distance matrix to save time
            if comparison_count >= count:
                # stack sample values
                comparison_array = np.stack((vector, comparison_vector), axis=0)
                # drop any comparison value with NA for any sample, note this is only done pairwise between
                # samples to get best distance comparison estimate
                comparison_array = np.ma.compress_cols(np.ma.masked_invalid(comparison_array))
                try:
                    if len(comparison_array[0]) and len(comparison_array[1]):
                        pairwise_distance[count].append(np.linalg.norm(comparison_array[0] - comparison_array[1]) /
                                                        len(comparison_array[1]))
                    else:
                        raise IndexError
                # if all values are missing the missing neighbor is global
                except IndexError:
                    if isinstance(global_neighbors, (np.ndarray, np.generic)):
                        pairwise_distance[count].append(global_neighbors[count][comparison_count])
                    else:
                        pairwise_distance[count].append(np.nan)
            else:
                pairwise_distance[count].append(0)
    # add upper and lower matrices
    pairwise_upper = np.asarray(pairwise_distance)
    pairwise_lower = np.transpose(pairwise_upper)
    # return list to ensure compatibility with downstream tools
    pairwise_distance = (pairwise_upper + pairwise_lower).tolist()
    return pairwise_distance
