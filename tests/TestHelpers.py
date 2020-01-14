import os
import numpy as np


test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'


def z_test_of_proportion(a_yes, a_no, b_yes, b_no):
    a_total = a_yes + a_no
    b_total = b_yes + b_no
    a_prop = a_yes / a_total
    b_prop = b_yes / b_total
    p_hat = (a_yes + b_yes) / (a_total + b_total)
    try:
        return (a_prop - b_prop) / np.sqrt(p_hat * (1 - p_hat) * (1 / a_total + 1 / b_total))
    except RuntimeWarning:
        return 0
