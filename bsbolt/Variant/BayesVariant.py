from typing import Tuple
import numpy as np

def generate_log_likelihood_matrix(error_rates: np.ndarray = np.array([0.005]), methylation_rate: float = 0.01, 
                                   conversion_rate: float = 0.99) -> np.ndarray:
    """
    Generates log likelihood matrix for Bayesian bisulfite sequencing variant caller. The likelihoods are dependent 
    on the expected methylation level for a particular site. The likelihood of observing one of 4 possible bases is a 
    function of the genotype ($n=10$ for a diploid organism) and the strand ($n=2$). The probability of observing 
    a cytosine / thymine for sense reads and guanine / adenine for anti-sense reads is dependent on the methylation 
    rate, the conversion rate, and the base call error rate. 

    Params:

    * *error_rates (np.ndarray)*: base call error probabilities
    * *methylation_rate (float)*: methylation rate for matrix generation
    * *conversion_rate (float)*: bisulfite conversion efficiency  

    Returns:

    * *llm (np.ndarray)*: $8\\times10$matrix of read base log likelihoods for sense and anti-sense 
                            bisulfite converted reads. *Rows (+A, +T, +C, +G, -A, -T, -C, -G), 
                            Columns (AA, AT, AC, AG, TT, TC, TG, CC, CG, GG)  
    """
    # use median error rate for likelihood construction
    er = np.median(error_rates) / 3
    cr = conversion_rate
    mr = methylation_rate

    # set convertable base probability C/G
    cbp = (1 - 3 * er) * (mr + (1 - cr) * (1 - mr))

    # set convertable base target probability T/A
    ctp = er + cr * (1 - mr)

    # set ref base probability
    rbp = 1 - er * 3

    # set values for genotype combinations
    cbp_er = (cbp + er) / 2
    ctp_er = (ctp + er) / 2
    ctp_rbp = (ctp + rbp) / 2
    er_rbp = (er + rbp) / 2

    llm = np.log(np.array([[rbp, er_rbp, er_rbp, er_rbp, er, er, er, er, er, er],
                            [er, er_rbp,ctp_er, er, rbp,ctp_rbp,er_rbp,ctp, ctp_er, er],
                            [er, er, cbp_er, er, er, cbp_er, er, cbp, cbp_er, er],
                            [er, er, er, er_rbp, er, er, er_rbp, er, er_rbp, rbp],
                            [rbp, er_rbp, er_rbp, ctp_rbp, er, er, ctp_er, er, ctp_er, ctp],
                            [er, er_rbp, er, er, rbp, er_rbp, er_rbp, er, er, er],
                            [er, er, er_rbp, er, er, er_rbp, er, rbp, er_rbp, er],
                            [er, er, er, cbp_er, er, er, cbp_er, er, cbp_er, cbp]]))
    return(llm)


def predict_bayes_genotype(base_calls: np.ndarray, log_likelihood_matrix: np.ndarray, 
                           priors: np.ndarray = np.ones((10))) -> Tuple[float, float, float, str]:
    """
    Provides genotype calls given an array of base counts at a site of interest and a generated log likelihood matrix.
    Weighted priors may also be provided, defaults to 1.0 for all genotypes. The genotype with the maximum probability
    is returned.

    Params:

    * *base_calls (np.ndarray)*: counts of base calls in order (+A, +T, +C, +G, -A, -T, -C, -G)
    * *log_likelihood_matrix (np.ndarray)*: $8\\times10$matrix of read base log likelihoods for sense and anti-sense
                            bisulfite converted reads. *Rows (+A, +T, +C, +G, -A, -T, -C, -G),
                            Columns (AA, AT, AC, AG, TT, TC, TG, CC, CG, GG)
    * *priors (np.ndarray)*: expected diploid genotype frequency, defaults to 1.0 for all genotypes

    Returns:

    * *genotype call: Tuple[float, float, float, str]*: genotype probability, 1 - genotype probability,
        genotype score (log_{10}(Max(genotype probability)/ Max -1(genotype probability )), genotype
    """
    # order of likelihood matrix
    genotype_order = ('AA', 'AT', 'AC', 'AG', 'TT', 'TC', 'TG', 'CC', 'CG', 'GG')
    # probability calculation
    genotype_likelihoods = np.exp(np.dot(base_calls, log_likelihood_matrix)) * priors
    probs = genotype_likelihoods / sum(genotype_likelihoods)
    # sort by probability
    genotype_probs = sorted(zip(probs, genotype_order), reverse=True)
    prob, genotype = genotype_probs[0]
    # score probability by comparing max to second best call
    score = np.log10(genotype_probs[0][0] / genotype_probs[1][0])
    return prob, 1 - prob, score, genotype
