import numpy as np
from bsbolt.Variant.BayesVariant import generate_log_likelihood_matrix, predict_bayes_genotype

class CallVariant:

    def __init__(self, error_rates: np.ndarray = np.array([0.001]), conversion_rate: float = 0.995,
                 lower_methylation_rate: float = 0.01, upper_methylation_rate: float = 0.9):
        self.lower_llm = generate_log_likelihood_matrix(error_rates=error_rates, 
                                                        methylation_rate=lower_methylation_rate, 
                                                        conversion_rate=conversion_rate)
        self.upper_llm = generate_log_likelihood_matrix(error_rates=error_rates, 
                                                        methylation_rate=upper_methylation_rate, 
                                                        conversion_rate=conversion_rate)
 
    def call_variant(self, base_calls: np.ndarray, priors: np.ndarray = np.ones(10), cg_site=False):
        lower_call = predict_bayes_genotype(base_calls, self.lower_llm, priors)
        if cg_site:
            pass
            # calls = lower_call, predict_bayes_genotype(base_calls, self.upper_llm, priors)
            # return calls[0] if calls[0][0] > calls[1][0] else calls[1]
        return lower_call
