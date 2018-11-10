import lss_likelihood._lss_likelihood as c
import numpy as np

class Likelihood:
    def __init__(self, data, model, cov):
        cov_inv = np.linalg.inv(cov)
        self._likelihood = c._likelihood_alloc(data._data, model._model, cov_inv)

    def chi2(self, params):
        return c._likelihood_chi2(self._likelihood, params)
