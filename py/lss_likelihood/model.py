import lss_likelihood._lss_likelihood as c
import numpy as np

from scipy.interpolate import interp1d

def exp_moment(a):
    return c._model_exp_moment(a)

class Kaiser:
    def __init__(self, k, P, k_eval):
        P_eval = interp1d(k, P, kind='cubic')(k_eval)
        self._model = c._model_kaiser_alloc(k_eval, P_eval)
        self._k = k_eval
        self.P = P_eval 

    def evaluate(self, b, f, sigma):
        P0 = np.empty_like(self._k)
        P2 = np.empty_like(self._k)
        P4 = np.empty_like(self._k)

        c._model_kaiser_evaluate(self._model, b, f, sigma,
                                 P0, P2, P4)

        d = {}
        d['b'] = b
        d['f'] = f
        d['sigma'] = sigma
        d['k'] = self._k
        d['P0'] = P0
        d['P2'] = P2
        d['P4'] = P4

        return d

