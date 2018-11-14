import lss_likelihood._lss_likelihood as c
import numpy as np

#from scipy.interpolate import interp1d

def exp_moment(a):
    return c._model_exp_moment(a)

class Kaiser:
    def __init__(self, k_min, k_max, dk, boxsize, k, P):
        self._model = c._model_kaiser_alloc(k_min, k_max, dk, boxsize,
                                            k, P)
        self.nbin = c._model_get_nbin(self._model)

    def nmodes(self):
        nmodes = np.empty(self.nbin)
        c._model_nmodes(self._model, nmodes)

        return nmodes

    def evaluate(self, b, f, sigma):
        P0 = np.empty(self.nbin)
        P2 = np.empty_like(P0)
        P4 = np.empty_like(P0)

        c._model_kaiser_evaluate(self._model, b, f, sigma,
                                 P0, P2, P4)

        d = {}
        d['b'] = b
        d['f'] = f
        d['sigma'] = sigma
        d['P0'] = P0
        d['P2'] = P2
        d['P4'] = P4

        return d

