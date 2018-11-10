import lss_likelihood._lss_likelihood as c

from scipy.interpolate import interp1d

class Kaiser:
    def __init__(self, k, P, k_eval):
        P_eval = interp1d(k, P, kind='cubic')(k_eval)
        self._model = c._model_kaiser_alloc(k_eval, P_eval)

    def evaluate(self, b, f, sigma):
        pass
