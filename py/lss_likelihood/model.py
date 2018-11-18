import lss_likelihood._lss_likelihood as c
import numpy as np

#from scipy.interpolate import interp1d

def exp_moment(a):
    return c._model_exp_moment(a)

class Model:
    def __init__(self, _model):
        self._model = _model
        self.nbin = c._model_get_nbin(self._model)

    def nmodes(self):
        nmodes = np.empty(self.nbin)
        c._model_nmodes(self._model, nmodes)

        return nmodes

    def evaluate(self, b, f, sigma):
        P0 = np.empty(self.nbin)
        P2 = np.empty_like(P0)
        P4 = np.empty_like(P0)

        c._model_evaluate(self._model, b, f, sigma, P0, P2, P4)

        d = {}
        d['b'] = b
        d['f'] = f
        d['sigma'] = sigma
        d['P0'] = P0
        d['P2'] = P2
        d['P4'] = P4

        return d


def Kaiser(k_min, dk, nbin, boxsize, k, P):
    """
    Create Kaiser model object
    
    Args:
      k_min (float)  : minimum edge of k bin [h/Mpc]
      dk (float)     : bin width [h/Mpc]
      nbin (int)     : number of k bins
      boxsize (float): boxsize of FFT grid
      k (array)      : k of real-space power spectrum P(k)
      P (array)      : Real-space power spectrum P(k)

    Returns:
      Model object
    """
    _model = c._model_kaiser_alloc(k_min, dk, nbin, boxsize, k, P)
    
    return Model(_model)

def Scoccimarro(k_min, dk, nbin, boxsize, k, Pdd, Pdt, Ptt):
    """
    Create Scoccimarro model object
    
    Args:
      k_min (float)  : minimum edge of k bin [h/Mpc]
      dk (float)     : bin width [h/Mpc]
      nbin (int)     : number of k bins
      boxsize (float): boxsize of FFT grid
      k (array)      : k of real-space power spectra
      Pdd (array)    : Real-space matter power spectrum P_delta_delta(k)
      Pdt (array)    : Real-space matter theta power P_delta_theta(k)
      Ptt (array)    : Real-space matter theta power P_theta_theta(k)

    Returns:
      Model object
    """
    _model = c._model_scoccimarro_alloc(k_min, dk, nbin, boxsize, k,
                                        Pdd, Pdt, Ptt)
    
    return Model(_model)


def Scoccimarro(k_min, dk, nbin, boxsize, k, Pdd, Pdt, Ptt):
    """
    Create Scoccimarro model object
    
    Args:
      k_min (float)  : minimum edge of k bin [h/Mpc]
      dk (float)     : bin width [h/Mpc]
      nbin (int)     : number of k bins
      boxsize (float): boxsize of FFT grid
      k (array)      : k of real-space power spectra
      Pdd (array)    : Real-space matter power spectrum P_delta_delta(k)
      Pdt (array)    : Real-space matter theta power P_delta_theta(k)
      Ptt (array)    : Real-space matter theta power P_theta_theta(k)

    Returns:
      Model object
    """
    _model = c._model_scoccimarro_alloc(k_min, dk, nbin, boxsize, k,
                                        Pdd, Pdt, Ptt)
    
    return Model(_model)


def Taruya(k_min, dk, nbin, boxsize, k, Pdd, Pdt, Ptt, AB):
    """
    Create Scoccimarro model object
    
    Args:
      k_min (float)  : minimum edge of k bin [h/Mpc]
      dk (float)     : bin width [h/Mpc]
      nbin (int)     : number of k bins
      boxsize (float): boxsize of FFT grid
      k (array)      : k of real-space power spectra
      Pdd (array)    : Real-space matter power spectrum P_delta_delta(k)
      Pdt (array)    : Real-space matter theta power P_delta_theta(k)
      Ptt (array)    : Real-space matter theta power P_theta_theta(k)
      AB (array)     : 15-column array with k, A11, A12, A22, A23, A33,
                       B111, B211, B112 + B121, B212 + B221, B312 + B321,
                       B122, B222, B322, B422

    Returns:
      Model object
    """

    if AB.ndim != 2:
        raise TypeError('AB must be a 2-dimensional array: shape = {}'.format(
                        AB.shape))
    elif AB.shape[1] != 15:
        raise TypeError('AB must be a 15-column array; shape = {}'.format(
                        AB.shape))

    _model = c._model_taruya_alloc(k_min, dk, nbin, boxsize, k,
                                   Pdd, Pdt, Ptt, AB)
    
    return Model(_model)

