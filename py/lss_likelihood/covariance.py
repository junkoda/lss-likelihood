import numpy as np
from scipy.special import legendre
from scipy.integrate import quad
from scipy.interpolate import interp1d

def integ():
    pass
        
def simple(k, P, k_eval, nmodes, b, f, nbar):
    """
    Compute simple 3x3 covariance matrix between multipole moments
    based on linear theory: Ps(k) = (b + f mu^2)^2 P(k).

    Args:
      k P (array): linear power spectrum
      k_eval (array): k of likelihood evaluation (k_data)
      nmodes (array): number of modes in the k bin
    """

    if len(k_eval) != len(nmodes):
        raise ValueError('Expected the same length for k_eval and nmodes: %d %d' % (len(k_eval), len(nmodes)))

    # linear power spectrum
    P = interp1d(k, P)(k_eval)
    
    Pl = [legendre(0), legendre(2), legendre(4)]
    cov = np.zeros((len(P), 3, 3))

    nbar = 1.0/nbar

    # nbar <- 1/nbar
    cov[:, 0, 0] = 1.0/nmodes*(
        P**2*(b**4 + 4*b**3*f/3 + 6*b**2*f**2/5 + 4*b*f**3/7 + f**4/9) + P*(2*b**2*nbar + 4*b*f*nbar/3 + 2*f**2*nbar/5) + nbar**2)

    cov[:, 0, 1] = 1.0/nmodes*(
        P**2*(8*b**3*f/3 + 24*b**2*f**2/7 + 40*b*f**3/21 + 40*f**4/99) + P*(8*b*f*nbar/3 + 8*f**2*nbar/7))

    cov[:, 1, 0] = cov[:, 0, 1]

    cov[:, 1, 1] = 1.0/nmodes*(
        P**2*(5*b**4 + 220*b**3*f/21 + 90*b**2*f**2/7 + 1700*b*f**3/231 + 2075*f**4/1287) + P*(10*b**2*nbar + 220*b*f*nbar/21 + 30*f**2*nbar/7) + 5*nbar**2)

    cov[:, 0, 2] = 1.0/nmodes*(
        P**2*(48*b**2*f**2/35 + 96*b*f**3/77 + 48*f**4/143) + 16*P*f**2*nbar/35)

    cov[:, 2, 0] = cov[:, 0, 2]

    cov[:, 1, 2] = 1.0/nmodes*(
        P**2*(48*b**3*f/7 + 816*b**2*f**2/77 + 6960*b*f**3/1001 + 240*f**4/143) + P*(48*b*f*nbar/7 + 272*f**2*nbar/77))

    cov[:, 2, 1] = cov[:, 1, 2]
    
    cov[:, 2, 2] = 1.0/nmodes*(
        P**2*(9*b**4 + 1404*b**3*f/77 + 104166*b**2*f**2/5005 + 11772*b*f**3/1001 + 6399*f**4/2431) + P*(18*b**2*nbar + 1404*b*f*nbar/77 + 34722*f**2*nbar/5005) + 9*nbar**2)
    
    return cov
            



    
    
