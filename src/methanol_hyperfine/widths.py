"""Voigt width parameterization used by L26."""
import numpy as np
from functools import lru_cache
from scipy.optimize import brentq
from scipy.special import voigt_profile
GAUSSIAN_FWHM=np.sqrt(8*np.log(2))

def voigt_fwhm(sigma, gamma):
    if sigma < 0 or gamma < 0 or sigma+gamma <= 0:
        raise ValueError('Widths must be nonnegative, with positive total width')
    if sigma == 0:
        return 2*gamma
    if gamma == 0:
        return GAUSSIAN_FWHM*sigma
    peak = voigt_profile(0., sigma, gamma)
    return 2*brentq(lambda x: voigt_profile(x, sigma, gamma)-peak/2,
                    0., 5*(sigma+gamma), xtol=2e-13)

@lru_cache(maxsize=8192)
def _unit_width(eta):
    return voigt_fwhm((1-eta)/GAUSSIAN_FWHM, eta/2)

def sigma_gamma(fwhm, eta):
    if fwhm <= 0 or not 0 <= eta <= 1:
        raise ValueError('Positive FWHM and shape coordinate in [0,1] required')
    scale = fwhm/_unit_width(float(eta))
    return scale*(1-eta)/GAUSSIAN_FWHM, scale*eta/2

def fwhm_eta(sigma, gamma):
    return voigt_fwhm(sigma, gamma), 2*gamma/(GAUSSIAN_FWHM*sigma+2*gamma)
