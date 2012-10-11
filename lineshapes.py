from math import sqrt, log
import numpy as np
from scipy.special import wofz

# TODO: Rewrite for SciPy curve_fit routine
def lorentzian(x, gamma):
    return 1/(np.pi * gamma * (1 + np.power(x/gamma, 2)))

# TODO: Rewrite for SciPy curve_fit routine
def gaussian(x, fwhm):
    sigma = fwhm/(2*sqrt(2*log(2)))
    return 1/(sigma*sqrt(2*np.pi)) * np.exp(-0.5 * np.power(x/sigma, 2))

# TODO: Rewrite for SciPy curve_fit routine
def voigt(x, fwhm_d, fwhm_a, offset=0.0, center=0.0):
    """ Returns a Voigt function using wofz.
    
    This function generates a Voigt function and returns it to the user.
    It is an "exact" calculation of the complex error function making it
    more expensive than the pseudo-Voigt.
    
    Keyword arguments:
    fwhm_d -- FWHM of the Gaussian portion
    fwhm_a -- FWHM of the Lorentzian portion
    """
    sigma = fwhm_d/(2*sqrt(2*log(2)))
    # TODO: Verify the factor of 0.5 for gamma!
    gamma = fwhm_a/2
    z = (x + offset + center + 1j*gamma)/(sigma*sqrt(2))
    return wofz(z).real/(sigma*sqrt(2*np.pi))
    
# TODO: Rewrite for SciPy curve_fit routine
def pseudo_voigt(gamma, eta):
    """Generates a pseudo-Voigt function.
    
    The pseudo-Voigt is a less expensive way of calculating the Voigt
    lineshape. It uses a weighted sum of the individual Gaussian and
    Lorentzian distributions.
    """
    sigma = 2*gamma*sqrt(2*log(2))
    G = gaussian(sigma)
    L = lorentzian(gamma)
    return (1 - eta) * G(x) + eta * L(x)