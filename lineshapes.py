from math import sqrt, log
import numpy as np
from scipy.special import wofz

def lorentzian(gamma):
    def func(x):
        return 1/(np.pi * gamma * (1 + np.power(x/gamma, 2)))
    return func

def gaussian(fwhm):
    sigma = fwhm/(2*sqrt(2*log(2)))
    def func(x):
        return 1/(sigma*sqrt(2*np.pi)) * np.exp(-0.5 * np.power(x/sigma, 2))
    return func

def voigt(fwhm_d, fwhm_p):
    """ Returns a Voigt function using wofz.
    
    This function generates a Voigt function and returns it to the user.
    It is an "exact" calculation of the complex error function making it
    more expensive than the pseudo-Voigt.
    
    Keyword arguments:
    sigma -- Gaussian linewidth
    gamma -- Lorentzian linewidth
    """
    sigma = fwhm_d/(2*sqrt(2*log(2)))
    gamma = fwhm_p/2
    def func(x):
        z = (x + 1j*gamma)/(sigma*sqrt(2))
        return wofz(z).real/(sigma*sqrt(2*np.pi))
    return func
    
def pseudo_voigt(gamma, eta):
    """Generates a pseudo-Voigt function.
    
    The pseudo-Voigt is a less expensive way of calculating the Voigt
    lineshape. It uses a weighted sum of the individual Gaussian and
    Lorentzian distributions.
    """
    sigma = 2*gamma*sqrt(2*log(2))
    G = gaussian(sigma)
    L = lorentzian(gamma)
    def func(x):
        return (1 - eta) * G(x) + eta * L(x)
    return func