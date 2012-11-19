""" Convenient lineshape definitions.

Several different lineshapes which are frequently found in spectroscopy.
These can be used to build more complex models.
"""

from math import sqrt, log
import numpy as np
from scipy.special import wofz

def lorentzian(x, gamma):
    """Evaluation of a lorentzian or cauchy distribution
    """
    return 1/(np.pi * gamma * (1 + np.power(x/gamma, 2)))

def gaussian(x, sigma):
    """ Evaluation of a gaussian distribution.
    """
    return 1/(sigma*sqrt(2*np.pi)) * np.exp(-0.5 * np.power((x+x0)/sigma, 2))

def voigt(x, sigma, gamma):
    """Evaluation of a Voigt profile using the complex error function.
    
    This function evaluates a Voigt function and returns it to the user.
    It is an "exact" calculation of the complex error function making it
    more expensive than the pseudo-Voigt.
    
    Keyword arguments:
    sigma -- Standard deviation of the Gaussian portion
    gamma -- Standard deviation of the Lorentzian portion
    """
    z = (x + 1j*gamma)/(sigma*sqrt(2))
    return wofz(z).real/(sigma*sqrt(2*np.pi))
    
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