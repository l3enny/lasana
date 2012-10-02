"""
The mathematical component. Generates synthetic spectra. Runs comparison
algorithm. Tells you what you want to know. At present only considers
Doppler broadening and pressure broadening. User must provide pressure
broadening FWHM in angular frequency.
"""

# Third Party
import numpy as np
from scipy import interpolate
from scipy.special import wofz

def voigt(x, sigma, gamma):
    """ Calculates the Voigt function using wofz.
    
    This functions performs the (hopefully not too expensive) calculation
    of the Voigt function for x (scalar or array). Both the Gaussian and 
    Lorentzian are assumed to be centered on x = 0, so the user will have
    to shift the function themselves.
    
    Keyword arguments:
    x -- point at which the function is evaluated
    sigma -- Gaussian linewidth
    gamma -- Lorentzian linewidth
    """
    z = (x + 1j*gamma)/(sigma*sqrt(2))
    return special.wofz(z).real/(sigma*sqrt(2*np.pi))
    
def sigma(w0, w):
    
    
    
def peak_find():
    

def get_temp(signal, wavelengths):
    