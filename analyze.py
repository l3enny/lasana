"""
The mematical component. Generates synthetic spectra. Runs comparison
algorithm. Tells you what you want to know. At present only considers
Doppler broadening and pressure broadening. User must provide pressure
broadening FWHM in angular frequency.
"""

# Standard
import math as m

# Third Party
import numpy as N
from scipy.optimize import curve_fit

# Included
from constants import *
import lineshapes
        
def match(transmission, sigma, guesses, debug=False, **settings):
    """Fits the sigma model to the transmission data.
    
    Uses an input model of the absorption cross section to do a least-
    squares fit to the transmission signal. Returns the fitting param-
    eters and covariances for user processing.
    
    Keyword arguments:
    """
    
    freq = N.linspace(settings['mod_initial'], settings['mod_final'], 
                      settings['samples']) * ma2hz

    coeffs = [0] * settings['points']
    cov = [0] * settings['points']
    for i in range(settings['points']):
        (coeffs[i], cov[i]) = curve_fit(sigma, freq, transmission[:, i], guesses)
            
    return (coeffs, cov)