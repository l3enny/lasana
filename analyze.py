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
import gui
import lineshapes
        
def match(transmission, sigma, guesses, debug=False, **settings):
    """Fits the sigma model to the transmission data.
    
    Uses an input model of the absorption cross section to do a least-
    squares fit to the transmission signal. Returns the fitting param-
    eters and covariances for user processing.
    
    Keyword arguments:
    """

    absorption = 1 - transmission + settings['voffset']
    
    freq = N.linspace(settings['mod_initial'], settings['mod_final'], 
                      settings['samples']) * ma2hz

    coeffs = [0] * settings['points']
    cov = [0] * settings['points']
    for i in range(settings['points']):
        (coeffs[i], cov[i]) = curve_fit(sigma, freq, absorption[:, i], guesses)
    if debug:
        check = [400, 800, 1200, 1600, 2000, 2400]
        for i in check:
            print "coeffs[%g]:" % i, coeffs[i]
            set = {'Matched':sigma(freq, *coeffs[i]),
                   'Measured':absorption[:, i]}
            gui.slice(**set)
            
    return (coeffs, cov)
    
def densities(transmission, sigma, params, debug=False, **settings):
    freq = N.linspace(settings['mod_initial'], settings['mod_final'], 
                      settings['samples']) * ma2hz
    zero = N.abs(freq).argmin()
    set = {'Peak Falloff':transmission[zero, :]}
    gui.slice(**set)
    
    n = N.zeros(settings['points'])
    nint = N.zeros(settings['points'])
    for i in range(settings['points']):
        n[i] = -m.log(abs(transmission[zero, i])) / (sigma(0.0, *params[i]))
        nint[i] = -N.mean(m.log(abs(transmission[:, i])) / (sigma(freq, *params[i])))
    set = {'Single Point':n, 'Averaged':nint}
    gui.slice(**set)
    return n