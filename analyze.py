"""
The mematical component. Generates synthetic spectra. Runs comparison
algorithm. Tells you what you want to know. At present only considers
Doppler broadening and pressure broadening. User must provide pressure
broadening FWHM in angular frequency.
"""

# Third Party
from scipy.optimize import curve_fit
        
def match(sigma, freq, measured, guesses, debug=False):
    """Fits the sigma model to the transmission data.
    
    Uses an input model of the absorption cross section to do a least-
    squares fit to the transmission signal. Returns the fitting param-
    eters and covariances for user processing. Thin wrapper over the SciPy
    curve_fit routine
    
    Keyword arguments:
    """
    return curve_fit(sigma, freq, measured, guesses)