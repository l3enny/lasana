"""
The mematical component. Generates synthetic spectra. Runs comparison
algorithm. Tells you what you want to know. At present only considers
Doppler broadening and pressure broadening. User must provide pressure
broadening FWHM in angular frequency.
"""

# Standard
import math as m
import matplotlib.pyplot as plt
import matplotlib.text as mtext

# Third Party
import numpy as N
from scipy import interpolate

# Included
from constants import *
import lineshapes
import misc
    
def sigma(transition, frequencies, Tg, debug=False, **settings):
    M = transition.M
    f_0 = (transition.dE)/h
    l = misc.vac2air(c/f_0)
    fwhm_g = N.sqrt((8*m.log(2)) * kB*Tg / (M*c**2)) * f_0
    fwhm_a = settings['pressure'] * torr2hz
    
    # Change this line to determine the lineshape that you desire
    function = lineshapes.voigt(fwhm_g, fwhm_a)
    
    profile = function(frequencies + settings['offset'])
    return profile
    # return final['g'] / initial['g'] * l**2 * A / 
    
def peak_find(debug=False):
    pass

def optimize(x, errors, debug=False):
    tck1 = interpolate.splrep(x, errors, s=0)
    der1 = interpolate.splev(x, tck1, der=1)
    tck2 = interpolate.splrep(x, der1, s=0)
    roots = interpolate.sproot(tck2)
    if len(roots) is not 1:
        print "Warning! Multiple minima, returning average."
        return N.mean(roots)
    return roots
    
def get_temp(Ti, Tf, transitions, signal, steps=50, debug=False, **settings):
    """ Finds the best-fit temperature for a given signal.
    
    Generates a set of synthetic absorption spectra for a given set of
    transitions and compares each one to the measured spectrum. Then
    creates a spline for the resulting surface errors and determines 
    the minimum error.
    
    Keyword arguments:
    Ti -- 
    """
    temperatures = N.linspace(Ti, Tf, steps)

def quickndirty(signal, transitions, debug=False, **settings):

    if debug:
        check = [250, 450, 700, 1000, 1500, 1999]

    Ti = 200
    Tf = 1500
    temperatures = N.linspace(Ti, Tf, 50)

    freq_shift = (settings['mod_initial']*ma2hz, settings['mod_final']*ma2hz)
    frequencies = N.linspace(freq_shift[0], freq_shift[1], settings['samples'])
    
    absorption_synthetic = N.zeros((len(temperatures), settings['samples']))
    
    for transition in transitions:
        for i in range(len(temperatures)):
            absorption_synthetic[i, :] += sigma(transition, frequencies,
                                                temperatures[i], **settings)
    
    # TODO: Eliminate normalization by using proper calculation of pathlength
    # absorption
    norms = N.array(N.max(absorption_synthetic, axis=1))
    absorption_synthetic = absorption_synthetic/norms[:, N.newaxis]
    
    absorption = 1 - signal

    # Initialize calculation arrays
    Tcalc = N.zeros((absorption.shape[1]))
    
    # iterate through each time point
    for i in range(settings['points']):
        # TODO: Associated with normalization, should not be necessary once
        # real metastable absorption is included
        peak = N.max(absorption[:,i])
        errors = N.abs(peak * absorption_synthetic - absorption[:, i])
        
        # Optimize via spline interpolation, slower, better?
        errors_collapsed = N.sum(errors, axis=1)
        Tcalc[i] = optimize(temperatures, errors_collapsed)
    
    if debug:
        plt.hold(True)
        pos = 230
        for i in check:
            match = sigma(transition, frequencies, Tcalc[i], **settings)
        
            pos = pos + 1
            plt.subplot(pos)
            plt.plot(frequencies, N.max(absorption[:, i]) * match/N.max(match), '-k')
            plt.plot(frequencies, absorption[:, i], '.r')
            mtext.Text(3e9, 0.9, 'T$_g$ = %g' % Tcalc[i])
            plt.axis([-4e9, 4e9, 0, 1])
            t = i*settings['dt']*1e6
            plt.title('Time = %g $\mu$s' % t)
            
            tmp_error = sum(abs(absorption[:, i] - N.max(absorption[:, i]) * match/N.max(match)))
            
            print 'Temperature (step %g):' % i, Tcalc[i]
            print 'Error (step %g):' % i, tmp_error
    
    plt.show()
    return Tcalc