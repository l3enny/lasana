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
    
def sigma(initial, final, frequencies, Tg, debug=False, **settings):
    A = 1.0216e7

    f_0 = (final['E'] - initial['E'])/h
    l = misc.vac2air(c/f_0)
    fwhm_g = N.sqrt((8*m.log(2))*kB*Tg/(M_He*c**2)) * f_0
    fwhm_a = settings['pressure'] * torr2hz
    
    # Change this line to determine the lineshape that you desire
    function = lineshapes.voigt(fwhm_g, fwhm_a)
    profile = function(frequencies + settings['offset'])
    return profile
    # return 2 * N.pi * m.log(2) * final['g'] / initial['g'] * l**2 * A / 
    
def peak_find(debug=False):
    pass

def get_temp(signal, frequencies, debug=False):
    pass
    
def quickndirty(signal, initial, final, debug=False, **settings):

    if debug:
        check = [250, 450, 700, 1000, 1500, 1999]

    Ti = 200
    Tf = 1500
    temperatures = N.linspace(Ti, Tf, Tf - Ti + 1)

    freq_shift = (settings['mod_initial']*ma2hz, settings['mod_final']*ma2hz)
    frequencies = N.linspace(freq_shift[0], freq_shift[1], settings['samples'])
    
    absorption_synthetic = N.zeros((len(temperatures), settings['samples']))
    
    for i in range(len(temperatures)):
        absorption_synthetic[i, :] = sigma(initial, final, frequencies,
                                           temperatures[i], **settings)
    
    # TODO: Eliminate normalization by using proper calculation of pathlength
    # absorption
    norms = N.array(N.max(absorption_synthetic, axis=1))
    absorption_synthetic = absorption_synthetic/norms[:, N.newaxis]
    
    absorption = 1 - signal

    # Initialize calculation arrays
    Tcalc = N.zeros((absorption.shape[1]))
    Terror = N.zeros(absorption.shape[1])
    
    # iterate through each time point
    for i in range(settings['points']):
        # TODO: Associated with normalization, should not be necessary once
        # real metastable absorption is included
        peak = N.max(absorption[:,i])
        errors = N.abs(peak * absorption_synthetic - absorption[:, i])
        
        errors_collapsed = N.sum(errors, axis=1)
        minimum = N.min(errors_collapsed)
        Terror[i] = minimum/N.sum(absorption[:,i])
        index = N.where(errors_collapsed == minimum)
        Tcalc[i] = temperatures[index[0][0]]
        
        if debug:
            if i == 0:
                position = 230
            if i in check:
                position = position + 1
                plt.subplot(position)
                print 'Temperature (step %g):' % i, Tcalc[i]
                print 'Error (step %g):' % i, Terror[i]
                plt.plot(frequencies, peak * absorption_synthetic[index[0][0], :], '-k')
                plt.plot(frequencies, absorption[:, i], '.r')
                mtext.Text(3e9, 0.9, 'T$_g$ = %g' % Tcalc[i])
                plt.axis([-4e9, 4e9, 0, 1])
                t = i*settings['dt']*1e6
                plt.title('Time = %g $\mu$s' % t)
    plt.show()
    return Tcalc