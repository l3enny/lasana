"""
The mathematical component. Generates synthetic spectra. Runs comparison
algorithm. Tells you what you want to know. At present only considers
Doppler broadening and pressure broadening. User must provide pressure
broadening FWHM in angular frequency.
"""

# Standard
import math
import matplotlib.pyplot as plt

# Third Party
import numpy as np
from scipy import interpolate

# Included
import lineshapes
    
def sigma(w0, w):
    pass
    
def peak_find():
    pass

def get_temp(signal, wavelengths):
    pass
    
def quickndirty(signal, wavelengths):
    torr = 7.0
    torr2hz = 25.6e6
    M = 4.002602 * 1.660538921e-27
    c = 299792458
    f_0 = c/1082.9091140e-9
    kB = 1.3806503e-23
    temperatures = np.linspace(200, 1400, 1200)
    adjust = 0.05e9
    wavelengths = wavelengths*1e9 + adjust
    fwhm_p = torr*torr2hz
    fwhm_d = np.sqrt((8*math.log(2))*kB*temperatures/(M*c**2)) * f_0
    synthetic_profiles = np.zeros((len(temperatures), len(wavelengths)))
    for i in range(len(temperatures)):
        function = lineshapes.voigt(fwhm_d[i], fwhm_p)
        # function = lineshapes.gaussian(fwhm_d[i])
        synthetic_profiles[i, :] = function(wavelengths)
    norms = np.array(np.max(synthetic_profiles, axis=1))
    synthetic_profiles = synthetic_profiles/norms[:, np.newaxis]
    
    transmission = 1 - signal
    # norms = np.max(transmission, axis=0)
    # transmission = transmission/norms
    # plt.contourf(transmission, levels=np.linspace(np.min(transmission), np.max(transmission), 100))
    # plt.show()

    calculated_temperatures = np.zeros((transmission.shape[1]))
    # iterate through each time point
    for i in range(transmission.shape[1]):
        peak = np.max(transmission[:,i])
        errors = np.abs(peak * synthetic_profiles - transmission[:, i])
        errors_collapsed = np.sum(errors, axis=1)
        minimum = np.min(errors_collapsed)
        index = np.where(errors_collapsed == minimum)
        calculated_temperatures[i] = temperatures[index[0][0]]
        if i == 365:
            print "Minimum Error =", minimum
            print "Temperature = ", temperatures[index[0][0]]
            plt.hold(True)
            plt.plot(wavelengths, peak * synthetic_profiles[index[0][0], :])
            plt.plot(wavelengths+adjust, transmission[:, i])
            plt.legend(['Synthetic'])
            plt.show()
        
    plt.plot(calculated_temperatures)
    plt.ylabel('Temperature (K)')
    plt.show()
    return calculated_temperatures