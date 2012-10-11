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
import scipy.optimize

# Included
from constants import *
import lineshapes
        
def get_temp(signal, transitions, debug=False, **settings):
    """ Finds the best-fit temperature for a given signal.
    
    Generates a set of synthetic absorption spectra for a given set of
    transitions and compares each one to the measured spectrum. Then
    creates a spline for the resulting surface errors and determines 
    the minimum error.
    
    Keyword arguments: 
    """
    
    def model(x, amp, center, drift, T):
        fwhm_a = settings['pressure'] * torr2hz
        f = 0
        V = lineshapes.voigt
        origin = transitions[0].f
        for t in transitions:
            fwhm_d = N.sqrt((8*m.log(2)) * kB*T / (t.M*c**2)) * t.f
            temp = amp * 0.5 * (V(x, fwhm_d, fwhm_a, drift, t.f - origin)
                                + V(x, fwhm_d, fwhm_a, -drift, t.f - origin))
            # temp = temp * 2 * N.pi * (t.gj/t.gi) * t.l**2 * t.A
            f += temp
        return f
   
    guesses = (5e7, -0.1e7, 0.1e7, 300)
    absorption = 1 - signal + settings['voffset']
    
    df = (settings['mod_initial']*ma2hz, settings['mod_final']*ma2hz)
    freq = N.linspace(df[0], df[1], settings['samples'])
    
    # import matplotlib.pyplot as plt
    # x = N.linspace(-4e9, 40e9, 1e3)
    # plt.plot(x, model(x, 1.0, 0.0, 0.0, 300))
    # plt.show()
    
    temperatures = N.zeros(settings['points'])
    covariances = [0] * settings['points']
    for i in range(settings['points']):
        try:
            p = scipy.optimize.curve_fit(model, freq, absorption[:, i],
                                         guesses)
        except ValueError:
            temperatures[i] = 0
            covariances[i] = None
        else:
            # print p
            temperatures[i] = p[0][3]
            covariances[i] = p[1]
    return temperatures

def quickndirty(signal, transitions, debug=False, **settings):

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
    
    absorption = 1 - signal + settings['voffset']

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
    
    baseline = N.mean(absorption[:, 0:200], axis=1)
    peak = N.max(baseline)
    errors = N.abs(peak * absorption_synthetic - baseline)
    errors_collapsed = N.sum(errors, axis=1)
    Tbase = optimize(temperatures, errors_collapsed)
    print "Baseline Temperature = %g (K)" % Tbase
    match = sigma(transition, frequencies, Tbase, **settings)
    
    plt.plot(1e-9*frequencies, baseline, '.r')
    plt.hold(True)
    plt.plot(1e-9*frequencies, peak * match/N.max(match), '-k')
    plt.axis([ma2hz*1e-9*settings['mod_initial'], ma2hz*1e-9*settings['mod_final'], 0, 1])
    plt.xlabel('Frequency Shift (GHz)')
    plt.ylabel('Absorption (a.u.)')
    plt.title('Baseline Temperature, T = %g' % Tbase)
    plt.show()
    
    if debug:
        # check = [220, 225, 230, 235, 240, 245]
        check = 1e-6 * N.array([69, 100, 150, 200, 350, 300])/settings['dt']
        check = check.astype(int)
        plt.hold(True)
        pos = 230
        for i in check:
            match = sigma(transition, frequencies, Tcalc[i], **settings)
            
            header = "Frequency Shifts (Hz), Calculated Signal) (a.u.), Measured Signal (a.u.)"
            data = N.column_stack((frequencies, match, absorption[:, i]))
            
            N.savetxt("t=%gus (T=%g).csv" % (i*settings['dt']*1e6, Tcalc[i]), data, delimiter=",")
            
            pos = pos + 1
            plt.subplot(pos)
            plt.plot(1e-9*frequencies, N.max(absorption[:, i]) * match/N.max(match), '-k')
            plt.plot(1e-9*frequencies, absorption[:, i], '.r')
            plt.axis([ma2hz*1e-9*settings['mod_initial'], ma2hz*1e-9*settings['mod_final'], 0, 1])
            t = i*settings['dt']*1e6
            plt.title('Time = %g $\mu$s' % t)
            
            tmp_error = sum(abs(absorption[:, i] - N.max(absorption[:, i]) * match/N.max(match)))
            
            print 'Temperature (step %g):' % i, Tcalc[i]
            print 'Error (step %g):' % i, tmp_error
        plt.show()
        
        # Plot temperature decay rate
        # plt.clf()
        # time = settings['dt'] * N.arange(settings['points'])
        # fit = lambda x, a, b, c: a*N.exp(-b*x) + c
        # start = 750
        # param = scipy.optimize.curve_fit(fit, time[start:], Tcalc[start:], p0=[Tcalc[start]-Tbase[0], 1e5*settings['dt'], Tbase[0]])
        # print "Estimated Baseline: %g (K)" % param[0][2]
        # efold = 1e6/param[0][1]
        # print "Estimate e-Folding time: %g (us)" % efold
        # plt.plot(time[start:], Tcalc[start:], '.r')
        # plt.plot(time[start:], fit(time[start:], *param[0]), '-k')
        # plt.xlabel('Time ($\mu$s)')
        # plt.ylabel('Temperature (K)')
        # plt.legend(('Calculated Temperatures', 'f(t) = %gexp(-%gt)+%g' % (param[0][0], param[0][1], param[0][2])))
        # plt.show()
    
    return Tcalc