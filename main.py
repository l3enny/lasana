"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

# Standard
import math as m

# Part of package
import analyze
from atoms import He
from constants import *
import gui
import lineshapes
import parse
import preprocess
import transition

# Third Party
import numpy as N

# Define transitions of interest to simulation
A = 1.0216e7
D0 = transition.Transition(He.II3S1(), He.II3P0(), A)
D1 = transition.Transition(He.II3S1(), He.II3P1(), A)
D2 = transition.Transition(He.II3S1(), He.II3P2(), A)
transitions = [D0, D1, D2]

# Pick directory of signal scan, then parse and load data
signal_dir = gui.pickdir('Pick the signal directory')
signal_settings = parse.config(signal_dir)
signal, signal_t = parse.data(signal_dir, **signal_settings)

# Pick directory of signal scan, then parse and load data
reference_dir = gui.pickdir('Pick the reference directory', dir=signal_dir)
reference_settings = parse.config(reference_dir)
reference, reference_t = parse.data(reference_dir, **signal_settings)

# Calculated transmission profiles with preprocessor
transmitted = preprocess.transmission(signal, reference)

# Define model used to evaluate the transition probabilities
guesses = [320, 1e13, 1e9]
def bimodal_voigt(x, T, amp, drift, center=0.0):
    # Pressure broadening/lorentzian part of profile
    fwhm_a = signal_settings['pressure'] * torr2hz
    gamma = fwhm_a/2
    f = 0
    V = lineshapes.voigt
    # TODO: Manual origin setting, assumes the 0.0 frequency shift
    # is the D0 peak. Should be a better way to do this.
    origin = D0.f
    for t in transitions:
        # Doppler broadening/gaussian part of the profile
        fwhm_d = N.sqrt((8*m.log(2)) * kB*T / (t.M*c**2)) * t.f
        sigma = fwhm_d/(2*m.sqrt(2*m.log(2)))
        temp = 0.5 * (V(x, sigma, gamma, t.f - origin + drift)
                      + V(x, sigma, gamma, t.f - origin - drift))
        temp = temp * 2 * N.pi * (t.gj/t.gi) * t.l**2 * t.A
        f += temp
    return amp * f

# Pass transmission profiles to analysis routine
(params, cov) = analyze.match(transmitted, bimodal_voigt, guesses, 
                              debug=False, **signal_settings)
                                    
temperatures = N.array([i[0] for i in params])
amplitudes = N.array([i[1] for i in params])
drifts = N.array([i[2] for i in params])

modded = [(temperatures[i], 1, drifts[i]) for i in range(signal_settings['points'])]
metastables = analyze.densities(transmitted, bimodal_voigt, modded, 
                                debug=False, **signal_settings)
                                   

import matplotlib.pyplot as plt
plt.plot(1e6*signal_t, temperatures, '-k')
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Temperature (K)')
plt.show()
plt.plot(1e6*signal_t, metastables, '-k')
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Line-Integrated Metastable Density (m$^{-2}$)')
plt.show()