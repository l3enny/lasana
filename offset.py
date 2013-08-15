"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

# Standard
from os import path, makedirs
import math as m

# Part of package
from atoms import He

import gui
import lineshapes
import parse
import preprocess
import transition

# Third Party
import numpy as N
from scipy.optimize import curve_fit
from scipy.constants import k, c

torr2hz = 25.6e6 # pressure broadening coefficient

# Define transitions to simulate
# TODO: Find better way to set up a database of transmission constants than
# the one used by crammer. SQLite perhaps? Is there an API for the NIST ASD?
A = 1.0216e7
D0 = transition.Transition(He.II3S1(), He.II3P0(), A)
D1 = transition.Transition(He.II3S1(), He.II3P1(), A)
D2 = transition.Transition(He.II3S1(), He.II3P2(), A)
transitions = [D0, D1, D2]

# Pick directory of signal scan, then parse and load data
target = gui.pickdir('Pick the data directory')
print "\nProcessing", target
settings = parse.config(target)
print "Loading data ..."
plasma, times, freq = parse.data(path.join(target, 'Plasma'), obsolete=False, **settings)
background = parse.data(path.join(target, 'Background'), obsolete=False, **settings)[0]

# Calculated transmission profiles with preprocessor
print "Running preprocessor ..."
transmitted = preprocess.transmission(plasma, background, **settings)

# Define model and some sensible estimates of the parameters
def model(x, T, amp, offset):
    # Pressure broadening/lorentzian part of profile
    f = 0
    V = lineshapes.voigt

    # Assumes origin is located at the center of the first listed transition
    origin = transitions[0].f
    for t in transitions:
        fwhm_a = t.A + settings['pressure'] * torr2hz
        gamma = fwhm_a/2
        # Doppler broadening/gaussian part of the profile
        fwhm_d = N.sqrt((8*m.log(2)) * k*T / (t.M*c**2)) * t.f
        sigma = fwhm_d/(2*m.sqrt(2*m.log(2)))
        temp = V(x + t.f - origin + offset, sigma, gamma)
        temp = temp * t.l**2 * t.A * (t.gj/t.gi) / (8 * N.pi)
        f += temp
    return N.exp(- amp * f)
guesses = [300, 5e15, -10e6]

# Find the offset for the "largest signal"
mn = N.std(transmitted, axis=0)
largest = N.where(max(mn) == mn)[0]
(p, co) = curve_fit(model, freq, transmitted[:, largest][:,0], guesses)

print "Largest:", largest
print "Offset is:", p[2], "Hz"
