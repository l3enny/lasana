"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

# Standard
from os import path, makedirs

# Part of package
import analyze
from atoms import He
import gui
import models
import parse
import preprocess
import transition

# Third Party
import numpy as N

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
plasma, times, freq = parse.data(path.join(target, 'Plasma'), **settings)
background = parse.data(path.join(target, 'Background'), **settings)[0]

# Calculated transmission profiles with preprocessor
print "Running preprocessor ..."
transmitted = preprocess.transmission(plasma, background, **settings)

# Define model and some sensible estimates of the parameters
model = models.voigt(transitions, settings['pressure'])
guesses = [300, 1e16]

# Pass transmission profiles to analysis routine
params = [0] * settings['points']
cov = [0] * settings['points']
print "Analyzing data ..."
for i in range(settings['points']):
    try:
        (params[i], cov[i]) = analyze.match(model, freq, transmitted[:, i],
                                            guesses)
    except RuntimeError:
        # Set values to zero in case of failure to find a match
        params[i] = N.zeros(len(guesses))
        cov[i] = N.zeros((len(guesses), len(guesses)))

# Everything below this is just data processing
temperatures = N.array([i[0] for i in params])
metastables = N.array([i[1] for i in params])

#time = N.array([0.35, 0.4, 0.5, 1.0, 1.5, 1.75])
#check = N.round(1e-6 * time / settings['dt']).astype(int)
check = [900, 1100, 1200, 1400, 1500, 1600]
N.savetxt('freq.csv', freq, delimiter=',')
for i in check:
    N.savetxt('m_step%d.csv' % i, transmitted[:, i], delimiter=',')
    N.savetxt('s_step%d.csv' % i, model(freq, *params[i]), delimiter=',')
