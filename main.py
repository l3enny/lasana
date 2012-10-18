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
import models
import parse
import preprocess
import transition

# Third Party
import numpy as N

# Define transitions to simulate
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
transmitted = preprocess.transmission(signal, reference, **signal_settings)

# Define model and some sensible estimates of the parameters
model = models.voigt(transitions, settings['pressure'])
guesses = [320, 1e13]#, 5e7]
    
# Pass transmission profiles to analysis routine
(params, cov) = analyze.match(transmitted, model, guesses, debug=True,
                              **signal_settings)
                                    
temperatures = N.array([i[0] for i in params])
metastables = N.array([i[1] for i in params])
if model is bimodal_voigt:
    drifts = c * N.array([abs(i[2]) for i in params])/D0.f

name = "fit_params.csv"
N.savetxt(name, params, delimiter=",")

import matplotlib.pyplot as plt
freq = N.linspace(signal_settings['mod_initial'], signal_settings['mod_final'], signal_settings['samples']) * ma2hz + signal_settings['offset']
time = N.array([62.4, 64.6, 65.2, 65.8, 66.4, 67.2])
check = 1e-6 * time / signal_settings['dt']
check = check.astype(int)
pos = 230
plt.hold(True)
for i in check:
    pos = pos + 1
    plt.subplot(pos)
    plt.plot(1e-9 * freq, transmitted[:, i], '.r')
    plt.plot(1e-9 * freq, model(freq, *params[i]), '-k')
    t = i * 1e6 * signal_settings['dt']
    plt.title('Time = %g $\mu$s' % t)
    plt.axis([1e-9 * N.min(freq), 1e-9 * N.max(freq), 0, 1])
plt.hold(False)
# plt.show()
plt.savefig("samples.pdf")
plt.savefig("samples.png")
plt.clf()

plt.plot(1e6*signal_t, temperatures, '-k')
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Temperature (K)')
plt.axis([0, 200, 0, 600])
# plt.show()
plt.savefig("temperatures.pdf")
plt.savefig("temperatures.png")
plt.clf()

plt.plot(1e6*signal_t, metastables, '-k')
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Line-Integrated Metastable Density (m$^{-2}$)')
plt.axis([0, 200, 0, 5e16])
# plt.show()
plt.savefig("metastables.pdf")
plt.savefig("metastables.png")
plt.clf()

if model is bimodal_voigt:
    plt.plot(1e6*signal_t, drifts, '-k')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Drift Velocities (m/s)')
    plt.axis([0, 200, 0, 500])
    # plt.show()
    plt.savefig("drifts.pdf")
    plt.savefig("drifts.png")