"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

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
signal_dir = gui.pickdir('Pick the signal directory')
signal_settings = parse.config(signal_dir)
signal, signal_t, freq = parse.data(signal_dir, **signal_settings)

# Pick directory of signal scan, then parse and load data
reference_dir = gui.pickdir('Pick the reference directory', dir=signal_dir)
reference_settings = parse.config(reference_dir)
reference, reference_t, freq = parse.data(reference_dir, **signal_settings)

# Calculated transmission profiles with preprocessor
transmitted = preprocess.transmission(signal, reference, **signal_settings)

# Define model and some sensible estimates of the parameters
model = models.bimodal_voigt(transitions, signal_settings['pressure'])
guesses = [320, 1e15, 5e7]

# Pass transmission profiles to analysis routine
params = [0] * signal_settings['points']
cov = [0] * signal_settings['points']
for i in range(signal_settings['points']):
    (params[i], cov[i]) = analyze.match(model, freq, transmitted[:, i], guesses)

temperatures = N.array([i[0] for i in params])
metastables = N.array([i[1] for i in params])
drifts = N.abs(N.array([i[2] for i in params]))

name = "fit_params.csv"
N.savetxt(name, params, delimiter=",")

import matplotlib.pyplot as plt
time = N.array([5, 10, 25, 45, 70, 100])
check = N.round(1e-6 * time / signal_settings['dt']).astype(int)
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

prepulse = 19e-6/signal_settings['dt']
baseline = N.mean(transmitted[:, :prepulse], axis=1)
param_base, cov_base = analyze.match(model, freq, baseline, guesses)
Tbase = param_base[0]

plt.plot(1e6*signal_t, temperatures, '-k')
plt.hold(True)
plt.plot([0, 200], [Tbase, Tbase], '--k')
plt.hold(False)
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Temperature (K)')
plt.axis([0, 200, 0, 600])
plt.legend(['Temperatures', 'Pre-pulse, T =%g' % Tbase])
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

plt.plot(1e6*signal_t, drifts, '-k')
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Drift-induced Frequency Shift(GHz)')
plt.axis([0, 200, 0, 1e9])
# plt.show()
plt.savefig("drifts.pdf")
plt.savefig("drifts.png")
plt.clf()