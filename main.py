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
guesses = [300, 1e16]#, 5e7]

# Pass transmission profiles to analysis routine
params = [0] * settings['points']
cov = [0] * settings['points']
print "Analyzing data ..."
for i in range(settings['points']):
    try:
        (params[i], cov[i]) = analyze.match(model, freq, transmitted[:, i],
                                            guesses)
        if cov[i] is N.inf:
            cov[i] = N.zeros((len(guesses), len(guesses)))
    except RuntimeError:
        # Set values to zero in case of failure to find a match
        params[i] = N.zeros(len(guesses))
        cov[i] = N.zeros((len(guesses), len(guesses)))


# Everything below this is just data processing
# TODO: Move to a separate file, automate loading of previous calculations
temperatures = N.array([i[0] for i in params])
temperatures_stdev = N.sqrt(N.array([i[0, 0] for i in cov]))
metastables = N.array([i[1] for i in params])
metastables_stdev = N.sqrt(N.array([i[1, 1] for i in cov]))

adir = path.join(target, "Analysis")
if not path.exists(adir):
    makedirs(adir)

output = N.array([temperatures, temperatures_stdev, metastables,
                 metastables_stdev])
with open(path.join(adir, "fit_params.csv"), mode="wb") as f:
    f.write("Temperatures,+-,Metastables,+-\n")
    N.savetxt(f, output.T, delimiter=",")

import matplotlib.pyplot as plt

# plt.plot([i[0][0] for i in cov])
# plt.show()

#time = N.array([0, 29.5, 50, 70.4, 75.6, 99])
time = N.array([0.35, 0.4, 0.5, 1.0, 1.5, 1.75])
check = N.round(1e-6 * time / settings['dt']).astype(int)
pos = 230
plt.hold(True)
for i in check:
    pos = pos + 1
    plt.subplot(pos)
    plt.plot(1e-9 * freq, transmitted[:, i], '.r')
    plt.plot(1e-9 * freq, model(freq, *params[i]), '-k')
    plt.hlines(1.0, 1e-9 * N.min(freq), 1e-9 * N.max(freq), colors='k', linestyles='dashed')
    t = i * 1e6 * settings['dt']
    plt.title('Time = %g $\mu$s' % t)
    plt.axis([1e-9 * N.min(freq), 1e-9 * N.max(freq), 0, 1.1])
plt.hold(False)
# plt.show()
plt.savefig(path.join(adir, r"samples.pdf"))
plt.savefig(path.join(adir, r"samples.png"))
plt.clf()

prepulse = 100
baseline = N.mean(transmitted[:, :prepulse], axis=1)
param_base, cov_base = analyze.match(model, freq, baseline, guesses)
Tbase = param_base[0]

plt.plot(1e6*times, temperatures, '-k')
plt.hold(True)
plt.plot([0, 200], [Tbase, Tbase], '--k')
plt.hold(False)
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Temperature (K)')
plt.axis([0, max(1e6*times), 0, 600])
plt.legend(['Temperatures', 'Pre-pulse, T =%g' % Tbase])
# plt.show()
plt.savefig(path.join(adir, r"temperatures.pdf"))
plt.savefig(path.join(adir, r"temperatures.png"))
plt.clf()

plt.plot(1e6*times, metastables, '-k')
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Line-Integrated Metastable Density (m$^{-2}$)')
plt.axis([0, max(1e6*times), 0, 5e16])
# plt.show()
plt.savefig(path.join(adir, r"metastables.pdf"))
plt.savefig(path.join(adir, r"metastables.png"))
plt.clf()

# plt.plot(1e6*times, drifts, '-k')
# plt.xlabel('Time ($\mu$s)')
# plt.ylabel('Drift-induced Frequency Shift(GHz)')
# plt.axis([0, 200, 0, 1e9])
# # plt.show()
# plt.savefig("drifts.pdf")
# plt.savefig("drifts.png")
# plt.clf()
