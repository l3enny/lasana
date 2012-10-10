"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

# Part of package
import analyze
from atoms import He
from constants import *
import gui
import parse
import preprocess
import transition

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

# Check to ensure that both scans were over the same range
# TODO: Move this out of the main body, doesn't make sense here
if (signal_settings['mod_initial'] != reference_settings['mod_initial']) \
   or (signal_settings['mod_final'] != reference_settings['mod_final']):
    raise ValueError('Signal and reference modulations must be identical')
    
reference, reference_t = parse.data(reference_dir, **signal_settings)

# Calculated transmission profiles with preprocessor
transmitted = preprocess.transmission(signal, reference)

# Pass transmission profiles to analysis routine
temperatures = analyze.quickndirty(transmitted, transitions, debug=True,
                                   **signal_settings)
                                   
import matplotlib.pyplot as plt
plt.plot(1e6*signal_t, temperatures, '-k')
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Temperature (K)')
plt.show()