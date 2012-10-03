"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

# Standard
import sys

# Part of package
import analyze
import gui
import lineshapes
import parse
import preprocess

# Third party
import numpy as np

signal_dir = gui.pickdir('Pick the signal directory')
signal_settings = parse.config(signal_dir)
signal, signal_shifts, signal_t = parse.data(signal_dir, **signal_settings)

reference_dir = gui.pickdir('Pick the reference directory', dir=signal_dir)
reference_settings = parse.config(reference_dir)
reference, reference_shifts, reference_t = parse.data(reference_dir, 
                                                      **reference_settings)
if (signal_shifts != reference_shifts).all():
    raise ValueError('Signal and reference modulations must be identical')

signal_corrected = preprocess.transmission(signal, reference)

temperatures = analyze.quickndirty(signal_corrected, signal_shifts)

# plots = {'t=20 us':signal_corrected[:,500]}
# gui.slice(**plots)

# gui.contour(signal_corrected)