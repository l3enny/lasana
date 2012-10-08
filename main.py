"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

# Standard
import sys

# Part of package
import analyze
from constants import *
import gui
import lineshapes
import misc
import parse
import preprocess

# Third party
import numpy as np

initial = {'E':19.81961363388*q, 'J':1, 'g':3}
final = {'E':20.96421789026*q, 'J':0, 'g':1}

signal_dir = gui.pickdir('Pick the signal directory')
signal_settings = parse.config(signal_dir)
signal, signal_t = parse.data(signal_dir, **signal_settings)

reference_dir = gui.pickdir('Pick the reference directory', dir=signal_dir)
reference_settings = parse.config(reference_dir)
reference, reference_t = parse.data(reference_dir, **reference_settings)

if (signal_settings['mod_initial'] != reference_settings['mod_initial']) \
   or (signal_settings['mod_final'] != reference_settings['mod_final']):
    raise ValueError('Signal and reference modulations must be identical')

transmitted = preprocess.transmission(signal, reference)

temperatures = analyze.quickndirty(transmitted, initial, final, debug=True, **signal_settings)