"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

# Standard

# Part of package
import parse
import preprocess

# Third party
import numpy as np

dir = parse.pickdir()
signal, wavelength, v, i = parse.grab(dir)