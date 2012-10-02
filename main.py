"""
Simple frontend to the Laser Absorption Spectroscopy ANalysis package.
Makes use of the modules 'parse.py', and 'analyze.py'.
Submodules notwithstanding.
"""

# Standard

# Part of package
import gui
import parse
import preprocess

# Third party
import numpy as np
from matplotlib.pyplot import *

dir = gui.pickdir()
settings = parse.config(dir)
signal, wavelength_range, t = parse.data(dir, **settings)

plots = {'t = 0':signal[:,0], 't=20 us':signal[:,500]}
gui.slice(**plots)