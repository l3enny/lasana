﻿"""
Processes the data from whatever format it was captured in and renders
it in a usable form to whoever's knocking.
"""

# Standard modules
import ConfigParser
import os
from os import path
import sys

# Third party
import numpy as np

# Part of package
import preprocess
    
def config(top_dir='.'):
    # Parse the settings file, check for proper formatting
    config = ConfigParser.RawConfigParser()
    try:
        read_return = config.read(path.join(top_dir, 'settings.cfg'))
    except ConfigParser.MissingSectionHeaderError:
        print "Settings file is not properly formatted, quitting."
        sys.exit(1)
    if read_return == []:
        raise ConfigParser.Error('No settings file detected.')
    
    try:
        samples = config.getint('Settings', 'Samples')
        points = config.getint('Settings', 'Points')
        mod_initial = config.getfloat('Settings', 'Initial Modulation')
        mod_final = config.getfloat('Settings', 'Final Modulation')
        averages = config.getint('Settings', 'Averages')
        dt = config.getfloat('Settings', 'dt')
    except ConfigParser.NoOptionError:
        print "Settings file is not properly formatted, quitting."
        sys.exit(1)
    return {'samples':samples, 'points':points, 'mod_initial':mod_initial,
            'mod_final':mod_final, 'averages':averages, 'dt':dt}
    
def data(top_dir='.', samples=1, points=1, averages=1, mod_initial=0,
         mod_final=0, dt=1e-6):
    """Custom parser for measurement files.
    
    Reads in the data related to a particular set of measurements. Assumes
    that the tuning of the diode has a 0.8 GHz/mA relation. Reads in
    voltage and current signals, but does not return them. Only returns
    the wavelength *range* and signals. The wavelength range is a 1D numpy
    array. The signals are a 2D numpy array where the first index is the
    sample index (essentially wavelength) and the second index is the time
    index.
    
    Keyword arguments:
    top_dir -- the reference directory for measurement files (default '.')
    """
   
    # Initialize storage arrays
    pd_b = np.zeros((samples, points))
    ref_b = np.zeros((samples, points))
    v_b = np.zeros((samples, points))
    i_b = np.zeros((samples, points))
    pd_s = np.zeros((samples, points))
    ref_s = np.zeros((samples, points))
    v_s = np.zeros((samples, points))
    i_s = np.zeros((samples, points))

    trace_dir = path.join(top_dir, 'Scopes')
    
    # Loop over all acquired files
    for i in range(samples):
        background = np.loadtxt(path.join(trace_dir, 'Background%g.dat' % i),
                                delimiter='\t', skiprows=2)
        signal = np.loadtxt(path.join(trace_dir, 'Signal%g.dat' % i),
                            delimiter='\t', skiprows=2)
    
        pd_b[i,:] = background[:,0]
        ref_b[i,:] = background[:,1]
        v_b[i,:] = background[:,2]
        i_b[i,:] = background[:,3]
    
        pd_s[i,:] = signal[:,0]
        ref_s[i,:] = signal[:,1]
        v_s[i,:] = signal[:,2]
        i_s[i,:] = signal[:,3]
    
    # Pass to preprocessor for background subtraction and baseline shift
    (pd, v, i) = preprocess.data(pd_b, ref_b, v_b, i_b, pd_s, ref_s,
                                    v_s, i_s)
                                        
    ma2ghz = 0.8
    wavelength_range = np.linspace(mod_initial*ma2ghz, mod_final*ma2ghz,
                                   samples)
    t = dt * np.arange(samples)
    return pd, wavelength_range, t