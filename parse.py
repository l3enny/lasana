"""
Processes the data from whatever format it was captured in and renders
it in a usable form to whoever's knocking.
"""

# Standard modules
import ConfigParser
import os
from os import path
import sys

# Third party
import numpy as N

# Part of package
from constants import *
import preprocess
    
def config(top_dir='.', debug=False):
    # Parse the settings file, check for proper formatting
    config = ConfigParser.RawConfigParser()
    read_return = config.read(path.join(top_dir, 'settings.cfg'))
    if read_return == []:
        raise ConfigParser.Error('No settings file detected.')
    
    samples = config.getint('Settings', 'Wavelength Samples')
    points = config.getint('Settings', 'Data Points')
    mod_initial = config.getfloat('Settings', 'Initial Modulation')
    mod_final = config.getfloat('Settings', 'Final Modulation')
    averages = config.getint('Settings', 'Averages')
    time_domain = config.getfloat('Settings', 'Time Domain')
    pressure = config.getfloat('Settings', 'Pressure')
    offset = config.getfloat('Settings', 'Offset')
    voffset = config.getfloat('Settings', 'Vertical Offset')
    return {'samples':samples,
            'points':points,
            'mod_initial':mod_initial,
            'mod_final':mod_final,
            'averages':averages,
            'dt':time_domain/points,
            'pressure':pressure,
            'offset':offset,
            'voffset':voffset}
    
def _load(dir, samples, points, debug=False):
    """Generate data arrays and load data from file.
    
    Private method for reading in data in the format determined by my exper-
    iment. Modify as necessary.
    """
    # Initialize storage arrays
    pd_r = N.zeros((samples, points))
    ref_r = N.zeros((samples, points))
    v_r = N.zeros((samples, points))
    i_r = N.zeros((samples, points))
    pd_s = N.zeros((samples, points))
    ref_s = N.zeros((samples, points))
    v_s = N.zeros((samples, points))
    i_s = N.zeros((samples, points))
    
    # Loop over all acquired files
    for i in range(samples):
        reference = N.loadtxt(path.join(dir, 'Reference%g.dat' % i),
                              delimiter='\t', usecols=(1,2,3,4))
        signal = N.loadtxt(path.join(dir, 'Signal%g.dat' % i),
                           delimiter='\t', usecols=(1,2,3,4))

        pd_r[i,:] = reference[:,0]
        ref_r[i,:] = reference[:,1]
        v_r[i,:] = reference[:,2]
        i_r[i,:] = reference[:,3]
    
        pd_s[i,:] = signal[:,0]
        ref_s[i,:] = signal[:,1]
        v_s[i,:] = signal[:,2]
        i_s[i,:] = signal[:,3]
    
    return (pd_s, pd_r), (ref_s, ref_r)#, (v_s, v_b), (i_s, i_b)
    
def data(dir='.', debug=False, **settings):
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
   
    trace_dir = path.join(dir)
    signal, reference = _load(trace_dir, settings['samples'],
                              settings['points'])
    time = settings['dt'] * N.arange(settings['points'])
    freq = N.linspace(settings['mod_initial'], settings['mod_final'],
                      settings['samples']) * ma2hz + settings['offset']
    return signal, time, freq