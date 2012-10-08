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
import numpy as np

# Part of package
import preprocess
    
def config(top_dir='.', debug=False):
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
        samples = config.getint('Settings', 'Wavelength Samples')
        points = config.getint('Settings', 'Data Points')
        mod_initial = config.getfloat('Settings', 'Initial Modulation')
        mod_final = config.getfloat('Settings', 'Final Modulation')
        averages = config.getint('Settings', 'Averages')
        time_domain = config.getfloat('Settings', 'Time Domain')
        pressure = config.getfloat('Settings', 'Pressure')
        offset = config.getfloat('Settings', 'Offset')
    except ConfigParser.NoOptionError:
        print "Settings file is not properly formatted, quitting."
        sys.exit(1)
    return {'samples':samples, 'points':points, 'mod_initial':mod_initial,
            'mod_final':mod_final, 'averages':averages, 'dt':time_domain/points, 'pressure':pressure, 'offset':offset}
    
def _load(dir, samples, points, debug=False):
    """Generate data arrays and load data from file.
    
    Private method for reading in data in the format determined by my exper-
    iment. Modify as necessary.
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
    
    # Loop over all acquired files
    for i in range(samples):
        background = np.loadtxt(path.join(dir, 'Background%g.dat' % i),
                                delimiter='\t', skiprows=2)
        signal = np.loadtxt(path.join(dir, 'Signal%g.dat' % i),
                            delimiter='\t', skiprows=2)
    
        pd_b[i,:] = background[:,0]
        ref_b[i,:] = background[:,1]
        v_b[i,:] = background[:,2]
        i_b[i,:] = background[:,3]
    
        pd_s[i,:] = signal[:,0]
        ref_s[i,:] = signal[:,1]
        v_s[i,:] = signal[:,2]
        i_s[i,:] = signal[:,3]
    
    return (pd_s, pd_b), (ref_s, ref_b)#, (v_s, v_b), (i_s, i_b)
    
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
   
    trace_dir = path.join(dir, 'Scopes')
    pd, ref = _load(trace_dir, settings['samples'], settings['points'])
    t = settings['dt'] * np.arange(settings['points'])
    return pd, t