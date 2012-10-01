"""
Processes the data from whatever format it was captured in and renders
it in a usable form to whoever's knocking.
"""

# Standard modules
import Tkinter, tkFileDialog
import os
from os import path

# Third party
import numpy as np

# Part of package
import preprocess


def pickdir():
    root = Tkinter.Tk()
    root.withdraw()
    top_dir = tkFileDialog.askdirectory(parent=root, initialdir="..", 
                                        title='Pick a directory')
    if top_dir is '':
        print 'No directory provided, quitting...'
        sys.exit(0)
    return top_dir
    

def grab(top_dir):
    trace_dir = path.join(top_dir, 'Scopes')
    
    samples = len(os.listdir(trace_dir))/2
    points = np.loadtxt(path.join(trace_dir, 'Background0.dat'),
                        delimiter='\t', skiprows=2).shape[0]
    
    pd_b = np.zeros((samples, points))
    ref_b = np.zeros((samples, points))
    v_b = np.zeros((samples, points))
    i_b = np.zeros((samples, points))
    pd_s = np.zeros((samples, points))
    ref_s = np.zeros((samples, points))
    v_s = np.zeros((samples, points))
    i_s = np.zeros((samples, points))

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
        
    (pd, v, i) = preprocess.convert(pd_b, ref_b, v_b, i_b, pd_s, ref_s,
                                    v_s, i_s)
                                        
    return pd, v, i