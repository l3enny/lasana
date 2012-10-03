"""
Pre-process information which is usually experiment-specific. Customize
to your liking (just make sure it jives with parse.py if you're using
it).
"""

# Included
import gui

# Third party
import numpy as np

def transmission(signal, reference): 
    unabsorbed = np.mean(reference[0] - reference[1], axis=1)
    
    signal_total = signal[0]
    signal_plasma = signal[1]
    signal_corrected = signal_total - signal_plasma
    
    baseline_adjusted = signal_corrected/unabsorbed[:, np.newaxis]
    
    return baseline_adjusted