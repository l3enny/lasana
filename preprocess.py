"""
Pre-process information which is usually experiment-specific. Customize
to your liking (just make sure it jives with parse.py if you're using
it).
"""

# Included
import gui

# Third party
import numpy as np

def transmission(signal, reference, debug=False):
    """Generates a proper transmission curve from input data.
    
    Method subtracts out the plasma-induced emissions (unrelated to the
    laser) and corrects for the baseline shift resulting from increased
    diode power.
    """
    # Determine 100% transmission signal and average over acquisition time
    unabsorbed = np.mean(reference[0] - reference[1], axis=1)
    
    # Correct transmission signal for normal plasma emissions
    signal_total = signal[0]
    signal_plasma = signal[1]
    transmitted = signal_total - signal_plasma
    
    # Normalize transmission signal to full transmission value
    baseline_adjusted = transmitted/unabsorbed[:, np.newaxis]
    
    return baseline_adjusted