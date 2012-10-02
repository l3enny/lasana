"""
Pre-process information which is usually experiment-specific. Customize
to your liking (just make sure it jives with parse.py if you're using
it).
"""

# Third party
import numpy as np

def data(pd_b, ref_b, v_b, i_b, pd_s, ref_s, v_s, i_s):
    pd = pd_s - pd_b # subtract out plasma-related background
    v = (v_b+v_s)/2 # both should be identical phenomena, average
    i = (i_b+i_s)/2 # both should be identical phenomena, average
    
    # Correct for power variations due to current modulation
    baseline_zero = np.mean(ref_b)   # Figure out zero intensity signal
    baseline_max = np.mean(ref_s[-1,:]) # Find max signal (assumes last trace
                                        # is the most intense!
    baseline_range = baseline_max - baseline_zero
    correction = baseline_range/np.mean(ref_s, axis=1)
    pd = pd * correction[:, np.newaxis]
    
    return pd, v, i