"""
Pre-process information which is usually experiment-specific. Customize
to your liking (just make sure it jives with parse.py if you're using
it).
"""

# Third party
import numpy as np

def convert(pd_b, ref_b, v_b, i_b, pd_s, ref_s, v_s, i_s):
    pd = pd_s - pd_b # subtract out plasma-related background
    v = (v_b+v_s)/2 # both should be identical phenomena, average
    i = (i_b+i_s)/2 # both should be identical phenomena, average
    
    baseline_zero = np.mean(ref_b)   # Figure out zero intensity signal
    baseline_min, baseline_max = (np.mean(ref_s[0,:]), np.mean(ref_s[-1,:]))
    baseline_min = baseline_min - baseline_zero
    baseline_max = baseline_max - baseline_zero
    baseline_nmin = baseline_min/baseline_max
    correction = np.linspace(baseline_nmin, 1, pd.shape[0])

    pd = pd * np.flipud(correction[:, np.newaxis])
    
    return pd, v, i