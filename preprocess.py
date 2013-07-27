"""
Pre-process information which is usually experiment-specific. Customize
to your liking (just make sure it jives with parse.py if you're using
it).
"""

# Third party
import numpy as np

def transmission(signal, background, debug=False, **settings):
    """Generates a proper transmission curve from input data.
    
    Method subtracts out the plasma-induced emissions (unrelated to the
    laser) and corrects for the baseline shift resulting from increased
    diode power.
    """
    # Determine 100% transmission signal and average over acquisition time
    unabsorbed = np.mean(background[0] - background[1], axis=1)
    
    # Correct transmission signal for normal plasma emissions
    signal_total = signal[0]
    signal_plasma = signal[1]
    transmitted = signal_total - signal_plasma + settings['voffset']
    
    # Normalize transmission signal to full transmission value
    baseline_adjusted = transmitted/unabsorbed[:, np.newaxis]
    
    return baseline_adjusted

def transmission2(signal, background, debug=False, **settings):
    signal_total = signal[0]
    signal_plasma = signal[1]
    transmitted = signal_total - signal_plasma + settings['voffset']

    low = np.mean(signal_total[0, :])
    high = np.mean(signal_total[-1, :])
    q = (high - low)/settings['samples'] * np.arange(0, settings['samples']) + \
        + (low + high)/2
    
    #import matplotlib.pyplot as plt
    #plt.plot(transmitted[:, 0])
    #plt.hold(True)
    #plt.plot(q)
    #plt.legend(('Measured', 'Fit'))
    #plt.show()

    return transmitted/q[:, np.newaxis]
