def vac2air(w):
    """ Converts wavelengths from vacuum to air.
    
    Input wavelength should be in meters.
    """
    w = w * 1e10
    return (w / (1.0 + 2.735182e-4 + 131.4182/w**2 + 2.76249e8/w**4)) * 1e-10