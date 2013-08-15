from scipy.constants import c, h

def vac2air(w):
    """ Converts wavelengths from vacuum to air.
    Input wavelength should be in meters.
    """
    w = w * 1e10
    return (w / (1.0 + 2.735182e-4 + 131.4182/w**2 +
            2.76249e8/w**4)) * 1e-10

class Transition(object):
    def __init__(self, state_i, state_f, A):
        self.A = A
        self.M = state_i.M
        
        self.Ei = state_i.E
        self.Ej = state_f.E
        self.dE = self.Ej - self.Ei
        
        self.f = self.dE/h
        self.l = vac2air(c/self.f)
        
        if self.dE <= 0:
            raise ValueError("Transition is not spontaneous.")
                
        self.gi = state_i.g
        self.gj = state_f.g
