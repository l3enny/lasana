import misc
from constants import *

class Transition(object):
    def __init__(self, state_i, state_f, A):
        self.A = A
        self.M = state_i.M
        
        self.Ei = state_i.E
        self.Ej = state_f.E
        self.dE = self.Ej - self.Ei
        
        self.f = self.dE/h
        self.l = misc.vac2air(c/self.f)
        
        if self.dE <= 0:
            raise ValueError("Transition is not spontaneous.")
                
        self.gi = state_i.g
        self.gj = state_f.g