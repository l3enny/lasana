from scipy.constants import e, atomic_mass

torr2hz = 25.6e6 # pressure broadening coefficient

class II3S1(object):
    M = 4.002602 * atomic_mass
    E = 19.81961363388 * e
    J = 1.
    g = 2 * J + 1
            
class II3P2(object):
    M = 4.002602 * atomic_mass
    E = 20.96408592885 * e
    J = 2.
    g = 2 * J + 1
    
class II3P1(object):
    M = 4.002602 * atomic_mass
    E = 20.96409540439 * e
    J = 1.
    g = 2 * J + 1
    
class II3P0(object):
    M = 4.002602 * atomic_mass
    E = 20.96421789026 * e
    J = 0.
    g = 2 * J + 1
