import math as m
import numpy as N
from constants import *
import lineshapes

def bimodal_voigt(transitions, p):
    def func(x, T, amp, drift):
        # Pressure broadening/lorentzian part of profile
        fwhm_a = p * torr2hz
        gamma = fwhm_a/2
        f = 0
        V = lineshapes.voigt
        # TODO: Manual origin setting, assumes the 0.0 frequency shift
        # is the D0 peak. Should be a better way to do this.
        origin = transitions[0].f
        for t in transitions:
            # Doppler broadening/gaussian part of the profile
            fwhm_d = N.sqrt((8*m.log(2)) * kB*T / (t.M*c**2)) * t.f
            sigma = fwhm_d/(2*m.sqrt(2*m.log(2)))
            temp = 0.5 * (V(x, sigma, gamma, t.f - origin + drift)
                          + V(x, sigma, gamma, t.f - origin - drift))
            temp = temp * t.l**2 * t.A * (t.gj/t.gi) / (8 * N.pi)
            f += temp
        return N.exp(- amp * f)
    return func

def voigt(transitions, p):
    def func(x, T, amp):
        # Pressure broadening/lorentzian part of profile
        fwhm_a = p * torr2hz
        gamma = fwhm_a/2
        f = 0
        V = lineshapes.voigt
        # TODO: Manual origin setting, assumes the 0.0 frequency shift
        # is the D0 peak. Should be a better way to do this.
        origin = transitions[0].f
        for t in transitions:
            # Doppler broadening/gaussian part of the profile
            fwhm_d = N.sqrt((8*m.log(2)) * kB*T / (t.M*c**2)) * t.f
            sigma = fwhm_d/(2*m.sqrt(2*m.log(2)))
            temp = V(x, sigma, gamma, t.f - origin)
            temp = temp * t.l**2 * t.A * (t.gj/t.gi) / (8 * N.pi)
            f += temp
        return N.exp(- amp * f)
    return func

def gaussian(transitions):
    def func(x, T, amp):
        f = 0
        G = lineshapes.gaussian
        # TODO: Manual origin setting, assumes the 0.0 frequency shift
        # is the D0 peak. Should be a better way to do this.
        origin = transitions[0].f
        for t in transitions:
            # Doppler broadening/gaussian part of the profile
            fwhm_d = N.sqrt((8*m.log(2)) * kB*T / (t.M*c**2)) * t.f
            sigma = fwhm_d/(2*m.sqrt(2*m.log(2)))
            temp = G(x, sigma, t.f - origin)
            temp = temp * t.A * t.l**2 * (t.gj/t.gi) / (8*N.pi)
            f += temp
        return N.exp(- amp * f)
    return func