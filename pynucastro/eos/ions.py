import numpy as np
from interface import EosState
from scipy.constants import physical_constants

from pynucastro.networks.rate_collection import Composition
from pynucastro.nucdata.nucleus import Nucleus

""" In this module we introduce the corrections of the free ions under the ideal gas approximation.
This approximation is valid in the low density, high temperature regime."""

class Ions():

    def __init__(self):
        pass

    def update(self, state):

        #The following units are in the SI system.
        h,_,_ = physical_constants['Planck constant']
        amu,_,_ = physical_constants['atomic mass constant']
        Na,_,_ = physical_constants['Avogadro constant']
        kB,_,_ = physical_constants['Boltzmann constant']

        #Converting the previous units to CGS:
        h *= 1.0e7
        amu *= 1.0e3
        kB *= 1.0e7

        #Let us write a set of definitions that would help us:
        rho_i = 1.0 / state.rho
        T_i = 1.0 / state.T

        xni = Na / state.abar
        sconst = 2.0 * np.pi * amu * kB / (h * h)

        p = xni * kB * state.T * state.rho
        dpdr =  xni * kB * state.T
        dpdT = state.rho * xni * kB

        e = 1.5 * p * rho_i
        dedr = (1.5 * dpdr - e) * rho_i
        dedT = 1.5 * rho_i * dpdT

        x = state.abar * state.abar * np.sqrt(state.abar) * rho_i / Na
        y = sconst * state.T
        y = y*y*np.sqrt(y) * x

        s = xni * kB * np.log(y) + (e + p * rho_i) * T_i
        dsdr = (dedr - p * rho_i * rho_i + rho_i * dpdr)*T_i - xni *rho_i * kB
        dsdT = (dpdT * rho_i + dedT)*T_i - (e + p*rho_i)*T_i*T_i + 1.5 * xni * kB

        state.p += p
        state.dpdr += dpdr
        state.dpdT += dpdT

        state.e += e
        state.dedr += dedr
        state.dedT += dedT

        state.s += s
        state.dsdr += dsdr
        state.dsdT += dsdT

if __name__ == "__main__":

    nuc_list = ["p", "he4",
               "c12", "c13",
               "n13", "n14", "n15",
               "o15"]

    nuc_iter = [Nucleus(nuc) for nuc in nuc_list]

    composition = Composition(nuc_iter)
    composition.set_solar_like()

    state = EosState(rho=1.0e5, T=1.0e4, composition=composition)
    ions = Ions()
    ions.update(state)

    print(f"The pressure               : {state.p}")
    print(f"Entropy                    : {state.s}")
    print(f"Energy                     : {state.e}")
    print(f"Pressure density derivative: {state.dpdr}")
    print(f"Chemical potential         : {state.eta}")
    print(f"Number Fraction            : {state.xn}")
