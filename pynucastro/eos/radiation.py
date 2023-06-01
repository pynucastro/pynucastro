import numpy as np
from interface import EosState
from scipy.constants import physical_constants

from pynucastro.networks.rate_collection import Composition
from pynucastro.nucdata.nucleus import Nucleus

class Radiation():

    def __init__(self):
        pass

    def update(self, state):

        #Let start defining the constant that will participate in the
        #calculation of the readiation background.

        sconst, _, _ = physical_constants["Stefan-Boltzmann constant"]
        c, _, _ = physical_constants["speed of light in vacuum"]

        #Nowm let us change the units of our measurements to the cgs
        #system

        sconst *= 1.0e3
        c *= 1.0e2

        a = 4.0 * sconst / (3.0 * c)

        rho_i = 1.0 / state.rho
        T_i = 1.0 / state.T

        p = a * state.T * state.T * state.T * state.T
        dp_d = 0.0
        dp_t = 4.0 * p * T_i

        e = 3.0 * p * rho_i
        de_d = -3.0 * p * rho_i * rho_i
        de_t = 3.0 * dp_t * rho_i

        s = (e + p * rho_i) * T_i
        ds_d = (de_d + dp_d * rho_i - p * rho_i * rho_i) * T_i
        ds_t = (de_t + dp_t * rho_i - s) * T_i

        state.p += p
        state.dpdr += dp_d
        state.dpdT += dp_t

        state.s += s
        state.dsdr += ds_d
        state.dsdT += ds_t

        state.e += e
        state.dedr += de_d
        state.dedT += de_t


if __name__ == "__main__":

    nuc_list = ["p", "he4",
               "c12", "c13",
               "n13", "n14", "n15",
               "o15"]

    #nuc_iter = [Nucleus(nuc) for nuc in nuc_list]

    composition = Composition(nuc_list)
    composition.set_solar_like()

    state = EosState(rho=1.0e5, T=1.0e4, composition=composition)
    light = Radiation()
    light.update(state)

    print(f"Pressure                   : {state.p}")
    print(f"Entropy                    : {state.s}")
    print(f"Energy                     : {state.e}")
    print(f"dp_dr                      : {state.dpdr}")
    print(f"Chemical potential         : {state.eta}")
    print(f"Number Fraction            : {state.xn}")
