import numpy as np
from interface import EosState
from scipy.constants import physical_constants

from pynucastro.networks.rate_collection import Composition
from pynucastro.nucdata.nucleus import Nucleus

class Coulomb():

    def __init__(self):
        pass

    def update(self, state):

        return 0

if __name__ == "__main__":

    nuc_list = ["p", "he4",
               "c12", "c13",
               "n13", "n14", "n15",
               "o15"]

    #nuc_iter = [Nucleus(nuc) for nuc in nuc_list]

    composition = Composition(nuc_list)
    composition.set_solar_like()

    state = EosState(rho=1.0e5, T=1.0e4, composition=composition)
    coulomb = Coulomb()
    coulomb.update(state)

    print(f"Pressure                   : {state.p}")
    print(f"Entropy                    : {state.s}")
    print(f"Energy                     : {state.e}")
    print(f"dp_dr                      : {state.dpdr}")
    print(f"Chemical potential         : {state.eta}")
    print(f"Number Fraction            : {state.xn}")
