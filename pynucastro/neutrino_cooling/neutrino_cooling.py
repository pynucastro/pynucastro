import matplotlib.pyplot as plt
import numpy as np
from sneut5 import sneut5


class NeutrinoCooling:
    """A class to provide an interface to explore neutrino cooling.
    This calls a specific implementation of the cooling that includes
    contributions from pairs, plasma, recombination, bremsstrahlung,
    and photoneutrinos."""

    def __init__(self):

        
