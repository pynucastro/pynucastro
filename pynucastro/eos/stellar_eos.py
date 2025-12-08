"""Classes a methods for a full stellar equation of state (ions +
electrons + radiation)

"""


import numpy as np
from scipy.special import zeta

from pynucastro.constants import constants

from .electron_eos import ElectronEOS
from .eos_components import EOSState, IdealGasEOS, RadiationEOS


class StellarEOS:
    pass
