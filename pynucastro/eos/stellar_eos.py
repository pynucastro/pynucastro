"""Classes a methods for a full stellar equation of state (ions +
electrons + radiation)

"""


import numpy as np
from scipy.special import zeta

from pynucastro.constants import constants

from .electron_eos import ElectronEOS
from .eos_components import EOSState, IdealGasEOS, RadiationEOS


class StellarEOS:
    """A full stellar equation of state with contributions from ions,
    electrons, and radiation.

    Parameters
    ----------
    electrons_are_degenerate : bool
        Do we treat electrons as a Fermi gas with arbitrary degeneracy
        and relativity? or just treat them as an ideal gas?
    include_positrons : bool
        If electrons are Fermi gas, do we consider both positrons and
        electrons?

    """

    def __init__(self, electrons_are_degenerate=True,
                 include_positrons=True):
        self.electrons_are_degenerate = electrons_are_degenerate
        self.include_positrons = include_positrons

    def pe_state(self, rho=None, T=None, comp=None):
        """Find the pressure and energy given density, temperature,
        and composition

        Parameters
        ----------
        rho : float
            Density (g/cm**3)
        T : float
            Temperature (K)
        comp : Composition
            Composition (abundances of each nucleus)

        Returns
        -------
        EOSState

        """

        ion_eos = None
        if self.electrons_are_degenerate:
            ion_eos = IdealGasEOS()
        else:
            ion_eos = IdealGasEOS(include_electrons=True)

        rad_eos = RadiationEOS()

        ele_eos = None
        if self.electrons_are_degenerate:
            ele_eos = ElectronEOS(include_positrons=self.include_positrons)

        # evaluate the full state
        ion_state = ion_eos.pe_state(rho, T, comp)
        rad_state = rad_eos.pe_state(rho, T, comp)
        ele_state, pos_state = ele_eos.pe_state(rho, T, comp)

_
