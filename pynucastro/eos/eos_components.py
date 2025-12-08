"""Classes a methods for a full stellar equation of state (ions +
electrons + radiation)

"""


import numpy as np
from scipy.special import zeta

from pynucastro.constants import constants


class EOSComponentState:
    """A container to hold a thermodynamic state for a single
    component of a full stellar plasma (i.e., ions, electrons, or
    radiation) that depends on density, temperature, and composition

    """

    def __init__(self,
                 eta=0.0,
                 n=0.0, p=0.0, e=0.0,
                 dn_drho=0.0, dn_dT=0.0,
                 dp_drho=0.0, dp_dT=0.0,
                 de_drho=0.0, de_dT=0.0):

        self.eta = eta

        self.n = n
        self.p = p
        self.e = e

        self.dn_drho = dn_drho
        self.dn_dT = dn_dT

        self.dp_drho = dp_drho
        self.dp_dT = dp_dT

        self.de_drho = de_drho
        self.de_dT = de_dT


class IdealGasEOS:
    """An idea gas equation of state for ions (and optionally
    electrons).

    Parameters
    ----------
    include_electrons : bool
        do we include electrons together with the ions, assuming
        full ionization?

    """

    def __init__(self, include_electrons=False):
        self.include_electrons = include_electrons

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
        EOSComponentState

        """

        mu = None
        abar = comp.abar

        if not self.include_electrons:
            mu = abar
        else:
            mu = 1.0 / (1.0 / abar + comp.Ye)

        n = rho / (mu * constants.m_u)
        p = rho * constants.k * T / (mu * constants.m_u)
        e = 1.5 * p / rho

        dn_drho = n / rho
        dn_dT = 0.0

        dp_drho = p / rho
        dp_dT = p / T

        de_drho = 0.0
        de_dT = e / T

        return EOSComponentState(eta=0.0,
                                 n=n, p=p, e=e,
                                 dn_drho=dn_drho, dn_dT=dn_dT,
                                 dp_drho=dp_drho, dp_dT=dp_dT,
                                 de_drho=de_drho, de_dT=de_dT)


class RadiationEOS:
    """A blackbody radiation equation of state"""

    def __init__(self):
        pass

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
        EOSComponentState

        """

        # composition is not used
        _ = comp

        n_const = 2 * np.pi * zeta(3) * (2.0 * constants.k / (constants.h * constants.c_light))**3
        n = n_const * T**3

        p = constants.a * T**4 / 3.0
        e = 3.0 * p / rho

        dn_drho = 0.0
        dn_dT = 3.0 * n_const * T**2

        dp_drho = 0.0
        dp_dT = 4.0 * p / T

        de_drho = -e / rho
        de_dT = 4.0 * e / T

        return EOSComponentState(eta=0.0,
                                 n=n, p=p, e=e,
                                 dn_drho=dn_drho, dn_dT=dn_dT,
                                 dp_drho=dp_drho, dp_dT=dp_dT,
                                 de_drho=de_drho, de_dT=de_dT)
