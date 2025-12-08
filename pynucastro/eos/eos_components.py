"""Classes a methods for a full stellar equation of state (ions +
electrons + radiation)

"""


import numpy as np
from scipy.special import zeta

from pynucastro.constants import constants


class EOSState:
    """A container to hold a thermodynamic state that depends on
    density, temperature, and composition

    """

    def __init__(self,
                 eta=0.0,
                 n=0.0, p=0.0, e=0.0,
                 dn_drho=0.0, dn_dT=0.0,
                 dp_drho=0.0, dp_dT=0.0,
                 de_drho=0.0, de_dT=0.0,
                 gamma1=0.0):

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

        self.gamma1 = gamma1


class IdealGas:
    """An idea gas equation of state for ions (and optionally
    electrons)."""

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
        EOSState

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

        gamma1 = 5./3.

        return EOSState(eta=0.0,
                        n=n, p=p, e=e,
                        dn_drho=dn_drho, dn_dT=dn_dT,
                        dp_drho=dp_drho, dp_dT=dp_dT,
                        de_drho=de_drho, de_dT=de_dT,
                        gamma1=gamma1)


class RadiationGas:
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
        EOSState

        """

        n = 2 * np.pi * zeta(3) * (2.0 * constants.k * T / (constants.h * constants.c_light))**3
        p = constants.a * T**4 / 3.0
        e = 3.0 * p / rho

        dn_drho = n / rho
        dn_dT = 0.0

        dp_drho = p / rho
        dp_dT = p / T

        de_drho = 0.0
        de_dT = e / T

        gamma1 = 4.0 / 3.0

        return EOSState(eta=0.0,
                        n=n, p=p, e=e,
                        dn_drho=dn_drho, dn_dT=dn_dT,
                        dp_drho=dp_drho, dp_dT=dp_dT,
                        de_drho=de_drho, de_dT=de_dT,
                        gamma1=gamma1)
