"""Classes a methods for a full stellar equation of state (ions +
electrons + radiation)

"""


from pynucastro.constants import constants

from .electron_eos import ElectronEOS


class EOSState:
    """A container to hold a thermodynamic state that depends on
    density, temperature, and composition

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


class IdealGas:
    """An idea gas equation of state for ions (and optionally
    electrons)."""

    def __init__(self, include_electrons=False):
        self.include_electrons = include_electrons

    def pe_state(self, rho=None, T=None, comp=None, *,
                 compute_derivs=True):
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

        return EOSState(eta=0.0,
                        n=n, p=p, e=e,
                        dn_drho=dn_drho, dn_dT=dn_dT,
                        dp_drho=dp_drho, dp_dT=dp_dT,
                        de_drho=de_drho, de_dT=de_dT)
