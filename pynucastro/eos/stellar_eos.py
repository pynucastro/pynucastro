"""Classes a methods for a full stellar equation of state (ions +
electrons + radiation)

"""


from .electron_eos import ElectronEOS
from .eos_components import IdealGasEOS, RadiationEOS


class EOSState:
    """The full thermodynamic state of a multi-component plasma (ions,
    electrons, and radiation)

    """

    def __init__(self,
                 eta=0.0,
                 n_ele=0.0, n_pos=0.0,
                 p=0.0, e=0.0,
                 dp_drho=0.0, dp_dT=0.0,
                 de_drho=0.0, de_dT=0.0,
                 gamma1=0.0):

        self.eta = eta

        self.n_ele = n_ele
        self.n_pos = n_pos

        self.p = p
        self.e = e

        self.dp_drho = dp_drho
        self.dp_dT = dp_dT

        self.de_drho = de_drho
        self.de_dT = de_dT

        self.gamma1 = gamma1


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

        p = ion_state.p + rad_state.p + ele_state.p + pos_state.p
        e = ion_state.e + rad_state.e + ele_state.e + pos_state.e

        dp_drho = ion_state.dp_drho + rad_state.dp_drho + ele_state.dp_drho + pos_state.dp_drho
        dp_dT = ion_state.dp_dT + rad_state.dp_dT + ele_state.dp_dT + pos_state.dp_dT

        de_drho = ion_state.de_drho + rad_state.de_drho + ele_state.de_drho + pos_state.de_drho
        de_dT = ion_state.de_dT + rad_state.de_dT + ele_state.de_dT + pos_state.de_dT

        dT_drho_s = (p / rho**2 - de_drho) / (de_dT)
        gamma1 = rho / p * (dp_drho + dp_dT * dT_drho_s)

        return EOSState(eta=ele_state.eta,
                        n_ele=ele_state.n, n_pos=pos_state.n,
                        p=p, e=e,
                        dp_drho=dp_drho, dp_dT=dp_dT,
                        de_drho=de_drho, de_dT=de_dT,
                        gamma1=gamma1)
