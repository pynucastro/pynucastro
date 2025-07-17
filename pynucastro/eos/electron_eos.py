from collections import namedtuple

import numpy as np
from scipy.optimize import brentq

from pynucastro.constants import constants

from .fermi_integrals import FermiIntegral

EOSState = namedtuple("EOSState", ["p_e", "e_e", "p_pos", "e_pos", "eta"])


class ElectronEOS:
    """A electron EOS that works for arbitrary degeneracy or
    relativity.  This works by performing the Fermi-Dirac integrals
    directly.  This assumes complete ionization.

    Parameters
    ----------
    include_positrons : bool
        consider both positrons and electrons.

    """

    def __init__(self, include_positrons=False):
        self.include_positrons = include_positrons

    def pe_state(self, rho, T, comp, *,
                 eta_guess_min=-500, eta_guess_max=1.e7):
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

        # compute the number density of electrons
        zbar = comp.zbar
        abar = comp.abar

        n_e = (zbar / abar) * constants.N_A * rho

        # our Fermi integrals will use a dimensionless temperature
        beta = constants.k * T / (constants.m_e * constants.c_light**2)

        coeff = 8 * np.pi * np.sqrt(2) * (constants.m_e * constants.c_light / constants.h)**3 * beta**1.5

        def n_e_fermi(eta):
            f12 = FermiIntegral(0.5, eta, beta)
            f12.evaluate(do_first_derivs=False, do_second_derivs=False)

            f32 = FermiIntegral(1.5, eta, beta)
            f32.evaluate(do_first_derivs=False, do_second_derivs=False)

            return coeff * (f12.F + beta * f32.F)

        def n_pos_fermi(eta):
            eta_pos = -eta - 2.0/beta

            f12 = FermiIntegral(0.5, eta_pos, beta)
            f12.evaluate(do_first_derivs=False, do_second_derivs=False)

            f32 = FermiIntegral(1.5, eta_pos, beta)
            f32.evaluate(do_first_derivs=False, do_second_derivs=False)

            return coeff * (f12.F + beta * f32.F)

        # compute the degeneracy parameter
        if self.include_positrons:
            eta = brentq(lambda eta: n_e - (n_e_fermi(eta) - n_pos_fermi(eta)),
                         eta_guess_min, eta_guess_max)
        else:
            eta = brentq(lambda eta: n_e - n_e_fermi(eta),
                         eta_guess_min, eta_guess_max)

        # for positrons
        eta_pos = -eta - 2.0/beta

        # compute the pressure and energy
        pcoeff = coeff * (2.0 / 3.0) * constants.m_e * constants.c_light**2 * beta
        ecoeff = coeff * constants.m_e * constants.c_light**2 * beta

        f32 = FermiIntegral(1.5, eta, beta)
        f32.evaluate(do_first_derivs=False, do_second_derivs=False)

        f52 = FermiIntegral(2.5, eta, beta)
        f52.evaluate(do_first_derivs=False, do_second_derivs=False)

        p_e = pcoeff * (f32.F + 0.5 * beta * f52.F)
        e_e = ecoeff * (f32.F + beta * f52.F) / rho

        p_pos = 0.0
        e_pos = 0.0

        if self.include_positrons:
            f32_pos = FermiIntegral(1.5, eta_pos, beta)
            f32_pos.evaluate(do_first_derivs=False, do_second_derivs=False)

            f52_pos = FermiIntegral(2.5, eta_pos, beta)
            f52_pos.evaluate(do_first_derivs=False, do_second_derivs=False)

            p_pos = pcoeff * (f32_pos.F + 0.5 * beta * f52_pos.F)
            e_pos = ecoeff * (f32_pos.F + beta * f52_pos.F) / rho

        return EOSState(eta=eta, p_e=p_e, e_e=e_e, p_pos=p_pos, e_pos=e_pos)
