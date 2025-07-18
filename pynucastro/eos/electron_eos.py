from collections import namedtuple

import numpy as np
from scipy.optimize import brentq

from pynucastro.constants import constants

from .fermi_integrals import FermiIntegral

EOSState = namedtuple("EOSState", ["n_e", "p_e", "e_e",
                                   "n_pos", "p_pos", "e_pos",
                                   "eta"])


class ElectronEOS:
    """A electron EOS that works for arbitrary degeneracy or
    relativity.  This works by performing the Fermi-Dirac integrals
    directly.  This assumes complete ionization.

    Parameters
    ----------
    include_positrons : bool
        consider both positrons and electrons.

    """

    def __init__(self, include_positrons=True):
        self.include_positrons = include_positrons

    def pe_state(self, rho, T, comp, *,
                 eta_guess_min=-100, eta_guess_max=1.e7):
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
        inv_compton_wavelength = constants.m_e * constants.c_light / constants.h
        rest_mass = constants.m_e * constants.c_light**2

        beta = constants.k * T / rest_mass

        coeff = 8 * np.pi * np.sqrt(2) * inv_compton_wavelength**3 * beta**1.5

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

        # compute the number density, pressure and energy
        pcoeff = coeff * (2.0 / 3.0) * rest_mass * beta
        ecoeff = coeff * rest_mass * beta

        f12 = FermiIntegral(0.5, eta, beta)
        f12.evaluate(do_first_derivs=False, do_second_derivs=False)

        f32 = FermiIntegral(1.5, eta, beta)
        f32.evaluate(do_first_derivs=False, do_second_derivs=False)

        f52 = FermiIntegral(2.5, eta, beta)
        f52.evaluate(do_first_derivs=False, do_second_derivs=False)

        n_e = coeff * (f12.F + beta * f32.F)
        p_e = pcoeff * (f32.F + 0.5 * beta * f52.F)
        e_e = ecoeff * (f32.F + beta * f52.F) / rho

        n_pos = 0.0
        p_pos = 0.0
        e_pos = 0.0

        if self.include_positrons:
            f12_pos = FermiIntegral(0.5, eta_pos, beta)
            f12_pos.evaluate(do_first_derivs=False, do_second_derivs=False)

            f32_pos = FermiIntegral(1.5, eta_pos, beta)
            f32_pos.evaluate(do_first_derivs=False, do_second_derivs=False)

            f52_pos = FermiIntegral(2.5, eta_pos, beta)
            f52_pos.evaluate(do_first_derivs=False, do_second_derivs=False)

            n_pos = coeff * (f12_pos.F + beta * f32_pos.F)
            p_pos = pcoeff * (f32_pos.F + 0.5 * beta * f52_pos.F)
            e_pos = ecoeff * (f32_pos.F + beta * f52_pos.F) / rho

        return EOSState(eta=eta,
                        n_e=n_e, p_e=p_e, e_e=e_e,
                        n_pos=n_pos, p_pos=p_pos, e_pos=e_pos)
