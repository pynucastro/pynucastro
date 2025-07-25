"""Classes and methods for managing an electron / positron equation of
state.

"""


from collections import namedtuple

import numpy as np
from scipy.optimize import brentq

from pynucastro.constants import constants

from .degeneracy_parameter_bounds import get_eta_bounds
from .fermi_integrals import FermiIntegral

EOSState = namedtuple("EOSState", ["n_e", "n_pos",
                                   "p_e", "p_pos",
                                   "e_e", "e_pos",
                                   "eta",
                                   "dne_drho", "dne_dT",
                                   "dnp_drho", "dnp_dT",
                                   "dpe_drho", "dpe_dT",
                                   "dpp_drho", "dpp_dT",
                                   "dee_drho", "dee_dT",
                                   "dep_drho", "dep_dT"])


class ElectronEOS:
    """An electron/positron EOS that works for arbitrary degeneracy or
    relativity.  This works by performing the Fermi-Dirac integrals
    directly.  This assumes complete ionization.

    Parameters
    ----------
    include_positrons : bool
        consider both positrons and electrons.

    """

    def __init__(self, include_positrons=True):
        self.include_positrons = include_positrons

    def pe_state(self, rho=None, T=None, comp=None, *,
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

        n_e_net = (zbar / abar) * constants.N_A * rho

        # our Fermi integrals will use a dimensionless temperature
        inv_compton_wavelength = constants.m_e * constants.c_light / constants.h
        rest_mass = constants.m_e * constants.c_light**2

        beta = constants.k * T / rest_mass

        coeff = 8 * np.pi * np.sqrt(2) * inv_compton_wavelength**3

        def n_e_fermi(eta):
            f12 = FermiIntegral(0.5, eta, beta)
            f12.evaluate(do_first_derivs=False, do_second_derivs=False)

            f32 = FermiIntegral(1.5, eta, beta)
            f32.evaluate(do_first_derivs=False, do_second_derivs=False)

            return coeff * beta**1.5 * (f12.F + beta * f32.F)

        def n_pos_fermi(eta):
            eta_pos = -eta - 2.0/beta

            f12 = FermiIntegral(0.5, eta_pos, beta)
            f12.evaluate(do_first_derivs=False, do_second_derivs=False)

            f32 = FermiIntegral(1.5, eta_pos, beta)
            f32.evaluate(do_first_derivs=False, do_second_derivs=False)

            return coeff * beta**1.5 * (f12.F + beta * f32.F)

        # compute the degeneracy parameter
        root_find_limits = get_eta_bounds(rho * comp.ye, T,
                                          include_positrons=self.include_positrons)
        if self.include_positrons:
            try:
                eta = brentq(lambda eta: n_e_net - (n_e_fermi(eta) - n_pos_fermi(eta)),
                             root_find_limits[0], root_find_limits[1], xtol=1.e-15)
            except ValueError:
                eta = brentq(lambda eta: n_e_net - (n_e_fermi(eta) - n_pos_fermi(eta)),
                             eta_guess_min, eta_guess_max, xtol=1.e-15)

        else:
            try:
                eta = brentq(lambda eta: n_e_net - n_e_fermi(eta),
                             root_find_limits[0], root_find_limits[1], xtol=1.e-15)
            except ValueError:
                eta = brentq(lambda eta: n_e_net - n_e_fermi(eta),
                             eta_guess_min, eta_guess_max, xtol=1.e-15)

        # for positrons
        eta_pos = -eta - 2.0/beta

        # compute the number density, pressure and energy
        pcoeff = coeff * (2.0 / 3.0) * rest_mass
        ecoeff = coeff * rest_mass

        f12 = FermiIntegral(0.5, eta, beta)
        f12.evaluate(do_first_derivs=True, do_second_derivs=False)

        f32 = FermiIntegral(1.5, eta, beta)
        f32.evaluate(do_first_derivs=True, do_second_derivs=False)

        f52 = FermiIntegral(2.5, eta, beta)
        f52.evaluate(do_first_derivs=True, do_second_derivs=False)

        n_e = coeff * beta**1.5 * (f12.F + beta * f32.F)
        p_e = pcoeff * beta**2.5 * (f32.F + 0.5 * beta * f52.F)

        # this is (energy / volume, what we usually write as rho e)
        E_e = ecoeff * beta**2.5 * (f32.F + beta * f52.F)

        n_pos = 0.0
        p_pos = 0.0
        E_pos = 0.0
        if self.include_positrons:
            f12_pos = FermiIntegral(0.5, eta_pos, beta)
            f12_pos.evaluate(do_first_derivs=True, do_second_derivs=False)

            f32_pos = FermiIntegral(1.5, eta_pos, beta)
            f32_pos.evaluate(do_first_derivs=True, do_second_derivs=False)

            f52_pos = FermiIntegral(2.5, eta_pos, beta)
            f52_pos.evaluate(do_first_derivs=True, do_second_derivs=False)

            n_pos = coeff * beta**1.5 * (f12_pos.F + beta * f32_pos.F)
            p_pos = pcoeff * beta**2.5 * (f32_pos.F + 0.5 * beta * f52_pos.F)
            E_pos = ecoeff * beta**2.5 * (f32_pos.F + beta * f52_pos.F) + 2 * rest_mass * n_pos

        # compute the derivatives of eta and beta with respect to
        # density and temperature
        dne_deta = coeff * beta**1.5 * (f12.dF_deta + beta * f32.dF_deta)
        dne_dbeta = 0.5 * coeff * np.sqrt(beta) * (3.0 * f12.F + 5.0 * beta * f32.F +
                                                   2 * beta * (f12.dF_dbeta + beta * f32.dF_dbeta))

        dnp_deta = 0.0
        dnp_dbeta = 0.0
        if self.include_positrons:
            dnp_deta = coeff * beta**1.5 * (-f12_pos.dF_deta - beta * f32_pos.dF_deta)
            dnp_dbeta = 0.5 * coeff / np.sqrt(beta) * (beta * (3.0 * f12_pos.F + 5.0 * beta * f32_pos.F) +
                                                       2.0 * beta**2 * (f12_pos.dF_dbeta +
                                                                        beta * f32_pos.dF_dbeta) +
                                                       4.0 * (f12_pos.dF_deta + beta * f32_pos.dF_deta))

        #dbeta_drho = 0.0
        dbeta_dT = constants.k / rest_mass

        deta_drho = constants.N_A * (zbar / abar) / (dne_deta - dnp_deta)
        deta_dT = -dbeta_dT * (dne_dbeta - dnp_dbeta) / (dne_deta - dnp_deta)

        # Compute partials of number density with density and temperature
        # For debugging
        dne_drho = dne_deta * deta_drho
        dne_dT = dne_deta * deta_dT + dne_dbeta * dbeta_dT

        dnp_drho = 0.0
        dnp_dT = 0.0
        if self.include_positrons:
            dnp_drho = dnp_deta * deta_drho
            dnp_dT = dnp_deta * deta_dT + dnp_dbeta * dbeta_dT

        # Compute partials of pressure with density and temperature via the chain rule
        dpe_deta = pcoeff * beta**2.5 * (f32.dF_deta + 0.5 * beta * f52.dF_deta)
        dpe_dbeta = 0.25 * pcoeff * beta**1.5 * (10.0 * f32.F + 7.0 * beta * f52.F +
                                                 4.0 * beta * (f32.dF_dbeta + 0.5 * beta * f52.dF_dbeta))

        dpe_drho = dpe_deta * deta_drho
        dpe_dT = dpe_deta * deta_dT + dpe_dbeta * dbeta_dT

        dpp_drho = 0.0
        dpp_dT = 0.0
        if self.include_positrons:
            dpp_deta = -pcoeff * beta**2.5 * (f32_pos.dF_deta + 0.5 * beta * f52_pos.dF_deta)
            dpp_dbeta = pcoeff * np.sqrt(beta) * (beta * (2.5 * f32_pos.F + 1.75 * beta * f52_pos.F) +
                                                  beta**2 * (f32_pos.dF_dbeta + 0.5 * beta * f52_pos.dF_dbeta) +
                                                  2.0 * (f32_pos.dF_deta + 0.5 * beta * f52_pos.dF_deta))

            dpp_drho = dpp_deta * deta_drho
            dpp_dT = dpp_deta * deta_dT + dpp_dbeta * dbeta_dT

        # Compute partials of energy with density and temperature
        dEe_deta = ecoeff * beta**2.5 * (f32.dF_deta + beta * f52.dF_deta)
        dEe_dbeta = 0.5 * ecoeff * beta**1.5 * (5.0 * f32.F + 7.0 * beta * f52.F +
                                                2.0 * beta * (f32.dF_dbeta + beta * f52.dF_dbeta))

        dEe_drho = dEe_deta * deta_drho   # dbeta_drho = 0
        dEe_dT = dEe_deta * deta_dT + dEe_dbeta * dbeta_dT

        dEp_drho = 0.0
        dEp_dT = 0.0
        if self.include_positrons:
            # these don't include the 2 m_e c**2 n_p part
            dEp_deta = -ecoeff * beta**2.5 * (f32_pos.dF_deta + beta * f52_pos.dF_deta)
            dEp_dbeta = ecoeff * np.sqrt(beta) * (beta * (2.5 * f32_pos.F + 3.5 * beta * f52_pos.F) +
                                                  beta**2 * (f32_pos.dF_dbeta + beta * f52_pos.dF_dbeta) +
                                                  2.0 * (f32_pos.dF_deta + beta * f52_pos.dF_deta))

            dEp_drho = dEp_deta * deta_drho + 2.0 * rest_mass * dnp_drho
            dEp_dT = dEp_deta * deta_dT + dEp_dbeta * dbeta_dT + 2.0 * rest_mass * dnp_dT

        # correct energy to be specific energy
        e_e = E_e / rho
        dee_drho = (dEe_drho - E_e/rho) / rho
        dee_dT = dEe_dT / rho

        e_pos = 0.0
        dep_drho = 0.0
        dep_dT = 0.0
        if self.include_positrons:
            e_pos = E_pos / rho
            dep_drho = (dEp_drho - E_pos/rho) / rho
            dep_dT = dEp_dT / rho

        return EOSState(eta=eta,
                        n_e=n_e, p_e=p_e, e_e=e_e,
                        n_pos=n_pos, p_pos=p_pos, e_pos=e_pos,
                        dne_drho=dne_drho, dne_dT=dne_dT,
                        dnp_drho=dnp_drho, dnp_dT=dnp_dT,
                        dpe_drho=dpe_drho, dpe_dT=dpe_dT,
                        dpp_drho=dpp_drho, dpp_dT=dpp_dT,
                        dee_drho=dee_drho, dee_dT=dee_dT,
                        dep_drho=dep_drho, dep_dT=dep_dT)
