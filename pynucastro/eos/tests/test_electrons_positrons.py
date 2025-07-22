# unit tests for rates
import numpy as np
from pytest import approx

from pynucastro.eos.difference_utils import fourth_order_rho, fourth_order_temp

from pynucastro import Composition
from pynucastro.constants import constants
from pynucastro.eos import ElectronEOS


def zero_temperature_eos(rho, comp):
    """This is the analytic expression for the T = 0 Fermi-Dirac
    electron gas.  See, e.g, Hansen, Kawaler, and Trimble Eq. 3.53 and
    following.

    """

    rest_mass = constants.m_e * constants.c_light**2
    inv_compton_wavelength = constants.m_e * constants.c_light / constants.h
    A = np.pi / 3.0 * inv_compton_wavelength**3 * rest_mass
    B = (8.0 * np.pi / 3.0) / constants.N_A * inv_compton_wavelength**3

    x = np.cbrt(rho * comp.ye / B)

    f = (x * (2*x**2 - 3) * np.sqrt(1.0 + x**2) + 3.0 * np.asinh(x))
    g = 8 * x**3 * (np.sqrt(1.0 + x**2) - 1.0) - f

    p = A * f
    e = A * g

    return p, e / rho


def ideal_gas(n_e, n_pos, T):
    """An electron ideal gas---the pressure is the same regardless if
    we are relativistic or non-relativistic."""

    p = (n_e + n_pos) * constants.k * T

    return p


class TestElectronPositronEOS:

    def test_low_temperature(self):
        # at low temperatures, we should be completely degenerate and
        # follow the zero-temperature expressions approximately.
        # we also don't expect there to be many positrons.

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        T = 1.e5

        for rho in [1.e4, 1.e5, 1.e7, 1.e9]:
            es = e.pe_state(rho, T, comp)
            p_zt, e_zt = zero_temperature_eos(rho, comp)

            assert es.p_e == approx(p_zt, rel=1.e-4)
            assert es.e_e == approx(e_zt, rel=1.e-4)

    def test_high_temperature(self):
        # at high temperatures, we should be an ideal electron gas.

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # we need to pick a set of thermodynamic conditions such that
        # eta << -1, in order for the ideal gas approximation to be
        # valid.  If we make T too high, eta actually increases
        # because of pair production (number density of fermions
        # increases)
        T = 1.e9

        for rho in [1.e-3, 1.e-2, 1.e-1]:
            es = e.pe_state(rho, T, comp)
            p_ideal = ideal_gas(es.n_e, es.n_pos, T)
            assert es.p_e + es.p_pos == approx(p_ideal, rel=1.e-3)

    def test_positrons(self):
        # at low densities and high temperatures (T ~ 1.e10 K), we
        # should have a lot of positrons

        e = ElectronEOS(include_positrons=True)
        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # too cold to have positrons (kT << 2 m_e c**2)

        rho = 1.e4
        T = 1.e6

        es = e.pe_state(rho, T, comp)

        assert es.n_pos == 0.0

        # positrons and electrons are almost equal

        rho = 1.e4
        T = 1.e10

        es = e.pe_state(rho, T, comp)

        assert es.n_pos == approx(es.n_e, rel=1.e-3)

        # For non-degenerate and kT < m_e c**2, we can get a Saha-like
        # relation (see Clayton Eq. 3-297).  This is *very*
        # approximate, so let's just check order of magnitude

        rho = 1.e3
        T = 1.e9

        n_e_0 = (comp.zbar / comp.abar) * constants.N_A * rho

        beta = constants.k * T / (constants.m_e * constants.c_light**2)

        n_1 = 1.0 / np.sqrt(2) * (constants.m_e * constants.k * T /
                                  np.pi / constants.hbar**2)**1.5 * np.exp(-1.0 / beta)

        n_pos_approx = -0.5 * n_e_0 + np.sqrt((0.5 * n_e_0)**2 + n_1**2)

        es = e.pe_state(rho, T, comp)

        assert es.n_pos == approx(n_pos_approx, rel=0.5)

    def test_ne_rho_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_rho = 1.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                drho = eps_rho * rho
                deriv = fourth_order_rho(e, (rho, T, comp), "n_e", drho)
                assert es.dne_drho == approx(deriv, rel=5.e-5)

    def test_ne_temp_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # use a relatively large eps because of how weakly beta enters when beta << 1

        eps = 1.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                es_r = e.pe_state(rho, T * (1.0 + eps), comp)

                dnedr_approx = (es_r.n_e - es.n_e) / (eps * T)
                if es.dne_dT == 0:
                    # no pair production
                    scale = es.n_e / T
                    assert es.dne_dT == approx(dnedr_approx / scale, abs=5.e-9)
                else:
                    assert es.dne_dT == approx(dnedr_approx, rel=1.e-4)
