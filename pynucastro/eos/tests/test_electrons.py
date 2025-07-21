# unit tests for rates
import numpy as np
from pytest import approx

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


def ideal_gas(rho, T, comp):
    """An electron ideal gas---the pressure is the same regardless if
    we are relativistic or non-relativistic."""

    p = rho * constants.k * T * comp.ye * constants.N_A

    return p


class TestElectronEOS:

    def test_low_temperature(self):
        # at low temperatures, we should be completely degenerate and
        # follow the zero-temperature expressions approximately.

        e = ElectronEOS(include_positrons=False)

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

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        T = 1.e10

        for rho in [1.e1, 1.e2, 1.e3]:
            es = e.pe_state(rho, T, comp)
            p_ideal = ideal_gas(rho, T, comp)

            assert es.p_e == approx(p_ideal, rel=1.e-4)

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

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps = 1.e-8

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                es_r = e.pe_state(rho * (1.0 + eps), T, comp)

                dnedr_approx = (es_r.n_e - es.n_e) / (eps * rho)
                assert es.dne_drho == approx(dnedr_approx)

    def test_ne_temp_derivs(self):

        # formally, dn_e/dT = 0 when we don't include positrons.  The
        # EOS gives this correctly.  So to compare with a finite
        # difference approximation, we need a scale to compare with,
        # so we'll compare to n_e/T.

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # use a relatively large eps because of how weakly beta enters when beta << 1

        eps = 1.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                es_r = e.pe_state(rho, T * (1.0 + eps), comp)

                dnedr_approx = (es_r.n_e - es.n_e) / (eps * T)
                scale = es.n_e / T
                assert es.dne_dT == approx(dnedr_approx / scale, abs=5.e-9)

    def test_pres_rho_derivs(self):

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_rho = 1.e-8

        for T in [1.e4, 1.e7, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:
                es = e.pe_state(rho, T, comp)
                es_r = e.pe_state(rho * (1.0 + eps_rho), T, comp)
                dpdr_approx = (es_r.p_e - es.p_e) / (eps_rho * rho)
                assert es.dpe_drho == approx(dpdr_approx)

    def test_pres_temp_derivs(self):

        # the temperature derivatives are very hard, since beta enters
        # into the integrand as sqrt(1 + 0.5 * x * beta).  For small
        # beta (beta <~ 1.e-5), we can have a hard time with
        # finite-differencing.  We use a relatively large epsilon
        # here, a fourth-order difference approximation, and also skip
        # very cool temperature (< 1.e6 K)

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_T = 5.e-4

        for T in [1.e6, 1.e7, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                dtemp = eps_T * T
                fvals = []
                for i in [-2, -1, 0, 1, 2]:
                    _es = e.pe_state(rho, T + i*dtemp, comp)
                    fvals.append(_es.p_e)

                deriv = (fvals[0] - 8.0 * fvals[1] + 8.0 * fvals[3] - fvals[4]) / (12 * dtemp)
                assert es.dpe_dT == approx(deriv, rel=1.e-3)
