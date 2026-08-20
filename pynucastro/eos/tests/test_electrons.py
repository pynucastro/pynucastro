# unit tests for rates

import sys

import numpy as np
import pytest
from pytest import approx

from pynucastro import Composition
from pynucastro.constants import constants
from pynucastro.eos import ElectronEOS
from pynucastro.eos.difference_utils import fourth_order_diff

#np.seterr(all="warn")


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
            es, _ = e.pe_state(rho, T, comp)
            p_zt, e_zt = zero_temperature_eos(rho, comp)

            assert es.p == approx(p_zt, rel=1.e-4)
            assert es.e == approx(e_zt, rel=1.e-4)

    def test_high_temperature(self):
        # at high temperatures, we should be an ideal electron gas.

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        T = 1.e10

        for rho in [1.e1, 1.e2, 1.e3]:
            es, _ = e.pe_state(rho, T, comp)
            p_ideal = ideal_gas(rho, T, comp)

            assert es.p == approx(p_ideal, rel=1.e-4)

    def test_ne_rho_derivs(self):

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps = 1.e-8

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es, _ = e.pe_state(rho, T, comp)
                es_eps, _ = e.pe_state(rho * (1.0 + eps), T, comp)

                dnedr_approx = (es_eps.n - es.n) / (eps * rho)
                assert es.dn_drho == approx(dnedr_approx)

    @pytest.mark.skipif(sys.platform == "darwin", reason="Macs give different roundoff when the value is ~ 0")
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

                es, _ = e.pe_state(rho, T, comp)
                es_eps, _ = e.pe_state(rho, T * (1.0 + eps), comp)

                dnedT_approx = (es_eps.n - es.n) / (eps * T)
                scale = es.n / T
                assert es.dn_dT == approx(dnedT_approx/scale, abs=5.e-9)

    def test_pres_rho_derivs(self):

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_rho = 1.e-8

        for T in [1.e4, 1.e7, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:
                es, _ = e.pe_state(rho, T, comp)
                es_eps, _ = e.pe_state(rho * (1.0 + eps_rho), T, comp)
                dpdr_approx = (es_eps.p - es.p) / (eps_rho * rho)
                assert es.dp_drho == approx(dpdr_approx)

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

                es, _ = e.pe_state(rho, T, comp)
                dtemp = eps_T * T
                deriv = fourth_order_diff(lambda _T: e.pe_state(rho, _T, comp)[0],  # pylint: disable=cell-var-from-loop
                                         T, dtemp, component="p")
                assert es.dp_dT == approx(deriv, rel=3.e-4)

    def test_e_rho_derivs(self):

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # we need to use a slightly larger eps for this because we
        # divide out the rho dependence in e
        eps_rho = 5.e-4

        for T in [1.e4, 1.e7, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es, _ = e.pe_state(rho, T, comp)
                drho = eps_rho * rho
                deriv = fourth_order_diff(lambda _rho: e.pe_state(_rho, T, comp)[0],  # pylint: disable=cell-var-from-loop
                                          rho, drho, component="e")
                assert es.de_drho == approx(deriv, rel=1.e-3)

    def test_e_temp_derivs(self):

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

                es, _ = e.pe_state(rho, T, comp)
                dtemp = eps_T * T
                deriv = fourth_order_diff(lambda _T: e.pe_state(rho, _T, comp)[0],  # pylint: disable=cell-var-from-loop
                                         T, dtemp, component="e")
                assert es.de_dT == approx(deriv, rel=5.e-4)

    def test_gamma_limits(self):

        e = ElectronEOS(include_positrons=False)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # at low density, we are either non-relativistic degenerate
        # or ideal gas, and should have rho e / p = 3/2 (gamma = 5/3)

        rho = 1.e-2
        T = 1.e6

        es, _ = e.pe_state(rho, T, comp)
        assert rho * es.e / es.p == approx(1.5, rel=5.e-4)
