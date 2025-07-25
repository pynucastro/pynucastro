# unit tests for rates

import numpy as np
from pytest import approx

from pynucastro import Composition
from pynucastro.constants import constants
from pynucastro.eos import ElectronEOS
from pynucastro.eos.difference_utils import fourth_order_diff, sixth_order_diff


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

            assert es.p_e == approx(p_zt, rel=1.e-5)
            assert es.e_e == approx(e_zt, rel=1.e-5)

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

        eps_rho = 5.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                drho = eps_rho * rho
                deriv = sixth_order_diff(lambda _rho: e.pe_state(_rho, T, comp),  # pylint: disable=cell-var-from-loop
                                         rho, drho, component="n_e")
                assert es.dne_drho == approx(deriv, rel=3.e-6)

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

    def test_np_rho_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_rho = 1.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                drho = eps_rho * rho
                deriv = fourth_order_diff(lambda _rho: e.pe_state(_rho, T, comp),  # pylint: disable=cell-var-from-loop
                                          rho, drho, component="n_pos")
                assert es.dnp_drho == approx(deriv, rel=1.e-5)

                # since n_e - n_pos = N_A (Z/A) rho, it should be the case that
                # dne/drho = dnp/drho + N_A (Z/A)

                assert es.dne_drho == approx(constants.N_A * comp.zbar / comp.abar + es.dnp_drho)

    def test_np_temp_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # use a relatively large eps because of how weakly beta enters when beta << 1

        eps_T = 1.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                dtemp = eps_T * T
                deriv = fourth_order_diff(lambda _T: e.pe_state(rho, _T, comp),  # pylint: disable=cell-var-from-loop
                                          T, dtemp, component="n_pos")
                if es.dnp_dT == 0 and es.n_pos != 0:
                    # no pair production
                    scale = es.n_pos / T
                    assert es.dnp_dT == approx(deriv / scale, abs=5.e-9)
                else:
                    assert es.dnp_dT == approx(deriv, rel=1.e-5)

                # since n_e - n_pos = N_A (Z/A) rho, it should be the case that
                # dne/T = dnp/T

                if es.dne_dT == 0.0:
                    # if dne_dT == 0, then just make sure that they are equal to
                    # roundoff error
                    scale = es.n_e / T
                    assert es.dne_dT == approx(es.dnp_dT/scale, abs=1.e-15)
                else:
                    assert es.dne_dT == approx(es.dnp_dT)

    def test_pe_rho_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_rho = 5.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                drho = eps_rho * rho
                deriv = fourth_order_diff(lambda _rho: e.pe_state(_rho, T, comp),  # pylint: disable=cell-var-from-loop
                                          rho, drho, component="p_e")
                assert es.dpe_drho == approx(deriv, rel=1.e-5)

    def test_pe_temp_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # use a relatively large eps because of how weakly beta enters when beta << 1

        eps_T = 1.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)

                # at high density and low T, the temperature derivative
                # is effectively zero, and differencing struggles
                if es.eta > 1.e4:
                    continue

                dtemp = eps_T * T
                deriv = sixth_order_diff(lambda _T: e.pe_state(rho, _T, comp),  # pylint: disable=cell-var-from-loop
                                         T, dtemp, component="p_e")
                assert es.dpe_dT == approx(deriv, rel=5.e-4)

    def test_pp_rho_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_rho = 5.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                drho = eps_rho * rho
                deriv = fourth_order_diff(lambda _rho: e.pe_state(_rho, T, comp),  # pylint: disable=cell-var-from-loop
                                          rho, drho, component="p_pos")
                assert es.dpp_drho == approx(deriv, rel=1.e-5)

    def test_pp_temp_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # use a relatively large eps because of how weakly beta enters when beta << 1

        eps_T = 1.e-4

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 5.e8, 1.e9]:

                es = e.pe_state(rho, T, comp)

                # at high density and low T, the temperature derivative
                # is effectively zero, and differencing struggles
                if es.eta > 1.e4:
                    continue

                dtemp = eps_T * T
                deriv = fourth_order_diff(lambda _T: e.pe_state(rho, _T, comp),  # pylint: disable=cell-var-from-loop
                                          T, dtemp, component="p_pos")
                assert es.dpp_dT == approx(deriv, rel=1.e-5)

    def test_ee_rho_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_rho = 5.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                drho = eps_rho * rho
                deriv = fourth_order_diff(lambda _rho: e.pe_state(_rho, T, comp),  # pylint: disable=cell-var-from-loop
                                          rho, drho, component="e_e")
                assert es.dee_drho == approx(deriv, rel=1.e-5)

    def test_ee_temp_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # use a relatively large eps because of how weakly beta enters when beta << 1

        eps_T = 1.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)

                # at high density and low T, the temperature derivative
                # is effectively zero, and differencing struggles
                if es.eta > 1.e4:
                    continue

                dtemp = eps_T * T
                deriv = sixth_order_diff(lambda _T: e.pe_state(rho, _T, comp),  # pylint: disable=cell-var-from-loop
                                         T, dtemp, component="e_e")
                assert es.dee_dT == approx(deriv, rel=5.e-4)

    def test_ep_rho_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        eps_rho = 5.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)
                drho = eps_rho * rho
                deriv = fourth_order_diff(lambda _rho: e.pe_state(_rho, T, comp),  # pylint: disable=cell-var-from-loop
                                          rho, drho, component="e_pos")
                assert es.dep_drho == approx(deriv, rel=1.e-5)

    def test_ep_temp_derivs(self):

        e = ElectronEOS(include_positrons=True)

        comp = Composition(["h1", "he4", "c12", "ne22"])
        comp.set_equal()

        # use a relatively large eps because of how weakly beta enters when beta << 1

        eps_T = 1.e-5

        for T in [1.e4, 1.e6, 1.e9]:
            for rho in [1.e-2, 1.e2, 1.e5, 1.e9]:

                es = e.pe_state(rho, T, comp)

                # at high density and low T, the temperature derivative
                # is effectively zero, and differencing struggles
                if es.eta > 1.e4:
                    continue

                dtemp = eps_T * T
                deriv = sixth_order_diff(lambda _T: e.pe_state(rho, _T, comp),  # pylint: disable=cell-var-from-loop
                                         T, dtemp, component="e_pos")
                assert es.dep_dT == approx(deriv, rel=5.e-4)
