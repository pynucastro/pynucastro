import pytest
from pytest import approx

import pynucastro as pyna


class TestStellarEOS:

    @pytest.fixture(scope="class")
    def eos(self):
        return pyna.StellarEOS(include_positrons=True)

    @pytest.fixture(scope="class")
    def comp(self):
        comp = pyna.Composition(["he4"])
        comp.set_equal()
        return comp

    def test_low_dens_high_T(self, eos, comp):

        # these comparison numbers come from Microphysics eos_cell /
        # Frank Timmes' EOS (with Coulomb corrections disabled)

        rho = 1.e-2
        T = 1.e9

        state = eos.pe_state(rho, T, comp)

        assert state.eta == approx(-5.929886677, rel=1.e-4)
        assert state.n_ele + state.n_pos == approx(1.085890355e+27, rel=1.e-4)

        assert state.p == approx(2.671956308e+21, rel=1.e-4)
        assert state.e == approx(8.720184159e+23, rel=1.e-4)

        assert state.dp_dT == approx(1.139214265e+13, rel=1.e-4)
        assert state.dp_drho == approx(2.078629692e+16, rel=1.e-4)

        assert state.de_dT == approx(3.944397464e+15, rel=1.e-4)
        assert state.de_drho == approx(-8.720186344e+25, rel=1.e-4)

        assert state.c_v == approx(3.944397464e+15, rel=1.e-4)
        assert state.c_p == approx(6.243580707e+22, rel=1.e-4)

        assert state.gamma1 == approx(1.231404775, rel=1.e-4)

    def test_high_dens_high_T(self, eos, comp):

        # these comparison numbers come from Microphysics eos_cell /
        # Frank Timmes' EOS (with Coulomb corrections disabled)

        rho = 1.e7
        T = 1.e9

        state = eos.pe_state(rho, T, comp)

        assert state.eta == approx(5.566355932, rel=1.e-4)
        assert state.n_ele + state.n_pos == approx(3.011070908e+30, rel=1.e-4)

        assert state.p == approx(1.148636161e+24, rel=1.e-4)
        assert state.e == approx(2.145687528e+17, rel=1.e-4)

        assert state.dp_dT == approx(3.934408418e+14, rel=1.e-4)
        assert state.dp_drho == approx(1.483248593e+17, rel=1.e-4)

        assert state.de_dT == approx(76651621.61, rel=1.e-4)
        assert state.de_drho == approx(7551953195, rel=1.e-4)

        assert state.c_v == approx(76651621.61, rel=1.e-4)
        assert state.c_p == approx(87087882.72, rel=1.e-4)

        assert state.gamma1 == approx(1.467127459, rel=1.e-4)

    def test_high_dens_med_T(self, eos, comp):

        # these comparison numbers come from Microphysics eos_cell /
        # Frank Timmes' EOS (with Coulomb corrections disabled)

        rho = 1.e7
        T = 1.e7

        state = eos.pe_state(rho, T, comp)

        assert state.eta == approx(589.4169538, rel=1.e-4)
        assert state.n_ele + state.n_pos == approx(3.011070896e+30, rel=1.e-4)

        assert state.p == approx(8.480873295e+23, rel=1.e-4)
        assert state.e == approx(1.607501131e+17, rel=1.e-4)

        assert state.dp_dT == approx(2.097958465e+14, rel=1.e-4)
        assert state.dp_drho == approx(1.228495432e+17, rel=1.e-4)

        assert state.de_dT == approx(31642880.04, rel=1.e-4)
        assert state.de_drho == approx(8459893711, rel=1.e-4)

        assert state.c_v == approx(31642880.04, rel=1.e-4)
        assert state.c_p == approx(31678707.84, rel=1.e-4)

        assert state.gamma1 == approx(1.450188394, rel=1.e-4)

    def test_high_dens_low_T(self, eos, comp):

        # these comparison numbers come from Microphysics eos_cell /
        # Frank Timmes' EOS (with Coulomb corrections disabled)

        rho = 1.e8
        T = 1.e5

        state = eos.pe_state(rho, T, comp)

        assert state.eta == approx(168930.2028, rel=1.e-4)
        assert state.n_ele + state.n_pos == approx(3.011070897e+31, rel=1.e-4)

        assert state.p == approx(2.150884839e+25, rel=1.e-4)
        assert state.e == approx(4.871994835e+17, rel=1.e-4)

        assert state.dp_dT == approx(2.078685394e+15, rel=1.e-4)
        assert state.dp_drho == approx(2.949191274e+17, rel=1.e-4)

        assert state.de_dT == approx(31181167.26, rel=1.e-4)
        assert state.de_drho == approx(2150864052, rel=1.e-4)

        assert state.c_v == approx(31181167.26, rel=1.e-4)
        assert state.c_p == approx(31181313.77, rel=1.e-4)

        assert state.gamma1 == approx(1.371159013, rel=1.e-4)
