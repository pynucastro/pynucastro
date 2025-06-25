import pytest
from pytest import approx

import pynucastro as pyna
from pynucastro.screening import (chugunov_2007, chugunov_2009, debye_huckel,
                                  make_plasma_state, make_screen_factors,
                                  potekhin_1998, screen5, screening_check)


class TestScreen:
    @pytest.fixture(scope="class")
    def nuclei(self):
        return [pyna.Nucleus("h1"),
                pyna.Nucleus("he4"),
                pyna.Nucleus("c12"),
                pyna.Nucleus("o16"),
                pyna.Nucleus("n14"),
                pyna.Nucleus("ca40")]

    @pytest.fixture(scope="class")
    def plasma_state(self, nuclei):
        temp = 1e6
        dens = 1e5
        comp = pyna.Composition(nuclei)
        comp.set_solar_like()
        return make_plasma_state(temp, dens, comp.get_molar())

    @pytest.fixture(scope="class")
    def scn_fac(self):
        c12 = pyna.Nucleus("c12")
        he4 = pyna.Nucleus("he4")
        return make_screen_factors(c12, he4)

    def test_plasma_state(self, plasma_state):
        assert plasma_state.temp == approx(1e6)
        assert plasma_state.dens == approx(1e5)
        assert plasma_state.qlam0z == approx(82.64519344765309)
        assert plasma_state.taufac == approx(14.162368047477718)
        assert plasma_state.aa == approx(10.00149809343336)
        assert plasma_state.abar == approx(1.2966614825934775)
        assert plasma_state.zbar == approx(1.1021622602044556)
        assert plasma_state.z2bar == approx(1.4036360549074394)
        assert plasma_state.n_e == approx(5.118819647768954e+28)
        assert plasma_state.gamma_e_fac == approx(10001498.09343337)

    def test_screen_factors(self, scn_fac):
        assert scn_fac.z1 == 6
        assert scn_fac.a1 == 12
        assert scn_fac.z2 == 2
        assert scn_fac.a2 == 4
        assert scn_fac.zs13 == approx(2.0)
        assert scn_fac.zhat == approx(9.013634402695844)
        assert scn_fac.zhat2 == approx(-1.0661692705686794)
        assert scn_fac.lzav == approx(0.675775180180274)
        assert scn_fac.aznut == approx(7.55952629936924)
        assert scn_fac.ztilde == approx(1.5385208213635064)

    def test_chugunov_2007(self, plasma_state, scn_fac):
        scor = chugunov_2007(plasma_state, scn_fac)
        assert scor == approx(7.785569477042635e+33)

    def test_chugunov_2009(self, plasma_state, scn_fac):
        scor = chugunov_2009(plasma_state, scn_fac)
        assert scor == approx(2.87983449091315e+33)

    def test_debye_huckel(self, plasma_state, scn_fac):
        scor = debye_huckel(plasma_state, scn_fac)
        assert scor == approx(1.9424263952412558e+130)

    def test_potekhin_1998(self, plasma_state, scn_fac):
        scor = potekhin_1998(plasma_state, scn_fac)
        assert scor == approx(1.0508243810383098e+36)

    def test_screen5(self, plasma_state, scn_fac):
        scor = screen5(plasma_state, scn_fac)
        assert scor == approx(4.049488384394272e+33)

    @pytest.mark.parametrize("screen_func", [chugunov_2007, chugunov_2009, potekhin_1998, screen5])
    def test_screening_check(self, plasma_state, scn_fac, screen_func):
        wrapped_screen_func = screening_check()(screen_func)
        scor_wrapped = wrapped_screen_func(plasma_state, scn_fac)
        scor = screen_func(plasma_state, scn_fac)
        assert scor_wrapped == approx(scor)
