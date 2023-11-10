import math

from pytest import approx

from pynucastro.constants import constants


class TestConstants:
    def test_constants(self):

        assert constants.m_u == approx(1.66054e-24, abs=0, rel=1.e-6)
        assert constants.m_u_MeV == approx(constants.m_u * constants.c_light**2 * constants.erg2MeV, abs=0, rel=1.e-6)
        assert constants.m_p == approx(1.67262192e-24, abs=0, rel=1.e-6)
        assert constants.m_p_MeV == approx(constants.m_p * constants.c_light**2 * constants.erg2MeV, abs=0, rel=1.e-6)
        assert constants.m_n == approx(1.674927498e-24, abs=0, rel=1.e-6)
        assert constants.m_n_MeV == approx(constants.m_n * constants.c_light**2 * constants.erg2MeV, abs=0, rel=1.e-6)
        assert constants.m_e == approx(9.1093837e-28, abs=0, rel=1.e-6)
        assert constants.m_e_MeV == approx(constants.m_e * constants.c_light**2 * constants.erg2MeV)

        assert constants.N_A == approx(6.02214e23)

        assert constants.k == approx(1.38064e-16, abs=0, rel=1.e-5)
        assert constants.k_MeV == approx(constants.k * constants.erg2MeV, abs=0, rel=1.e-5)

        assert constants.h == approx(6.62607e-27, abs=0, rel=1.e-6)
        assert constants.hbar == approx(6.62607e-27 / 2 / math.pi, abs=0, rel=1.e-6)

        assert constants.q_e == approx(4.803e-10, abs=0, rel=1.e-4)
        assert constants.c_light == approx(2.997924e10)

        assert constants.erg2MeV == approx(1.0 / 1.60218e-6, abs=0, rel=1.e-4)
        assert constants.MeV2erg == approx(1.60218e-6, abs=0, rel=1.e-4)
