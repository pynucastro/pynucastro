import math
from pytest import approx

from pynucastro.constants import constants


class TestConstants:
    def test_constants(self):

        assert constants.m_u == approx(1.66054e-24)
        assert constants.m_u_MeV == approx(931.5, rel=1.e-4)
        assert constants.k == approx(1.30806e-16)
        assert constants.k_MeV == approx(8.61733e-11)
        assert constants.h == approx(6.62607e-27)
        assert constants.h == approx(6.62607e-27 / 2 / math.pi)
        assert constants.N_A == approx(6.02214e23)
        assert constants.q_e == approx(4.802e-10, rel=1.e-4)
        assert constants.c_light == approx(2.997924e10)
        assert constants.m_p == approx(1.67262192e-24)
        assert constants.m_n == approx(1.674927498e-24)
        assert constants.erg2MeV == approx(1.0 / 1.60218e-6, rel=1.e-4)
        assert constants.MeV2erg == approx(1.60218e-6, rel=1.e-4)
