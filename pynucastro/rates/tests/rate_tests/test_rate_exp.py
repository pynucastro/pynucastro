# test Rate.get_rate_exponent

import pytest
from pytest import approx

from pynucastro.rates import ReacLibRate, SingleSet


class TestRateExponent:

    @pytest.fixture(scope="class")
    @classmethod
    def sl_rate(cls, starlib_library):
        return starlib_library.get_rate_by_name("n14(p,g)o15")

    @pytest.fixture(scope="class")
    @classmethod
    def rl_rate(cls, reaclib_library):
        return reaclib_library.get_rate_by_name("n14(p,g)o15")

    def test_get_rate_exp(self, sl_rate, rl_rate):

        T0 = 3.e7

        assert sl_rate.get_rate_exponent(T0) == approx(15.5036928097374,
                                                       rel=1.e-10, abs=1.e-100)

        assert rl_rate.get_rate_exponent(T0) == approx(15.601858116021104,
                                                       rel=1.e-10, abs=1.e-100)

    def test_underflow(self):

        rate = ReacLibRate(reactants=["c12"],
                           products=["he4", "he4", "he4"],
                           sets=[SingleSet([-800., 0., 0., 0., 0., 0., 2.],
                                           "test  ")])

        # this underflows because of the first coefficient, exp(-800)
        assert rate.eval(1.e9) == 0.0

        # but our derivative should recover the T**2 form from
        # our single set (the last coefficient)
        assert rate.get_rate_exponent(1.e9) == approx(2.0)
