import numpy as np
import pytest
from pytest import approx

from pynucastro import Nucleus
from pynucastro.rates import TemperatureTabularRate


def rate_fcn(T):
    # we can't just do a straight powerlaw, since that is a
    # straight line in log T - log rate, which is what we
    # interpolate on.  So add a cutoff
    return 1.e-2 * (T / 1.e8)**12 * (101 - T / 1.e8)**10


class TestTemperatureTabular:

    # we'll just makeup a fake power law to see
    # if the interpolation recovers the correct values

    @pytest.fixture(scope="class")
    def myrate(self):
        Ts = np.logspace(np.log10(1.e7), np.log10(1.e10), 75)
        rate = rate_fcn(Ts)

        # it's a fake rate, so we can make up any
        # nuclei involved
        reactants = [Nucleus("c12"), Nucleus("p")]
        products = [Nucleus("n13")]

        return TemperatureTabularRate(np.log(Ts/1.e9),
                                      np.log(rate),
                                      reactants=reactants,
                                      products=products)

    def test_eval(self, myrate):

        for T in [3.12e7, 9.81e7, 2.56e8, 4.7e8, 1.2e9, 5.7e9]:
            assert myrate.eval(T) == approx(rate_fcn(T),
                                            rel=1.e-4, abs=1.e-100)

    def test_ydot_string(self, myrate):
        ydot_str = "rho*Y[jp]*Y[jc12]*rate_eval.C12_p_to_N13_temptab"
        assert myrate.ydot_string_py() == ydot_str
