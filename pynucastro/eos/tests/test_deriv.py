# unit tests for rates
import numpy as np
from pytest import approx

import pynucastro.eos.difference_utils as du


def f1(x):
    """this is the test function from Ridders' paper"""
    return np.exp(x) / (np.sin(x) - x**2)


def f1prime(x):
    """analytic derivative for f1"""
    return np.exp(x) * (-x**2 + 2*x - np.sqrt(2) * np.cos(x + np.pi/4)) / (x**2 - np.sin(x))**2


class TestDerivs:

    def test_fourth(self):

        x0 = 1
        h = 1.e-3
        actual_deriv = f1prime(x0)

        deriv = du.fourth_order_diff(f1, x0, h)
        assert deriv == approx(actual_deriv)
        err = np.abs(deriv - actual_deriv)

        h /= 2
        deriv2 = du.fourth_order_diff(f1, x0, h)
        assert deriv == approx(actual_deriv)
        err2 = np.abs(deriv2 - actual_deriv)

        assert err / err2 == approx(2**4, rel=1.e-3, abs=1.e-3)

    def test_sixth(self):

        x0 = 1
        h = 1.e-3
        actual_deriv = f1prime(x0)

        deriv = du.sixth_order_diff(f1, x0, h)
        assert deriv == approx(actual_deriv)
        err = np.abs(deriv - actual_deriv)

        h /= 2
        deriv2 = du.sixth_order_diff(f1, x0, h)
        assert deriv == approx(actual_deriv)
        err2 = np.abs(deriv2 - actual_deriv)

        # this seems to converge better than h**6, so just
        # check for improvement

        assert err2 < err

    def test_eighth(self):

        x0 = 1

        # we seem to need a larger h for good convergence
        h = 0.002
        actual_deriv = f1prime(x0)

        deriv = du.eighth_order_diff(f1, x0, h)
        assert deriv == approx(actual_deriv)
        err = np.abs(deriv - actual_deriv)

        h /= 2
        deriv2 = du.eighth_order_diff(f1, x0, h)
        assert deriv == approx(actual_deriv)
        err2 = np.abs(deriv2 - actual_deriv)

        # this seems to converge better than h**6, so just
        # check for improvement

        assert err2 < err

    def test_ridders(self):

        x0 = 1

        h = 0.1
        actual_deriv = f1prime(x0)

        deriv, _ = du.adaptive_diff(f1, x0, h)
        assert deriv == approx(actual_deriv, rel=1.e-9, abs=1.e-11)
