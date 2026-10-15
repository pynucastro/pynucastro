import dual

import numpy as np

from pytest import approx


class TestDual:

    def test_addition(self):

        # here, both of these seed the derivative to 1
        x = dual.Dual(2.0)
        y = dual.Dual(5.5)

        f = x + y
        assert f.value == 7.5
        assert f.deriv == 2.0

        g = x + 2
        assert g.value == 4
        assert g.deriv == 1

    def test_subtraction(self):

        x = dual.Dual(2.0)
        y = dual.Dual(5.5)

        f = x - y
        assert f.value == -3.5
        assert f.deriv == 0.0

        g = x - 2
        assert g.value == 0
        assert g.deriv == 1

    def test_multiply(self):

        x = dual.Dual(2.0)
        y = dual.Dual(5.5)

        assert (x * y).value == 11.0
        assert (x * y).deriv == 7.5

    def test_func1(self):

        x = dual.Dual(3.0)

        f = 2 * x**2 + 3 * x + 4

        assert f.value == 31
        assert f.deriv == 15

    def test_func2(self):

        x = dual.Dual(-1.0)

        g = (x**2 + 10) / dual.exp(x) + 5

        assert g.value == approx(5.0 + 11.0 * np.exp(1.0))
        assert g.deriv == approx(-13.0 * np.exp(1.0))

    def test_func3(self):

        x = dual.Dual(2.0)

        h = dual.sqrt(x**3.5 + 1)

        assert h.value == approx(3.509089411654363)
        assert h.deriv == approx(2.821100797185029)

    def test_func4(self):

        x = dual.Dual(4.0)

        q = dual.cbrt(x * dual.exp(2.0 * x) + 1)

        assert q.value == approx(22.8463813891676)
        assert q.deriv == approx(17.13333491422937)
