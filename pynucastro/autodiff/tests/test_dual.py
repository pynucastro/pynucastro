"""Regression tests for Python and compiled scalar differentiation."""

import operator
from collections import namedtuple

import numpy as np
import pytest
from numba import njit, typeof, types
from pytest import approx

from pynucastro.autodiff import dual


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


@pytest.mark.parametrize("op", [operator.add, operator.sub, operator.mul,
                               operator.truediv, operator.pow])
@pytest.mark.parametrize("a,b", [(dual.Dual(2.0, 3.0), dual.Dual(4.0, -2.0)),
                                 (dual.Dual(2.0, 3.0), 4.0),
                                 (dual.Dual(2.0, 3.0), 4)])
def test_compiled_arithmetic(op, a, b):
    def evaluate(a, b):
        return op(a, b)

    expected = evaluate(a, b)
    result = njit(evaluate)(a, b)
    assert result.value == approx(expected.value)
    assert result.deriv == approx(expected.deriv)


@pytest.mark.parametrize("op", [operator.add, operator.sub, operator.mul, operator.truediv])
@pytest.mark.parametrize("scalar", [3, 3.0])
def test_reflected_arithmetic(op, scalar):
    def evaluate(x):
        return op(scalar, x)

    x = dual.Dual(2.0, 5.0)
    expected = evaluate(x)
    result = njit(evaluate)(x)
    assert result.value == approx(expected.value)
    assert result.deriv == approx(expected.deriv)


@pytest.mark.parametrize("op", [operator.lt, operator.le, operator.gt,
                               operator.ge, operator.eq, operator.ne])
@pytest.mark.parametrize("a,b", [(dual.Dual(2.0, 3.0), dual.Dual(2.0, -2.0)),
                                 (dual.Dual(2.0), 3.0), (3.0, dual.Dual(2.0))])
def test_comparisons(op, a, b):
    def evaluate(a, b):
        return op(a, b)

    assert njit(evaluate)(a, b) == evaluate(a, b)


@pytest.mark.parametrize("func", [min, max])
@pytest.mark.parametrize("a,b", [(dual.Dual(2.0, 3.0), dual.Dual(2.0, -2.0)),
                                 (dual.Dual(2.0), 3.0), (3.0, dual.Dual(2.0)),
                                 (dual.Dual(3.0), 3.0), (3.0, dual.Dual(3.0))])
def test_clipping(func, a, b):
    def evaluate(a, b):
        return func(a, b)

    expected = evaluate(a, b)
    result = njit(evaluate)(a, b)
    assert result.value == (expected.value if isinstance(expected, dual.Dual) else expected)
    assert result.deriv == (expected.deriv if isinstance(expected, dual.Dual) else 0.0)


@pytest.mark.parametrize("func", [dual.exp, dual.log, dual.log1p, dual.sqrt,
                                 dual.cbrt, dual.sin, dual.cos, dual.arctan])
def test_math_derivatives(func):
    x = 1.3
    step = 1.e-5
    result = func(dual.Dual(x, 2.0))
    assert result.value == approx(func(x))
    assert result.deriv == approx(2.0 * (func(x + step) - func(x - step)) / (2 * step), rel=1.e-8)


@pytest.mark.parametrize("x", [-2.0, 0.0, 2.0])
def test_unary(x):
    @njit
    def evaluate(x):
        return -x, abs(x)

    negative, magnitude = evaluate(dual.Dual(x, 3.0))
    assert negative.value == -x
    assert negative.deriv == -3.0
    assert magnitude.value == abs(x)
    assert magnitude.deriv == 3.0 * np.sign(x)


def test_polynomial():
    @njit
    def evaluate(x):
        return x * x + 3.0 * x - 1.0

    result = evaluate(dual.Dual(1.3))
    assert result.value == approx(evaluate(1.3))
    assert result.deriv == approx(5.6)


def test_constructor_and_identity():
    @njit
    def construct(x):
        return dual.Dual(x)

    result = construct(2.0)
    assert result.deriv == 1.0
    assert isinstance(typeof(result), types.BaseNamedTuple)
    # pylint: disable=protected-access
    assert not dual._is_dual(typeof(namedtuple("Other", "value deriv")(2.0, 1.0)))
    with pytest.raises(AttributeError):
        result.value = 4.0


def test_power_edges():
    @njit
    def evaluate(x, y):
        return x**y

    result = evaluate(dual.Dual(0.0), 0)
    assert result.value == 1.0
    assert result.deriv == 0.0
    with pytest.raises(ValueError, match="positive base"):
        evaluate(dual.Dual(-2.0), dual.Dual(3.0))
