"""Two-component seeds in Python and compiled forward-mode differentiation."""

import numpy as np
import pytest
from numba import njit

from pynucastro.autodiff import dual
from pynucastro.autodiff.derivative_pair import DerivativePair


@pytest.mark.parametrize("compiled", [False, True])
def test_two_variables(compiled):
    """Independent seeds recover both analytic partial derivatives."""
    def evaluate(x, y):
        return dual.exp(x * y) + x / y - 3.0 / x + x**y - y + 2 * x

    func = njit(evaluate) if compiled else evaluate
    x, y = 1.3, 2.1
    result = func(dual.Dual(x, DerivativePair(1.0, 0.0)),
                  dual.Dual(y, DerivativePair(0.0, 1.0)))
    assert result.value == pytest.approx(evaluate(x, y))
    assert result.deriv.d1 == pytest.approx(y * np.exp(x * y) + 1 / y + 3 / x**2 + y * x**(y - 1) + 2)
    assert result.deriv.d2 == pytest.approx(x * np.exp(x * y) - x / y**2 + x**y * np.log(x) - 1)


@pytest.mark.parametrize("func", [dual.exp, dual.log, dual.log1p, dual.sqrt,
                                  dual.cbrt, dual.cos, dual.sin, dual.arctan])
def test_math_helpers(func):
    """Both components follow the chain rule for every math helper."""
    @njit
    def evaluate(x):
        return func(x)

    step = 1.e-5
    expected = (func(1.3 + step) - func(1.3 - step)) / (2 * step)

    # we can think of x here as depending on two separate quantities, x(u, v).
    # here we seed things such that ∂x/∂u = 2 and ∂x/∂v = -3
    x = dual.Dual(1.3, DerivativePair(2.0, -3.0))

    for result in (func(x), evaluate(x)):
        # this will be ∂f/∂x ∂x/∂u
        assert result.deriv.d1 == pytest.approx(2 * expected, rel=1.e-8)
        # this will be ∂f/∂x ∂x/∂v
        assert result.deriv.d2 == pytest.approx(-3 * expected, rel=1.e-8)


def test_constant_results():
    """Expressions independent of x have zero derivatives despite nonzero seeds."""
    @njit
    def evaluate(x):
        return x**0, dual.constant_like(x, 7.0)

    x = dual.Dual(2.0, DerivativePair(2.0, -3.0))
    for func in (evaluate.py_func, evaluate):
        for result, value in zip(func(x), (1.0, 7.0)):
            assert result.value == value
            assert result.deriv.d1 == 0.0
            assert result.deriv.d2 == 0.0


def test_clipped_constant_branches():
    """Selecting a scalar bound gives a locally constant result."""
    @njit
    def evaluate(x):
        return min(x, -1.0), max(1.0, x)

    x = dual.Dual(0.0, DerivativePair(2.0, -3.0))
    # Near x = 0, min selects -1 and max selects 1, independently of x.

    # Python returns the scalar bounds; we use the .py_func attribute
    # to get the uncompiled version
    assert evaluate.py_func(x) == (-1.0, 1.0)

    # Numba preserves the dual seed shape.
    for result, value in zip(evaluate(x), (-1.0, 1.0)):
        assert result.value == value
        assert result.deriv.d1 == 0.0
        assert result.deriv.d2 == 0.0


def test_abs_at_zero():
    """At its nondifferentiable point, abs uses our zero-derivative convention."""
    @njit
    def evaluate(x):
        return abs(x)

    x = dual.Dual(0.0, DerivativePair(2.0, -3.0))
    for func in (evaluate.py_func, evaluate):
        result = func(x)
        assert result.value == 0.0
        assert result.deriv.d1 == 0.0
        assert result.deriv.d2 == 0.0


def test_negation():
    """Negation reverses both nonzero seeds; it is not a constant expression."""
    @njit
    def evaluate(x):
        return -x

    x = dual.Dual(1.3, DerivativePair(2.0, -3.0))
    for func in (evaluate.py_func, evaluate):
        result = func(x)
        assert result.value == -1.3
        assert result.deriv.d1 == -2.0
        assert result.deriv.d2 == 3.0


@pytest.mark.parametrize("seed", [2, DerivativePair(2, -3)])
def test_clipping_integer_seeds(seed):
    """Clipping promotes integer seeds consistently with scalar bounds."""
    @njit
    def evaluate(x):
        return min(x, 1.0), max(x, 3.0)

    for result in evaluate(dual.Dual(2.0, seed)):
        if isinstance(seed, DerivativePair):
            assert result.deriv.d1 == 0.0
            assert result.deriv.d2 == 0.0
        else:
            assert result.deriv == 0.0
