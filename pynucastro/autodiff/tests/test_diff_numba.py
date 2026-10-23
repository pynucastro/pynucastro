"""Exercise dual arithmetic through Numba's nopython dispatch."""

import operator

import numpy as np
import pytest
from numba import njit

from pynucastro.autodiff import dual


@pytest.mark.parametrize("name, derivative", [
    ("exp", np.exp(2.0)),
    ("log", 0.5),
    ("log1p", 1.0 / 3.0),
    ("sqrt", 0.5 / np.sqrt(2.0)),
    ("cbrt", 1.0 / (3.0 * np.cbrt(2.0)**2)),
    ("sin", np.cos(2.0)),
    ("cos", -np.sin(2.0)),
    ("arctan", 0.2),
])
def test_intrinsic(name, derivative):
    """Check values, non-unit seeds, and scalar specializations."""
    func = getattr(dual, name)

    @njit
    def evaluate(x):
        return func(x)

    result = evaluate(dual.Dual(2.0, 3.0))
    expected = getattr(np, name)(2.0)
    assert result.value == pytest.approx(expected)
    assert result.deriv == pytest.approx(3.0 * derivative)
    assert evaluate(2.0) == pytest.approx(expected)


@pytest.mark.parametrize("op, expected", [
    (operator.add, ((6.0, 3.0), (6.0, 3.0), (6.0, 8.0))),
    (operator.sub, ((-2.0, 3.0), (2.0, -3.0), (-2.0, -2.0))),
    (operator.mul, ((8.0, 12.0), (8.0, 12.0), (8.0, 22.0))),
    (operator.truediv, ((0.5, 0.75), (2.0, -3.0), (0.5, 0.125))),
])
@pytest.mark.parametrize("scalar", [4, 4.0])
def test_arithmetic(op, expected, scalar):
    """Check both scalar operand orders and independently seeded duals."""
    @njit
    def evaluate(a, b):
        return op(a, b)

    x = dual.Dual(2.0, 3.0)
    y = dual.Dual(4.0, 5.0)
    for operands, (value, deriv) in zip(((x, scalar), (scalar, x), (x, y)), expected):
        result = evaluate(*operands)
        assert result.value == pytest.approx(value)
        assert result.deriv == pytest.approx(deriv)


@pytest.mark.parametrize("op", [operator.lt, operator.le, operator.gt,
                               operator.ge, operator.eq, operator.ne])
@pytest.mark.parametrize("value", [1, 2, 3, 1.0, 2.0, 3.0])
def test_comparison(op, value):
    """Compiled comparisons use values for both operand orders and ties."""
    @njit
    def evaluate(a, b):
        return op(a, b)

    x = dual.Dual(2.0, 3.0)
    y = dual.Dual(value, -7.0)
    assert evaluate(x, value) == op(2.0, value)
    assert evaluate(value, x) == op(value, 2.0)
    assert evaluate(x, y) == op(2.0, value)


def test_powers():
    """Propagate constant and variable exponents, including zero powers."""
    @njit
    def evaluate(x, y):
        return x**y

    x = dual.Dual(2.0, 3.0)
    result = evaluate(x, dual.Dual(4.0, 5.0))
    assert result.value == 16.0
    assert result.deriv == pytest.approx(16.0 * (5.0 * np.log(2.0) + 6.0))
    result = evaluate(x, 4.0)
    assert result.value == 16.0
    assert result.deriv == 96.0
    result = evaluate(dual.Dual(0.0, 3.0), 0)
    assert result.value == 1.0
    assert result.deriv == 0.0
    with pytest.raises(ValueError, match="positive base"):
        evaluate(dual.Dual(-1.0), dual.Dual(2.0))


@pytest.mark.parametrize("func", [min, max])
@pytest.mark.parametrize("bound", [1, 2, 3, 1.0, 2.0, 3.0])
def test_clipping(func, bound):
    """Clipping preserves the chosen seed, with zero for scalar bounds."""
    @njit
    def evaluate(a, b):
        return func(a, b)

    x = dual.Dual(2.0, 3.0)
    y = dual.Dual(bound, -7.0)
    for a, b in ((x, bound), (bound, x), (x, y)):
        av = a.value if isinstance(a, dual.Dual) else a
        bv = b.value if isinstance(b, dual.Dual) else b
        # The first operand wins ties, including when it is a scalar.
        selected = a if func(av, bv) == av else b
        expected_deriv = selected.deriv if isinstance(selected, dual.Dual) else 0.0
        result = evaluate(a, b)
        assert result.value == func(av, bv)
        assert result.deriv == expected_deriv


def test_temperature_dependence():
    """A weak-screening-like enhancement has df/dT = -3 h exp(h)/(2T)."""
    @njit
    def enhancement(temp):
        h = 2.0 * dual.sqrt(1.0 / temp) / temp
        return dual.exp(max(h, 0.0))

    temp = 4.0
    h = 2.0 * temp**(-1.5)
    result = enhancement(dual.Dual(temp, 1.0))
    assert result.value == pytest.approx(np.exp(h))
    assert result.deriv == pytest.approx(-1.5 * h * np.exp(h) / temp)
