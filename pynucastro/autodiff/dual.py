"""Forward-mode derivatives of scalar values, usable from Python and Numba.

Seed the independent variable with ``Dual(value, 1.0)`` and constants with
``Dual(value, 0.0)``. Use this module's math helpers in place of NumPy ufuncs
when an argument may be dual. CPU Numba overloads support scalar-left
arithmetic, value comparisons, and two-argument ``min``/``max``. Compiled
``min``/``max`` return a Dual if either operand is dual, with zero derivative
when a scalar bound is selected. Ties select the first operand, as in Python.
Use DerivativePair seeds to propagate two derivatives simultaneously. All dual
operands in an expression must use the same seed kind. Other branches must
likewise keep dual-valued variables consistently typed.
"""

import operator
from collections import namedtuple

import numpy as np

from pynucastro.autodiff.derivative_pair import zero_derivative
from pynucastro.numba_util import overload, types


class Dual(namedtuple("DualTuple", ("value", "deriv"), defaults=(1.0,))):
    """Carry a quantity that can propagate its derivative.

    A scalar seed carries one derivative; a DerivativePair seed carries
    two independent derivatives. Instances are immutable namedtuples,
    lowered by Numba to native tuple values.

    Parameters
    ----------
    value : float
        Scalar value.
    deriv : float or DerivativePair, optional
        Derivative seed, defaulting to one. Use the same seed kind for all
        dual operands in an expression.

    """

    __slots__ = ()
    # NumPy scalars must defer to our reflected arithmetic, not build arrays.
    __array_priority__ = 1000
    __hash__ = None

    def __abs__(self):
        if self.value > 0:
            return Dual(self.value, self.deriv)
        if self.value < 0:
            return Dual(-self.value, -self.deriv)
        return Dual(abs(self.value), zero_derivative(self.deriv))

    def __lt__(self, other):
        if isinstance(other, Dual):
            return self.value < other.value
        return self.value < other

    def __le__(self, other):
        if isinstance(other, Dual):
            return self.value <= other.value
        return self.value <= other

    def __gt__(self, other):
        if isinstance(other, Dual):
            return self.value > other.value
        return self.value > other

    def __ge__(self, other):
        if isinstance(other, Dual):
            return self.value >= other.value
        return self.value >= other

    def __eq__(self, other):
        if isinstance(other, Dual):
            return self.value == other.value
        return self.value == other

    def __ne__(self, other):
        return not self == other

    def __add__(self, other):
        if isinstance(other, Dual):
            return Dual(self.value + other.value, self.deriv + other.deriv)
        return Dual(self.value + other, self.deriv)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Dual):
            return Dual(self.value - other.value, self.deriv - other.deriv)
        return Dual(self.value - other, self.deriv)

    def __rsub__(self, other):
        if isinstance(other, Dual):
            return other.__sub__(self)
        return Dual(other - self.value, -self.deriv)

    def __neg__(self):
        return Dual(-self.value, -self.deriv)

    def __mul__(self, other):
        if isinstance(other, Dual):
            # Product rule, using the original operand values.
            return Dual(self.value * other.value,
                        self.value * other.deriv + other.value * self.deriv)
        return Dual(self.value * other, self.deriv * other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Dual):
            # Quotient rule, using the original operand values.
            return Dual(self.value / other.value,
                        -self.value / other.value**2 * other.deriv + self.deriv / other.value)
        return Dual(self.value / other, self.deriv / other)

    def __rtruediv__(self, other):
        if isinstance(other, Dual):
            return other.__truediv__(self)
        return Dual(other / self.value, -other / self.value**2 * self.deriv)

    def __pow__(self, y):
        if isinstance(y, Dual):
            if self.value <= 0:
                raise ValueError("Dual exponents require a positive base")
            value = self.value**y.value
            deriv = value * (y.deriv * np.log(self.value) +
                             y.value * self.deriv / self.value)
            return Dual(value, deriv)

        if y == 0:
            return Dual(1.0, zero_derivative(self.deriv))
        return Dual(self.value**y,
                    y * self.value**(y - 1.0) * self.deriv)


def exp(x):
    """Compute exp(x)."""
    if isinstance(x, Dual):
        value = np.exp(x.value)
        return Dual(value, value * x.deriv)
    return np.exp(x)


@overload(exp)
def _ol_exp(x):
    if _is_dual(x):
        def dual_impl(x):
            value = np.exp(x.value)
            return Dual(value, value * x.deriv)
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x):
            return np.exp(x)
        return scalar_impl
    return None


def log(x):
    """Compute log(x), the natural logarithm."""
    if isinstance(x, Dual):
        value = np.log(x.value)
        return Dual(value, x.deriv / x.value)
    return np.log(x)


@overload(log)
def _ol_log(x):
    if _is_dual(x):
        def dual_impl(x):
            value = np.log(x.value)
            return Dual(value, x.deriv / x.value)
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x):
            return np.log(x)
        return scalar_impl
    return None


def log1p(x):
    """Compute log(1 + x), using the natural logarithm."""
    if isinstance(x, Dual):
        value = np.log1p(x.value)
        return Dual(value, x.deriv / (1.0 + x.value))
    return np.log1p(x)


@overload(log1p)
def _ol_log1p(x):
    if _is_dual(x):
        def dual_impl(x):
            value = np.log1p(x.value)
            return Dual(value, x.deriv / (1.0 + x.value))
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x):
            return np.log1p(x)
        return scalar_impl
    return None


def sqrt(x):
    """Compute sqrt(x)."""
    if isinstance(x, Dual):
        value = np.sqrt(x.value)
        return Dual(value, 0.5 * x.deriv / value)
    return np.sqrt(x)


@overload(sqrt)
def _ol_sqrt(x):
    if _is_dual(x):
        def dual_impl(x):
            value = np.sqrt(x.value)
            return Dual(value, 0.5 * x.deriv / value)
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x):
            return np.sqrt(x)
        return scalar_impl
    return None


def cbrt(x):
    """Compute cbrt(x)."""
    if isinstance(x, Dual):
        value = np.cbrt(x.value)
        return Dual(value, (1.0 / 3.0) * x.deriv / value**2)
    return np.cbrt(x)


@overload(cbrt)
def _ol_cbrt(x):
    if _is_dual(x):
        def dual_impl(x):
            value = np.cbrt(x.value)
            return Dual(value, 1.0 / 3.0 * x.deriv / value ** 2)
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x):
            return np.cbrt(x)
        return scalar_impl
    return None


def cos(x):
    """Compute cos(x), with x in radians."""
    if isinstance(x, Dual):
        value = np.cos(x.value)
        return Dual(value, -np.sin(x.value) * x.deriv)
    return np.cos(x)


@overload(cos)
def _ol_cos(x):
    if _is_dual(x):
        def dual_impl(x):
            value = np.cos(x.value)
            return Dual(value, -np.sin(x.value) * x.deriv)
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x):
            return np.cos(x)
        return scalar_impl
    return None


def sin(x):
    """Compute sin(x), with x in radians."""
    if isinstance(x, Dual):
        value = np.sin(x.value)
        return Dual(value, np.cos(x.value) * x.deriv)
    return np.sin(x)


@overload(sin)
def _ol_sin(x):
    if _is_dual(x):
        def dual_impl(x):
            value = np.sin(x.value)
            return Dual(value, np.cos(x.value) * x.deriv)
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x):
            return np.sin(x)
        return scalar_impl
    return None


def arctan(x):
    """Compute arctan(x)."""
    if isinstance(x, Dual):
        value = np.arctan(x.value)
        return Dual(value, x.deriv / (1.0 + x.value**2))
    return np.arctan(x)


@overload(arctan)
def _ol_arctan(x):
    if _is_dual(x):
        def dual_impl(x):
            value = np.arctan(x.value)
            return Dual(value, x.deriv / (1.0 + x.value ** 2))
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x):
            return np.arctan(x)
        return scalar_impl
    return None


_REAL_TYPES = (types.Integer, types.Float, types.Boolean)


def _is_dual(t):
    """Recognize only this namedtuple class, not unrelated namedtuples."""
    return isinstance(t, types.BaseNamedTuple) and t.instance_class is Dual


@overload(operator.add, target="cpu")
def _ol_add(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_impl(a, b):
            return Dual(a.value + b.value, a.deriv + b.deriv)
        return dual_impl
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def scalar_impl(a, b):
            return Dual(a.value + b, a.deriv)
        return scalar_impl
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def reflected_impl(a, b):
            return Dual(b.value + a, b.deriv)
        return reflected_impl
    return None


@overload(operator.sub, target="cpu")
def _ol_sub(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_impl(a, b):
            return Dual(a.value - b.value, a.deriv - b.deriv)
        return dual_impl
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def scalar_impl(a, b):
            return Dual(a.value - b, a.deriv)
        return scalar_impl
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def reflected_impl(a, b):
            return Dual(a - b.value, -b.deriv)
        return reflected_impl
    return None


@overload(operator.mul, target="cpu")
def _ol_mul(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_impl(a, b):
            return Dual(a.value * b.value, a.value * b.deriv + b.value * a.deriv)
        return dual_impl
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def scalar_impl(a, b):
            return Dual(a.value * b, a.deriv * b)
        return scalar_impl
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def reflected_impl(a, b):
            return Dual(b.value * a, b.deriv * a)
        return reflected_impl
    return None


@overload(operator.truediv, target="cpu")
def _ol_truediv(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_impl(a, b):
            return Dual(a.value / b.value, -a.value / b.value ** 2 * b.deriv + a.deriv / b.value)
        return dual_impl
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def scalar_impl(a, b):
            return Dual(a.value / b, a.deriv / b)
        return scalar_impl
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def reflected_impl(a, b):
            return Dual(a / b.value, -a / b.value ** 2 * b.deriv)
        return reflected_impl
    return None


@overload(operator.pow, target="cpu")
def _ol_pow(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_impl(a, b):
            if a.value <= 0:
                raise ValueError("Dual exponents require a positive base")
            value = a.value ** b.value
            deriv = value * (b.deriv * np.log(a.value) + b.value * a.deriv / a.value)
            return Dual(value, deriv)
        return dual_impl
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def scalar_impl(a, b):
            if b == 0:
                return Dual(1.0, zero_derivative(a.deriv))
            return Dual(a.value ** b, b * a.value ** (b - 1.0) * a.deriv)
        return scalar_impl
    return None


@overload(operator.neg, target="cpu")
def _ol_neg(x):
    if _is_dual(x):
        def impl(x):
            return Dual(-x.value, -x.deriv)
        return impl
    return None


@overload(abs, target="cpu")
def _ol_abs(x):
    if _is_dual(x):
        def impl(x):
            if x.value > 0:
                return Dual(x.value, x.deriv)
            if x.value < 0:
                return Dual(-x.value, -x.deriv)
            return Dual(abs(x.value), zero_derivative(x.deriv))
        return impl
    return None


# Override tuple comparisons: only values, not derivative seeds, are compared.
@overload(operator.lt, target="cpu")
def _ol_lt(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_dual(a, b):
            return a.value < b.value
        return dual_dual
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def dual_scalar(a, b):
            return a.value < b
        return dual_scalar
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def scalar_dual(a, b):
            return a < b.value
        return scalar_dual
    return None


@overload(operator.le, target="cpu")
def _ol_le(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_dual(a, b):
            return a.value <= b.value
        return dual_dual
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def dual_scalar(a, b):
            return a.value <= b
        return dual_scalar
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def scalar_dual(a, b):
            return a <= b.value
        return scalar_dual
    return None


@overload(operator.gt, target="cpu")
def _ol_gt(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_dual(a, b):
            return a.value > b.value
        return dual_dual
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def dual_scalar(a, b):
            return a.value > b
        return dual_scalar
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def scalar_dual(a, b):
            return a > b.value
        return scalar_dual
    return None


@overload(operator.ge, target="cpu")
def _ol_ge(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_dual(a, b):
            return a.value >= b.value
        return dual_dual
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def dual_scalar(a, b):
            return a.value >= b
        return dual_scalar
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def scalar_dual(a, b):
            return a >= b.value
        return scalar_dual
    return None


@overload(operator.eq, target="cpu")
def _ol_eq(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_dual(a, b):
            return a.value == b.value
        return dual_dual
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def dual_scalar(a, b):
            return a.value == b
        return dual_scalar
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def scalar_dual(a, b):
            return a == b.value
        return scalar_dual
    return None


@overload(operator.ne, target="cpu")
def _ol_ne(a, b):
    if _is_dual(a) and _is_dual(b):
        def dual_dual(a, b):
            return a.value != b.value
        return dual_dual
    if _is_dual(a) and isinstance(b, _REAL_TYPES):
        def dual_scalar(a, b):
            return a.value != b
        return dual_scalar
    if isinstance(a, _REAL_TYPES) and _is_dual(b):
        def scalar_dual(a, b):
            return a != b.value
        return scalar_dual
    return None


def _register_minmax(func, compare):
    """Keep both branches of binary clipping dual-valued."""
    @overload(func)
    def minmax(a, b):
        if _is_dual(a) and _is_dual(b):
            def dual_dual(a, b):
                if compare(b.value, a.value):
                    return Dual(float(b.value), b.deriv * 1.0)
                return Dual(float(a.value), a.deriv * 1.0)
            return dual_dual
        if _is_dual(a) and isinstance(b, _REAL_TYPES):
            def dual_scalar(a, b):
                if compare(b, a.value):
                    return Dual(float(b), zero_derivative(a.deriv))
                return Dual(float(a.value), a.deriv * 1.0)
            return dual_scalar
        if isinstance(a, _REAL_TYPES) and _is_dual(b):
            def scalar_dual(a, b):
                if compare(b.value, a):
                    return Dual(float(b.value), b.deriv * 1.0)
                return Dual(float(a), zero_derivative(b.deriv))
            return scalar_dual
        return None


_register_minmax(min, operator.lt)
_register_minmax(max, operator.gt)


def constant_like(x, value):
    """Give a constant the same scalar or dual representation as x.

    Parameters
    ----------
    x : float or Dual
        Determines the result type and derivative seed shape.
    value : float
        Constant value, whose derivatives are zero.

    Returns
    -------
    float or Dual
        The constant with zero derivatives when x is dual.

    """
    if isinstance(x, Dual):
        return Dual(float(value), zero_derivative(x.deriv))
    return float(value)


@overload(constant_like)
def _ol_constant_like(x, value):
    if _is_dual(x):
        def dual_impl(x, value):
            return Dual(float(value), zero_derivative(x.deriv))
        return dual_impl
    if isinstance(x, _REAL_TYPES):
        def scalar_impl(x, value):  # pylint: disable=unused-argument
            return float(value)
        return scalar_impl
    return None
