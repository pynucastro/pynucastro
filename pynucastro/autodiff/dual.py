"""Scalar forward-mode derivatives, usable from Python and Numba.

Seed the independent variable with ``Dual(value, 1.0)`` and constants with
``Dual(value, 0.0)``. Use this module's math helpers in place of NumPy ufuncs
when an argument may be dual. CPU Numba overloads support scalar-left
arithmetic, value comparisons, and two-argument ``min``/``max``. Compiled
``min``/``max`` return a Dual if either operand is dual, with zero derivative
when a scalar bound is selected. Ties select the first operand, as in Python.
Other branches must likewise keep dual-valued variables consistently typed.
"""

import operator

import numpy as np
from numba import float64, njit, types
from numba.experimental import jitclass
from numba.extending import overload


@jitclass([("value", float64), ("deriv", float64)])
class Dual:  # noqa: PLW1641 -- mutable duals are intentionally unhashable
    """Carry a quantity that can propagate its derivative.

    We carry a scalar derivative, so that means we are always
    differentiating with respect to the same quantity.

    """

    def __init__(self, x, deriv=1):
        self.value = x
        self.deriv = deriv

    def __abs__(self):
        if self.value > 0:
            return Dual(self.value, self.deriv)
        if self.value < 0:
            return Dual(-self.value, -self.deriv)
        return Dual(abs(self.value), 0.0)

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
            return Dual(1.0, 0.0)
        return Dual(self.value**y,
                    y * self.value**(y - 1.0) * self.deriv)


@njit
def exp(x):
    """Compute exp(x)."""
    if isinstance(x, Dual):
        value = np.exp(x.value)
        return Dual(value, value * x.deriv)
    return np.exp(x)


@njit
def log(x):
    """Compute log(x), the natural logarithm."""
    if isinstance(x, Dual):
        value = np.log(x.value)
        return Dual(value, x.deriv / x.value)
    return np.log(x)


@njit
def log1p(x):
    """Compute log(1 + x), using the natural logarithm."""
    if isinstance(x, Dual):
        value = np.log1p(x.value)
        return Dual(value, x.deriv / (1.0 + x.value))
    return np.log1p(x)


@njit
def sqrt(x):
    """Compute sqrt(x)."""
    if isinstance(x, Dual):
        value = np.sqrt(x.value)
        return Dual(value, 0.5 * x.deriv / value)
    return np.sqrt(x)


@njit
def cbrt(x):
    """Compute cbrt(x)."""
    if isinstance(x, Dual):
        value = np.cbrt(x.value)
        return Dual(value, (1.0 / 3.0) * x.deriv / value**2)
    return np.cbrt(x)


@njit
def cos(x):
    """Compute cos(x), with x in radians."""
    if isinstance(x, Dual):
        value = np.cos(x.value)
        return Dual(value, -np.sin(x.value) * x.deriv)
    return np.cos(x)


@njit
def sin(x):
    """Compute sin(x), with x in radians."""
    if isinstance(x, Dual):
        value = np.sin(x.value)
        return Dual(value, np.cos(x.value) * x.deriv)
    return np.sin(x)


@njit
def arctan(x):
    """Compute arctan(x)."""
    if isinstance(x, Dual):
        value = np.arctan(x.value)
        return Dual(value, x.deriv / (1.0 + x.value**2))
    return np.arctan(x)


# jitclass adds class_type dynamically.
_DUAL_TYPE = Dual.class_type.instance_type  # pylint: disable=no-member
_REAL_TYPES = (types.Integer, types.Float, types.Boolean)


# Call the reflected methods directly: converting the scalar to Dual here
# can recursively compile the dual/dual operator before type pruning.
@overload(operator.add)
def _scalar_add(a, b):
    if isinstance(a, _REAL_TYPES) and b == _DUAL_TYPE:
        def impl(a, b):
            return b.__radd__(a)  # pylint: disable=unnecessary-dunder-call
        return impl
    return None


@overload(operator.sub)
def _scalar_sub(a, b):
    if isinstance(a, _REAL_TYPES) and b == _DUAL_TYPE:
        def impl(a, b):
            return b.__rsub__(a)  # pylint: disable=unnecessary-dunder-call
        return impl
    return None


@overload(operator.mul)
def _scalar_mul(a, b):
    if isinstance(a, _REAL_TYPES) and b == _DUAL_TYPE:
        def impl(a, b):
            return b.__rmul__(a)  # pylint: disable=unnecessary-dunder-call
        return impl
    return None


@overload(operator.truediv)
def _scalar_div(a, b):
    if isinstance(a, _REAL_TYPES) and b == _DUAL_TYPE:
        def impl(a, b):
            return b.__rtruediv__(a)  # pylint: disable=unnecessary-dunder-call
        return impl
    return None


def _register_comparison(op):
    """Compare values independently of the derivative seeds."""
    # CPU specificity gives these overloads priority over Numba's generic
    # jitclass comparisons, which assume both operands have jit_methods.
    @overload(op, target="cpu")
    def compare(a, b):
        if a == _DUAL_TYPE and b == _DUAL_TYPE:
            def dual_dual(a, b):
                return op(a.value, b.value)
            return dual_dual
        if a == _DUAL_TYPE and isinstance(b, _REAL_TYPES):
            def dual_scalar(a, b):
                return op(a.value, b)
            return dual_scalar
        if isinstance(a, _REAL_TYPES) and b == _DUAL_TYPE:
            def scalar_dual(a, b):
                return op(a, b.value)
            return scalar_dual
        return None


for _op in (operator.lt, operator.le, operator.gt, operator.ge,
            operator.eq, operator.ne):
    _register_comparison(_op)


def _register_minmax(func, compare):
    """Keep both branches of binary clipping dual-valued."""
    @overload(func)
    def minmax(a, b):
        if a == _DUAL_TYPE and b == _DUAL_TYPE:
            def dual_dual(a, b):
                if compare(b.value, a.value):
                    return b
                return a
            return dual_dual
        if a == _DUAL_TYPE and isinstance(b, _REAL_TYPES):
            def dual_scalar(a, b):
                if compare(b, a.value):
                    return Dual(b, 0.0)
                return a
            return dual_scalar
        if isinstance(a, _REAL_TYPES) and b == _DUAL_TYPE:
            def scalar_dual(a, b):
                if compare(b.value, a):
                    return b
                return Dual(a, 0.0)
            return scalar_dual
        return None


_register_minmax(min, operator.lt)
_register_minmax(max, operator.gt)
