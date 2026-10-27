"""Two independent forward-mode derivative components, stored without arrays."""

import operator
from collections import namedtuple

from pynucastro.numba_util import overload, types


class DerivativePair(namedtuple("DerivativePairTuple", ("d1", "d2"))):
    """Store derivatives with respect to two independent variables.

    Parameters
    ----------
    d1 : float
        Derivative with respect to the first variable.
    d2 : float
        Derivative with respect to the second variable.

    """

    __slots__ = ()
    # NumPy scalars must defer to our reflected arithmetic, not build arrays.
    __array_priority__ = 1000

    def __add__(self, other):
        return DerivativePair(self.d1 + other.d1, self.d2 + other.d2)

    def __sub__(self, other):
        return DerivativePair(self.d1 - other.d1, self.d2 - other.d2)

    def __neg__(self):
        return DerivativePair(-self.d1, -self.d2)

    def __mul__(self, other):
        return DerivativePair(self.d1 * other, self.d2 * other)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        return DerivativePair(self.d1 / other, self.d2 / other)


def _is_pair(t):
    return isinstance(t, types.BaseNamedTuple) and t.instance_class is DerivativePair


_REAL_TYPES = (types.Integer, types.Float, types.Boolean)


@overload(operator.add, target="cpu")
def _ol_add(a, b):
    if _is_pair(a) and _is_pair(b):
        def impl(a, b):
            return DerivativePair(a.d1 + b.d1, a.d2 + b.d2)
        return impl
    return None


@overload(operator.sub, target="cpu")
def _ol_sub(a, b):
    if _is_pair(a) and _is_pair(b):
        def impl(a, b):
            return DerivativePair(a.d1 - b.d1, a.d2 - b.d2)
        return impl
    return None


@overload(operator.neg, target="cpu")
def _ol_neg(x):
    if _is_pair(x):
        def impl(x):
            return DerivativePair(-x.d1, -x.d2)
        return impl
    return None


@overload(operator.mul, target="cpu")
def _ol_mul(a, b):
    if _is_pair(a) and isinstance(b, _REAL_TYPES):
        def impl(a, b):
            return DerivativePair(a.d1 * b, a.d2 * b)
        return impl
    if isinstance(a, _REAL_TYPES) and _is_pair(b):
        def reflected_impl(a, b):
            return DerivativePair(a * b.d1, a * b.d2)
        return reflected_impl
    return None


@overload(operator.truediv, target="cpu")
def _ol_div(a, b):
    if _is_pair(a) and isinstance(b, _REAL_TYPES):
        def impl(a, b):
            return DerivativePair(a.d1 / b, a.d2 / b)
        return impl
    return None


def zero_derivative(deriv):
    """Return a zero seed of the same kind as deriv."""
    if isinstance(deriv, DerivativePair):
        return DerivativePair(0.0, 0.0)
    return 0.0


@overload(zero_derivative)
def _ol_zero_derivative(deriv):
    if _is_pair(deriv):
        def pair_impl(deriv):  # pylint: disable=unused-argument
            return DerivativePair(0.0, 0.0)
        return pair_impl
    if isinstance(deriv, _REAL_TYPES):
        def scalar_impl(deriv):  # pylint: disable=unused-argument
            return 0.0
        return scalar_impl
    return None
