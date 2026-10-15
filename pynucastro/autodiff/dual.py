from numba import float64
from numba.experimental import jitclass

import numpy as np


@jitclass([("value", float64), ("deriv", float64)])
class Dual:
    """Carry a quantity that can propagate its derivatve.

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
            deriv = value * (y.deriv * np.log(self.value)
                             + y.value * self.deriv / self.value)
            return Dual(value, deriv)

        if y == 0:
            return Dual(1.0, 0.0)
        return Dual(self.value**y,
                    y * self.value**(y - 1.0) * self.deriv)


def exp(x):
    """Compute exp(x)."""
    if isinstance(x, Dual):
        value = np.exp(x.value)
        return Dual(value, value * x.deriv)
    return np.exp(x)


def log(x):
    """Compute log(x), the natural logarithm."""
    if isinstance(x, Dual):
        value = np.log(x.value)
        return Dual(value, x.deriv / x.value)
    return np.log(x)


def log1p(x):
    """Compute log(1 + x), using the natural logarithm."""
    if isinstance(x, Dual):
        value = np.log1p(x.value)
        return Dual(value, x.deriv / (1.0 + x.value))
    return np.log1p(x)


def sqrt(x):
    """Compute sqrt(x)."""
    if isinstance(x, Dual):
        value = np.sqrt(x.value)
        return Dual(value, 0.5 * x.deriv / value)
    return np.sqrt(x)


def cbrt(x):
    """Compute cbrt(x)."""
    if isinstance(x, Dual):
        value = np.cbrt(x.value)
        return Dual(value, (1.0 / 3.0) * x.deriv / value**2)
    return np.cbrt(x)


def cos(x):
    """Compute cos(x), with x in radians."""
    if isinstance(x, Dual):
        value = np.cos(x.value)
        return Dual(value, -np.sin(x.value) * x.deriv)
    return np.cos(x)


def sin(x):
    """Compute sin(x), with x in radians."""
    if isinstance(x, Dual):
        value = np.sin(x.value)
        return Dual(value, np.cos(x.value) * x.deriv)
    return np.sin(x)


def arctan(x):
    """Compute arctan(x)."""
    if isinstance(x, Dual):
        value = np.arctan(x.value)
        return Dual(value, x.deriv / (1.0 + x.value**2))
    return np.arctan(x)
