Scalar automatic differentiation
================================

``pynucastro.autodiff.dual`` provides experimental forward-mode automatic
differentiation with one derivative per value. Seed the independent variable
with ``Dual(value, 1.0)`` (the default derivative is one), and constants with
``Dual(value, 0.0)``. For example::

    from numba import njit
    from pynucastro.autodiff.dual import Dual

    @njit
    def polynomial(x):
        return x * x + 3.0 * x - 1.0

    result = polynomial(Dual(1.3))
    # result.value = 4.59; result.deriv = 5.6

``Dual`` is an immutable namedtuple with fields ``value`` and ``deriv``.
Numba operator overloads support dual/dual arithmetic and real scalars on
either side of addition, subtraction, multiplication, and division.
Powers support a dual base and a scalar or dual exponent; a dual exponent
requires a positive base. Compiled operations use a native tuple structure
rather than a heap-allocated jitclass. Python arithmetic uses the same
derivative formulas.

Use the module's ``exp``, ``log``, ``log1p``, ``sqrt``, ``cbrt``, ``sin``,
``cos``, and ``arctan`` helpers when an argument may be dual. Direct NumPy
ufunc calls do not support dual arguments. Comparisons use only the value.
The derivative of ``abs`` is defined as zero at zero.

Compiled two-argument ``min`` and ``max`` return a dual with a floating-point
value, even when a scalar bound wins; that bound has zero derivative
components with the same seed kind as the dual operand.
Ties select the first operand. Python's built-in ``min`` and ``max`` return
the selected operand itself, which may be a scalar. In other compiled
branches, keep the return types compatible explicitly.

For a temperature seed of one, ``deriv`` measures the output change per
kelvin when temperature is expressed in kelvin.

Two independent derivatives
---------------------------

Use ``DerivativePair`` seeds to differentiate with respect to two variables
in one evaluation. Both components are scalar fields in a namedtuple, so
this does not require a NumPy derivative array::

    from pynucastro.autodiff.derivative_pair import DerivativePair

    x = Dual(2.0, DerivativePair(1.0, 0.0))
    y = Dual(3.0, DerivativePair(0.0, 1.0))
    result = x * y
    # result.value = 6; result.deriv.d1 = 3; result.deriv.d2 = 2

All dual operands in an expression must use the same kind of seed.
Constants may be ordinary real scalars. ``constant_like(x, value)`` creates
a constant with the same scalar or dual representation as ``x``, including
zero derivative components; this is useful for compiled clipping branches.

Screening composition derivatives
---------------------------------

All five screening prescriptions support derivatives with respect to the
composition moments

.. math::

   Y_e = \sum_i Z_i Y_i, \qquad S_2 = \sum_i Z_i^2 Y_i,

at fixed temperature and density. Here :math:`Y_i = X_i/A_i` is the molar
abundance used by the network. These are unconstrained abundance derivatives:
individual abundances are varied without renormalizing the composition.

Construct a derivative-enabled state using arrays of abundances and charges::

    import numpy as np
    from pynucastro.autodiff import dual
    from pynucastro.screening import (
        make_plasma_state_with_derivatives, make_screen_factors, chugunov_2009,
    )

    Y = np.array([0.125, 1.0 / 24.0])  # He4 and C12, equal mass fractions
    Z = np.array([2, 6])
    state = make_plasma_state_with_derivatives(1.e8, 1.e6, Y, Z)
    pair = make_screen_factors("he4", "c12")
    h = chugunov_2009(state, pair)

    dh_dY = Z * h.deriv.d1 + Z**2 * h.deriv.d2
    factor = dual.exp(h)
    dfactor_dY = Z * factor.deriv.d1 + Z**2 * factor.deriv.d2

The existing screening routines return the **logarithm** of the enhancement,
:math:`h = \log f`, so ``h.deriv.d1`` is :math:`\partial h/\partial Y_e`
and ``h.deriv.d2`` is :math:`\partial h/\partial S_2`. The species derivatives
follow from

.. math::

   \frac{\partial h}{\partial Y_i}
   = Z_i\frac{\partial h}{\partial Y_e}
     + Z_i^2\frac{\partial h}{\partial S_2},
   \qquad
   \frac{\partial f}{\partial Y_i} = f\frac{\partial h}{\partial Y_i}.

The state is reusable across ion pairs, and array-based state construction
and screening can run together inside ``@njit`` functions. For a molar-abundance
dictionary, use ``make_plasma_state(T, rho, molar, compute_derivatives=True)``.
The ordinary ``PlasmaState`` constructor and default ``make_plasma_state``
behavior remain float-only. The derivative-enabled state is an immutable
namedtuple with the same field names and dual-valued composition factors.

State construction holds the total molar abundance fixed when seeding the
two moments; that total cancels out of the supported screening expressions.
Temperature and density are not differentiated by this constructor.
Clipped quantities have zero derivatives on their constant branches, and
comparisons select branches using primal values. Derivatives at a branch
boundary follow the selected branch, rather than a finite difference across
that boundary.

This exposes screening derivatives for use in a composition Jacobian;
generated network Jacobians do not yet include these new terms automatically.
