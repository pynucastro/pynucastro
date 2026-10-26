Screening composition derivatives
=================================

Derivatives of the screening functions with respect to composition
are performed using the :py:mod:`autodiff <pynucastro.autodiff>` module.

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
