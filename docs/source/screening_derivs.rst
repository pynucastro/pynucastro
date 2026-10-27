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

Below we show an example of computing the derivatives of the screening term.
We'll consider a plasma consisting only of ${}^4\mathrm{He}$,
${}^{12}\mathrm{C}$, and ${}^{16}\mathrm{O}$, and look at the
screening factor for
${}^{12}\mathrm{C}(\alpha,\gamma){}^{16}\mathrm{O}$.

Plasma state
------------

First we setup a version of the plasma state that can be
differentiated, using :py:func:`make_plasma_state
<pynucastro.screening.screen.make_plasma_state>`.

.. important::

   When computing derivatives, the plasma state is stored as a
   ``namedtuple``.  It has the same field names as the float-only
   :py:obj:`PlasmaState <pynucastro.screening.screen.PlasmaState>`.

Here's the construction of the plasma state:

.. code:: python

    from pynucastro.screening import (make_plasma_state,
                                      make_screen_factors,
                                      chugunov_2009)

    comp = pyna.Composition(["he4", "c12", "o16"], init="uniform")
    T = 1.e8  # K
    rho = 1.e6  # g/cm^3

    pstate = make_plasma_state(T, rho, comp, compute_derivatives=True)

.. note::

   This will create a :py:obj:`Dual <pynucastro.autodiff.dual.Dual>` seeded with a
   :py:obj:`DerivativePair <pynucastro.autodiff.derivative_pair.DerivativePair>`,
   with ``d1`` corresponding to the derivative with respect to $Y_e$ and ``d2``
   corresponding to the derivative with respect to $S_2$.

.. tip::

   The plasma state, ``pstate``, only needs to be computed once, and can then
   be used for all screening pairs, as long as density, temperature, and composition
   are unchanged.

.. important::

   When working inside a ``@njit`` function, the explicit Numpy array
   interface, :py:func:`make_plasma_state_with_derivatives
   <pynucastro.screening.screen.make_plasma_state_with_derivatives>`
   should be used.  This higher-level ``make_plasma_state`` interface
   is only intended to be used outside of an ``@njit`` function.

Screening pair
--------------

Next we create the :py:obj:`ScreenFactors <pynucastro.screening.screen.ScreenFactors>`,
which describes the screening pair.  Note that this is fixed, so there are no derivatives
with respect to the nuclei that are in the pair.

.. code:: python

   pair = make_screen_factors("he4", "c12")


Screening factor
----------------

Finally, we can compute the screening factor and the composition derivatives.

.. note::

   The screening routines return the **logarithm** of the enhancement,
   :math:`h = \log f`, so ``h.deriv.d1`` is :math:`\partial h/\partial Y_e`
   and ``h.deriv.d2`` is :math:`\partial h/\partial S_2`.

   The species derivatives follow from

   .. math::

      \frac{\partial h}{\partial Y_i}
      &= Z_i\frac{\partial h}{\partial Y_e}
        + Z_i^2\frac{\partial h}{\partial S_2} \\
      \frac{\partial f}{\partial Y_i} &= f\frac{\partial h}{\partial Y_i}.


.. code:: python

    h = chugunov_2009(state, pair)

    dh_dY = Z * h.deriv.d1 + Z**2 * h.deriv.d2
    factor = dual.exp(h)
    dfactor_dY = Z * factor.deriv.d1 + Z**2 * factor.deriv.d2


.. note::

   Although intermediate plasma quantities depend on the mean
   molecular weight, $\bar{A}$, this dependence cancels in the
   supported screening expressions.  We therefore hold $\bar{A}$ fixed
   when seeding derivatives with respect to $Y_e$ and $S_2$.

