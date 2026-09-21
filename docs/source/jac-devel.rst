********************
Form of the Jacobian
********************

The exported networks include the Jacobian needed to use implicit integration
methods.  When energy evolution is included, the full form of the Jacobian looks like:

.. math::

   {\bf J} = \left ( \begin{array}{ccc|c}
                     \ddots & \vdots & \iddots & \vdots \\
                     \cdots & \partial \dot{Y}_i/\partial Y_j & \cdots & \partial \dot{Y}_i / \partial e  \\
                     \iddots & \vdots & \ddots & \vdots \\
                     \hline
                     \cdots & \partial \dot{e} /\partial Y_j & \cdots & \partial \dot{e}/\partial e
                     \end{array} \right )

Here we describe how and where each of these contributions are computed.

:math:`\partial \dot{Y}_i/\partial Y_j`
=======================================

For most rates, this is simply the derivative of the composition factors that appear explicitly
in the reactive flux.  There are a few instances where the composition is internal to the rate
itself that we also need to capture.

Basic rates
-----------

For a 2-body strong-mediated reaction rate, the reactive flux has the form:

$$F_{AB} = \rho Y(A) Y(B) \lambda_{AB}$$

and the contributions to :math:`\partial \dot{Y}_i/\partial Y_j` would be $\partial F_{AB}/\partial Y_j$,

.. math::

   \frac{\partial F_{AB}}{\partial Y_j} = \rho (Y(B) \delta_{Aj} + Y(A) \delta_{Bj}) \lambda_{AB}

For weak rates, the explicit composition dependence of the parent nucleus is also computed this way.  E.g.,
for a decay with a parent nucleus $P$, the flux is:

$$F_{P,\mathrm{weak}} = Y(P) \lambda_{P,\mathrm{weak}}$$

and the contribution to the Jacobian is:

.. math::

   \frac{\partial F_{P,\mathrm{weak}}}{\partial Y_j} = \delta_{Pj} \lambda_{P,\mathrm{weak}}


This is computed as:

* ``RateCollection`` : directly via
  :py:meth:`Rate.eval_jacobian_term <pynucastro.rates.rate.Rate.eval_jacobian_term>`
* ``PythonNetwork`` : as a string via
  :py:meth:`Rate.jacobian_string_py <pynucastro.rates.rate.Rate.jacobian_string_py>`
* ``AmrexAstroCxxNetwork`` / ``SimpleCxxNetwork`` : symbolically using the SymPy methods
  in :py:meth:`SympyRates.jacobian_term_symbol <pynucastro.networks.sympy_network_support.SympyRates.jacobian_term_symbol>`



``ApproximateRate``
-------------------

For some rate approximations, the composition appears explicitly in
the effective rate, $\lambda$.  In this case, we need to also need to
compute $\partial\lambda/\partial Y_j$.  An example of such a rate is the
``ApproximateRate`` for $(nn,\gamma)$.

In a rate class, we set ``Rate.rate_comp_dependence = True`` to indicate that
we need to compute this derivative.

The rate class itself will then compute the explicit $\partial\lambda/\partial Y_j$ term
and store it in the python ``RateEval`` class or the C++ ``rate_derivs_t`` struct.

Status of this term:

* ``RateCollection`` : not currently included
* ``PythonNetwork`` : included in the string returned via
  :py:meth:`Rate.jacobian_string_py <pynucastro.rates.rate.Rate.jacobian_string_py>`
* ``AmrexAstroCxxNetwork`` / ``SimpleCxxNetwork`` : stored in
  ``rate_derivs_t`` and used symbolically using the SymPy methods in
  :py:meth:`SympyRates.jacobian_term_symbol
  <pynucastro.networks.sympy_network_support.SympyRates.jacobian_term_symbol>`



Weak-tabulated rates
--------------------


