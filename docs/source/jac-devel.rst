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

.. note::

   For a :py:obj:`RateCollection
   <pynucastro.networks.rate_collection.RateCollection>`, the Jacobian
   is only used for visualization, via :py:meth:`plot_jacobian
   <pynucastro.networks.rate_collection.RateCollection.plot_jacobian>`
   and for the stiffness assessment in :py:meth:`spectral_radius
   <pynucastro.networks.rate_collection.RateCollection.spectral_radius>`.

.. note::

   For a :py:obj:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>`,
   temperature is integrated instead of specific internal energy, $e$.

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
:py:obj:`ApproximateRate <pynucastro.rates.approximate_rates.ApproximateRate>` for $(nn,\gamma)$.

For this type of rate, the flux is:

$$F_{AB} = \rho Y(A) Y(B) \lambda_{AB}(Y)$$

and the contribution to the Jacobian is then:

.. math::

   \frac{\partial F_{AB}}{\partial Y_j} = \rho (Y(B) \delta_{Aj} + Y(A) \delta_{Bj}) \lambda_{AB} +
          \rho Y(A) Y(B) \frac{\partial \lambda_{AB}}{\partial Y_j}


In a rate class, we set ``Rate.rate_comp_dependence = True`` to indicate that
we need to compute this derivative.

The rate class itself will then compute the explicit $\partial\lambda/\partial Y_j$ term
and store it in the python ``RateEval`` class or the C++ ``rate_derivs_t`` struct.

Status of this term:

* ``RateCollection`` : not currently included
* ``PythonNetwork`` : (to be done) included in the string returned via
  :py:meth:`Rate.jacobian_string_py <pynucastro.rates.rate.Rate.jacobian_string_py>`
* ``AmrexAstroCxxNetwork`` / ``SimpleCxxNetwork`` : stored in
  ``rate_derivs_t`` in the ``ApproximateRate`` evaluation and used
  symbolically using the SymPy methods in
  :py:meth:`SympyRates.jacobian_term_symbol
  <pynucastro.networks.sympy_network_support.SympyRates.jacobian_term_symbol>`



Weak-tabulated rates
--------------------

:py:obj:`TabularWeakRate <pynucastro.rates.tabular_rate.TabularWeakRate>` store the rate
data as functions of $T$ and $\rho Y_e$, where

$$Y_e = \sum_k Z_k Y_k$$

is the electron fraction.  The flux for a tabulated weak-rate decay is:

$$F_{P,\mathrm{weak}} = Y(P) \lambda_{P,\mathrm{weak}}(T, \rho Y_e)$$

this means that in addition to the derivative with respect to the explicit $Y(P)$ composition term,
we also need to account for the $Y_e$ in the tabulation.  For the parent $P$ and child $C$ of the decay,
we need to accumulate the contributions due to $Y_e$.  We do this in an array ``dweak_rates_dYe``
that is part of ``RateEval`` or ``rate_derivs_t``:

.. math::

   \begin{split}
   \texttt{dweak\_rates\_dYe}(P)&\mathrel{-}=
       Y_p\,\rho\,\frac{\partial\lambda}{\partial(\rho Y_e)},\\
   \texttt{dweak\_rates\_dYe}(C)&\mathrel{+}=
   Y_p\,\rho\,\frac{\partial\lambda}{\partial(\rho Y_e)}
   \end{split}

After all of the contributions are accumulated, they are added to every species column:

.. math::

   J_{ij}\mathrel{+}=\texttt{dweak\_rates\_dYe}(i)\,Z_j

.. important::

   This contribution affects all species, not just the parent and child.  As a result,
   the Jacobian with weak rates in it will not be sparse.

Status of this term:

* ``RateCollection`` : not currently included
* ``PythonNetwork`` : not currently include
* ``AmrexAstroCxxNetwork`` / ``SimpleCxxNetwork`` : stored in
  ``rate_derivs_t`` in the ``TabularWeakRate`` evaluation and explicitly
  added to the Jacobian during the final construction of the Jacobian
  in the template C++ code.


Screening
---------

The screening function applied to rates is a function of composition, through the plasma state.
This contribution enters as $Y_e$ and $\overline{Z^2}$:

.. math::

   Y_e &= \sum_k Z_k Y_k \\
   \overline{Z^2} &= \sum_k Z_k^2 Y_k

This means that for each rate flux (where we now explicitly include the screening factor, $f$:

$$F_{AB} = \rho Y(A) Y(B) f_{AB} \lambda_{AB}$$

the contribution to the Jacobian would be:

.. math::

   \frac{\partial F_{AB}}{\partial Y_j} = \rho (Y(B) \delta_{Aj} + Y(A) \delta_{Bj}) f_{AB} \lambda_{AB}
          + \rho Y(A) Y(B) \frac{\partial f_{AB}}{\partial Y_j} \lambda_{AB}

with

.. math::

   \frac{\partial f_{AB}}{\partial Y_j} = \frac{\partial f_{AB}}{\partial Y_e} Z_j
      + \frac{\partial f_{AB}}{\partial \overline{Z^2}} Z_j^2



.. note::

   This contribution affects all composition Jacobian elements, which
   means that the Jacobian is not sparse.


.. important::

   This contribution is not currently implemented.
