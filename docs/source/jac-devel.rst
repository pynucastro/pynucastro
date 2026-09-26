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

which we note as Regions I through IV:

.. math::

   {\bf J} = \left ( \begin{array}{c|c}
                   I & III \\
                  \hline
                  II & IV \end{array} \right )

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

   Furthermore, we presently require a numerical (finite-difference) Jacobian when
   evolving energy in a ``PythonNetwork``.

.. note::

   A :py:obj:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>`
   does not include the energy evolution terms, so its Jacobian is species-only.


Region I: :math:`\partial \dot{Y}_i/\partial Y_j`
=================================================

For most rates, this is simply the derivative of the composition factors that appear explicitly
in the reactive flux.  There are a few instances where the composition is internal to the rate
itself that we also need to capture.

Basic rates
-----------

For a 2-body strong-mediated reaction rate (with distinct reactants $A$ and $B$), the reactive flux has the form:

$$F_{AB} = \rho Y(A) Y(B) \lambda_{AB}$$

and the contributions to :math:`\partial \dot{Y}_i/\partial Y_j` would be $\partial F_{AB}/\partial Y_j$,

.. math::

   \frac{\partial F_{AB}}{\partial Y_j} = \rho (Y(B) \delta_{Aj} + Y(A) \delta_{Bj}) \lambda_{AB}

For weak rates, the explicit composition dependence of the parent nucleus is also computed this way.  E.g.,
for a decay with a parent nucleus $P$, the flux is:

$$F_{P,\mathrm{weak}} = Y(P) \lambda_{P,\mathrm{weak}}$$

The composition derivative is:

.. math::

   \frac{\partial F_{P,\mathrm{weak}}}{\partial Y_j} = \delta_{Pj} \lambda_{P,\mathrm{weak}}

and the contribution to the Jacobian is:

.. math::

   J_{ij}\mathrel{+}=(n_i^{\mathrm{products}}-n_i^{\mathrm{reactants}})
     \frac{\partial F}{\partial Y_j}.

where $n_i$ are the stoichiometric coefficients for species $i$.


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


In a rate class, we set ``Rate.rate_comp_dependence`` to a list of the
nuclei that the rate ($\lambda$) depends on internally.  This
indicates that we need to compute this derivative.

The rate class itself will then compute the explicit $\partial\lambda/\partial Y_j$ term
and store it in the python ``RateEval`` class or the C++ ``rate_derivs_t`` struct.

This contribution was added in `pynucastro PR #1537 <https://github.com/pynucastro/pynucastro/pull/1537>`_.

Status of this term:

* ``RateCollection`` : not currently included
* ``PythonNetwork`` : (to be done) included in the string returned via
  :py:meth:`Rate.jacobian_string_py <pynucastro.rates.rate.Rate.jacobian_string_py>`
* ``AmrexAstroCxxNetwork`` / ``SimpleCxxNetwork`` : stored in
  ``rate_derivs_t`` in the ``ApproximateRate`` evaluation and used
  symbolically using the SymPy methods in
  :py:meth:`BaseCxxNetwork.compose_jacobian
  <pynucastro.networks.base_cxx_network.BaseCxxNetwork.compose_jacobian>`



Weak-tabulated rates
--------------------

:py:obj:`TabularWeakRate <pynucastro.rates.tabular_rate.TabularWeakRate>` store the rate
data as functions of $T$ and $\rho Y_e$, where

$$Y_e = \sum_k Z_k Y_k$$

is the electron fraction.  The flux for a tabulated weak-rate decay is:

$$F_{P,\mathrm{weak}} = Y(P) \lambda_{P,\mathrm{weak}}(T, \rho Y_e)$$

this means that in addition to the derivative with respect to the explicit $Y(P)$ composition term,
we also need to account for the $Y_e$ in the tabulation.  For the parent $P$ and child $C$ of the decay,
we need to accumulate the contributions due to $Y_e$.  We do this in an array ``dweak_ydot_dYe``
that is part of ``RateEval`` or ``rate_derivs_t``:

.. math::

   \begin{split}
   \texttt{dweak\_ydot\_dYe}(P)&\mathrel{-}=
       Y_p\,\rho\,\frac{\partial\lambda}{\partial(\rho Y_e)},\\
   \texttt{dweak\_ydot\_dYe}(C)&\mathrel{+}=
   Y_p\,\rho\,\frac{\partial\lambda}{\partial(\rho Y_e)}
   \end{split}

After all of the contributions are accumulated, they are added to every species column:

.. math::

   J_{ij}\mathrel{+}=\texttt{dweak\_ydot\_dYe}(i)\,Z_j

.. important::

   This contribution affects all charged species, not just the parent and child.  As a result,
   the Jacobian with weak rates in it will not be nearly as sparse.

This contribution was added in `pynucastro PR #1539 <https://github.com/pynucastro/pynucastro/pull/1539>`_.

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
This contribution enters as $Y_e$ and $\zeta$:

.. math::

   Y_e &= \sum_k Z_k Y_k \\
   \zeta &= \sum_k Z_k^2 Y_k

.. note::

   $\zeta$ is related to the quantity $\overline{Z^2}$ that the ``plasma_state_t`` calls ``z2bar`` as:

   $$\overline{Z^2} = \bar{A} \zeta$$

This means that for each rate flux (where we now explicitly include the screening factor, $f$):

$$F_{AB} = \rho Y(A) Y(B) f_{AB} \lambda_{AB}$$

the contribution to the Jacobian would be:

.. math::

   \frac{\partial F_{AB}}{\partial Y_j} = \rho (Y(B) \delta_{Aj} + Y(A) \delta_{Bj}) f_{AB} \lambda_{AB}
          + \rho Y(A) Y(B) \frac{\partial f_{AB}}{\partial Y_j} \lambda_{AB}

with

.. math::

   \frac{\partial f_{AB}}{\partial Y_j} = \frac{\partial f_{AB}}{\partial Y_e} Z_j
      + \frac{\partial f_{AB}}{\partial \zeta} Z_j^2



.. note::

   This contribution affects all composition Jacobian elements for charged species, which
   means that the Jacobian is no longer really sparse.


.. important::

   This contribution is not currently implemented.


Region II: :math:`\partial \dot{e}/\partial Y_j`
================================================

The energy evolution equation is:

$$\frac{de}{dt} = \epsilon_\mathrm{nuc} + \epsilon_{\nu,\mathrm{weak}} - \epsilon_{\nu,\mathrm{thermal}}$$

.. important::

   The neutrino contribution in $\epsilon_{\nu,\mathrm{weak}}$ is constructed to be negative in
   the table interpolation routines.  In general, $\epsilon_{\nu,\mathrm{weak}}$ can also include
   a positive $\gamma$ contribution, hence the $+$ in this evolution equation.
   This term will represents an energy loss when neutrinos dominate.

Each of these energy terms depends on temperature and composition.


Binding energy contribution
---------------------------

The dominant contribution to energy comes from the mass difference between reactants and products:

$$\epsilon_\mathrm{nuc} = -N_A \sum_i \frac{\partial Y_i}{\partial t} m_i c^2$$

Differentiating with respect to $Y_j$, we get the contribution to each
column of the $\partial \dot{e}/\partial Y_j$ row:

$$\frac{\partial \epsilon_\mathrm{nuc}}{\partial Y_j} = -N_A \sum_i \frac{\partial \dot{Y}_i}{\partial Y_j} m_i c^2$$

We already computed the elements $\partial \dot{Y}_i/\partial Y_j$ in
Region I, so we can interpret this sum as simply summing down a
column, weighting by $m_i c^2$.

Status of this term:

* ``RateCollection`` : N/A (energy not considered)
* ``PythonNetwork`` : N/A (numerical Jacobian is used with self-heating networks)
* ``SimpleCxxNetwork`` : N/A (energy not considered)
* ``AmrexAstroCxxNetwork`` : computed directly in
  the template C++ code.


Weak rate energy contribution
-----------------------------

For a decay $P \rightarrow C$, there is an energy release:

$$\epsilon_{\nu.\mathrm{weak}} = N_A \, Y(P)\, (\dot{e}_\nu + \dot{e}_\gamma)$$

where $\dot{e}_\nu$ is the neutrino energy release rate (erg/s) from the rate table,
and $\dot{e}_\gamma$ is the gamma energy release rate (erg/s) (note: most rate tabulations
do not provide this).

There are 2 contributions we need to account for here.

Explicit Y dependence
^^^^^^^^^^^^^^^^^^^^^

First we account for the explicit dependence on the
parent nucleus.  We accumulate this in ``rate_derivs_t.denuc_weak_dY``:

.. math::

   \texttt{denuc\_weak\_dY}(P) \mathrel{+}=
        N_A\, (\dot{e}_\nu + \dot{e}_\gamma)

.. important::

   $\dot{e}_\nu$ is negative, so this quantity represents an energy loss.

This term is then added as $\partial
\epsilon_{\nu,\mathrm{weak}}/\partial Y_j$ to the respective species
column in the energy row.

Electron fraction dependence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The second term reflects the fact that $\dot{e}_\nu = \dot{e}_\nu(T, \rho Y_e)$,
so we have an additional species contribution through $Y_e$.  We accumulate this
in ``rate_derivs_t.denuc_weak_dYe``:

.. math::

   \texttt{denuc\_weak\_dYe} \mathrel{+}=
        N_A\, \rho Y(P) \, \frac{\partial \dot{e}_\nu}{\partial \rho Y_e}

This term is then added as $\texttt{denuc\_weak\_dYe} \, Z_j$
in the same place as the term above.

This contribution was added in `pynucastro PR #1539 <https://github.com/pynucastro/pynucastro/pull/1539>`_.

Status of these terms
^^^^^^^^^^^^^^^^^^^^^

* ``RateCollection`` : N/A (energy not considered)
* ``PythonNetwork`` : N/A (numerical Jacobian is used with self-heating networks)
* ``SimpleCxxNetwork`` : N/A (energy not considered)
* ``AmrexAstroCxxNetwork`` : stored in
  ``rate_derivs_t`` in the ``TabularWeakRate`` evaluation and explicitly
  added to the Jacobian during the final construction of the Jacobian
  in the template C++ code.

Thermal neutrinos
-----------------

The thermal neutrino loss rate, $\epsilon_{\nu,\mathrm{thermal}}$, is a positive
energy loss per unit mass (erg/g/s).  In ``AmrexAstroCxxNetwork``, its composition
dependence enters through the mean mass number and mean charge:

.. math::

   \bar{A} &= \left (\sum_k Y_k \right )^{-1} \\
   \bar{Z} &= \bar{A} \sum_k Z_k Y_k = \bar{A} Y_e .

At fixed $T$ and $\rho$, their composition derivatives are:

.. math::

   \frac{\partial \bar{A}}{\partial Y_j} &= -\bar{A}^2 \\
   \frac{\partial \bar{Z}}{\partial Y_j} &= \bar{A}(Z_j - \bar{Z}).

Since the thermal neutrino loss is subtracted in the energy equation, its
contribution to the energy row is:

.. math::

   J_{e j} \mathrel{+}= -\frac{\partial \epsilon_{\nu,\mathrm{thermal}}}{\partial Y_j}
     = -\left[-\bar{A}^2
         \frac{\partial \epsilon_{\nu,\mathrm{thermal}}}{\partial \bar{A}}
       + \bar{A}(Z_j - \bar{Z})
         \frac{\partial \epsilon_{\nu,\mathrm{thermal}}}{\partial \bar{Z}}\right].

In ``actual_jac`` in the ``actual_rhs.H`` template, the call to
``neutrino_cooling<1>`` returns these derivatives as ``dsnuda`` and ``dsnudz``.
This term is then computed and added to the Jacobian.


Region III: :math:`\partial \dot{Y}_i/\partial e`
=================================================

All rate objects can compute their temperature derivative, so we can compute
the temperature derivative of the flux of each rate, e.g., for our 
2-body strong-mediated reaction rate:

$$F_{AB} = \rho Y(A) Y(B) \lambda_{AB}$$

we have:

$$\frac{\partial F_{AB}}{\partial T} = \rho Y(A) Y(B) \frac{\partial \lambda_{AB}}{\partial T}$$

likewise, for tabulate weak rates, we can compute the derivative with
respect to temperature by differentiating the interpolant.

We take advantage of the fact that each rate's flux contributing to $\partial Y_i/\partial t$ is linear in
$\lambda$, and simply construct the algebraic form of $\partial Y_i/\partial t$ using
$\partial \lambda/\partial T$ instead of $\lambda$ by calling the ``rhs_nuc`` function
in ``actual_rhs.H``.  This gives us $\partial \dot{Y}_i/\partial T$.


Screening
---------

There is an additional contribution from the temperature derivative of screening.  Taking into
account screening $f$, the flux is:

$$F_{AB} = \rho Y(A) Y(B) f_{AB} \lambda_{AB}$$

and the total temperature derivative

.. math::

   \frac{\partial F_{AB}}{\partial T} = \rho Y(A) Y(B) \left [ \frac{\partial f_{AB}}{\partial T} \lambda_{AB} + f_{AB} \frac{\partial \lambda_{AB}}{\partial T} \right ]

We store the quantity in $[ \ldots ]$ in ``rate_derivs_t.dscreened_rates_dT`` when we evaluate the rates.


Status of these terms
---------------------

* ``RateCollection`` : N/A (energy not considered)
* ``PythonNetwork`` : N/A (numerical Jacobian is used with self-heating networks)
* ``SimpleCxxNetwork`` : N/A (energy not considered)
* ``AmrexAstroCxxNetwork`` : computed directly
  in the template C++ code using the rate derivatives with respect to $T$.



As energy derivative
--------------------

We convert this to an energy derivative as:

.. math::

   \frac{\partial F_{AB}}{\partial e} = \frac{1}{c_v} \frac{\partial F_{AB}}{\partial T}




Region IV: :math:`\partial \dot{e}/\partial e`
==============================================

Binding energy contribution
---------------------------

Starting with

$$\epsilon_\mathrm{nuc} = -N_A \sum_i \frac{\partial Y_i}{\partial t} m_i c^2$$

and differentiating with respect to temperature, we have:

$$\frac{\partial \epsilon_\mathrm{nuc}}{\partial T} = -N_A \sum_i \frac{\partial \dot{Y}_i}{\partial T} m_i c^2$$

we have $\frac{\partial \dot{Y}_i}{\partial T}$ from Region III.  So
we can just compute $\frac{\partial \epsilon_\mathrm{nuc}}{\partial T}$
from these.

This is computed as:

* ``RateCollection`` : N/A (energy not considered)
* ``PythonNetwork`` : N/A (numerical Jacobian is used with self-heating networks)
* ``SimpleCxxNetwork`` : N/A (energy not considered)
* ``AmrexAstroCxxNetwork`` : directly in the C++ template.


Weak-rate neutrino contribution
-------------------------------

The weak rate energy is:

$$\epsilon_{\nu.\mathrm{weak}} = N_A \, Y(P)\, (\dot{e}_\nu + \dot{e}_\gamma)$$

We can differentiate this with respect to temperature (ignoring $\dot{e}_\gamma$:

$$\frac{\partial\epsilon_{\nu.\mathrm{weak}}}{\partial T} = N_A \, Y(P)\, \frac{\partial \dot{e}_\nu}{\partial T}$$

we accumulate this contribution in ``rate_derivs_t.denuc_weak_dT`` when we evaluate the
tabular rates as:

.. math::

   \texttt{denuc\_weak\_dT} \mathrel{+}=
       N_A\, Y(P)\, \frac{\partial\dot{e}_\nu}{\partial T}

and then add it to the $\partial \dot{e} /\partial T$ term in the Jacobian function.

This is computed as:

* ``RateCollection`` : N/A (energy not considered)
* ``PythonNetwork`` : N/A (numerical Jacobian is used with self-heating networks)
* ``SimpleCxxNetwork`` : N/A (energy not considered)
* ``AmrexAstroCxxNetwork`` : directly in the C++ template.



Thermal neutrinos
-----------------

At fixed density and composition, the thermal neutrino contribution to the
temperature derivative of the energy RHS is:

.. math::

   \frac{\partial \dot{e}}{\partial T} \mathrel{+}=
       -\frac{\partial \epsilon_{\nu,\mathrm{thermal}}}{\partial T}.

The same ``neutrino_cooling<1>`` call in the ``AmrexAstroCxxNetwork`` template
returns this loss-rate derivative as ``dsneutdt``.  In ``actual_jac``, we subtract
it from ``jac_e_T``.


As energy derivative
--------------------

We convert this to an energy derivative as:

.. math::

   \frac{\partial \dot{e}}{\partial e} = \frac{1}{c_v} \frac{\partial  \dot{e}}{\partial T}


Final conversion to energy
==========================

There is one last part of the conversion from $T$ to $e$.  If we take the derivative with respect
to species, with $e$ and $\rho$ held constant, then $T$ will change.  We need to take into account
how this affects the reactions.

This is handled by the integration wrappers in `AMReX Microphysics <https://github.com/amrex-astro/Microphysics>`_.
