Reaction Rate Types
===================

pynucastro can understand reaction rates in a variety of formats.
The basic :func:`Rate <pynucastro.rates.rate.Rate>` class serves
as the base class for all the different rate types and provides
the core functionality.

Here we describe the most commonly used rates.

ReacLib
-------

The `JINA ReacLib database <https://reaclib.jinaweb.org/>`_ provides
the core thermonuclear rates.  These are managed in pynucastro through
the :func:`ReacLibRate <pynucastro.rates.reaclib_rate.ReacLibRate>` class.

A ReacLib rate is composed of several sets (managed as
:func:`SingleSet <pynucastro.rates.reaclib_rate.SingleSet>` objects, each
taking the form:

.. math::

   \lambda = \exp{\left (a_0 + \sum_{i=1}^5  a_i T_9^{(2i-5)/3}  + a_6 \log T_9\right )}

where :math:`\lambda` is of the form:

.. math::

   \lambda = \left \{ \begin{array}{cc} \tau^{-1}  & \mbox{unary reaction} \\
                                     N_A \langle \sigma v\rangle & \mbox{binary reaction} \\
                                     N_A^2 \langle \sigma v\rangle & \mbox{trinary reaction}
                   \end{array} \right .


Chapters
^^^^^^^^

The rates are divided into chapters, based on how many nuclei are on each side of the reaction:

========  ====================================================
chapter    reaction format
--------  ----------------------------------------------------
1         :math:`e_1 \rightarrow e_2`
2         :math:`e_1 \rightarrow e_2 + e_3`
3         :math:`e_1 \rightarrow e_2 + e_3 + e_4`
4         :math:`e_1 + e_2 \rightarrow e_3`
5         :math:`e_1 + e_2 \rightarrow e_3 + e_4`
6         :math:`e_1 + e_2 \rightarrow e_3 + e_4 + e_5`
7         :math:`e_1 + e_2 \rightarrow e_3 + e_4 + e_5 + e_6`
8         :math:`e_1 + e_2 + e_3 \rightarrow e_4`
9         :math:`e_1 + e_2 + e_3 \rightarrow e_4 + e_5`
10        :math:`e_1 + e_2 + e_3 + e_4 \rightarrow e_5 + e_6`
11        :math:`e_1 \rightarrow e_2 + e_3 + e_4 + e_5`
========  ====================================================

Labels
^^^^^^

The ReacLib database lists the source / reference of each rate with a 6 character string.

* The first 4 characters are the label that gives the source of the rate, according to:
  https://reaclib.jinaweb.org/labels.php

  This is stored in both ``SingleSet`` and ``ReacLibRate`` as the
  ``.label`` attribute.

* The next character is ``n`` for a non-resonant set, ``r`` for a
  resonance, or ``w`` to indicate that it is a weak rate.

  Different sets can have either `n` or `r` (and there can be multiple
  resonances).

  * For a ``SingleSet``, ``SingleSet.resonant`` is ``True`` if it is a
    resonance

  * For the combined ``ReacLibRate`` object, we combine all of the
    resonances and the non-resonance set, and change the flag to ``c``
    for _combined_.

  If the weak flag, ``w`` is set, then ``ReacLibRate.weak`` will be ``True``.

* The 6th character indicates it the rate is a derived reverse rate,
  by the presence of a `v`.  This is stored in both ``SingleSet`` and
  ``ReacLibRate`` as the ``.reverse`` attribute.


Weak rates
^^^^^^^^^^

Electron capture capture rates are those identified with a label of either ``ec`` or ``bec``.  These
have ``ReacLibRate.weak_type`` set to ``electron_capture``.  For these rates, and only these rates,
we include :math:`Y_e` in the overall rate (multiplying density).

.. note::

   There is some ambiguity as to whether :math:`Y_e` should be included for ``bec`` rates.


The only other place ``weak_type`` is used for ``ReacLibRate`` is to set the print representation.

.. note::

   ReacLib only provides a fit to the temperature dependence of a
   rate, so for electron-captures, it may not be a very good
   approximation, and tabulated electron capture rates should be used
   instead.


Reverse rates
^^^^^^^^^^^^^

As noted above, rates with the ``v`` flag in the label are reverse
rates that were computed from the forward rate via detailed balance.
These rates do not include the corrections from partition functions,
and therefore, should not be used directly.  Instead, the
:func:`DerivedRate <pynucastro.rates.derived_rate.DerivedRate>` functionality in pynucastro can be used to redo the
detailed balance including the effects of partition functions.

.. note::

   In ReacLib, *reverse* does not always mean :math:`Q < 0`.  Sometimes the
   rate with :math:`Q < 0` is easier to measure experimentally, and so
   that is measured and then the :math:`Q > 0` rate is computed via
   detailed balance.

Rate evaluation functions
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``ReacLibRate`` class has functions :func:`function_string_py
<pynucastro.rates.reaclib_rate.ReacLibRate.function_string_py>` and
:func:`function_string_cxx
<pynucastro.rates.reaclib_rate.ReacLibRate.function_string_cxx>` to write out
the python and C++ code needed to evaluate the T-dependent portion of
the reaction rate (basically what is encoded in the ReacLib database).

For the C++ version, templating is used to allow for the derivative
with respect to temperature to also be computed.

ydot term
^^^^^^^^^

The ``ReacLibRate`` class knows how to output the contribution to the
molar fraction evolution (:math:`\dot{Y}`) as a python expression (for
C++, this is handled separately via SymPy in the ``network`` module).
This is handled by ``ReacLibRate.ydot_string_py()``.

For a unary reaction involving nucleus :math:`A`, it takes the form:

.. math::

   \dot{Y} = \frac{Y(A)}{\tau}

for a binary reaction, :math:`A + B`, it takes the form:

.. math::

   \dot{Y} = \rho Y(A) Y(B) \frac{N_A \langle \sigma v \rangle}{1 + \delta_{AB}}

where the :math:`1 + \delta_{AB}` factor is stored in the rate as ``Rate.prefactor``.

and for a trinary reaction, it is:

.. math::

   \dot{Y} = \rho^2 Y(A)^{n_A} Y(B)^{n_B} Y(C)^{n_C} \frac{N_A^2 \langle \sigma v \rangle}{n_A! n_B! n_C!}

where :math:`n_A` is the number of nucleus :math:`A` in the reaction.

.. note::

   The rate class does not include the stoichiometric factors -- that
   is the responsibility of the network.


Similarly,  :func:`jacobian_string_py <pynucastro.rates.rate.Rate.jacobian_string_py>`
outputs the contribution to the Jacobian for this rate.


Tabulated Rates
---------------

For electron captures and beta-decays (which are of the form
:math:`\rm{A \rightarrow B}`), we use tabulated rates.  These are
two-dimensional tables, in terms of :math:`T` and :math:`\rho Y_e`,
and managed by the :func:`TabularRate <pynucastro.rates.tabular_rate.TabularRate>` class.

.. note::

   If positron captures and decays are available in the data source,
   then these are combined with the appropriate electron counterpart
   into a single rate.

A tabular rate is described by a single file for each reaction
(i.e., beta-decays and electron-captures are in separate files).

The data reading and interpolation are managed by the
:py:obj:`TabularRate <pynucastro.rates.tabular_rate.TabularRate>`
class.

ydot term
^^^^^^^^^

The form of the reaction :math:`A \rightarrow B` is:

.. math::

   \dot{Y}_A = -Y(A) \lambda

where :math:`\lambda` is the rate returned from the table.


Table format
^^^^^^^^^^^^

Each rate table has a header (with lines starting with ``!``), followed
by the data.  An example can be seen as:
`suzuki-23na-23ne_electroncapture.dat <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/library/tabular/suzuki/suzuki-23na-23ne_electroncapture.dat>`_ in
``pynucastro/library/tabular/suzuki``

.. important::

   The first line header needs to include a line giving the reaction.
   It should take the form, e.g.,

   ::

      !65zn -> 65cu, e- capture

   where ``!`` is the comment character, ``65zn`` is the reactant, and ``65cu`` is
   the product.  Any additional information is ignored.

   An alternate form including the spins works as well, but the above
   it preferred.

The columns of the tables (and
units) are:

* $\log_{10} (\rho Y_e)$: electron density in $\mathrm{g~cm^{-3}}$

* $\log_{10} T$: temperature in $\mathrm{K}$

* $\mu$: chemical potential in erg

* $\Delta Q$: threshold energy in erg

* $V_s$: Coulomb potential at the origin in erg

* $\log_{10} (\lambda)$: electron capture or $\beta$-decay rate in $\mathrm{s^{-1}}$.
  For some rates, this is ($e^-$-capture and $e^+$-decay) or ($\beta$-decay + $e^+$-capture)

* $\log_{10} (\epsilon_\nu)$: neutrino energy loss in $\mathrm{erg~s^{-1}}$

* $\log_{10} (\epsilon_\gamma)$: gamma energy loss in $\mathrm{erg~s^{-1}}$

and the data is ordered with ``rhoY`` varying the slowest (i.e., for a
given ``rhoY`` we loop over all of the temperatures).

.. note::

   The number of temperature and $\rho Y_e$ points used in the tabulation
   is inferred from the data file itself.


Rate evaluation functions
^^^^^^^^^^^^^^^^^^^^^^^^^

Analogous to ``ReacLibRate``, ``TabularRate`` provides functions to
evaluate the rate and output the python code.  The function
:func:`function_string_py
<pynucastro.rates.tabular_rate.TabularRate.function_string_py>`
outputs the python code for managing the interpolation of the data.
For C++ networks, this interpolation is handled directly by the network class.
