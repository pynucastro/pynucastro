Data Sources and Rate Types
===========================

pynucastro makes use of data sources from several different nuclear
physics compilations.  Here we explain how these are used in
pynucastro.

ReacLib
-------

The `JINA ReacLib database <https://reaclib.jinaweb.org/>`_ provides
the core thermonuclear rates.  These are managed in pynucastro through
the :func:`ReacLibRate <pynucastro.rates.rate.ReacLibRate>` class.

A ReacLib rate is composed of several sets (managed as
:func:`SingleSet <pynucastro.rates.rate.SingleSet>` objects, each
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
:func:`DerivedRate <pynucastro.rates.rate.DerivedRate>` functionality in pynucastro can be used to redo the
detailed balance including the effects of partition functions.

.. note::

   In ReacLib, _reverse_ does not always mean :math:`Q < 0`.  Sometimes the
   rate with :math:`Q < 0` is easier to measure experimentally, and so
   that is measured and then the :math:`Q > 0` rate is computed via
   detailed balance.

Rate evaluation functions
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``ReacLibRate`` class has functions :func:`function_string_py
<pynucastro.rates.rate.ReacLibRate.function_string_py>` and
:func:`function_string_cxx
<pynucastro.rates.rate.ReacLibRate.function_string_cxx>` to write out
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


Similarly, ``ReacLib.jacobian_string_py()`` outputs the contribution to the Jacobian for this rate.


Tabulated Rates
---------------

For electron captures and beta-decays, we used tabulated rates.  These are two-dimensional tables,
in terms of :math:`T` and :math:`\rho Y_e`.

A tabular rate is described by 2 files.  The first file mimics the
ReacLib header, with a chapter indicated as ``t`` and gives the name
of the table and the number of columns, density, and temperature
points.

The second file is the table itself.  The columns of the tables (and
units) are:

::

   !rhoY        T            mu           dQ           Vs           e-cap-rate   nu-energy-loss gamma-energy
   !g/cm^3      K            erg          erg          erg          1/s          erg/s        erg/s

and the data is ordered with ``rhoY`` varying the slowest (i.e., for a
given ``rhoY`` we loop over all of the temperatures).

pynucastro uses linear interpolation to return the rate given the
temperature and electron density.

The form of the reaction :math:`A \rightarrow B`

.. math::

   \dot{Y} = Y(A) \lambda

where :math:`\lambda` is the rate returned from the table.

