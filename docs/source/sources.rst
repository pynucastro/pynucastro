Third Party Data
================

pynucastro incorporates the following publicly-available
third-party data. Links to this data as well as citations to the
relevant publications are as follows.

Reaction rates
--------------

Nuclear reaction rates from JINA Reaclib
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reaction rate parameterizations in `pynucastro/library <https://github.com/pynucastro/pynucastro/tree/main/pynucastro/library>`_
were obtained from the `JINA Reaclib database <https://reaclib.jinaweb.org/>`_, :cite:t:`ReacLib`.

.. _tabulated_rate_sources:

Tabulated weak nuclear reaction rates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The weak rates come from several different sources, each of which focuses on
a particular range of masses.  There is overlap in the mass range between
these sources, so pynucastro provides a way to give priority to particular
sources.

.. warning::

   In the C++ networks, no error is produced if you try to evaluate a rate outside of the table
   limits.  Instead, an extrapolation is done using the data at the edge of the table.

The available sources are:

* :cite:t:`ffn` (covers $A = 21$ to $60$):

  The density (g/cm$^3$) and temperature (K) ranges of the rates are:

  * $1 \le \log_{10}(\rho Y_e) \le 11$
  * $7 \le \log_{10}(T) \le 11$

* :cite:t:`oda:1994` (covers $A = 17$ to $39$):

  The density (g/cm$^3$) and temperature (K) ranges of the rates are:

  * $1 \le \log_{10}(\rho Y_e) \le 11$
  * $7 \le \log_{10}(T) \le 10.48$

* :cite:t:`langanke:2001` (covers $A = 45$ to $65$):

  The data tables are provided as part of that journal link.

  The density (g/cm$^3$) and temperature (K) ranges of the rates are:

  * $1 \le \log_{10}(\rho Y_e) \le 11$
  * $7 \le \log_{10}(T) \le 11$

  .. note::

     These rates should be preferred to the FFN rates where there is
     overlap.

* :cite:t:`suzuki:2016` (covers $A = 17$ to $28$):

  These data tables were obtained from
  `<https://www.phys.chs.nihon-u.ac.jp/suzuki/data2/link.html>`_.

  The density (g/cm$^3$) and temperature (K) ranges of the rates are:

  * $7 \le \log_{10}(\rho Y_e) \le 11$
  * $7 \le \log_{10}(T) \le 9.65$

  .. note::

     The paper :cite:`suzuki:2016` says that the rates are evaluated
     are in the range $8 \le \log_{10}(\rho Y_e) \le 11$, but the
     tables provided have the lower limit as $\log_{10}(\rho Y_e) =
     7$.



Physical constants
------------------

We use the `scipy.constants <https://docs.scipy.org/doc/scipy/reference/constants.html>`_ module
from SciPy to get all the physical constants.  This in turn gets the constants from the CODATA
recommended values (for SciPy versions 1.14 and earlier, this is CODATA 2018, :cite:t:`codata:2018`;
for SciPy 1.15 they are CODATA 2022).


Nuclei properties
-----------------

We get the basic nuclear properties from the Nubase 2020 evaluation :cite:t:`nubase:2020`.  This
is available online at `Nuclear Data Services <https://www-nds.iaea.org/amdc/>`_.
We are currently using the file `nubase_4.mas20.txt <https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt>`_.

In particular, we get the mass excesses, $\Delta m$, and spins from there.  We then compute
mass of the nucleus as:

.. math::

   m = \Delta M + A m_u

and the binding energies from the mass excesses as:

.. math::

   B = Z m_H + N m_n - (A m_u + \Delta m)

where $m_H$ is the mass of the hydrogen atom, computed from the mass
excess of ``1H`` listed in the table.  This is consistent with the
discussion in section 2 of the AME 2020 paper :cite:`ame2020_1`, and
these numbers match the binding energies computed in the AME tables to
the uncertainty in the nuclear masses.

Binding energies are also computed and tablulated in the AME mass
evaluation (see `AME2020 mass table
<https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt>`_).  But note
that the Nubase evaluation seems to more closely follow the "rounded"
version of the table `AME2020 rounded mass
table <https://www-nds.iaea.org/amdc/ame2020/massround.mas20.txt>`_.
The rounding procedure is discussed in Table I on the `AME 2020 paper
II <https://iopscience.iop.org/article/10.1088/1674-1137/abddaf>`_ (also
see the `Nubase2020
paper <https://iopscience.iop.org/article/10.1088/1674-1137/abddae>`_,
Table I).

Partition functions
-------------------

We use the tabulated partition functions from the following sources:

* :cite:t:`rauscher:1997`

* :cite:t:`rauscher:2003`
