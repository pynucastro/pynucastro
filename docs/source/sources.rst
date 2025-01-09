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

Tabulated weak nuclear reaction rates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* For nuclei with $A = 17$ to $28$ we use the weak rates from
  :cite:t:`suzuki:2016`.

  The data tables were obtained from `<https://www.phys.chs.nihon-u.ac.jp/suzuki/data2/link.html>`_.

* For nuclei with $A = 45$ to $65$ we use the weak rates from
  :cite:t:`langanke:2001`.


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
