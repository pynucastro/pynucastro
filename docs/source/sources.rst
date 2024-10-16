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
were obtained from the `JINA Reaclib database <https://reaclib.jinaweb.org/>`_.

* :cite:t:`ReacLib`

Tabulated weak nuclear reaction rates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* For nuclei with $A = 17$ to $28$ we use the weak rates from
  :cite:t:`suzuki:2016`.

  The data tables were obtained from `<https://www.phys.chs.nihon-u.ac.jp/suzuki/data2/link.html>`_.

* For nuclei with $A = 45$ to $65$ we use the weak rates from
  :cite:t:`langanke:2001`.


Physical constants
------------------

We use the [scipy.constants](https://docs.scipy.org/doc/scipy/reference/constants.html) module
from SciPy to get all the physical constants.  This in turn gets the constants from the CODATA
recommended values (currently CODATA 18)

* :cite:t:`codata:2018`


Nuclei properties
-----------------

We get the basic nuclear properties from the Nubase 2020 evaluation.  This
is available online at `Nuclear Data Services <https://www-nds.iaea.org/amdc/>`_.
We are currently using the file `nubase_4.mas20.txt <https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt>`_.

* :cite:t:`nubase:2020`

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


Partition functions
-------------------

We use the tabulated partition functions from the following sources:

* :cite:t:`rauscher:1997`

* :cite:t:`rauscher:2003`
