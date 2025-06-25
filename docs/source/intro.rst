Overview of pynucastro
======================

pynucastro is a set of python interfaces to nuclear reaction rates. It
is meant for both interactive exploration of rates (through Jupyter
notebooks) and to create reaction networks for use in simulation
codes.

The preferred way of importing pynucastro is as follows:

.. code-block:: python

   import pynucastro as pyna


The main classes are:

* :func:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>`: This is a single
  nucleus.  It knows its proton number, ``Z``, neutron number, ``N``,
  weight, ``A``, and binding energy, ``nucbind``, as well as
  T-dependent partition function and ground state spin.

* :func:`Rate <pynucastro.rates.rate.Rate>`: This is a single rate.
  It knows the reactants and products and has methods that allow you
  to evaluate it at a specified temperature and plot its temperature
  dependence.  A `Rate` also knows how to the generate code (C++ and
  python) needed to evaluate it.

  There are a few special rates derived from `Rate`:

  * :func:`ReacLibRate <pynucastro.rates.reaclib_rate.ReacLibRate>`: This is a rate in the
    JINA ReacLib format, with the temperature dependence specified by an interpolant
    with 7 different coefficients.

  * :func:`TabularRate <pynucastro.rates.tabular_rate.TabularRate>`: This is a
    rate that is tabulated in terms of :math:`(T, \rho Y_e)`.  This is
    how the weak rate (electron captures and beta-decays) are stored.
    Interpolation is used to find the rate at any thermodynamic state.

  * :func:`ApproximateRate <pynucastro.rates.approximate_rates.ApproximateRate>`:
    An approximate rate assumes equilibration of intermediate nuclei to create
    an approximation for a rate sequence.  Currently, there are two
    approximations that can be made:

    * grouping $A(\alpha, \gamma)B$ and $A(\alpha,p)X(p,\gamma)B$ into
      a single effective rate, assuming equilibrium of $p$ and $X$.

    * converting $A(n,\gamma)X(n,\gamma)B$ into $A(nn,\gamma)B$
      by assuming equilibrium of $X$.

  * :func:`ModifiedRate <pynucastro.rates.modified_rate.ModifiedRate>`:
    A container for a single rate that allows for different stoichiometry
    or products.

  * :func:`DerivedRate <pynucastro.rates.derived_rate.DerivedRate>`: A
    derived rate uses detailed balance to recompute a reverse rate from the forward rate.

* :func:`RatePair <pynucastro.rates.rate.RatePair>`: For a single nuclear process,
  this holds the corresponding forward and reverse rates.

* :func:`Library <pynucastro.rates.library.Library>`: This is a collection of
  rates (for example, the entire ReacLib library).  It provides methods
  for filtering out rates based on different sets of rules.

  There are two important subclasses:

  * :func:`ReacLibLibrary <pynucastro.rates.library.ReacLibLibrary>`: The
    entire ReacLib rate library (> 80,000 rates)

  * :func:`TabularLibrary <pynucastro.rates.library.TabularLibrary>`: A
    `Library` containing all known tabular weak rates.

* :func:`Composition
  <pynucastro.networks.rate_collection.Composition>`: This is a
  collection of nuclei and their mass fractions.  A ``Composition`` is
  used when evaluating the full rates in a network.

* :func:`RateCollection
  <pynucastro.networks.rate_collection.RateCollection>`: This is the
  most basic form of a network.  It is a collection of rates and
  nuclei, that knows about the connectivity of the nuclei through
  different reactions.  This acts as the base class for different
  reaction networks.  A ``RateCollection`` has methods to evaluate the
  rates and make a plot of the links between rates.

  There are a few important subclasses:

  * :func:`NSENetwork
    <pynucastro.networks.nse_network.NSENetwork>`: This allows
    a user to find the nuclear statistical equilibrium state
    of a collection of nuclei.

  * :func:`PythonNetwork
    <pynucastro.networks.python_network.PythonNetwork>`: This is a
    collection of rates with functions that know how to write python
    code to express the righthand side of the system of ODEs.

  * :func:`SimpleCxxNetwork
    <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>`:
    This is a simple C++ network that provides functions for
    computing the righthand side and Jacobian of a network.
    Not all pynucastro features are supported in this network.

  * :func:`FortranNetwork
    <pynucastro.networks.fortran_network.FortranNetwork>`:
    A network that provides Fortran wrappers to ``SimpleCxxNetwork``.

  * :func:`AmrexAstroCxxNetwork
    <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>`:
    This is a C++ network of the form needed by the `AMReX
    Astrophysics Microphysics
    <https://github.com/AMReX-Astro/Microphysics>`_ library used by
    the Castro and MAESTROeX simulation codes.


Usage
-----

There are two modes of usage for pynucastro.  

* Within a Jupyter notebook, one can evaluate the rates and
  interactively visualize a network and see the flow between nuclei as
  connections colored by the rate linking them.

* You can use pynucastro to write the righthand side routine for the
  system of ODEs that must be integrated to evolve a reaction network.
  A reaction network takes the form:

  .. math::

     \frac{dY_i}{dt} = - \sum_{j,k} Y_i Y_j \lambda_{i(j,k)l} + \sum_{j,k} Y_l Y_k \lambda_{l(j,k)i}

  where the :math:`\lambda`'s are the rates of destruction and creation
  of species i, represented by the molar fraction :math:`Y_i` (see,
  e.g., :cite:t:`timmes:1999`).  pynucastro
  will create the righthand sides of this system of ODEs (as python or
  C++ code) from the list of rates you provide. One can use this to
  add reaction networks to existing simulation codes, for example, the
  `MAESTROeX <https://amrex-astro.github.io/MAESTROeX/>`_ and `Castro
  <https://amrex-astro.github.io/Castro/>`_ codes.


