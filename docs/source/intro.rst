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

* :func:`Nucleus <pynucastro.rates.rate.Nucleus>`: This is a single
  nucleus.  It knows its proton number, ``Z``, neutron number, ``N``,
  weight, ``A``, and binding energy, ``nucbind``, as well as
  T-dependent partition function and ground state spin.

* :func:`Rate <pynucastro.rates.rate.Rate>`: This is a single rate.
  It knows the reactants and products and has methods that allow you
  to evaluate it at a specified temperature and plot its temperature
  dependence.  A `Rate` also knows how to the generate code (C++ and
  python) needed to evaluate it.

  There are a few special rates derived from `Rate`:

  * :func:`ReacLibRate <pynucastro.rates.rate.ReacLibRate>`: This is a rate in the
    JINA ReacLib format, with the temperature dependence specified by an interpolant
    with 7 different coefficients.

  * :func:`TabularRate <pynucastro.rates.rate.TabularRate>`: This is a
    rate that is tabulated in terms of :math:`(T, \rho Y_e)`.  This is
    how the weak rate (electron captures and beta-decays) are stored.
    Interpolation is used to find the rate at any thermodynamic state.

  * :func:`ApproximateRate <pynucastro.rates.rate.ApproximateRate>`:
    An approximate rate groups together :math:`A(\alpha, \gamma)B` and
    :math:`A(\alpha,p)X(p,\gamma)B` into a single effective rate, assuming
    equilibrium of :math:`p` and :math:`X`.

  * :func:`DerivedRate <pynucastro.rates.rate.DerivedRate>`: A
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

  There are two important subclasses:

  * :func:`PythonNetwork
    <pynucastro.networks.python_network.PythonNetwork>`: This is a
    collection of rates with functions that know how to write python
    code to express the righthand side of the system of ODEs.

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
  e.g., `Timmes 1999
  <http://adsabs.harvard.edu/abs/1999ApJS..124..241T>`_).  pynucastro
  will create the righthand sides of this system of ODEs (as python or
  C++ code) from the list of rates you provide. One can use this to
  add reaction networks to existing simulation codes, for example, the
  `MAESTROeX <https://amrex-astro.github.io/MAESTROeX/>`_ and `Castro
  <https://amrex-astro.github.io/Castro/>`_ codes.


Data sources
------------

pynucastro can currently read rates from:

* `JINA Reaclib <https://reaclib.jinaweb.org/>`_

* Electron-capture / :math:`\beta`-decay rates from Suzuki et al. 2016

