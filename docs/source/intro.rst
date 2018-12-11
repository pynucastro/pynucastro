Overview of pynucastro
==========================

pynucastro is a set of python interfaces to nuclear reaction rates
(currently via the JINA reaclib database). It is meant for both
interactive exploration of rates (through Jupyter notebooks) and to
create reaction networks for use in simulation codes.

The preferred way of importing pynucastro is as follows:

.. code-block:: python

   import pynucastro as pyrl


The main classes are:

* :func:`Nucleus <pynucastro.rates.rate.Nucleus>`: This is a single
  nucleus.  It knows its proton number, Z, neutron number, N, and
  weight, A.

* :func:`Rate <pynucastro.rates.rate.Rate>`: This is a single rate.  It
  knows the reactants and products and has methods that allow you to
  evaluate it at a specified temperature and plot its temperature
  dependence.

* :func:`Composition
  <pynucastro.networks.rate_collection.Composition>`: This is a
  collection of nuclei and their mass fractions.  A ``Composition`` is
  used when evaluating the full rates in a network.

* :func:`RateCollection <pynucastro.networks.rate_collection.RateCollection>`:
  This is a group of rates.  This acts as the base class for different
  reaction networks.  A ``RateCollection`` has methods to evaluate the
  rates and make a plot of the links between rates.

* :func:`PythonNetwork
  <pynucastro.networks.python_network.PythonNetwork>`: This is a
  collection of rates with functions that know how to write python
  code to express the righthand side of the system of ODEs.

Usage
-----

There are two modes of usage for pynucastro.  Within a Jupyter
notebook, one can evaluate the rates and interactively visualize a
network and see the flow between nuclei as connections colored by the
rate linking them.  

The other usage mode is to use pynucastro to write the righthand
side routine for the system of ODEs that must be integrated to evolve
a reaction network.  A reaction network takes the form:

.. math::

   \frac{dY_i}{dt} = - \sum_{j,k} Y_i Y_j \lambda_{i(j,k)l} + \sum_{j,k} Y_l Y_k \lambda_{l(j,k)i}

where the :math:`\lambda`'s are the rates of destruction and creation
of species i, represented by the molar fraction :math:`Y_i` (see,
e.g., `Timmes 1999
<http://adsabs.harvard.edu/abs/1999ApJS..124..241T>`_).  pynucastro
will create the righthand sides of this system of ODEs (as python or
Fortran code) from the list of rates you provide. One can use this to
add reaction networks to existing simulation codes, for example, the
`Maestro <https://amrex-astro.github.io/MAESTRO/>`_ and `Castro
<https://amrex-astro.github.io/Castro/>`_ codes.
