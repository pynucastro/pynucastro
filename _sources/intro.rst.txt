Overview of ``pynucastro``
==========================

``pynucastro`` is a set of python interfaces to the JINA reaclib
nuclear reaction rate database.  It is meant for both interactive
exploration of rates (through Jupyter notebooks) and to create
reaction networks for use in simulation codes.

The preferred way of importing ``pynucastro`` is as:

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
  reaction networks.  A `RateCollection` has methods to evaluate the
  rates and make a plot of the links between rates.

* :func:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>`:

Usage
-----

An example notebook showing the basic use of these classes is here:

:ref:`pynucastro-examples.ipynb`
