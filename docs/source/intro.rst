Overview of ``pyreaclib``
=========================

``pyreaclib`` is a set of python interfaces to the JINA reaclib
nuclear reaction rate database.  It is meant for both interactive
exploration of rates (through Jupyter notebooks) and to create
reaction networks for use in simulation codes.

The preferred way of importing ``pyreaclib`` is as:

.. code-block:: python

   import pyreaclib as pyrl


The main classes are:

* :func:`Nucleus <pyreaclib.rates.rate.Nucleus>`

  This is a single nucleus.  It knows its proton number, Z, neutron
  number, N, and weight, A.

* :func:`Rate <pyreaclib.rates.rate.Rate>`

  This is a single rate.  It knows the reactants and products and has
  methods that allow you to evaluate it at a specified temperature and
  plot its temperature dependence.

* :func:`Composition <pyreaclib.networks.rate_collection.Composition>`

* :func:`RateCollection <pyreaclib.networks.rate_collection.RateCollection>`

  This is a group of rates.  This acts as the base class for different
  reaction networks.  A `RateCollection` has methods to evaluate the
  rates and make a plot of the links between rates.

* :func:`PythonNetwork <pyreaclib.networks.python_network.PythonNetwork>`
