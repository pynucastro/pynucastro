Simple C++ network
==================

A ``SimpleCxxNetwork`` is a very basic C++ network.  Currently, the following features
are not supported:

* tabular rates
* screening
* partition functions
* NSE
* temperature evolution
* plasma neutrino losses

A simple driver, ``main.cpp`` and a ``GNUmakefile`` will be generated
that will simply evaluate the righthand side and Jacobian, and the
energy release for a single thermodynamic state.  The generated code
is intended to be used to interface with a hydrodynamics code.  That
code will be responsible for providing the equation of state and
adding any desired temperature / energy evolution to the network.

.. note::

   A C++17 compiler is required

The test driver can be built by simply doing:

.. prompt:: bash

   make

