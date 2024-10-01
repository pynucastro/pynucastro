Simple C++ network
==================

A ``SimpleCxxNetwork`` is a very basic C++ network.


.. important::

   Currently, the following features are not supported:

   * tabular rates
   * screening
   * partition functions
   * NSE
   * temperature evolution
   * plasma neutrino losses

A simple C++ network can be created as:

.. code:: python

   import pynucastro as pyna

   rl = pyna.ReacLibLibrary()
   lib = rl.linking_nuclei(["he4", "c12", "o16"])

   net = pyna.SimpleCxxNetwork(libraries=[lib])
   net.write_network()

The generated code is intended to be used to interface with a
hydrodynamics code.  That application code will be responsible for
providing the equation of state and adding any desired temperature /
energy evolution to the network.

.. note::

   A C++17 compiler is required


This will output the following files:

* ``actual_network_data.cpp`` : this contains the initialization
  function, ``actual_network_init()`` and the definitions of the mass
  and binding energy arrays.

* ``actual_network.H`` : this provides enums to index the rates and a
  vector of strings that give the rate names.

* ``actual_rhs.H`` : this provides the righthand side function and Jacobian.

* ``amrex_bridge.H`` : a header that defines some of the basic types that we
  use to store information.  It is derived from the AMReX library
  since we reuse some of the ``AmrexAstroCxxNetwork`` code to create a
  ``SimpleNetwork``.

* ``burn_type.H`` : this is a simple struct, ``burn_t``, that holds
  thermodynamic data.  An application code can add more members to the
  struct, as needed, to store additional data.

* ``fundamental_constants.H`` : this provides the fundamental constants
  needed throughout the network.

* ``GNUmakefile`` : a GNU makefile to build the test program.

* ``main.cpp`` : a simple driver that simply evaluates the righthand side and Jacobian
  for a single thermodynamic state and compute the energy release.

* ``network_properties.H`` : a header providing the properties of the nuclei.

* ``reaclib_rates.H`` : the functions that evaluate the ReacLib reaction rates.

* ``tfactors.H`` : a struct that stores the various temperature powers needed
  to compute reaction rates.


The test driver can be built by simply doing:

.. prompt:: bash

   make

