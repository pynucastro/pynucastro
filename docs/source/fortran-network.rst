Fortran Networks
================

A :class:`FortranNetwork <pynucastro.networks.fortran_network.FortranNetwork>` provides Fortran wrappers to a :class:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>`,
in particular, the righthand side and Jacobian functions.  All features
supported by ``SimpleCxxNetwork`` are supported by ``FortranNetwork``.

A Fortran network can be created as:

.. code:: python

   import pynucastro as pyna

   rl = pyna.ReacLibLibrary()
   lib = rl.linking_nuclei(["he4", "c12", "o16"])

   net = pyna.FortranNetwork(libraries=[lib])
   net.write_network()

The generated code is intended to be used to interface with a
hydrodynamics code.  That application code will be responsible for
providing the equation of state and adding any desired temperature /
energy evolution to the network.

.. note::

   A C++17 compiler and Fortran 2003 compiler is required.

In addition to the files output from a ``SimpleCxxNetwork``, the
following additional files are output:

* ``fortran_interface.f90`` : this contains a module ``pynet`` that
  provides the nuclear properties and wrappers (using ``iso_c_binding``)
  for the righthand side function (``rhs_f``) and Jacobian (``jac_f``).

* ``test.f90`` : a simple test driver that outputs the righthand side
  and Jacobian for a single thermodynamic state.

* ``wrapper.cpp`` : a set of C++ functions that convert the Fortran inputs
  into the data structures needed by the ``SimpleCxxNetwork`` functions.

The test driver can be built by simply doing:

.. prompt:: bash

   make

