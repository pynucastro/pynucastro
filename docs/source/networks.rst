Creating Networks
=================

pynucastro can output the righthand side functions for ODE integration
in Python, pure C++, or in a form for the `AMReX-Astrophysics Microphysics <https://github.com/amrex-astro/Microphysics>`_ package in C++.

To create the python or C++ code for a network, you start by building
the library representing the rates that you want and then you construct
the appropriate network class:

* :class:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>` for a pure python network.

* :class:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>` for a basic, pure C++ network.

* :class:`AmrexAstroCxxNetwork <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>` for an AMReX-Astro C++ network.

Usually the steps are something like:

* Create a working directory for your network (this does not have to
  be in the ``pynucastro/`` directory tree).

* Write a short python script to generate your network,
  e.g. ``mynet.py``.  As an example, for a network
  with just He4, C12, and O16, we would do:

  .. code:: python

     import pynucastro as pyna

     rl = pyna.ReacLibLibrary()
     lib = rl.linking_nuclei(["he4", "c12", "o16"])

     net = pyna.PythonNetwork(libraries=[lib])
     net.write_network()

  and we can replace ``PythonNetwork`` with one of the other network types to get the 
  output in a different format.

* Run your python script

  .. prompt:: bash

     python mynet.py

Python network
--------------

A ``PythonNetwork`` will output a python module with all the data needed for constructing
the righthand side and Jacobian of the ODE system.  To integrate a Python
network using the SciPy integration routines, customize
`examples/CNO/burn.py <https://github.com/pynucastro/pynucastro/blob/main/examples/CNO/burn.py>`_ to initialize and run your network using the
right hand side module you generated above.


Simple C++ network
------------------

A simple C++ network is a very basic C++ network.  Currently, the following features
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


AMReX-Astro Microphysics network
--------------------------------

An ``AmrexAstroCxxNetwork`` is intended to be used with an AMReX-Astro
code like `Castro <https://github.com/amrex-astro/Castro>`_ or `MAESTROeX <https://github.com/amrex-astro/MAESTROeX>`_.
All pynucastro features are supported and these networks can run on GPUs in the AMReX Astro framework.

The `examples/triple-alpha/triple-alpha-cxx.py
<https://github.com/pynucastro/pynucastro/blob/main/examples/triple-alpha/triple-alpha-cxx.py>`_
example builds the right hand side, Jacobian, and helper C++ modules
to copy into the ``networks/`` subdirectory of the AMReX-Astro
Microphysics repository.

No additional customization is required after running the steps for
all networks above.
