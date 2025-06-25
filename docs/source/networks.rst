Exporting Networks
==================

pynucastro can output the righthand side functions and Jacobian for
ODE integration in Python, pure C++, or in a form for the
`AMReX-Astrophysics Microphysics
<https://github.com/amrex-astro/Microphysics>`_ package in C++.

To create the python or C++ code for a network, you start by building
the library representing the rates that you want and then you construct
the appropriate network class:

* :class:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>` for a pure python network.

* :class:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>` for a basic, pure C++ network.

* :class:`FortranNetwork <pynucastro.networks.fortran_network.FortranNetwork>` a set of Fortran interfaces for a ``SimpleCxxNetwork``.

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
     net.write_network("triple_alpha.py")

  and we can replace ``PythonNetwork`` with one of the other network types to get the 
  output in a different format.

* Run your python script

  .. prompt:: bash

     python mynet.py
