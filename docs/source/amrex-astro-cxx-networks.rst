*********************************
AMReX Astro Microphysics Networks
*********************************

pynucastro can generate a C++ network that works in the
AMReX-Astro `Microphysics
<https://github.com/amrex-astro/Microphysics>`_.  This is done through
the ``AmrexAstroCxxNetwork`` class.  A simple 
C++ 3-alpha network that works with ``Microphysics`` can be created via:

.. code:: python

   import pynucastro as pyna

   reaclib_library = pyna.ReacLibLibrary()

   mylib = reaclib_library.linking_nuclei(["he4", "c12", "o16"])

   net = pyna.AmrexAstroCxxNetwork(libraries=[mylib], rate_params=[r])
   net.write_network()

The C++ network uses a set of template files in
``pynucastro/templates/amrexastro-cxx-microphysics/`` which have
placeholders of the form ``<func>(num)``.  Here,

* ``func`` is a tag that corresponds to a function that will be called
  to write the output in this location

* ``num`` is the number of indentations for that line of code

In ``BaseCxxNetwork``, there is a dictionary of function tags and the corresponding function
called ``ftags``.  Any network class that derives from this can add to this list to specialize
the output.

When the network is written, each template file is read and output line by line, injecting
the code output by a function if the line has a function tag in it.

For the generation of the righthand side of the network itself, we use
SymPy to build up an expression for each term and then use SymPy's
``cxxcode`` function to translate it into C++ code.


This will directly write out the C++ code into a collection of headers
and source files.  These are:

* ``actual_network.H``

  This header simply stores the number of rates, and separately the number of ReacLib and
  tabular rates, as well as defines the arrays that hold the rate data.

* ``actual_rhs.H``

  This defines the function to compute the energy release from the
  network as well as evaluate the rates and construct the "ydot" and
  Jacobian terms.  Function tags are used to compute the screening
  factors, compute the tabular rates, and build the ydot and Jacobian,
  and include any plasma neutrino losses.

* ``reaclib_rates.H``

  This computes the ReacLib reaction rates, with a function provided
  for each rate.

* ``table_rates.H``

  This manages reading in tabular rates and interpolating the data.
  It has function tags to define how many tables
  there are as well as declare the memory for storing the tables.

  .. warning::

     We do not check if the thermodynamic conditions for evaluating the
     rate fall outside of the tabulation.  Instead we simply extrapolate
     from the edge of the table.  See :ref:`tabulated_rate_sources` for
     the bounds of the data.

* ``tfactors.H``

  This stores the ``struct`` that holds the different temperature factors
  used in the reaction rates.

There are 2 C++ files that are essentially used to define the global arrays.

* ``actual_network_data.cpp``

* ``table_rates_data.cpp``

Finally, there are a few meta-data files:

* ``pynucastro.net``

  This lists the nuclei in the network and their properties in a format that
  Microphysics requires.

* ``Make.package``

  This is a stub for the build system that list the files needed to build
  the network.

* ``_parameters``

  This lists the runtime parameters that affect the network.

To test the library, you can use the unit tests in ``Microphysics/unit_test/``, for example,
``test_react``.


