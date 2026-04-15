*********************************
AMReX Astro Microphysics Networks
*********************************

pynucastro can generate a C++ network that works in the
AMReX-Astro `Microphysics
<https://github.com/amrex-astro/Microphysics>`_.  This is done through
the :class:`AmrexAstroCxxNetwork <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>` class.  A simple 
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

In :class:`BaseCxxNetwork <pynucastro.networks.base_cxx_network.BaseCxxNetwork>`, there is a dictionary of function tags and the corresponding function
called ``ftags``.  Any network class that derives from this can add to this list to specialize
the output.

When the network is written, each template file is read and output line by line, injecting
the code output by a function if the line has a function tag in it.

For the generation of the righthand side of the network itself, we use
SymPy to build up an expression for each term and then use SymPy's
`cxxcode <https://docs.sympy.org/latest/modules/printing.html#sympy.printing.codeprinter.cxxcode>`_ function to translate it into C++ code.  This is managed by the :class:`SympyRates <pynucastro.networks.sympy_network_support.SympyRates>` class.


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

* ``derived_rates.H``

  This computes any rates that are derived via detailed balance,
  with a function provided for each rate.

* ``interp_tools.H``

  This contains interpolation functions used by tabulated rates.

* ``partition_functions.H``

  This holds the partition function data (inline) and the functions
  to interpolate it.

* ``reaclib_rates.H``

  This computes the ReacLib reaction rates, with a function provided
  for each rate.

* ``table_rates.H``

  This manages ``TabularRate`` rates and interpolating the
  data.  The rate data itself is stored inline in the header file.

  .. warning::

     We do not check if the thermodynamic conditions for evaluating the
     rate fall outside of the tabulation.  Instead we simply extrapolate
     from the edge of the table.  See :ref:`tabulated_rate_sources` for
     the bounds of the data.

* ``temperature_table_rates``

  This manages ``TemperatureTabularRate`` rates and interpolating
  their data.  The rate data itself is stored inline in the header
  file.

* ``tfactors.H``

  This stores the ``struct`` that holds the different temperature factors
  used in the reaction rates.

Finally, there are a few meta-data files:

* ``pynucastro.net``

  This lists the nuclei in the network and their properties in a format that
  Microphysics requires.

* ``Make.package``

  This is a stub for the build system that list the files needed to build
  the network.

* ``_parameters``

  This lists the runtime parameters that affect the network.

Using the network
=================

To use the network, there are 2 options.

* You can use the unit tests in ``Microphysics/unit_test/``, for example,
  ``burn_cell`` will do a single-zone burn.

* You can add the argument ``standalone_build=True`` to the
  ``write_network()`` function to output a basic ``main.cpp`` and ``GNUmakefile``.
  This still requires a copy of ``Microphysics``, and expect the environment
  variable ``MICROPHYSICS_HOME`` to point to its root directory.  You can then
  build simply as:

  .. prompt:: bash

     make

  and run as:

  .. prompt:: bash

     ./main3d.gnu.ex testing.density=1.e6 testing.temperature=1.e9

  where the density and temperature are set via runtime parameters as shown.
  To control the screening, use the ``SCREEN_METHOD`` build parameter from
  Microphysics.  For example, to disable screening, we can do:

  .. prompt:: bash

     make SCREEN_METHOD=null
