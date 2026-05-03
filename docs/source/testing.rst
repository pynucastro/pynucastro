******************************
Testing and Comparing Backends
******************************

If we create a network with a set of rates, we should get the same
"ydots" ($dY_k/dt$) values regardless of the backend (to roundoff, and maybe
interpolation, error)

To facilitate comparing the backends, we have a :py:obj:`NetworkCompare <pynucastro.networks.network_compare.NetworkCompare>`
helper class that takes a library and evaluates the ydots for different backends.
Currently it can work with:

* ``RateCollection`` / inline python networks: this uses
  the built-in ``eval`` methods to evaluate each of the rates.

* ``PythonNetwork`` modules: this means writing out the network as a
  ``.py`` file, importing it, and then using the functions in the
  module to evaluate the rates.

* ``AmrexAstroCxxNetwork`` : this uses the ``standalone_build`` option
  to ``write_network()`` to output a driver and makefile.  The output
  is then parsed to get the ydot values.

* ``SimpleCxxNetwork`` : we build a simple test driver and parse the output
  to get the ydot values.

Tests can be run with or without screening.  Once the ``NetworkCompare``
object is setup, the comparison can be run for a single density and temperature
using the :py:func:`evaluate <pynucastro.networks.network_compare.NetworkCompare.evaluate>` method.  The data are then stored in the object.


Accessing comparison data
=========================

The $dY/dt$ data for each network is stored as a ``dict`` keyed by
:py:class:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>` as
``ydots_py_inline``, ``ydots_py_module``, ``ydots_amrex``, and ``ydots_cxx``.

The individual rate data (just the $N_A\langle \sigma v \rangle$ or
equivalent) is stored as a ``dict`` keyed by :py:class:`Rate
<pynucastro.rates.rate.Rate>` as ``rates_py_inline``,
``rates_py_module``, ``rates_amrex``, and ``rates_cxx``.

A summary of the comparison (including errors) can be printed using
:py:func:`print_summary <pynucastro.networks.network_compare.NetworkCompare.print_summary>`.

Current unit tests comparing networks
=====================================

``NetworkCompare`` is used in the following unit tests in
``pynucastro/networks/tests/comparing_nets_tests``:

* ``test_compare_big_net.py`` : this creates a network consisting of
  an $\alpha$-chain and an iron-group, with $(\alpha,p)(p,\gamma)$ and
  $(nn,\gamma)$ rate approximations, derived reverse rates, modified
  rates, and tabular weak rates.  It is modeled off of the `AMReX
  Astro Microphysics ase-iron network
  <https://amrex-astro.github.io/Microphysics/docs/networks.html#ase-iron>`_.
  Note that simple-C++ networks are not currently included.

* ``test_compare_cxx_and_python.py`` : this compares just strong
  reaction rates, but includes derived reverse rates.

* ``test_compare_cxx_and_python_screened.py`` : a simpler version of
  ``test_compare_cxx_and_python.py`` (no derived reverse rates), but
  with screening enabled.

* ``test_compare_temperature_tabular.py`` : this compares the
  evaluation of a ``TemperatureTabularRate``.  Note that simple-C++
  networks are not currently included.

As more features are ported to ``SimpleCxxNetwork``, we should extend
the testing to compare with python rates.

