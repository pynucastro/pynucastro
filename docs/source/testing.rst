******************************
Testing and Comparing Backends
******************************

If we create a network with a set of rates, we should get the same
"ydots" ($dY_k/dt$) values regardless of the backend (to roundoff, and maybe
interpolation, error)

To facilitate comparing the backends, we have a :py:obj:`NetworkCompare <pynucastro.networks.network_compare.NetworkCompare>`
helper class that takes a library and evaluates the ydots for different backends.
Currently it can work with:

* :py:obj:`RateCollection <pynucastro.networks.rate_collection.RateCollection>` / inline python networks: this uses
  the built-in ``eval`` methods to evaluate each of the rates.

* :py:obj:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>` modules: this means writing out the network as a
  ``.py`` file, importing it, and then using the functions in the
  module to evaluate the rates.

* :py:obj:`AmrexAstroCxxNetwork <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>` : this uses the ``standalone_build`` option
  to ``write_network()`` to output a driver and makefile.  The output
  is then parsed to get the ydot values.

* :py:obj:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>` : we build a simple test driver and parse the output
  to get the ydot values.

Tests can be run with or without screening.  Once the ``NetworkCompare``
object is setup, the comparison can be run for a single density and temperature
using the :py:meth:`evaluate <pynucastro.networks.network_compare.NetworkCompare.evaluate>` method.  The data are then stored in the object.


Accessing comparison data
=========================

* ODE righthand side: The $\partial Y/\partial t$ data for each network is stored as a ``dict`` keyed by
  :py:class:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>` as
  ``ydots_py_inline``, ``ydots_py_module``, ``ydots_amrex``, and ``ydots_cxx``.

* Jacobian: The composition terms of the Jacobian, $\partial \dot{Y}_i/\partial Y_j$,
  are stored in a dict keyed by a tuple ``(nuc_i, nuc_j)``, where ``nuc_i`` and
  ``nuc_j`` are ``Nucleus`` objects representing the row and column nucleus
  respectively.

* Raw rates : The individual rate data (just the $N_A\langle \sigma v
  \rangle$ or equivalent) is stored as a ``dict`` keyed by
  :py:class:`Rate <pynucastro.rates.rate.Rate>` as ``rates_py_inline``,
  ``rates_py_module``, ``rates_amrex``, and ``rates_cxx``.

* Energy : The nuclear energy from mass changes is stored as a scalar,
  ``enuc_py_inline``, ``enuc_py_module``, ``enuc_amrex``, and ``enuc_cxx``.

* Weak rate neutrino losses : The neutrino energy loss from weak rates
  is stored as a scalar, ``enu_weak_py_inline``,
  ``enu_weak_py_module``, ``enu_weak_amrex``, and ``enu_weak_cxx``.


Performing the comparison
=========================

The :py:meth:`compare_results
<pynucastro.networks.network_compare.NetworkCompare.compare_results>`
method manages the comparison.  It accepts a relative and absolute
tolerance and raises a ``ValueError`` if networks do not agree.

A summary of the comparison (including errors) can be printed using
:py:meth:`print_summary <pynucastro.networks.network_compare.NetworkCompare.print_summary>`.

Current unit tests comparing networks
=====================================

``NetworkCompare`` is used in the following unit tests in
``pynucastro/networks/tests/comparing_nets_tests``:

* `test_compare_big_net.py <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/networks/tests/comparing_nets_tests/test_compare_bignet.py>`_ : this creates a network consisting of
  an $\alpha$-chain and an iron-group, with $(\alpha,p)(p,\gamma)$ and
  $(nn,\gamma)$ rate approximations, derived reverse rates, modified
  rates, and tabular weak rates.  It is modeled off of the `AMReX
  Astro Microphysics ase-iron network
  <https://amrex-astro.github.io/Microphysics/docs/networks.html#ase-iron>`_.

* `test_compare_branched.py <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/networks/tests/comparing_nets_tests/test_compare_branched.py>`_ : this compares networks that use
  :py:obj:`BranchedRate <pynucastro.rates.branched_rate.BranchedRate>`
  for a reduced CNO cycle.

* `test_compare_co_approx.py <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/networks/tests/comparing_nets_tests/test_compare_co_approx.py>`_ : this compares the approximate
  rates for carbon and oxygen burning (as discussed in :doc:`co-approximations`).

* `test_compare_cxx_and_python.py <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/networks/tests/comparing_nets_tests/test_compare_cxx_and_python.py>`_ : this compares just strong
  reaction rates, but includes derived reverse rates.

* `test_compare_cxx_and_python_screened.py <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/networks/tests/comparing_nets_tests/test_compare_cxx_and_python_screened.py>`_ : a simpler version of
  ``test_compare_cxx_and_python.py`` (no derived reverse rates), but
  with screening enabled.

* `test_compare_temp_tabular_starlib.py <https://github.com/pynucastro/pynucastro/blob/main/pynucastro/networks/tests/comparing_nets_tests/test_compare_temp_tabular_starlib.py>`_ : this compares the
  evaluation of a :py:obj:`TemperatureTabularRate <pynucastro.rates.temperature_tabular_rate.TemperatureTabularRate>` and several
  :py:obj:`StarLibRate <pynucastro.rates.starlib_rate.StarLibRate>`.  Note that simple-C++ networks are not currently
  included.

As more features are ported to ``SimpleCxxNetwork``, we should extend
the testing to compare with python rates.

