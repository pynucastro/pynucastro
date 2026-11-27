Simple C++ network
==================

A :class:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>` is a very basic C++ network.  It supports
ReacLib rates, approximate rates, and screening (via
the method of :cite:t:`chugunov:2007`).

.. important::

   Currently, the following features are not supported:

   * an NSE solver
   * thermal / plasma neutrino losses
   * StarLib and ``TemperatureTabular`` rates

A simple C++ network can be created as:

.. code:: python

   import pynucastro as pyna

   rl = pyna.ReacLibLibrary()
   lib = rl.linking_nuclei(["he4", "c12", "o16"])

   net = pyna.SimpleCxxNetwork(libraries=[lib])
   net.write_network()

.. note::

   The ``SimpleCxxNetwork`` outputs the righthand side function
   (providing $dY/dt$ and the energy loss due to weak rate neutrinos)
   and Jacobian.  It is meant to be used in an application code that
   provides its own time integrator.

   There is no energy/temperature evolution, but the application code
   can augment the set of equations being integrated with an energy
   equation as needed.

.. note::

   A C++20 compiler is required

.. important::

   Derivatives of the screening factor with respect to composition
   are not currently implemented in the Jacobian.


This will output the following files:

* ``actual_network.H`` : this provides the nuclei masses, enums to
  index the rates, and a vector of strings that give the rate names.

* ``actual_rhs.H`` : this provides the righthand side function and Jacobian.

* ``amrex_bridge.H`` : a header that defines some of the basic types that we
  use to store information.  It is derived from the AMReX library
  since we reuse some of the ``AmrexAstroCxxNetwork`` code to create a
  ``SimpleNetwork``.

* ``autodiff.H`` : dual numbers used to compute screening composition derivatives.

* ``approximate_rates.H`` : the functions that evaluate ``AproximateRate`` rates.

* ``burn_type.H`` : this is a simple struct, ``burn_t``, that holds
  thermodynamic data.  An application code can add more members to the
  struct, as needed, to store additional data.

* ``derived_rates.H`` : the functions that evaluate rates recomputed via
  detailed balance.

* ``fundamental_constants.H`` : this provides the fundamental constants
  needed throughout the network.

* ``GNUmakefile`` : a GNU makefile to build the test program.

* ``interp_tools.H`` : functions used for interpolating partition functions
  and rates.

* ``main.cpp`` : a simple driver that simply evaluates the righthand side and Jacobian
  for a single thermodynamic state and compute the energy release.

* ``modified_rates.H`` : the functions that evaluate ``ModifiedRate`` rates.

* ``network_properties.H`` : a header providing the properties of the nuclei.

* ``partition_functions.H`` : data and functions used for interpolating
  partition functions.

* ``rate_type.H`` : the core ``struct`` used for holding rate evaluations

* ``reaclib_rates.H`` : the functions that evaluate the ReacLib reaction rates.

* ``screen_data.H`` : this defines some functions used in precomputing common
  screening factors.

* ``screen.H`` : the actual implementation of the Chugunov 2007 screening.

* ``table_rates.H`` : the data and functions for interpolating tabular rates.

* ``tfactors.H`` : a struct that stores the various temperature powers needed
  to compute reaction rates.


Screening derivatives
---------------------

With ``SCREENING`` enabled, evaluating rates into ``rate_derivs_t`` also
computes ``dlog_screen_dYe`` and ``dlog_screen_dS2`` for each screening pair.
The independent moments are $Y_e = \sum_i Z_i Y_i$ and
$S_2 = \sum_i Z_i^2 Y_i$; temperature, density, and the reacting nuclei
are held fixed.  The Chugunov 2007 prescription depends only on $Y_e$, so
its $S_2$ partial is zero.  Ordinary ``rate_t`` evaluations use scalar arithmetic.

``screen.H`` provides ``plasma_state_t<screening_dual_t>`` and
``fill_plasma_state`` to seed these moments.  Calling ``actual_log_screen``
with this state returns a dual number whose ``value`` is the logarithm of
the screening enhancement and whose ``deriv.d1`` and ``deriv.d2`` are the
$Y_e$ and $S_2$ partials.  ``actual_screen`` returns the enhancement itself
and propagates the derivatives through its exponential and clipping.

These partials are stored for use in constructing composition derivatives;
the generated Jacobian does not yet include the screening contribution.


Test driver
-----------

A test driver can be built by simply doing:

.. prompt:: bash

   make

and run as:

.. prompt:: bash

   ./main

This will simply evaluate the righthand side and Jacobian at a single
thermodynamic state and print out the output, including the $dY/dt$
terms, the individual rates evaluations (without density or
composition terms) and the weak-rate neutrino energy loss term.  It is
intended to show how to work with the interfaces provided by
pynucastro.

To change the density and temperature, you can pass values as arguments
(both must be specified if you are overriding the defaults), e.g.:

.. prompt:: bash

   ./main 100 1.e7

would use $\rho = 100~\mathrm{g~cm^{-3}}$ and $T = 10^7~\mathrm{K}$.

By default screening is included.  To disable screening, compile as:

.. prompt:: bash

   make DISABLE_SCREENING=TRUE
