"""pynucastro is a python library that allows for the creation and
interactive exploration of nuclear reaction networks for astrophysical
environments.  It has methods for selecting rates from various
sources, approximating rates, solving for nuclear statistical
equilibrium, and exporting networks for simulation codes.

Several different submodules provide the data structures needed to
interpret rates and build networks.

``constants``
-------------

:py:mod:`constants <pynucastro.constants>` provides the physical
constants used throughout pynucastro.

``eos``
-------

:py:mod:`eos <pynucastro.eos>` provides core routines needed to build
an equation of state describing stellar matter.

``networks``
------------

:py:mod:`networks <pynucastro.networks>` provides classes and
functions for organizing collection of rates, including:

* :py:obj:`Composition
  <pynucastro.networks.rate_collection.Composition>` : this is a
  container of :py:obj:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>`
  objects along with associated mass fractions.  This is used for
  evaluating rates when doing interactive exploration with pynucastro.

* :py:obj:`RateCollection
  <pynucastro.networks.rate_collection.RateCollection>` : this is
  a collection of nuclei and the reaction rates that link them
  together.  It serves as the most basic network and provides
  methods to evaluate and visualize the rates and more.

* :py:obj:`Explorer <pynucastro.networks.rate_collection.Explorer>` :
  this an explorer takes a ``RateCollection`` and ``Composition`` and
  produces an interactive visualization that can be explored in
  Jupyter notebooks.

* :py:obj:`NSENetwork <pynucastro.networks.nse_network.NSENetwork>` :
  this extends the ``RateCollection`` to allow for solving for the
  nuclear statistical equilibrium state of a collection of nuclei.

and the network classes intended for outputting the righthand side
of the ODE system for use in python or other applications:

* :py:obj:`PythonNetwork
  <pynucastro.networks.python_network.PythonNetwork>` : this extends
  the ``RateCollection`` to enable output of python code that can be used
  with ODE integrators to solve a network.

* :py:obj:`BaseCxxNetwork
  <pynucastro.networks.base_cxx_network.BaseCxxNetwork>` : this extends
  the ``RateCollection`` to enable output of C++ code.  It is not
  intended to be used directly, but rather serves as the base class for
  the other C++ network types.

* :py:obj:`AmrexAstroCxxNetwork
  <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>` :
  this extends the ``RateCollection`` to enable output of C++ code
  that can be used to add a network to the hydrodynamics codes Castro
  and MAESTROeX via the `AMReX-Astro Microphysics
  <https://amrex-astro.github.io/Microphysics>`_ repository.

* :py:obj:`SimpleCxxNetwork
  <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>` : this
  extends the ``RateCollection`` to enable output of basic C++ code that
  can be as the starting point for incorporating a pynucastro network
  into a non-AMReX simulation code.

* :py:obj:`FortranNetwork
  <pynucastro.networks.fortran_network.FortranNetwork>` : this is a
  set of wrappers around ``SimpleCxxNetwork`` that provide a Fortran
  interface.

``neutrino_cooling``
--------------------

:py:mod:`neutrino_cooling <pynucastro.neutrino_cooling>` provides
functions for modeling thermal neutrino losses (photo-, plasma-,
pair-, and Bremsstralung neutrinos).  Included are:

* :py:obj:`NeutrinoCooling
  <pynucastro.neutrino_cooling.neutrino_cooling.NeutrinoCooling>` : a
  simply wrapper class that allows for the visualization of the
  neutrino cooling terms.


``nucdata``
-----------

:py:mod:`nucdata <pynucastro.nucdata>` provides nuclear properties.
While there are a large number of classes here, most interaction
is done through :py:obj:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>`.
The nuclear data is derived from the Atomic Mass Evaluations
2020



``rates``
---------

:py:mod:`rates <pynucastro.rates>` provides classes and functions for interpreting individual
reaction rates.  The main core rate classes are:

* :py:obj:`Rate <pynucastro.rates.rate.Rate>` : the base class for all rates.

* :py:obj:`ReacLibRate <pynucastro.rates.reaclib_rate.ReacLibRate>` : a
  temperature-depended rate from the ReacLib library, including
  methods to visualize and evaluate it.

* :py:obj:`TabularRate <pynucastro.rates.tabular_rate.TabularRate>` :
  a tabulated (temperature and electron density) weak reaction rate),
  with methods to visualize and evaluate it.

* :py:obj:`Tfactors <pynucastro.rates.rate.Tfactors>` : a container
  class that holds common temperature factors.

and collections of rates are a library, described in the following
classes:

* :py:obj:`Library <pynucastro.rates.library.Library>` : a container
  that stores multiple rates and allows for filtering based on a set
  of rules.

* :py:obj:`RateFilter <pynucastro.rates.library.RateFilter>` : a class
  implementing search constraints to look up a desired rate or group
  of rates from a ``Library``.

with some specialized library collections available as:

* :py:obj:`ReacLibLibrary <pynucastro.rates.library.ReacLibLibrary>` :
  read in the entire collection of ReacLib rates.

* :py:obj:`TabularLibrary <pynucastro.rates.library.TabularLibrary>` :
  read in all known tabular reaction rates.

``reduction``
-------------

:py:mod:`reduction <pynucastro.reduction>` provides tools for directed
relation graph with error propagation and sensitivity analysis
reduction of reaction networks.

``screening``
-------------

:py:mod:`screening <pynucastro.screening>` provides electron screening
routines for modifying the reaction rates.

"""

from ._version import version

__version__ = version


import pynucastro.screening
from pynucastro.eos import FermiIntegral
from pynucastro.networks import (AmrexAstroCxxNetwork, BaseCxxNetwork,
                                 Composition, Explorer, FortranNetwork,
                                 NSENetwork, PythonNetwork, RateCollection,
                                 SimpleCxxNetwork, StarKillerCxxNetwork,
                                 SympyRates, network_helper)
from pynucastro.nucdata import Nucleus, get_all_nuclei, get_nuclei_in_range
from pynucastro.rates import (ApproximateRate, DerivedRate, FFNLibrary,
                              LangankeLibrary, Library, ModifiedRate,
                              OdaLibrary, PruetFullerLibrary, Rate, RateFilter,
                              ReacLibLibrary, SuzukiLibrary, TabularLibrary,
                              Tfactors, load_rate)
from pynucastro.reduction import drgep, sens_analysis
from pynucastro.screening import make_plasma_state, make_screen_factors
