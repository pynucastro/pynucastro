"""pynucastro is a python module that interprets the nuclear reaction rates,
including those cataloged by the JINA ReacLib project:

https://groups.nscl.msu.edu/jina/reaclib/db/

It provides both interactive access to the rates, for use in Jupyter
notebooks as well as methods for writing python and C++ nuclear
reaction networks, including the the righthand side and Jacobian
routines.

Several different packages provide the data structures needed to
interpret rates and build networks.

library
-------

library stores Reaclib rate files in Reaclib 1 or 2 formats as well as the
directory `library/tabular` containing tabulated rates. pynucastro
will search these directories for rate files as well as the current
working directory.

Adding Reaclib rate files simply requires downloading them from the
JINA ReacLib project and placing them in the library folder. They must
be in the Reaclib 1 format.

nucdata
-------

nucdata provides:

* Nucleus : a single nucleus, with a descriptive name and its
  properties.

nucdata provides tables of binding energy per nucleon in MeV,
partition function and the number of spin states for nuclides specified
by their number of neutrons N and atomic number Z.

The data for these tables is derived from the Atomic Mass Evaluations
2012, 2016 and 2020. By default, pynucastro uses Atomic Mass Evaluation
2016, and Nucleus Spin Evaluation 2020. Scripts for reading the Atomic
Mass Evaluation tables and generating binding energy tables for pynucastro
are provided in `pynucastro/nucdata/AtomicMassEvaluation`.

networks
--------

networks provides classes and functions for organizing collection of
rates, including:

* Composition : this is a container of Nucleus objects along with
  associated mass fractions.  This is used for evaluating rates
  when doing interactive exploration with pynucastro.

* RateCollection : this is simply a collection of rates.  It will
  determine all of the links between the rates and allow you to
  visualize the dependencies.

* Explorer : an explorer takes a rate composition and composition
  and produces an interactive visualization that can be explored
  in Jupyter notebooks.

* PythonNetwork : this extends the RateCollection to enable output
  of python code that can be used with ODE integrators to solve
  a network.

* BaseCxxNetwork : this extends the RateCollection to enable
  output of C++ code.

* AmrexAstroCxxNetwork : this extends the RateCollection to enable
  output of C++ code that can be used to add a network to the
  hydrodynamics codes Castro and MAESTROeX via the AMReX-Astro
  Microphysics repository.


rates
-----

rates provides classes and functions for interpreting individual
reaction rates, including:

* Rate : a single Reaclib rate, with methods for plotting and
  evaluating it.

* RatePair : a pair of rates representing the corresponding forward
  and reverse rates

* Tfactors : this is a simple container class that holds the various
  temperature powers needed to evaluate a rate.

* RateLibrary: a container class for Reaclib rates that provides a
  high level interface for reading Reaclib-formatted files containing
  one or more rates.

* RateFilter: a class implementing search constraints to look up a
  desired rate or group of rates from a RateLibrary.

screening
---------

screening provides python screening routines for the rates.

templates
---------

templates contains subdirectories for generating AmrexAstroCxxNetwork
C++ files implementing the ODE right hand side, jacobian, and
integration driver routines. pynucastro processes these template files
by replacing tags of the form `<tag>` with generated code specific to
a given choice of reaction rates.

"""

from ._version import version

__version__ = version


import pynucastro.screening
from pynucastro.networks import (AmrexAstroCxxNetwork, BaseCxxNetwork,
                                 Composition, Explorer, PythonNetwork,
                                 RateCollection, StarKillerCxxNetwork,
                                 SympyRates)
from pynucastro.nucdata import Nucleus
from pynucastro.rates import (ApproximateRate, DerivedRate, Library, Rate,
                              RateFilter, ReacLibLibrary, TabularLibrary,
                              Tfactors, list_known_rates, load_rate)
from pynucastro.screening import make_plasma_state, make_screen_factors
