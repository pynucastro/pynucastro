"""pynucastro is a python module that interprets the nuclear reaction rates,
including those cataloged by the JINA ReacLib project:

https://groups.nscl.msu.edu/jina/reaclib/db/

It provides both interactive access to the rates, for use in Jupyter
notebooks as well as methods for writing python and Fortran nuclear
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

nucdata provides tables of binding energy per nucleon in MeV for
nuclides specified by their number of neutrons N and atomic
number Z.

The data for these tables is derived from the Atomic Mass Evaluations
2012 and 2016. By default, pynucastro uses Atomic Mass Evaluation
2016. Scripts for reading the Atomic Mass Evaluation tables and
generating binding energy tables for pynucastro are provided in
`pynucastro/nucdata/AtomicMassEvaluation`.

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

* BaseFortranNetwork : this extends the RateCollection to enable
  output of Fortran code that can be used to integrate the network
  with the included VODE package for ODE integration.

* StarKillerNetwork : this extends the RateCollection to enable output
  of Fortran code that can be used to add a network to the
  hydrodynamics codes Castro and Maestro via the StarKiller
  Microphysics repository.


rates
-----

rates provides classes and functions for interpreting individual
reaction rates, including:

* Nucleus : a single nucleus, with a descriptive name and its
  properties.

* Rate : a single Reaclib rate, with methods for plotting and
  evaluating it.

* Tfactors : this is a simple container class that holds the various
  temperature powers needed to evaluate a rate.

* RateLibrary: a container class for Reaclib rates that provides a
  high level interface for reading Reaclib-formatted files containing
  one or more rates.

* RateFilter: a class implementing search constraints to look up a
  desired rate or group of rates from a RateLibrary.

screening
---------

screening provides routines used by the BaseFortranNetwork to screen
Reaclib reaction rates in the weak, intermediate, and strong
regimes. Tabulated rates are not screened.

The Fortran module in `pynucastro/screening` is only used for the
standalone Fortran network. StarKiller Microphysics networks also use
rate screening, but they use the screening module in the StarKiller
Microphysics repository.

templates
---------

templates contains subdirectories for generating BaseFortranNetwork
and StarKillerNetwork Fortran modules implementing the ODE right hand
side, jacobian, and integration driver routines. pynucastro processes
these template files by replacing tags of the form `<tag>` with
generated code specific to a given choice of reaction rates.

"""

__version__ = "1.1.0"

from pynucastro.networks import \
    RateCollection, \
    Composition, \
    Explorer, \
    PythonNetwork, \
    BaseFortranNetwork, \
    StarKillerNetwork

from pynucastro.rates import \
    Tfactors, \
    Nucleus, \
    Rate, \
    list_known_rates
