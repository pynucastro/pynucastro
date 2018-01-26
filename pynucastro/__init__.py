"""pynucastro is a python module that interprets the nuclear reaction rates,
including those cataloged by the JINA ReacLib project:

https://groups.nscl.msu.edu/jina/reaclib/db/

It provides both interactive access to the rates, for use in Jupyter
notebooks as well as methods for writing python and Fortran nuclear
reaction networks, including the the righthand side and Jacobian
routines.

Several different packages provide the data structures needed to
interpret rates and build networks

amemass
-------

amemass provides the 2012 atomic mass evaluation.  The main interfaces
is the AME2012 class, which reads in the data table and provides
methods to get nuclide properties.


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

"""

__version__ = "1.0"

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

