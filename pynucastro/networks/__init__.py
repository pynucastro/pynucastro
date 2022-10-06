"""The pynucastro modules that support the creation of networks.
There are several main submodules here:

:meth:`rate_collection <pynucastro.networks.rate_collection>`: this is
simply a collection of rates that knows about the links connecting
nuclei.  This is used as the base for the different classes the write
code to output networks for integration.

:meth:`python_network <pynucastro.networks.python_network>`: the
support routines to generate a full, integrable network in python.

:meth:`base_cxx_network <pynucastro.networks.base_cxx_network>`:
the support routines to generate a standalone integrable network in
pure C++.

:meth:`anrexastro_cxx_network <pynucastro.networks.amrexastro_cxx_network>`:
the support routines to generate a C++ network that can be incorporated
into the AMReX-Astro Microphysics routines supported by astrophysical
hydrodynamics codes.

"""

#__all__ = ["python_network", "rate_collection", "sympy_network_support"]

from .amrexastro_cxx_network import AmrexAstroCxxNetwork
from .base_cxx_network import BaseCxxNetwork
from .python_network import PythonNetwork
from .rate_collection import Composition, Explorer, RateCollection
from .sympy_network_support import SympyRates

StarKillerCxxNetwork = AmrexAstroCxxNetwork
