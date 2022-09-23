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

:meth:`starkiller_network <pynucastro.networks.starkiller_network>`:
the support routines to generate a network that can be incorporated
into the StarKiller microphysics routines supported by astrophysical
hydrodynamics codes.

:meth:`starkiller_cxx_network <pynucastro.networks.starkiller_cxx_network>`:
the support routines to generate a C++ network that can be incorporated
into the StarKiller microphysics routines supported by astrophysical
hydrodynamics codes.

"""

#__all__ = ["python_network", "rate_collection", "starkiller_network", "sympy_network_support"]

from .rate_collection import RateCollection, Composition, Explorer
from .python_network import PythonNetwork
from .sympy_network_support import SympyRates
from .base_cxx_network import BaseCxxNetwork
from .starkiller_cxx_network import StarKillerCxxNetwork
