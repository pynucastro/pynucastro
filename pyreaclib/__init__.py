"""
pyreaclib is a python module that interprets the nuclear reaction rates
cataloged by the JINA ReacLib project:

https://groups.nscl.msu.edu/jina/reaclib/db/

It provides both interactive access to the rates, for use in Jupyter
notebooks as well as methods for writing python and Fortran nuclear
reaction networks, including the the righthand side and Jacobian
routines.
"""

__version__ = "1.0"

from pyreaclib.networks import \
    RateCollection, \
    Composition, \
    Explorer, \
    PythonNetwork, \
    BaseFortranNetwork, \
    BoxLibNetwork, \
    SundialsNetwork
