"""The core routines needed to read and interpret nuclear reaction rates"""

from pynucastro.nucdata import UnsupportedNucleus

from .approximate_rates import ApproximateRate, create_double_neutron_capture
from .derived_rate import DerivedRate
from .files import RateFileError, _find_rate_file
from .known_duplicates import find_duplicate_rates, is_allowed_dupe
from .library import (FFNLibrary, LangankeLibrary, Library, RateFilter,
                      ReacLibLibrary, SuzukiLibrary, TabularLibrary,
                      list_known_rates)
from .modified_rate import ModifiedRate
from .rate import BaryonConservationError, Rate, RatePair, Tfactors
from .reaclib_rate import ReacLibRate, SingleSet
from .tabular_rate import TableIndex, TableInterpolator, TabularRate


def load_rate(rfile=None):
    """Try to load a rate of any type.

    :raises: :class:`.RateFileError`, :class:`.UnsupportedNucleus`
    """

    try:
        rate = TabularRate(rfile=rfile)
    except (RateFileError, UnsupportedNucleus):
        rate = ReacLibRate(rfile=rfile)

    return rate
