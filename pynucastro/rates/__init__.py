"""The core routines needed to read and interpret nuclear reaction rates"""

from pynucastro.nucdata import UnsupportedNucleus

from .approximate_rates import ApproximateRate, create_double_neutron_capture
from .derived_rate import DerivedRate
from .files import RateFileError, _find_rate_file
from .known_duplicates import find_duplicate_rates, is_allowed_dupe
from .library import (FFNLibrary, LangankeLibrary, Library, OdaLibrary,
                      PruetFullerLibrary, RateFilter, ReacLibLibrary,
                      SuzukiLibrary, TabularLibrary, list_known_rates)
from .modified_rate import ModifiedRate
from .rate import BaryonConservationError, Rate, RatePair, Tfactors
from .reaclib_rate import ReacLibRate, SingleSet
from .tabular_rate import TableIndex, TableInterpolator, TabularRate
from .temperature_tabular_rate import TemperatureTabularRate


def load_rate(rfile=None):
    """Try to load a rate of any type.

    Parameters
    ----------
    rfile : str
        the name of a file containing the reaction rate parameterization.

    Returns
    -------
    Rate

    Raises
    ------
    :py:obj:`pynucastro.rates.files.RateFileError`,
    :py:obj:`pynucastro.nucdata.nucleus.UnsupportedNucleus`

    """

    try:
        rate = TabularRate(rfile=rfile)
    except (AttributeError, RateFileError, UnsupportedNucleus):
        rate = ReacLibRate(rfile=rfile)

    return rate
