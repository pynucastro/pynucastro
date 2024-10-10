"""The core routines needed to read and interpret nuclear reaction rates"""

__all__ = ["approximate_rates", "rate", "library"]

from .approximate_rates import ApproximateRate
from .known_duplicates import find_duplicate_rates, is_allowed_dupe
from .library import (LangankeLibrary, Library, RateFilter, ReacLibLibrary,
                      SuzukiLibrary, TabularLibrary, list_known_rates)
from .rate import (DerivedRate, Rate, RateFileError, RatePair, ReacLibRate,
                   SingleSet, TableIndex, TableInterpolator, TabularRate,
                   Tfactors, _find_rate_file, load_rate)
