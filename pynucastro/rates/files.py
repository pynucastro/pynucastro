"""Functions for finding and working with the rate data directories.

"""
from pathlib import Path

_pynucastro_dir = Path(__file__).parents[1]
_pynucastro_rates_dir = _pynucastro_dir/"library"
_pynucastro_tabular_dir = _pynucastro_rates_dir/"tabular"
_pynucastro_suzuki_dir = _pynucastro_tabular_dir/"suzuki"
_pynucastro_langanke_dir = _pynucastro_tabular_dir/"langanke"
_pynucastro_ffn_dir = _pynucastro_tabular_dir/"ffn"
_dirs = [
    _pynucastro_dir, _pynucastro_rates_dir, _pynucastro_tabular_dir,
    _pynucastro_suzuki_dir, _pynucastro_langanke_dir, _pynucastro_ffn_dir
]


def get_rates_dir():
    """Return the top-level directory containing rate files

    Returns
    -------
    pathlib.Path

    """

    return _pynucastro_rates_dir


def get_tabular_dir():
    """Return the top-level directory containing the tabulated rates

    Returns
    -------
    pathlib.Path

    """

    return _pynucastro_tabular_dir


class RateFileError(Exception):
    """An error occurred while trying to read a Rate from a file."""


def _find_rate_file(ratename):
    """locate the Reaclib or tabular rate or library file given its name.  Return
    None if the file cannot be located, otherwise return its path."""

    # check to see if the rate file is in the working dir,
    # is already the full path, or is in _dirs

    for path in ("", *_dirs):
        x = Path(path, ratename).resolve()
        if x.is_file():
            return x

    # notify user we can't find the file
    raise RateFileError(f'File {ratename!r} not found in the working directory, {_pynucastro_rates_dir}, or {_pynucastro_tabular_dir}')
