import pytest

import pynucastro as pyna


@pytest.fixture(scope="package")
def reaclib_library():
    return pyna.ReacLibLibrary()


@pytest.fixture(scope="package")
def tabular_library():
    return pyna.TabularLibrary()
