import pynucastro as pyna
import pytest


@pytest.fixture(scope="package")
def reaclib_library():
    return pyna.ReacLibLibrary()
