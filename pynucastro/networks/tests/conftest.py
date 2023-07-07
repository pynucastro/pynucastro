import pytest

import pynucastro as pyna

# enable assert rewriting for the helper functions in ./helpers.py
pytest.register_assert_rewrite('pynucastro.networks.tests.helpers')


@pytest.fixture(scope="package")
def reaclib_library():
    return pyna.ReacLibLibrary()


@pytest.fixture(scope="package")
def tabular_library():
    return pyna.TabularLibrary()
