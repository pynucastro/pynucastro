# unit tests for rates
import shutil

import pytest

from pynucastro import networks
from pynucastro.networks.tests.helpers import compare_network_files


class TestSimpleCxxNetwork:
    # pylint: disable=protected-access
    @pytest.fixture(scope="class")
    def fn(self):
        files = ["c12-c12a-ne20-cf88",
                 "c12-c12n-mg23-cf88",
                 "c12-c12p-na23-cf88",
                 "c12-ag-o16-nac2",
                 "na23--ne23-toki",
                 "ne23--na23-toki",
                 "n--p-wc12"]

        fn = networks.SimpleCxxNetwork(files)
        return fn

    def test_write_network(self, fn):
        """ test the write_network function"""
        test_path = "_test_simple_cxx/"
        # subdirectory of pynucastro/networks/tests/
        reference_path = "_simple_cxx_reference/"
        # files that will be ignored if present in the generated directory
        skip_files = []

        # remove any previously generated files
        shutil.rmtree(test_path, ignore_errors=True)
        fn.write_network(odir=test_path)
        compare_network_files(test_path, reference_path, skip_files)

        # clean up generated files if the test passed
        shutil.rmtree(test_path)
