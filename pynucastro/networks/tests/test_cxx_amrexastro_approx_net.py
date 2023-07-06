# unit tests for rates
import os.path

import pytest

import pynucastro as pyna
from pynucastro.networks.tests.helpers import compare_network_files


class TestAmrexAstroCxxNetwork:
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library):

        mylib = reaclib_library.linking_nuclei(["mg24", "al27", "si28", "p31", "s32", "he4", "p"])
        net = pyna.AmrexAstroCxxNetwork(libraries=[mylib])
        net.make_ap_pg_approx()
        net.remove_nuclei(["al27", "p31"])
        fn = net
        return fn

    def test_write_network(self, fn):
        """ test the write_network function"""
        test_path = "_test_cxx_approx/"
        reference_path = "_amrexastro_cxx_approx_reference/"
        # files that will be ignored if present in the generated directory
        skip_files = []

        base_path = os.path.relpath(os.path.dirname(__file__))
        reference_path = os.path.join(base_path, reference_path)

        fn.write_network(odir=test_path)
        errors = compare_network_files(test_path, reference_path, skip_files)

        assert not errors, f"files don't match: {' '.join(errors)}"
