# unit tests for rates
import shutil

import pytest

import pynucastro as pyna


class TestAmrexAstroCxxNetwork:
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library):

        mylib = reaclib_library.linking_nuclei(["mg24", "al27", "si28", "p31", "s32", "he4", "p"])
        net = pyna.AmrexAstroCxxNetwork(libraries=[mylib])
        net.make_ap_pg_approx()
        net.remove_nuclei(["al27", "p31"])
        fn = net
        return fn

    def test_write_network(self, fn, compare_network_files):
        """ test the write_network function"""
        test_path = "_test_cxx_approx/"
        # subdirectory of pynucastro/networks/tests/
        reference_path = "_amrexastro_cxx_approx_reference/"
        # files that will be ignored if present in the generated directory
        skip_files = []

        # remove any previously generated files
        shutil.rmtree(test_path, ignore_errors=True)
        fn.write_network(odir=test_path)
        compare_network_files(test_path, reference_path, skip_files)

        # clean up generated files if the test passed
        shutil.rmtree(test_path)
