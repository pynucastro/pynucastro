# unit test for C++ network with derived rates using partition functions
import shutil

import pytest

import pynucastro as pyna


class TestAmrexAstroCxxNetwork:
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library):
        # based on the partition_test network in Microphysics
        nuclei = ["p", "he4", "fe52", "ni56", "co55"]
        lib = reaclib_library.linking_nuclei(nuclist=nuclei,
                                             with_reverse=False)

        fwd_rates = []
        derived = []
        for r in lib.get_rates():
            try:
                d = pyna.DerivedRate(rate=r, compute_Q=False, use_pf=True)
            except ValueError:
                continue
            fwd_rates.append(r)
            derived.append(d)

        fwd_rates_lib = pyna.Library(rates=fwd_rates)
        der_rates_lib = pyna.Library(rates=derived)
        full_lib = fwd_rates_lib + der_rates_lib
        return pyna.AmrexAstroCxxNetwork(libraries=[full_lib])

    def test_write_network(self, fn, compare_network_files):
        """ test the write_network function"""
        test_path = "_test_cxx_derived/"
        # subdirectory of pynucastro/networks/tests/
        reference_path = "_amrexastro_cxx_derived_reference/"
        # files that will be ignored if present in the generated directory
        skip_files = ["pynucastro-info.txt"]

        # remove any previously generated files
        shutil.rmtree(test_path, ignore_errors=True)
        fn.write_network(odir=test_path)
        compare_network_files(test_path, reference_path, skip_files)

        # clean up generated files if the test passed
        shutil.rmtree(test_path)
