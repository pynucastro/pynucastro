# unit tests for rates
import shutil

import pytest

from pynucastro import networks
from pynucastro.rates import DerivedRate


class TestSimpleCxxNetwork:
    @pytest.fixture(scope="class")
    def fn(self, reaclib_library):
        reaclib_rate_names = ["c12(c12,a)ne20",
                              "c12(c12,n)mg23",
                              "c12(c12,p)na23",
                              "c12(a,g)o16",
                              "ne20(a,g)mg24",
                              "mg23(n,g)mg24",
                              "n(,)p"]
        reaclib_rates = reaclib_library.get_rate_by_name(reaclib_rate_names)

        rev_rates = []
        for r in reaclib_rates:
            if not (r.derived_from_inverse or r.weak):
                # Don't use partition function as its not supported for now
                rev_rates.append(DerivedRate(source_rate=r, use_pf=False,
                                             use_unreliable_spins=False))

        rates = reaclib_rates + rev_rates
        fn = networks.SimpleCxxNetwork(rates=rates)
        return fn

    def test_write_network(self, fn, compare_network_files):
        """ test the write_network function"""
        test_path = "_test_simple_cxx/"
        # subdirectory of pynucastro/networks/tests/
        reference_path = "_simple_cxx_reference/"
        # files that will be ignored if present in the generated directory
        skip_files = ["pynucastro-info.txt"]

        # remove any previously generated files
        shutil.rmtree(test_path, ignore_errors=True)
        fn.write_network(odir=test_path)
        compare_network_files(test_path, reference_path, skip_files)

        # clean up generated files if the test passed
        shutil.rmtree(test_path)
