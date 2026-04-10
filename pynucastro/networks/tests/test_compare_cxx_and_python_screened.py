# this creates the same network as SimpleCxxNetwork and PythonNetwork.
# we then compare the ydots from the C++ net, the python network
# written to a module, and the python network evaluating as a
# RateCollection.  Here we use the Chugunov 2007 screening.

import sys
from pathlib import Path

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
        nuc = ["p", "he4", "c12", "o16", "ne20", "na23", "mg24"]
        lib = reaclib_library.linking_nuclei(nuc)
        return lib

    @pytest.mark.skipif(sys.platform.startswith("win"), reason="Does not run on Windows")
    def test_compare(self, lib):

        test_path = Path("_test_compare_screened/")

        # thermodynamic conditions
        # we set the composition to be uniform for all tests
        rho = 2.e8
        T = 1.e9

        nc = NetworkCompare(lib,
                            use_screening=True,
                            include_simple_cxx=True,
                            python_module_name="screened_cxx_py_compare.py",
                            cxx_test_path=test_path)
        nc.evaluate(rho=rho, T=T)

        # compare the simple C++ net to the python inline version

        for nuc in nc.ydots_cxx:
            assert nc.ydots_cxx[nuc] == approx(nc.ydots_py_inline[nuc],
                                               rel=1.e-6, abs=1.e-14)

        # compare the simple C++ net to the python module version

        for nuc in nc.ydots_cxx:
            assert nc.ydots_cxx[nuc] == approx(nc.ydots_py_module[nuc],
                                               rel=1.e-6, abs=1.e-14)

        # compare the python inline and module versions (shouldn't
        # really be needed)

        for nuc in nc.ydots_py_inline:
            assert nc.ydots_py_inline[nuc] == approx(nc.ydots_py_module[nuc],
                                                     rel=1.e-11, abs=1.e-14)
