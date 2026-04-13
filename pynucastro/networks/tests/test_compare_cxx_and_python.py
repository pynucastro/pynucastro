# this creates the same network as SimpleCxxNetwork and PythonNetwork.
# we then compare the ydots from the C++ net, the python network
# written to a module, and the python network evaluating as a
# RateCollection.  Note: screening is not considered.

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

    @pytest.fixture(scope="class")
    def nc(self, lib):
        cxx_test_path = Path("_test_compare_cxx/")
        amrex_test_path = Path("_test_compare_amrex/")

        nc = NetworkCompare(lib,
                            include_amrex=True,
                            include_simple_cxx=True,
                            python_module_name="basic_cxx_py_compare.py",
                            amrex_test_path=amrex_test_path,
                            cxx_test_path=cxx_test_path)
        return nc

    @pytest.mark.skipif(sys.platform.startswith("win"), reason="Does not run on Windows")
    def test_compare(self, nc):

        # thermodynamic conditions
        rho = 2.e8
        T = 1.e9

        nc.evaluate(rho=rho, T=T)

        # compare the simple C++ net to the python inline version

        for nuc in nc.ydots_cxx:
            assert nc.ydots_cxx[nuc] == approx(nc.ydots_py_inline[nuc],
                                               rel=1.e-11, abs=1.e-14)

        # compare the AMReX C++ net to the python inline version

        for nuc in nc.ydots_amrex:
            assert nc.ydots_amrex[nuc] == approx(nc.ydots_py_inline[nuc],
                                                 rel=1.e-11, abs=1.e-14)

        # compare the python module version to the python inline version

        for nuc in nc.ydots_py_module:
            assert nc.ydots_py_module[nuc] == approx(nc.ydots_py_inline[nuc],
                                                     rel=1.e-11, abs=1.e-14)

        # other comparisons -- these should not be needed, since we already
        # compared everything to the python inline version

        # compare the simple C++ net to the python module version

        for nuc in nc.ydots_cxx:
            assert nc.ydots_cxx[nuc] == approx(nc.ydots_py_module[nuc],
                                               rel=1.e-11, abs=1.e-14)

        # compare the simple C++ net to the AMReX C++ net

        for nuc in nc.ydots_cxx:
            assert nc.ydots_cxx[nuc] == approx(nc.ydots_amrex[nuc],
                                               rel=1.e-11, abs=1.e-14)

    @pytest.mark.skipif(sys.platform.startswith("win"), reason="Does not run on Windows")
    def test_compare2(self, nc):

        # thermodynamic conditions
        rho = 2.e7
        T = 4.e9

        nc.evaluate(rho=rho, T=T)

        # compare the simple C++ net to the python inline version

        for nuc in nc.ydots_cxx:
            assert nc.ydots_cxx[nuc] == approx(nc.ydots_py_inline[nuc],
                                               rel=1.e-11, abs=1.e-14)

        # compare the simple C++ net to the python module version

        for nuc in nc.ydots_cxx:
            assert nc.ydots_cxx[nuc] == approx(nc.ydots_py_module[nuc],
                                               rel=1.e-11, abs=1.e-14)

        # compare the python inline and module versions (shouldn't
        # really be needed)

        for nuc in nc.ydots_py_inline:
            assert nc.ydots_py_inline[nuc] == approx(nc.ydots_py_module[nuc],
                                                     rel=1.e-11, abs=1.e-14)
