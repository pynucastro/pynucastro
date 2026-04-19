# this creates a network with a TemperatureTabularRate and several
# StarLib rates compares to ensure that the RateCollection
# PythonNetwork, and C++ versions all give the same ydots

import sys
from pathlib import Path

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare
from pynucastro.rates.alternate_rates import IliadisO16pgF17


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library, starlib_library):
        nuc = ["p", "n15", "o16"]
        lib = reaclib_library.linking_nuclei(nuc, with_reverse=False)

        r = IliadisO16pgF17()
        lib.add_rate(r)

        r2 = starlib_library.get_rate_by_name("c12(a,g)o16")
        lib.add_rate(r2)

        r3 = starlib_library.get_rate_by_name("c12(p,g)n13")
        lib.add_rate(r3)

        r4 = starlib_library.get_rate_by_name("n13(a,p)o16")
        lib.add_rate(r4)

        return lib

    @pytest.fixture(scope="class")
    def nc(self, lib):
        amrex_test_path = Path("_test_tt_amrex/")

        nc = NetworkCompare(lib,
                            include_amrex=True,
                            include_simple_cxx=False,
                            python_module_name="net_tt.py",
                            amrex_test_path=amrex_test_path)
        return nc

    @pytest.mark.skipif(sys.platform == "darwin" or sys.platform.startswith("win"),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare(self, nc):

        # thermodynamic conditions
        # we set the composition to be uniform for all tests
        rho = 1.e5
        T = 2.e8

        nc.evaluate(rho=rho, T=T)

        # compare the inline and module python versions

        # compare the AMReX and python module nets to the python
        # inline version

        for other in [nc.ydots_amrex, nc.ydots_py_module]:
            for nuc in nc.ydots_py_inline:
                assert other[nuc] == approx(nc.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    @pytest.mark.skipif(sys.platform == "darwin" or sys.platform.startswith("win"),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare2(self, nc):

        # thermodynamic conditions
        # we set the composition to be uniform for all tests
        rho = 2.1e6
        T = 3.95e8

        nc.evaluate(rho=rho, T=T)

        # compare the inline and module python versions

        # compare the AMReX and python module nets to the python
        # inline version

        for other in [nc.ydots_amrex, nc.ydots_py_module]:
            for nuc in nc.ydots_py_inline:
                assert other[nuc] == approx(nc.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)
