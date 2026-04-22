# this creates the same network as in Python (inline / module) and C++
# (AMReX and simple-C++).  We then compare the ydots across each.
# Note: screening is not considered.

import sys
import warnings
from pathlib import Path

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare
from pynucastro.rates.derived_rate import DerivedRate


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
        nuc = ["p", "he4", "c12", "o16", "ne20", "na23", "mg24"]
        lib = reaclib_library.linking_nuclei(nuc)
        rates_to_derive = lib.backward().get_rates()
        for r in rates_to_derive:
            fr = lib.get_rate_by_nuclei(r.products, r.reactants)
            if fr:
                lib.remove_rate(r)
                d = DerivedRate(source_rate=fr, use_pf=True, use_unreliable_spins=True)
                lib.add_rate(d)
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

    @pytest.mark.skipif(sys.platform == "darwin" or sys.platform.startswith("win"),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare(self, nc):

        # thermodynamic conditions
        rho = 2.e8
        T = 1.e9

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            nc.evaluate(rho=rho, T=T)

        # compare the simple C++, AMReX, and python module nets to the
        # python inline version

        for other in [nc.ydots_cxx, nc.ydots_amrex, nc.ydots_py_module]:
            for nuc in nc.ydots_py_inline:
                assert other[nuc] == approx(nc.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    @pytest.mark.skipif(sys.platform == "darwin" or sys.platform.startswith("win"),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare2(self, nc):

        # thermodynamic conditions
        rho = 2.e7
        T = 4.e9

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            nc.evaluate(rho=rho, T=T)

        # compare the simple C++, AMReX, and python module nets to the
        # python inline version

        for other in [nc.ydots_cxx, nc.ydots_amrex, nc.ydots_py_module]:
            for nuc in nc.ydots_py_inline:
                assert other[nuc] == approx(nc.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)
