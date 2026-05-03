# this creates the same network as SimpleCxxNetwork and PythonNetwork.
# we then compare the ydots from the C++ net, the python network
# written to a module, and the python network evaluating as a
# RateCollection.  Here we use the Chugunov 2007 screening.

import sys
import warnings
from pathlib import Path

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare


def _skip_build():
    return sys.platform == "darwin" or sys.platform.startswith("win")


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
        nuc = ["p", "he4", "c12", "o16", "ne20", "na23", "mg24"]
        lib = reaclib_library.linking_nuclei(nuc)
        return lib

    @pytest.fixture(scope="class")
    def nc(self, lib):
        cxx_test_path = Path("_test_compare_cxx_screened/")
        amrex_test_path = Path("_test_compare_amrex_screened/")

        nc = NetworkCompare(lib,
                            use_screening=True,
                            include_amrex=True,
                            include_simple_cxx=True,
                            python_module_name="screened_cxx_py_compare.py",
                            amrex_test_path=amrex_test_path,
                            cxx_test_path=cxx_test_path)
        return nc

    @pytest.fixture(scope="class")
    def eval_cond(self, nc):
        # thermodynamic conditions
        # we set the composition to be uniform for all tests
        rho = 2.e8
        T = 1.e9

        if not _skip_build():
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                nc.evaluate(rho=rho, T=T)

        return nc

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_ydots(self, eval_cond):

        # compare the simple C++, AMReX, and python module nets to the
        # python inline version

        for other in [eval_cond.ydots_cxx, eval_cond.ydots_amrex, eval_cond.ydots_py_module]:
            for nuc in eval_cond.ydots_py_inline:
                assert other[nuc] == approx(eval_cond.ydots_py_inline[nuc],
                                            rel=1.e-6, abs=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_rates(self, eval_cond):

        # compare the simple C++, AMReX, and python module nets to the
        # python inline version

        for other in [eval_cond.rates_cxx, eval_cond.rates_amrex, eval_cond.rates_py_module]:
            for nuc in eval_cond.rates_py_inline:
                assert other[nuc] == approx(eval_cond.rates_py_inline[nuc],
                                            rel=1.e-6, abs=1.e-30)
