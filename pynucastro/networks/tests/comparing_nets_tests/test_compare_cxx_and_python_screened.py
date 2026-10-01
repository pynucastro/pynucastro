# this creates the same network as SimpleCxxNetwork and PythonNetwork.
# we then compare the ydots from the C++ net, the python network
# written to a module, and the python network evaluating as a
# RateCollection.  Here we use the Chugunov 2007 screening.

import sys
import warnings
from pathlib import Path

import pytest

from pynucastro.networks.network_compare import NetworkCompare


def _skip_build():
    return sys.platform == "darwin" or sys.platform.startswith("win")


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    @classmethod
    def lib(cls, reaclib_library):
        nuc = ["p", "he4", "c12", "o16", "ne20", "na23", "mg24"]
        lib = reaclib_library.linking_nuclei(nuc)
        return lib

    @pytest.fixture(scope="class")
    @classmethod
    def nc(cls, lib):
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

    @pytest.fixture(scope="class",
                    params=[(2.e8, 1.e9)],
                    ids=["rho2e8-T1e9"])
    @classmethod
    def eval_cond(cls, nc, request):
        # thermodynamic conditions come from the fixture
        # we group them as (rho, T)
        rho, T = request.param

        if not _skip_build():
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                nc.evaluate(rho=rho, T=T)

        return nc

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_ydots(self, eval_cond):
        eval_cond.compare_results(quantity="ydots",
                                  rtol=1.e-6, atol=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_jac(self, eval_cond):
        eval_cond.compare_results(quantity="jac",
                                  rtol=1.e-6, atol=1.e-80)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_rates(self, eval_cond):
        eval_cond.compare_results(quantity="rates",
                                  rtol=1.e-6, atol=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_energy(self, eval_cond):

        # we use a relaxed tolerance here because of differences
        # in constants in simple C++ nets (N_A)
        eval_cond.compare_results(quantity="enuc",
                                  rtol=1.e-6, atol=1.e-30)
        eval_cond.compare_results(quantity="enu_weak",
                                  rtol=1.e-6, atol=1.e-30)
