# this creates a network with a TemperatureTabularRate and several
# StarLib rates compares to ensure that the RateCollection
# PythonNetwork, and C++ versions all give the same ydots

import sys
import warnings
from pathlib import Path

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare
from pynucastro.rates.alternate_rates import IliadisO16pgF17


def _skip_build():
    return sys.platform == "darwin" or sys.platform.startswith("win")


class TestNetworkCompare:

    @pytest.fixture(scope="class")
    @classmethod
    def lib(cls, reaclib_library, starlib_library):
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
    @classmethod
    def nc(cls, lib):
        amrex_test_path = Path("_test_tt_amrex/")

        nc = NetworkCompare(lib,
                            include_amrex=True,
                            include_simple_cxx=False,
                            python_module_name="net_tt.py",
                            amrex_test_path=amrex_test_path)
        return nc

    @pytest.fixture(scope="class",
                    params=[(2.e8, 1.e9), (2.e7, 4.e9), (2.e6, 1.e8)],
                    ids=["rho2e8-T1e9", "rho2e7-T4e9", "rho1e6-T1e8"])
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
                                  rtol=1.e-11, atol=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_jac(self, eval_cond):
        eval_cond.compare_results(quantity="jac",
                                  rtol=1.e-11, atol=1.e-80)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_rates(self, eval_cond):
        eval_cond.compare_results(quantity="rates",
                                  rtol=1.e-11, atol=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_energy(self, eval_cond):

        # we use a relaxed tolerance here because of differences
        # in constants in simple C++ nets (N_A)
        eval_cond.compare_results(quantity="enuc",
                                  rtol=1.e-7, atol=1.e-30)
        eval_cond.compare_results(quantity="enu_weak",
                                  rtol=1.e-7, atol=1.e-30)
