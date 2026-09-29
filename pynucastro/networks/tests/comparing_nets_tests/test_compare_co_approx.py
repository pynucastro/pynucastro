# test the comparison of the C+C, C+O, and O+O approximation
# across network types

import sys
import warnings
from pathlib import Path

import pytest

from pynucastro.networks.network_compare import NetworkCompare
from pynucastro.rates.aprox_family_rates import make_CO_approx_rates
from pynucastro.rates.library import Library


def _skip_build():
    return sys.platform == "darwin" or sys.platform.startswith("win")


class TestNetworkCompare:

    # pylint: disable=duplicate-code
    @pytest.fixture(scope="class")
    @classmethod
    def lib(cls, reaclib_library):
        crates = make_CO_approx_rates(reaclib_library.get_rates(), "C")
        corates = make_CO_approx_rates(reaclib_library.get_rates(), "CO")
        orates = make_CO_approx_rates(reaclib_library.get_rates(), "O")

        c12ag = reaclib_library.get_rate_by_name("c12(a,g)o16")
        c12ag_reverse = reaclib_library.get_rate_by_name("o16(g,a)c12")
        o16ag = reaclib_library.get_rate_by_name("o16(a,g)ne20")
        o16ag_reverse = reaclib_library.get_rate_by_name("ne20(g,a)o16")
        other_rates = [c12ag, c12ag_reverse, o16ag, o16ag_reverse]

        lib = Library(rates=crates+corates+orates+other_rates)
        return lib

    @pytest.fixture(scope="class")
    @classmethod
    def nc(cls, lib):
        cxx_test_path = Path("_test_compare_coapprox_cxx/")
        amrex_test_path = Path("_test_compare_coapprox_amrex/")

        nc = NetworkCompare(lib,
                            include_amrex=True,
                            include_simple_cxx=True,
                            python_module_name="coapprox_compare.py",
                            amrex_test_path=amrex_test_path,
                            cxx_test_path=cxx_test_path)
        return nc

    @pytest.fixture(scope="class",
                    params=[(2.e6, 1.e9), (2.e9, 4.e9)],
                    ids=["rho2e6-T1e9", "rho2e9-T4e9"])
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

    # pylint: enable=duplicate-code
