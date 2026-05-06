# test the comparison of the C+C, C+O, and O+O approximation
# across network types

import sys
import warnings
from pathlib import Path

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare
from pynucastro.rates.aprox_family_rates import make_CO_approx_rates
from pynucastro.rates.library import Library


def _skip_build():
    return sys.platform == "darwin" or sys.platform.startswith("win")


class TestNetworkCompare:

    # pylint: disable=duplicate-code
    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
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
    def nc(self, lib):
        cxx_test_path = Path("_test_compare_coapprox_cxx/")
        amrex_test_path = Path("_test_compare_coapprox_amrex/")

        nc = NetworkCompare(lib,
                            include_amrex=True,
                            include_simple_cxx=True,
                            python_module_name="coapprox_compare.py",
                            amrex_test_path=amrex_test_path,
                            cxx_test_path=cxx_test_path)
        return nc

    @pytest.fixture(scope="class")
    def eval_cond1(self, nc):
        # thermodynamic conditions
        rho = 2.e6
        T = 1.e9

        if not _skip_build():
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                nc.evaluate(rho=rho, T=T)

        return nc

    @pytest.fixture(scope="class")
    def eval_cond2(self, nc):
        # thermodynamic conditions
        rho = 2.e9
        T = 4.e9

        if not _skip_build():
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                nc.evaluate(rho=rho, T=T)

        return nc

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_ydots(self, eval_cond1):

        # compare the simple C++, AMReX, and python module nets to the
        # python inline version

        for other in [eval_cond1.ydots_cxx, eval_cond1.ydots_amrex, eval_cond1.ydots_py_module]:
            for nuc in eval_cond1.ydots_py_inline:
                assert other[nuc] == approx(eval_cond1.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_rates(self, eval_cond1):

        # compare the simple C++, AMReX, and python module nets to the
        # python inline version

        for other in [eval_cond1.rates_cxx, eval_cond1.rates_amrex, eval_cond1.rates_py_module]:
            for nuc in eval_cond1.rates_py_inline:
                assert other[nuc] == approx(eval_cond1.rates_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_ydots2(self, eval_cond2):

        # compare the simple C++, AMReX, and python module nets to the
        # python inline version

        for other in [eval_cond2.ydots_cxx, eval_cond2.ydots_amrex, eval_cond2.ydots_py_module]:
            for nuc in eval_cond2.ydots_py_inline:
                assert other[nuc] == approx(eval_cond2.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_rates2(self, eval_cond2):

        # compare the simple C++, AMReX, and python module nets to the
        # python inline version

        for other in [eval_cond2.rates_cxx, eval_cond2.rates_amrex, eval_cond2.rates_py_module]:
            for nuc in eval_cond2.rates_py_inline:
                assert other[nuc] == approx(eval_cond2.rates_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    # pylint: enable=duplicate-code
