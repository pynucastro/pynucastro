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


def _skip_build():
    return sys.platform == "darwin" or sys.platform.startswith("win")


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

    @pytest.fixture(scope="class")
    def eval_cond1(self, nc):
        # thermodynamic conditions
        rho = 2.e8
        T = 1.e9

        if not _skip_build():
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                nc.evaluate(rho=rho, T=T)

        return nc

    @pytest.fixture(scope="class")
    def eval_cond2(self, nc):
        # thermodynamic conditions
        rho = 2.e7
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
