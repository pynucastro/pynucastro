# this creates a network with a TemperatureTabularRate and compares to
# ensure that both the RateCollection version and the PythonNetwork
# written to a module give the same ydots

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
    def lib(self, reaclib_library):
        nuc = ["p", "n15", "o16"]
        lib = reaclib_library.linking_nuclei(nuc, with_reverse=False)
        r = IliadisO16pgF17()
        lib.add_rate(r)
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
        rho = 2.1e6
        T = 3.95e8

        if not _skip_build():
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                nc.evaluate(rho=rho, T=T)

        return nc

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_ydots(self, eval_cond1):

        # compare the AMReX and python module nets to the python
        # inline version

        for other in [eval_cond1.ydots_amrex, eval_cond1.ydots_py_module]:
            for nuc in eval_cond1.ydots_py_inline:
                assert other[nuc] == approx(eval_cond1.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_rates(self, eval_cond1):

        # compare the AMReX and python module nets to the python
        # inline version

        for other in [eval_cond1.rates_amrex, eval_cond1.rates_py_module]:
            for nuc in eval_cond1.rates_py_inline:
                assert other[nuc] == approx(eval_cond1.rates_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_ydots2(self, eval_cond2):

        # compare the AMReX and python module nets to the python
        # inline version

        for other in [eval_cond2.ydots_amrex, eval_cond2.ydots_py_module]:
            for nuc in eval_cond2.ydots_py_inline:
                assert other[nuc] == approx(eval_cond2.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)

    @pytest.mark.skipif(_skip_build(),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare_rates2(self, eval_cond2):

        # compare the AMReX and python module nets to the python
        # inline version

        for other in [eval_cond2.rates_amrex, eval_cond2.rates_py_module]:
            for nuc in eval_cond2.rates_py_inline:
                assert other[nuc] == approx(eval_cond2.rates_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)
