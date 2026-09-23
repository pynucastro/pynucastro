# test the comparison of the BranchedRate (and ModifiedRate) for
# the reduced CNO network across network types

import sys
import warnings
from pathlib import Path

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare
from pynucastro.nucdata.nucleus import Nucleus
from pynucastro.rates.branched_rate import BranchedRate
from pynucastro.rates.library import Library
from pynucastro.rates.modified_rate import ModifiedRate


def _skip_build():
    return sys.platform == "darwin" or sys.platform.startswith("win")


class TestNetworkCompare:

    # pylint: disable=duplicate-code
    @pytest.fixture(scope="class")
    @classmethod
    def lib(cls, reaclib_library):

        rc12pg = reaclib_library.get_rate_by_name("c12(p,g)n13")
        rn14pg = reaclib_library.get_rate_by_name("n14(p,g)o15")
        rn15pa = reaclib_library.get_rate_by_name("n15(p,a)c12")
        rn15pg = reaclib_library.get_rate_by_name("n15(p,g)o16")
        ro16pg = reaclib_library.get_rate_by_name("o16(p,g)f17")

        # C12 + 2p -> N14
        rc12_2p_n14 = ModifiedRate(rc12pg,
                                   new_products=[Nucleus("n14")],
                                   stoichiometry={Nucleus("p"): 2})

        # N14 + 2p branches
        rn14_2p_c12 = BranchedRate(rn14pg,
                                   primary_branch=rn15pa,
                                   other_branch=rn15pg,
                                   stoichiometry={Nucleus("p"): 2},
                                   description="N14(p,g)O15(,e+nu)N15(p,a)C12")

        rn14_2p_o16 = BranchedRate(rn14pg,
                                   primary_branch=rn15pg,
                                   other_branch=rn15pa,
                                   stoichiometry={Nucleus("p"): 2},
                                   description="N14(p,g)O15(,e+nu)N15(p,g)O16")

        # O16 + 2p -> N14 + He4
        ro16_2p_n14_a = ModifiedRate(ro16pg,
                                     new_products=[Nucleus("n14"),
                                                   Nucleus("he4")],
                                     stoichiometry={Nucleus("p"): 2})

        lib = Library(rates=[rc12_2p_n14, rn14_2p_c12, rn14_2p_o16, ro16_2p_n14_a])
        return lib

    @pytest.fixture(scope="class")
    @classmethod
    def nc(cls, lib):
        cxx_test_path = Path("_test_compare_branched_cxx/")
        amrex_test_path = Path("_test_compare_branched_amrex/")

        nc = NetworkCompare(lib,
                            include_amrex=True,
                            include_simple_cxx=True,
                            python_module_name="branched_compare.py",
                            amrex_test_path=amrex_test_path,
                            cxx_test_path=cxx_test_path)
        return nc

    @pytest.fixture(scope="class")
    @classmethod
    def eval_cond1(cls, nc):
        # thermodynamic conditions
        rho = 200
        T = 2.e7

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
    def test_compare_energy(self, eval_cond1):

        # we use a relaxed tolerance here because of differences
        # in constants in simple C++ nets (N_A)
        for other in [eval_cond1.enuc_cxx, eval_cond1.enuc_amrex, eval_cond1.enuc_py_module]:
            assert other == approx(eval_cond1.enuc_py_inline,
                                   rel=1.e-7, abs=1.e-30)
        for other in [eval_cond1.enu_weak_cxx, eval_cond1.enu_weak_amrex, eval_cond1.enu_weak_py_module]:
            assert other == approx(eval_cond1.enu_weak_py_inline,
                                   rel=1.e-7, abs=1.e-30)

    # pylint: enable=duplicate-code
