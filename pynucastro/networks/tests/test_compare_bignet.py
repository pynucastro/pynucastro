# This creates the Microphysics ase-iron network and compares python
# inline, module, and AMReX ydots.  No screening is used.  At the moment,
# SimpleCxxNetwork is not included.

import warnings

import sys
from pathlib import Path

import pytest
from pytest import approx

from pynucastro.networks.network_compare import NetworkCompare
from pynucastro.rates.modified_rate import ModifiedRate
from pynucastro.rates.derived_rate import DerivedRate
from pynucastro.rates.library import TabularLibrary
from pynucastro.rates.library import Library

class TestNetworkCompare:

    @pytest.fixture(scope="class")
    def lib(self, reaclib_library):
        all_reactants = ["p",
                         "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                         "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                         "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                         "n13", "na23"]
        lib = reaclib_library.linking_nuclei(all_reactants, print_warning=False)

        # in this list, we have the reactants, the actual reactants,
        # and modified products that we will use instead
        other_rates = [("c12(c12,n)mg23", "mg24"),
                       ("o16(o16,n)s31", "s32"),
                       ("o16(c12,n)si27", "si28")]

        for r, mp in other_rates:
            _r = reaclib_library.get_rate_by_name(r)
            forward_rate = ModifiedRate(_r, new_products=[mp])
            derived_rate = DerivedRate(forward_rate, use_pf=True)
            lib += Library(rates=[forward_rate, derived_rate])

        # C12+Ne20 and reverse
        # (a,g) links between Na23 and Al27
        # (a,g) links between Al27 and P31

        rates_to_remove = ["p31(p,c12)ne20",
                           "si28(a,c12)ne20",
                           "ne20(c12,p)p31",
                           "ne20(c12,a)si28",
                           "na23(a,g)al27",
                           "al27(g,a)na23",
                           "al27(a,g)p31",
                           "p31(g,a)al27"]

        for r in rates_to_remove:
            _r = lib.get_rate_by_name(r)
            lib.remove_rate(_r)

        # iron group
        iron_peak = ["n", "p", "he4",
                     "mn51",
                     "fe52", "fe53", "fe54", "fe55", "fe56",
                     "co55", "co56", "co57",
                     "ni56", "ni57", "ni58"]
        lib += reaclib_library.linking_nuclei(iron_peak, print_warning=False)
        weak_lib = TabularLibrary(ordering=["ffn", "langanke", "oda"])
        iron_weak_lib = weak_lib.linking_nuclei(set(iron_peak + all_reactants),
                                                print_warning=False)
        lib += iron_weak_lib

        rates_to_derive = lib.backward().get_rates()

        # now for each of those derived rates, look to see if the pair exists

        for r in rates_to_derive:
            fr = lib.get_rate_by_nuclei(r.products, r.reactants)
            if fr:
                lib.remove_rate(r)
                d = DerivedRate(fr, use_pf=True)
                lib.add_rate(d)

        lib.eliminate_duplicates(rate_type_preference="tabular")

        return lib

    def check_lib(self, lib):
        assert lib.num_rates == 156

    @pytest.fixture(scope="class")
    def nc(self, lib):
        cxx_test_path = Path("_test_compare_bignet_cxx/")
        amrex_test_path = Path("_test_compare_bignet_amrex/")

        nc = NetworkCompare(lib,
                            include_amrex=True,
                            include_simple_cxx=False,
                            python_module_name="basic_cxx_py_bignet_compare.py",
                            amrex_test_path=amrex_test_path,
                            cxx_test_path=cxx_test_path)
        return nc

    @pytest.mark.skipif(sys.platform == "darwin" or sys.platform.startswith("win"),
                        reason="We do not build C++ on Mac or Windows")
    def test_compare(self, nc):

        # thermodynamic conditions
        rho = 2.e8
        T = 1.e9

        # filter the partition function warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            nc.evaluate(rho=rho, T=T)

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
        rho = 2.e7
        T = 4.e9

        # filter the partition function warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            nc.evaluate(rho=rho, T=T)

        # compare the AMReX and python module nets to the python
        # inline version

        for other in [nc.ydots_amrex, nc.ydots_py_module]:
            for nuc in nc.ydots_py_inline:
                assert other[nuc] == approx(nc.ydots_py_inline[nuc],
                                            rel=1.e-11, abs=1.e-30)
