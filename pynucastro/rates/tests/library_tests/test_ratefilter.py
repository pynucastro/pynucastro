# unit tests for rates

import pytest

import pynucastro as pyna


class TestRateFilter:
    @pytest.fixture(scope="class")
    def library(self, reaclib_library):
        """ this is run before each test """

        rate_names = ["c12(p,g)n13",
                      "c13(p,g)n14",
                      "n13(p,g)o14",
                      "n14(p,g)o15",
                      "n15(p,a)c12",
                      "o14(,)n14",
                      "o15(,)n15",
                      "o14(a,p)f17",
                      "f17(,)o17",
                      "f17(p,g)ne18",
                      "ne18(,)f18",
                      "f18(,)o18",
                      "n15(p,g)o16",
                      "o16(p,g)f17",
                      "o17(p,g)f18",
                      "f18(p,a)o15"]
        rates = reaclib_library.get_rate_by_name(rate_names)

        return pyna.Library(rates=rates)

    def test_inexact_filter(self, library):
        filt = pyna.RateFilter(reactants=['c12'], exact=False)
        newlib = library.filter(filt)

        rates = newlib.get_rates()

        assert len(rates) == 1
        assert str(rates[0]) == "C12 + p ⟶ N13 + 𝛾"

    def test_custom(self, library):

        # filter out all the rates with fluorine

        filt = pyna.RateFilter(filter_function=lambda r: len([q for q in r.reactants + r.products if q.Z == 9]))
        newlib = library.filter(filt)

        assert len(newlib.get_rates()) == 8

    def test_exact(self, library):

        filt = pyna.RateFilter(reactants=["n15", "p"])
        newlib = library.filter(filt)

        rates = newlib.get_rates()

        assert len(rates) == 2
        assert str(rates[0]) == "N15 + p ⟶ He4 + C12"
        assert str(rates[1]) == "N15 + p ⟶ O16 + 𝛾"

    def test_endpoint(self, reaclib_library):

        filt = pyna.RateFilter(endpoint="zn60")
        newlib = reaclib_library.filter(filt)

        assert newlib.num_rates == 7761
        assert newlib.heaviest() == pyna.Nucleus("s60")
