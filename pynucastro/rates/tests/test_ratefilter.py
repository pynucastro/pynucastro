# unit tests for rates

import pytest

import pynucastro as pyna


class TestRateFilter:
    @pytest.fixture(scope="class")
    def library(self):
        """ this is run before each test """

        files = ["c12-pg-n13-ls09",
                 "c13-pg-n14-nacr",
                 "n13-pg-o14-lg06",
                 "n14-pg-o15-im05",
                 "n15-pa-c12-nacr",
                 "o14--n14-wc12",
                 "o15--n15-wc12",
                 "o14-ap-f17-Ha96c",
                 "f17--o17-wc12",
                 "f17-pg-ne18-cb09",
                 "ne18--f18-wc12",
                 "f18--o18-wc12",
                 "n15-pg-o16-li10",
                 "o16-pg-f17-ia08",
                 "o17-pg-f18-il10",
                 "f18-pa-o15-il10"]

        rates = []
        for f in files:
            rates.append(pyna.load_rate(f))

        return pyna.Library(rates=rates)

    def test_inexact_filter(self, library):
        filter = pyna.RateFilter(reactants=['c12'], exact=False)
        newlib = library.filter(filter)

        rates = newlib.get_rates()

        assert len(rates) == 1
        assert str(rates[0]) == "C12 + p ⟶ N13 + 𝛾"

    def test_custom(self, library):

        # filter out all the rates with fluorine

        filter = pyna.RateFilter(filter_function=lambda r: len([q for q in r.reactants + r.products if q.Z == 9]))
        newlib = library.filter(filter)

        assert len(newlib.get_rates()) == 8

    def test_exact(self, library):

        filter = pyna.RateFilter(reactants=["n15", "p"])
        newlib = library.filter(filter)

        rates = newlib.get_rates()

        assert len(rates) == 2
        assert str(rates[0]) == "N15 + p ⟶ He4 + C12"
        assert str(rates[1]) == "N15 + p ⟶ O16 + 𝛾"
