"""Tests for the Rauscher library."""

from collections import Counter

import pytest

import pynucastro as pyna
from pynucastro.rates import RauscherRate


class TestRauscherLibrary:
    @pytest.fixture(scope="class")
    @classmethod
    def library(cls):
        return pyna.RauscherLibrary()

    def test_number_of_rates(self, library):
        assert library.num_rates == 9577

        counts = Counter(
            filename
            for rate in library.get_rates()
            for filename in rate.rfile_name
        )
        assert counts == {
            "fits_p.dat": 2755,
            "fits_a.dat": 1913,
            "fits_n.dat": 4909,
        }

    @pytest.mark.parametrize(
        "reaction, reactants, products",
        [
            ("ne20(p,g)na21", ["ne20", "p"], ["na21"]),
            ("ne20(a,g)mg24", ["ne20", "he4"], ["mg24"]),
            ("ne20(n,g)ne21", ["ne20", "n"], ["ne21"]),
            ("na23(p,a)ne20", ["na23", "p"], ["ne20", "he4"]),
        ],
    )
    def test_reaction_nuclei(self, library, reaction, reactants, products):
        rate = library.get_rate_by_name(reaction)

        assert isinstance(rate, RauscherRate)
        assert Counter(rate.reactants) == Counter(
            pyna.Nucleus(n) for n in reactants
        )
        assert Counter(rate.products) == Counter(
            pyna.Nucleus(n) for n in products
        )

    def test_ne20_evaluation(self, library):
        rate = library.get_rate_by_name("ne20(p,g)na21")

        assert rate.Q == pytest.approx(2.4320)
        assert rate.eval(2.0e9) == pytest.approx(
            1095.3513216021156, rel=1.0e-12
        )

    def test_derived_reverse(self, library):
        backward = library.derived_backward(use_pf=True)
        source = library.get_rate_by_name("ne20(p,g)na21")
        reverse = backward.get_rate_by_name("na21(g,p)ne20")

        assert isinstance(reverse, pyna.DerivedRate)
        assert reverse.source_rate is source
        assert reverse.reactants == source.products
        assert reverse.products == source.reactants
        assert reverse.derived_from_inverse
        assert all(
            isinstance(rate, pyna.DerivedRate) and rate.use_pf
            for rate in backward.get_rates()
        )
        assert reverse.eval(2.0e9) == pytest.approx(
            8.848159827195546e6, rel=1.0e-10
        )
