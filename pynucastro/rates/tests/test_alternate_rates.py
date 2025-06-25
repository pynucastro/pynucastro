# unit tests for rates

from pytest import approx

from pynucastro import Nucleus
from pynucastro.rates.alternate_rates import DeBoerC12agO16


class TestAlternateRates:
    def test_deboer(self):
        r = DeBoerC12agO16()

        assert r.reactants[0] == Nucleus("c12")
        assert r.reactants[1] == Nucleus("he4")

        assert r.products[0] == Nucleus("o16")

        assert r.eval(3.e8) == approx(3.667497144762534e-12)
