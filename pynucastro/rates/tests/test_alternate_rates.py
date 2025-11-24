# unit tests for rates

from pytest import approx

from pynucastro import Nucleus
from pynucastro.rates.alternate_rates import DeBoerC12agO16, IliadisO16pgF17


class TestAlternateRates:
    def test_deboer(self):
        r = DeBoerC12agO16()

        assert r.reactants[0] == Nucleus("c12")
        assert r.reactants[1] == Nucleus("he4")

        assert r.products[0] == Nucleus("o16")

        assert r.eval(3.e8) == approx(3.667497144762534e-12)


class TestIliadisO16:
    def test_o16pg(self, reaclib_library):
        # the new rate should agree with the reaclib rate
        # to within 5%

        o16pg_rl = reaclib_library.get_rate_by_name("o16(p,g)f17")
        new_rate = IliadisO16pgF17()

        for T in [4.2e7, 5.81e7, 1.02e8, 4.55e8, 1.23e9]:
            assert new_rate.eval(T) == approx(o16pg_rl.eval(T),
                                              rel=0.055, abs=1.e-50)
