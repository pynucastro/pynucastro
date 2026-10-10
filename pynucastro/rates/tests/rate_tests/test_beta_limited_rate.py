# test the implementation of BetaLimitedRate

import pytest

from pynucastro.nucdata import Composition
from pynucastro.rates import BetaLimitedRate


class TestBetaLimitedRate:
    @pytest.fixture(scope="class")
    @classmethod
    def raw_rates(cls, reaclib_library):
        rn14pg = reaclib_library.get_rate_by_name("n14(p,g)o15")
        ro15betap = reaclib_library.get_rate_by_name("o15(,)n15")
        ro14betap = reaclib_library.get_rate_by_name("o14(,)n14")

        return {"n14(p,g)o15": rn14pg,
                "o15(,)n15": ro15betap,
                "o14(,)n14": ro14betap}

    @pytest.fixture(scope="class")
    @classmethod
    def beta_limited_rate(cls, raw_rates):

        rn14pg_betalim = BetaLimitedRate(raw_rates["n14(p,g)o15"],
                                         [raw_rates["o15(,)n15"], raw_rates["o14(,)n14"]],
                                         limiter_nucleus="p",
                                         description="N14(p,g)O15 limited by O14(β+) and O15(β+)")

        return rn14pg_betalim

    def test_rates(self, raw_rates, beta_limited_rate):

        b = beta_limited_rate
        assert b.underlying_rate == raw_rates["n14(p,g)o15"]
        assert raw_rates["o15(,)n15"] in b.beta_limiting_rates
        assert raw_rates["o14(,)n14"] in b.beta_limiting_rates

    def test_get_child_rates(self, beta_limited_rate):

        assert len(beta_limited_rate.get_child_rates()) == 3

    def test_limiting(self, raw_rates, beta_limited_rate):

        b = beta_limited_rate
        r0 = raw_rates["n14(p,g)o15"]

        # at low temperatures, our beta-limited rate should just give
        # us the underlying rate
        rho = 100
        T = 2.e7
        comp = Composition(["p", "he4", "c12", "n14", "o16"], init="solar")

        assert b.eval(T, rho=rho, comp=comp) == r0.eval(T)

        # at high-temperatures, we should be beta-limited, the rate
        # will be independent of T

        T1 = 2.e8
        T2 = 5.e8

        assert b.eval(T1, rho=rho, comp=comp) == b.eval(T2, rho=rho, comp=comp)
