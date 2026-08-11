# test the implementation of BranchedRate

import pytest
from pytest import approx

from pynucastro.nucdata import Nucleus
from pynucastro.rates import BranchedRate


class TestBranchedRate:
    @pytest.fixture(scope="class")
    @classmethod
    def raw_rates(cls, reaclib_library):
        rn14pg = reaclib_library.get_rate_by_name("n14(p,g)o15")
        rn15pa = reaclib_library.get_rate_by_name("n15(p,a)c12")
        rn15pg = reaclib_library.get_rate_by_name("n15(p,g)o16")

        return {"n14(p,g)o15": rn14pg,
                "n15(p,a)c12": rn15pa,
                "n15(p,g)o16": rn15pg}

    @pytest.fixture(scope="class")
    @classmethod
    def branched_rates(cls, raw_rates):
        rn14_2p_c12 = BranchedRate(raw_rates["n14(p,g)o15"],
                                   primary_branch=raw_rates["n15(p,a)c12"],
                                   other_branch=raw_rates["n15(p,g)o16"],
                                   stoichiometry={Nucleus("p"): 2})

        rn14_2p_o16 = BranchedRate(raw_rates["n14(p,g)o15"],
                                   primary_branch=raw_rates["n15(p,g)o16"],
                                   other_branch=raw_rates["n15(p,a)c12"],
                                   stoichiometry={Nucleus("p"): 2})

        return {"n14(pp,ae+nu)c12": rn14_2p_c12,
                "n14(pp,e+nu)o16": rn14_2p_o16}

    def test_rates(self, raw_rates, branched_rates):

        b1 = branched_rates["n14(pp,ae+nu)c12"]
        assert b1.underlying_rate == raw_rates["n14(p,g)o15"]
        assert b1.primary_branch == raw_rates["n15(p,a)c12"]
        assert b1.other_branch == raw_rates["n15(p,g)o16"]

        b2 = branched_rates["n14(pp,e+nu)o16"]
        assert b2.primary_branch == b1.other_branch
        assert b2.other_branch == b1.primary_branch

    def test_products(self, branched_rates):

        b1 = branched_rates["n14(pp,ae+nu)c12"]
        assert len(b1.products) == 2
        assert Nucleus("he4") in b1.products
        assert Nucleus("c12") in b1.products

        b2 = branched_rates["n14(pp,e+nu)o16"]
        assert len(b2.products) == 1
        assert Nucleus("o16") in b2.products

    def test_branching_ratio(self, raw_rates, branched_rates):

        T = 2.e7
        f = raw_rates["n15(p,a)c12"].eval(T) / (raw_rates["n15(p,a)c12"].eval(T) +
                                                raw_rates["n15(p,g)o16"].eval(T))

        r0 = raw_rates["n14(p,g)o15"].eval(T)

        b1 = branched_rates["n14(pp,ae+nu)c12"]
        assert b1.eval(T) == approx(r0 * f, rel=1.e-10, abs=1.e-99)

        b2 = branched_rates["n14(pp,e+nu)o16"]
        assert b2.eval(T) == approx(r0 * (1.0 - f), rel=1.e-10, abs=1.e-99)

        assert b1.eval(T) + b2.eval(T) == approx(raw_rates["n14(p,g)o15"].eval(T))

    def test_ydot_form(self, branched_rates):

        # we mainly want to make sure that the Y[jp] term is linear
        # and not squared (from the stoichiometry)

        b1 = branched_rates["n14(pp,ae+nu)c12"]
        assert b1.ydot_string_py() == "rho*Y[jp]*Y[jn14]*rate_eval.p_N14_to_He4_C12_branched"

        b2 = branched_rates["n14(pp,e+nu)o16"]
        assert b2.ydot_string_py() == "rho*Y[jp]*Y[jn14]*rate_eval.p_N14_to_O16_branched"
