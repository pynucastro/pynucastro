# unit tests for rates

import pytest
from pytest import approx

import pynucastro as pyna


class TestTfactors:
    @pytest.fixture(scope="class")
    def rc(self, reaclib_library):
        mylib = reaclib_library.linking_nuclei(["mg24", "al27", "si28", "he4", "p"])
        return pyna.RateCollection(libraries=[mylib])

    @pytest.fixture(scope="class")
    def rp(self, rc):
        return rc.get_rate("he4_mg24__si28")

    @pytest.fixture(scope="class")
    def rs(self, rc):
        return [rc.get_rate("he4_mg24__p_al27"), rc.get_rate("p_al27__si28")]

    @pytest.fixture(scope="class")
    def ar(self, rc, rp, rs):
        rp_reverse = rc.get_rate("si28__he4_mg24")
        rs_reverse = [rc.get_rate("si28__p_al27"), rc.get_rate("p_al27__he4_mg24")]

        # approximate Mg24(a,g)Si28 together with Mg24(a,p)Al27(p,g)Si28

        return pyna.ApproximateRate(primary_rate=rp, secondary_rates=rs,
                                    primary_reverse=rp_reverse, secondary_reverse=rs_reverse)

    def test_label(self, ar):
        assert ar.fname == "mg24_he4__si28__approx"

    def test_low_temp(self, ar, rp, rs):
        # at low temperatures, the approximate (a,g) should be ~ (a,g) + (a,p)

        T = 5.e8
        assert ar.eval(T) == approx(rp.eval(T) + rs[0].eval(T), 1.e-6)

        T = 1.e9
        assert ar.eval(T) == approx(rp.eval(T) + rs[0].eval(T), 0.1)

    def test_child_rates(self, ar):

        cr = ar.get_child_rates()
        assert cr[0].fname == "he4_mg24__si28"
        assert cr[1].fname == "he4_mg24__p_al27"
        assert cr[2].fname == "p_al27__si28"
        assert cr[3].fname == "si28__he4_mg24"
        assert cr[4].fname == "si28__p_al27"
        assert cr[5].fname == "p_al27__he4_mg24"

        assert len(cr) == 6
