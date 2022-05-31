# unit tests for rates
import math

import pynucastro as pyna
from pytest import approx


class TestTfactors:
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """
        pass

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        pass

    def setup_method(self):
        """ this is run before each test """
        reaclib_library = pyna.ReacLibLibrary()
        mylib = reaclib_library.linking_nuclei(["mg24", "al27", "si28", "he4", "p"])
        rc = pyna.RateCollection(libraries=[mylib])

        self.rp = rc.get_rate("he4_mg24__si28")
        self.rs = [rc.get_rate("he4_mg24__p_al27"), rc.get_rate("p_al27__si28")]

        self.rp_reverse = rc.get_rate("si28__he4_mg24")
        self.rs_reverse = [rc.get_rate("si28__p_al27"), rc.get_rate("p_al27__he4_mg24")]

        # approximate Mg24(a,g)Si28 together with Mg24(a,p)Al27(p,g)Si28

        self.ar = pyna.ApproximateRate(primary_rate=self.rp, secondary_rates=self.rs,
                          primary_reverse=self.rp_reverse, secondary_reverse=self.rs_reverse)

    def teardown_method(self):
        """ this is run after each test """

    def test_label(self):
        assert self.ar.fname == "mg24_he4__si28"

    def test_low_temp(self):
        # at low temperatures, the approximate (a,g) should be ~ (a,g) + (a,p)

        T = 5.e8
        assert self.ar.eval(T) == approx(self.rp.eval(T) + self.rs[0].eval(T), 1.e-6)

        T = 1.e9
        assert self.ar.eval(T) == approx(self.rp.eval(T) + self.rs[0].eval(T), 0.1)
