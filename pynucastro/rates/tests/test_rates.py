# unit tests for rates
import math

import pytest
from pytest import approx

from pynucastro import rates
from pynucastro.nucdata import Nucleus


class TestTfactors:
    @pytest.fixture(scope="class")
    def tf(self):
        return rates.Tfactors(2.e9)

    def test_tfactors(self, tf):
        assert tf.T9 == approx(2.0)
        assert tf.T9i == approx(0.5)
        assert tf.T913i == approx(0.5**(1./3.))
        assert tf.T913 == approx(2.0**(1./3.))
        assert tf.T953 == approx(2.0**(5./3.))
        assert tf.lnT9 == approx(math.log(2.0))


class TestRate:

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

        # chapter-1
        self.rate1 = rates.load_rate("o15--n15-wc12")

        # chapter-2
        self.rate2 = rates.load_rate("t-gn-d-nk06")

        # chapter-3
        self.rate3 = rates.load_rate("he6-gnn-he4-cf88")

        # chapter-4
        self.rate4 = rates.load_rate("c12-ag-o16-nac2")

        # chapter-5
        self.rate5 = rates.load_rate("n15-pa-c12-nacr")

        # chapter-6
        self.rate6 = rates.load_rate("he3-he3pp-he4-nacr")

        # chapter-7
        self.rate7 = rates.load_rate("li7-tnna-he4-mafo")

        # chapter-8
        self.rate8 = rates.load_rate("he4-aag-c12-fy05")

        # chapter-8, historical format (same rate as chapter-9)
        self.rate8_hist = rates.load_rate("he4-pphe3-he3-nacr-historical")

        # chapter-9
        self.rate9 = rates.load_rate("he4-pphe3-he3-nacr")

        # chapter-10
        self.rate10 = rates.load_rate("he4-npahe3-li7-mafo")

        # chapter-11
        self.rate11 = rates.load_rate("b17-nnn-c14-wc12")

        self.n = Nucleus("n")

        self.p = Nucleus("p")
        self.h1 = Nucleus("H1")
        self.d = Nucleus("d")
        self.h3 = Nucleus("H3")

        self.he3 = Nucleus("He3")
        self.he4 = Nucleus("He4")
        self.he6 = Nucleus("He6")

        self.li7 = Nucleus("Li7")

        self.b17 = Nucleus("B17")

        self.c12 = Nucleus("C12")
        self.c14 = Nucleus("C14")

        self.n15 = Nucleus("N15")

        self.o15 = Nucleus("O15")
        self.o16 = Nucleus("O16")

        self.ni56 = Nucleus("Ni56")
        self.u238 = Nucleus("U238")
        self.he4_also = Nucleus("he4")

    def teardown_method(self):
        """ this is run after each test """
        pass

    def test_reactants(self):

        # o15--n15-wc12
        assert self.rate1.reactants[0] == self.o15
        assert len(self.rate1.reactants) == 1

        # t-gn-d-nk06
        assert self.rate2.reactants[0] == self.h3
        assert len(self.rate2.reactants) == 1

        # he6-gnn-he4-cf88
        assert self.rate3.reactants[0] == self.he6
        assert len(self.rate3.reactants) == 1

        # c12-ag-o16-nac2
        assert self.rate4.reactants[0] == self.he4
        assert self.rate4.reactants[1] == self.c12
        assert len(self.rate4.reactants) == 2

        # n15-pa-c12-nacr
        assert self.rate5.reactants[0] == self.h1
        assert self.rate5.reactants[1] == self.n15
        assert len(self.rate5.reactants) == 2

        # he3-he3pp-he4-nacr
        assert self.rate6.reactants[0] == self.he3
        assert self.rate6.reactants[1] == self.he3
        assert len(self.rate6.reactants) == 2

        # li7-tnna-he4-mafo
        assert self.rate7.reactants[0] == self.h3
        assert self.rate7.reactants[1] == self.li7
        assert len(self.rate7.reactants) == 2

        # he4-aag-c12-fy05
        assert self.rate8.reactants[0] == self.he4
        assert self.rate8.reactants[1] == self.he4
        assert self.rate8.reactants[2] == self.he4
        assert len(self.rate8.reactants) == 3

        # he4-pphe3-he3-nacr-historical
        assert self.rate8_hist.reactants[0] == self.p
        assert self.rate8_hist.reactants[1] == self.h1
        assert self.rate8_hist.reactants[2] == self.he4
        assert len(self.rate8_hist.reactants) == 3

        # he4-pphe3-he3-nacr
        assert self.rate9.reactants[0] == self.p
        assert self.rate9.reactants[1] == self.h1
        assert self.rate9.reactants[2] == self.he4
        assert len(self.rate9.reactants) == 3

        # he4-npahe3-li7-mafo
        assert self.rate10.reactants[0] == self.n
        assert self.rate10.reactants[1] == self.h1
        assert self.rate10.reactants[2] == self.he4
        assert self.rate10.reactants[3] == self.he4
        assert len(self.rate10.reactants) == 4

        # b17-nnn-c14-wc12
        assert self.rate11.reactants[0] == self.b17
        assert len(self.rate11.reactants) == 1

    def test_products(self):
        assert self.rate4.products[0] == self.o16
        assert self.rate8.products[0] == self.c12
        assert len(self.rate8.products) == 1

        # he4-pphe3-he3-nacr-historical
        assert self.rate8_hist.products[0] == self.he3
        assert self.rate8_hist.products[1] == self.he3
        assert len(self.rate8_hist.products) == 2

        # he4-pphe3-he3-nacr
        assert self.rate9.products[0] == self.he3
        assert self.rate9.products[1] == self.he3
        assert len(self.rate9.products) == 2

    def test_prefactor(self):
        assert self.rate4.prefactor == 1.0
        assert self.rate8.prefactor == approx(0.16666666)

    def test_rate_exponent(self):
        assert self.rate8.get_rate_exponent(1.e8) == approx(40.9106396)

    def test_eval(self):
        assert self.rate8.eval(1.e8) == approx(2.0403192412842946e-24)

    def test_comparison(self):
        assert self.rate1 > self.rate2
        assert self.rate1 > self.rate4
        assert self.rate8 > self.rate9

    def test_weak(self):
        assert self.rate1.weak
        assert not self.rate2.weak

    def test_screen(self):
        assert not self.rate1.ion_screen
        assert self.rate4.ion_screen == [Nucleus("he4"), Nucleus("c12")]
        assert self.rate8.ion_screen == 3*[Nucleus("he4")]

    def test_heaviest_lightest(self):
        assert self.rate4.heaviest() == Nucleus("o16")
        assert self.rate4.lightest() == Nucleus("he4")
        assert self.rate2.lightest() == Nucleus("n")
        assert self.rate2.heaviest() == Nucleus("t")


class TestDerivedRate:

    def a_a_ag_c12(self, reaclib_library):
        """
        Here we test the inverse rate, computed by the use of detailed balance
        of a:

        A + B -> C

        reaction type.
        """

        a_a_ag_c12 = reaclib_library.get_rate('he4 + he4 + he4 --> c12 <fy05_reaclib__>')
        c12_ga_a_a_reaclib = reaclib_library.get_rate('c12 --> he4 + he4 + he4 <fy05_reaclib__reverse>')
        c12_ga_a_a_derived = rates.DerivedRate(rate=a_a_ag_c12, compute_Q=False, use_pf=False)

        assert c12_ga_a_a_reaclib.eval(T=2.0e9) == approx(c12_ga_a_a_derived.eval(T=2.0e9), rel=1.7e-5)

    def test_a_a_ag_c12_with_pf(self, reaclib_library):
        """
        This function test the correct rate value if we take in consideration the partition
        functions on the range 1.0e9 to 100.0e9
        """

        a_a_ag_c12 = reaclib_library.get_rate('he4 + he4 + he4 --> c12 <fy05_reaclib__>')
        c12_ga_a_a_derived = rates.DerivedRate(rate=a_a_ag_c12, compute_Q=False, use_pf=True)

        assert c12_ga_a_a_derived.eval(T=2.0e9) == approx(2.8953989705969484e-07)

    def test_a_a_ag_c12_with_Q(self, reaclib_library):
        """
        This function test the correct rate value if we take in consideration the
        exact values of atomic nuclear weight in order to compute the Q capture value
        of the reaction rate.
        """

        a_a_ag_c12 = reaclib_library.get_rate('he4 + he4 + he4 --> c12 <fy05_reaclib__>')
        c12_ga_a_a_derived = rates.DerivedRate(rate=a_a_ag_c12, compute_Q=True, use_pf=False)

        assert c12_ga_a_a_derived.eval(T=2.0e9) == approx(2.899433744446781e-07)


class TestWeakRates:
    @pytest.fixture(scope="class")
    def rate1(self):
        return rates.TabularRate("o18--f18-toki")

    @pytest.fixture(scope="class")
    def rate2(self):
        return rates.TabularRate("na22--ne22-toki")

    def test_reactants(self, rate1, rate2):

        assert len(rate1.reactants) == 1 and len(rate1.products) == 1
        assert rate1.products[0] == Nucleus("f18")
        assert rate1.reactants[0] == Nucleus("o18")
        assert rate1.eval(1.e10, 1.e7) == approx(3.990249e-11)

        assert len(rate2.reactants) == 1 and len(rate2.products) == 1
        assert rate2.products[0] == Nucleus("ne22")
        assert rate2.reactants[0] == Nucleus("na22")
        assert rate2.eval(1.e9, 1.e6) == approx(1.387075e-05)


class TestModify:
    @pytest.fixture(scope="function")
    def rate(self):
        return rates.load_rate("c12-c12n-mg23-cf88")

    def test_modify(self, rate):

        rate.modify_products("mg24")

        assert rate.Q == approx(13.93356)
        assert rate.products == [Nucleus("mg24")]
        assert rate.modified
