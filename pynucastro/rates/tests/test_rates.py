# unit tests for rates
import math

import pynucastro.rates as rates
from pytest import approx


class TestTfactors(object):
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
        self.tf = rates.Tfactors(2.e9)

    def teardown_method(self):
        """ this is run after each test """
        self.tf = None

    def test_tfactors(self):
        assert self.tf.T9 == approx(2.0)
        assert self.tf.T9i == approx(0.5)
        assert self.tf.T913i == approx(0.5**(1./3.))
        assert self.tf.T913 == approx(2.0**(1./3.))
        assert self.tf.T953 == approx(2.0**(5./3.))
        assert self.tf.lnT9 == approx(math.log(2.0))


class TestNucleus(object):
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

        self.p = rates.Nucleus("p")
        self.h1 = rates.Nucleus("H1")
        self.d = rates.Nucleus("d")
        self.he4 = rates.Nucleus("He4")
        self.c12 = rates.Nucleus("C12")
        self.o16 = rates.Nucleus("O16")
        self.ni56 = rates.Nucleus("Ni56")
        self.u238 = rates.Nucleus("U238")
        self.he4_also = rates.Nucleus("he4")

    def teardown_method(self):
        """ this is run after each test """
        pass

    def test_atomic_weights(self):
        assert self.d.A == 2
        assert self.he4.A == 4
        assert self.c12.A == 12
        assert self.ni56.A == 56
        assert self.u238.A == 238

    def test_atomic_numbers(self):
        assert self.d.Z == 1
        assert self.he4.Z == 2
        assert self.c12.Z == 6
        assert self.ni56.Z == 28
        assert self.u238.Z == 92

    def test_comparisons(self):
        assert self.p == self.h1
        assert self.d != self.he4
        assert self.d < self.he4
        assert self.ni56 > self.c12
        assert self.he4_also == self.he4


class TestRate(object):

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
        self.rate1 = rates.Rate("o15--n15-wc12")

        # chapter-2
        self.rate2 = rates.Rate("t-gn-d-nk06")

        # chapter-3
        self.rate3 = rates.Rate("he6-gnn-he4-cf88")

        # chapter-4
        self.rate4 = rates.Rate("c12-ag-o16-nac2")

        # chapter-5
        self.rate5 = rates.Rate("n15-pa-c12-nacr")

        # chapter-6
        self.rate6 = rates.Rate("he3-he3pp-he4-nacr")

        # chapter-7
        self.rate7 = rates.Rate("li7-tnna-he4-mafo")

        # chapter-8
        self.rate8 = rates.Rate("he4-aag-c12-fy05")

        # chapter-9
        self.rate9 = rates.Rate("he4-pphe3-he3-nacr")

        # chapter-10
        self.rate10 = rates.Rate("he4-npahe3-li7-mafo")

        # chapter-11
        self.rate11 = rates.Rate("b17-nnn-c14-wc12")

        self.n = rates.Nucleus("n")

        self.p = rates.Nucleus("p")
        self.h1 = rates.Nucleus("H1")
        self.d = rates.Nucleus("d")
        self.h3 = rates.Nucleus("H3")

        self.he3 = rates.Nucleus("He3")
        self.he4 = rates.Nucleus("He4")
        self.he6 = rates.Nucleus("He6")

        self.li7 = rates.Nucleus("Li7")

        self.b17 = rates.Nucleus("B17")

        self.c12 = rates.Nucleus("C12")
        self.c14 = rates.Nucleus("C14")

        self.n15 = rates.Nucleus("N15")

        self.o15 = rates.Nucleus("O15")
        self.o16 = rates.Nucleus("O16")

        self.ni56 = rates.Nucleus("Ni56")
        self.u238 = rates.Nucleus("U238")
        self.he4_also = rates.Nucleus("he4")

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

    def test_prefactor(self):
        assert self.rate4.prefactor == 1.0
        assert self.rate8.prefactor == approx(0.16666666)

    def test_rate_exponent(self):
        assert self.rate8.get_rate_exponent(1.e8) == approx(40.9106396)

    def test_eval(self):
        assert self.rate8.eval(1.e8) == approx(2.0403192412842946e-24)
