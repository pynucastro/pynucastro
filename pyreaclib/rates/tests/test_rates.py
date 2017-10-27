# unit tests for rates
import math
import os

import pyreaclib.rates as rates
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
        self.d = rates.Nucleus("d")
        self.p = rates.Nucleus("p")
        self.h1 = rates.Nucleus("H1")
        self.he4 = rates.Nucleus("He4")
        self.c12 = rates.Nucleus("C12")
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
        base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__))).split("pyreaclib")[0]
        cls.library_path = os.path.join(base_dir, "pyreaclib/pyreaclib/library/")

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        pass

    def setup_method(self):
        """ this is run before each test """
        self.rate1 = rates.Rate(os.path.join(self.library_path, "c12-ag-o16-nac2"))
        self.rate2 = rates.Rate(os.path.join(self.library_path, "he4-aag-c12-fy05"))

        self.he4 = rates.Nucleus("he4")
        self.c12 = rates.Nucleus("c12")
        self.o16 = rates.Nucleus("o16")

    def teardown_method(self):
        """ this is run after each test """
        pass

    def test_reactants(self):
        assert self.rate1.reactants[0] == self.he4
        assert self.rate1.reactants[1] == self.c12

        assert self.rate2.reactants[0] == self.he4
        assert self.rate2.reactants[1] == self.he4
        assert self.rate2.reactants[2] == self.he4

    def test_products(self):
        assert self.rate1.products[0] == self.o16
        assert self.rate2.products[0] == self.c12

    def test_prefactor(self):
        assert self.rate1.prefactor == 1.0
        assert self.rate2.prefactor == approx(0.16666666)

    def test_rate_exponent(self):
        assert self.rate2.get_rate_exponent(1.e8) == approx(40.9106396)

    def test_eval(self):
        assert self.rate2.eval(1.e8) == approx(2.0403192412842946e-24)

