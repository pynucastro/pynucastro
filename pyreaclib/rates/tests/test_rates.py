# unit tests for rates
import math

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
