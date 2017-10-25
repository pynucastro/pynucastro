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
