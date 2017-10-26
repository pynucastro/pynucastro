# unit tests for rates
import math

import pyreaclib.networks as networks
import pyreaclib.rates as rates

from pytest import approx

class TestComposition(object):
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
        self.nuclei = [rates.Nucleus("h1"),
                       rates.Nucleus("he4"),
                       rates.Nucleus("c12"),
                       rates.Nucleus("o16"),
                       rates.Nucleus("n14"),
                       rates.Nucleus("ca40")]

        self.comp = networks.Composition(self.nuclei)

    def teardown_method(self):
        """ this is run after each test """
        self.tf = None

    def test_solar(self):
        self.comp.set_solar_like()

        sum = 0.0
        for k in self.comp.X:
            sum += self.comp.X[k]

        assert sum == approx(1.0)
        assert self.comp.X[rates.Nucleus("h1")] == approx(0.7)
