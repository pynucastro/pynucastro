# unit tests for AME database

import os

import pyreaclib.amemass as amemass

class TestAME(object):
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
        self.ame = amemass.ame2012.AME2012()

    def teardown_method(self):
        """ this is run after each test """
        self.ame = None

    def test_get(self):
        nuc = self.ame.get_nuclide(n=0, z=1, a=1)
        assert nuc.z == 1
        assert nuc.a == 1
        assert nuc.n == 0
        assert nuc.mexcess == 7.28897059

