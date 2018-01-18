# unit tests for Binding Energy database taken from AME 2016.

import os

from pyreaclib.nucdata import BindingTable

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
        self.bintable = BindingTable()

    def teardown_method(self):
        """ this is run after each test """
        self.bintable = None

    def test_get(self):
        nuc = self.bintable.get_nuclide(n=0, z=1)
        assert nuc.z == 1
        assert nuc.n == 0
        assert nuc.mexcess == 7.28897061

