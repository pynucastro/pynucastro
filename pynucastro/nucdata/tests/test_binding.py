# unit tests for Binding Energy database taken from AME 2016.
from pynucastro.nucdata import BindingTable


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
        nuc = self.bintable.get_nuclide(n=1, z=1)
        assert nuc.z == 1
        assert nuc.n == 1
        assert nuc.nucbind == 1.112283

        nuc = self.bintable.get_nuclide(n=5, z=6)
        assert nuc.z == 6
        assert nuc.n == 5
        assert nuc.nucbind == 6.676456

        nuc = self.bintable.get_nuclide(n=17, z=23)
        assert nuc.z == 23
        assert nuc.n == 17
        assert nuc.nucbind == 7.317

        nuc = self.bintable.get_nuclide(n=90, z=78)
        assert nuc.z == 78
        assert nuc.n == 90
        assert nuc.nucbind == 7.773605
