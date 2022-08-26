# unit tests for Binding Energy database taken from AME 2016.
from pynucastro.nucdata import BindingTable
import pytest


class TestAME:
    @pytest.fixture(scope="class")
    def bintable(self):
        return BindingTable()

    def test_get(self, bintable):
        nuc = bintable.get_nuclide(n=1, z=1)
        assert nuc.z == 1
        assert nuc.n == 1
        assert nuc.nucbind == 1.112283

        nuc = bintable.get_nuclide(n=5, z=6)
        assert nuc.z == 6
        assert nuc.n == 5
        assert nuc.nucbind == 6.676456

        nuc = bintable.get_nuclide(n=17, z=23)
        assert nuc.z == 23
        assert nuc.n == 17
        assert nuc.nucbind == 7.317

        nuc = bintable.get_nuclide(n=90, z=78)
        assert nuc.z == 78
        assert nuc.n == 90
        assert nuc.nucbind == 7.773605
