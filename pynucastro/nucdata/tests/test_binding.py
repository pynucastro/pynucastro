# unit tests for Binding Energy database taken from AME 2016.
from pytest import approx

from pynucastro.nucdata import BindingTable


class TestAME:
    def test_get(self):
        bintable = BindingTable()
        assert bintable.get_binding_energy(n=1, z=1) == approx(1.112283)
        assert bintable.get_binding_energy(n=5, z=6) == approx(6.676456)
        assert bintable.get_binding_energy(n=17, z=23) == approx(7.317)
        assert bintable.get_binding_energy(n=90, z=78) == approx(7.773605)
