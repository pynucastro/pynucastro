# unit tests for Binding Energy database taken from AME 2016.
from pytest import approx

from pynucastro import Nucleus


class TestAME:
    def test_get(self):

        d = Nucleus("d")
        assert d.nucbind == approx(1.1122831344998758)

        c11 = Nucleus("c11")
        assert c11.nucbind == approx(6.676456080363675)

        v40 = Nucleus("v40")
        assert v40.nucbind == approx(7.309718554299979)

        pt168 = Nucleus("pt168")
        assert pt168.nucbind == approx(7.77362126185714)

        c40 = Nucleus("c40")
        assert c40.nucbind is None
