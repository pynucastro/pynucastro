# unit tests for Binding Energy database taken from AME 2016.
from pytest import approx

from pynucastro import Nucleus


class TestAME:
    def test_get(self):

        d = Nucleus("d")
        assert d.nucbind == approx(1.112289947499903)

        c11 = Nucleus("c11")
        assert c11.nucbind == approx(6.676463512727423)

        v40 = Nucleus("v40")
        assert v40.nucbind == approx(7.309726389249954)

        pt168 = Nucleus("pt168")
        assert pt168.nucbind == approx(7.773627588214454)

        c40 = Nucleus("c40")
        assert c40.nucbind is None
