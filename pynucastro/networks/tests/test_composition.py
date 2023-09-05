# unit tests for networks
import pytest
from pytest import approx

from pynucastro import networks
from pynucastro.nucdata import Nucleus


class TestComposition:
    @pytest.fixture(scope="class")
    def nuclei(self):
        return [Nucleus("h1"),
                Nucleus("he4"),
                Nucleus("c12"),
                Nucleus("o16"),
                Nucleus("n14"),
                Nucleus("ca40")]

    @pytest.fixture(scope="class")
    def comp(self, nuclei):
        return networks.Composition(nuclei)

    def test_solar(self, comp):
        comp.set_solar_like()

        xsum = sum(comp.X.values())

        assert xsum == approx(1.0)
        assert comp.X[Nucleus("h1")] == approx(0.7)

    def test_set_all(self, nuclei, comp):
        val = 1.0/len(nuclei)
        comp.set_all(1.0/len(nuclei))
        for n in nuclei:
            assert comp.X[n] == val

    def test_set_nuc(self, nuclei, comp):
        n = nuclei[0]
        comp.set_nuc(n.raw, 0.55)
        assert comp.X[n] == 0.55

    def test_get_molar(self, comp):
        comp.set_solar_like(Z=0.02)
        molar = comp.get_molar()
        assert molar[Nucleus("he4")] == approx((0.3-0.02)/4.0)

    def test_set_equal(self, nuclei, comp):
        comp.set_equal()
        assert comp.X[nuclei[0]] == approx(1.0 / len(nuclei))


class TestCompBinning:
    @pytest.fixture(scope="class")
    def nuclei(self):
        return [Nucleus("p"),
                Nucleus("c12"),
                Nucleus("c13"),
                Nucleus("he4"),
                Nucleus("o14"),
                Nucleus("n14"),
                Nucleus("cr48"),
                Nucleus("fe52"),
                Nucleus("fe56"),
                Nucleus("cr56")]

    @pytest.fixture(scope="class")
    def comp(self, nuclei):
        c = networks.Composition(nuclei)
        c.set_equal()
        return c

    def test_bin_as(self, nuclei, comp):

        new_nuclei = [Nucleus("he4"), Nucleus("c12"), Nucleus("ni56")]

        c_new = comp.bin_as(new_nuclei)

        assert c_new.get_sum_X() == approx(comp.get_sum_X())
