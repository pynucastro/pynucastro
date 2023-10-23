# unit tests for networks
import pytest
from pytest import approx

import pynucastro as pyna
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


class TestCompositionVars:
    @pytest.fixture(scope="class")
    def nuclei(self):
        return [Nucleus("h1"),
                Nucleus("he4"),
                Nucleus("c12"),
                Nucleus("fe56")]

    @pytest.fixture(scope="class")
    def comp(self, nuclei):
        c = networks.Composition(nuclei)
        c.set_equal()
        return c

    def test_abar(self, comp):

        # 1/Abar = sum_k X_k / A_k
        Abar_inv = sum(comp.X[n] / n.A for n in comp.X)
        assert 1.0 / Abar_inv == approx(comp.eval_abar())

    def test_zbar(self, comp):

        # Zbar = Abar sum_k Z_k X_k / A_k
        Abar_inv = sum(comp.X[n] / n.A for n in comp.X)
        Zbar = 1.0 / Abar_inv * sum(comp.X[n] * n.Z / n.A for n in comp.X)
        assert Zbar == approx(comp.eval_zbar())

    def test_ye(self, comp):

        # Ye = sum_k Z_k X_k / A_k
        Ye = sum(comp.X[n] * n.Z / n.A for n in comp.X)
        assert Ye == approx(comp.eval_ye())


class TestCompBinning:
    """this example has several cases where there are multiple matches
    with the same A, so we exercise the secondary check on Z"""

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
                Nucleus("cr56"),
                Nucleus("co56"),
                Nucleus("zn60")]

    @pytest.fixture(scope="class")
    def comp(self, nuclei):
        c = networks.Composition(nuclei)
        c.set_equal()
        return c

    def test_bin_as(self, nuclei, comp, capsys):

        new_nuclei = [Nucleus("he4"),
                      Nucleus("c12"),
                      Nucleus("ti44"),
                      Nucleus("ni56"),
                      Nucleus("fe56")]

        c_new = comp.bin_as(new_nuclei, verbose=True)

        captured = capsys.readouterr()

        output = """storing p as He4
storing C12 as C12
storing C13 as C12
storing He4 as He4
storing O14 as C12
storing N14 as C12
storing Cr48 as Ti44
storing Fe52 as Ti44
storing Fe56 as Fe56
storing Cr56 as Fe56
storing Co56 as Fe56
storing Zn60 as Ni56
"""

        assert captured.out == output

        assert c_new.get_sum_X() == approx(comp.get_sum_X())

        orig_X = 1.0 / len(nuclei)

        # we should have placed p and He4 into He4
        assert c_new.X[Nucleus("he4")] == approx(2.0 * orig_X)

        # we should have placed C12, C13, N14, O14 into C12
        assert c_new.X[Nucleus("c12")] == approx(4.0 * orig_X)

        # we should have placed Cr48 and Fe52 into Ti44
        assert c_new.X[Nucleus("ti44")] == approx(2.0 * orig_X)

        # we should have placed Cr56, Co56, and Fe56 into Fe56
        assert c_new.X[Nucleus("fe56")] == approx(3.0 * orig_X)

        # we should have placed Zn60 into Ni56
        assert c_new.X[Nucleus("ni56")] == approx(orig_X)


class TestCompBinning2:
    """a more extreme example -- we'll be into a composition where
    none of the original nuclei are present"""
    @pytest.fixture(scope="class")
    def nuclei(self):
        return [Nucleus("d"),
                Nucleus("he3"),
                Nucleus("he4"),
                Nucleus("he5"),
                Nucleus("c12"),
                Nucleus("o14"),
                Nucleus("o15"),
                Nucleus("o16"),
                Nucleus("o17"),
                Nucleus("o18")]

    @pytest.fixture(scope="class")
    def comp(self, nuclei):
        c = networks.Composition(nuclei)
        c.set_equal()
        return c

    def test_bin_as(self, nuclei, comp):

        new_nuclei = [Nucleus("p"), Nucleus("n14"), Nucleus("f18")]
        c_new = comp.bin_as(new_nuclei)

        assert c_new.get_sum_X() == approx(comp.get_sum_X())

        orig_X = 1.0 / len(nuclei)

        # we should have placed d, He3. He4. He5, and C12 into p
        assert c_new.X[Nucleus("p")] == approx(5.0 * orig_X)

        # we should have placed O14, O15, O16, and O17 into N14
        assert c_new.X[Nucleus("n14")] == approx(4.0 * orig_X)

        # we should have placed O18 into F18
        assert c_new.X[Nucleus("f18")] == approx(orig_X)


class TestCompBinning3:
    """an example where we exclude Ni56 from the binning."""
    @pytest.fixture(scope="class")
    def nuclei(self):
        nuc_list = pyna.get_nuclei_in_range(26, 26, 52, 58)
        nuc_list += pyna.get_nuclei_in_range(28, 28, 56, 58)
        return nuc_list

    @pytest.fixture(scope="class")
    def comp(self, nuclei):
        c = networks.Composition(nuclei)
        c.set_equal()
        return c

    def test_bin_as(self, nuclei, comp):
        """first bin without exclusion"""

        new_nuclei = [Nucleus("fe52"), Nucleus("fe54"), Nucleus("ni56")]

        c_new = comp.bin_as(new_nuclei)

        assert c_new.get_sum_X() == approx(comp.get_sum_X())

        orig_X = 1.0 / len(nuclei)

        # we should have placed fe52, fe53 into fe52
        assert c_new.X[Nucleus("fe52")] == approx(2.0 * orig_X)

        # we should have placed fe54, fe55 into fe54
        assert c_new.X[Nucleus("fe54")] == approx(2.0 * orig_X)

        # everything else should be ni56
        assert c_new.X[Nucleus("ni56")] == approx(6.0 * orig_X)

    def test_bin_as_exclude(self, nuclei, comp):
        """exclude Ni56"""

        new_nuclei = [Nucleus("fe52"), Nucleus("fe54"), Nucleus("ni56")]

        c_new = comp.bin_as(new_nuclei, exclude=[Nucleus("ni56")])

        assert c_new.get_sum_X() == approx(comp.get_sum_X())

        orig_X = 1.0 / len(nuclei)

        # we should have placed fe52, fe53 into fe52
        assert c_new.X[Nucleus("fe52")] == approx(2.0 * orig_X)

        # we should have placed fe54, fe55, fe56, fe57, fe58, ni57, ni58 into fe54
        assert c_new.X[Nucleus("fe54")] == approx(7.0 * orig_X)

        # only ni56 should be ni56
        assert c_new.X[Nucleus("ni56")] == approx(orig_X)
