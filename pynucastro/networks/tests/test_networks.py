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

        sum = 0.0
        for k in comp.X:
            sum += comp.X[k]

        assert sum == approx(1.0)
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


class TestRateCollection:
    @pytest.fixture(scope="class")
    def rc(self):
        files = ["c12-pg-n13-ls09",
                 "c13-pg-n14-nacr",
                 "n13--c13-wc12",
                 "n13-pg-o14-lg06",
                 "n14-pg-o15-im05",
                 "n15-pa-c12-nacr",
                 "o14--n14-wc12",
                 "o15--n15-wc12"]
        return networks.RateCollection(files)

    def test_nuclei(self, rc):
        nuc = rc.get_nuclei()
        assert nuc == [Nucleus("p"),
                       Nucleus("he4"),
                       Nucleus("c12"),
                       Nucleus("c13"),
                       Nucleus("n13"),
                       Nucleus("n14"),
                       Nucleus("n15"),
                       Nucleus("o14"),
                       Nucleus("o15")]

    def test_eval(self, rc):
        c = networks.Composition(rc.unique_nuclei)
        c.set_solar_like()

        rates = {"c12 + p --> n13 <ls09_reaclib__>": 4.3825344233265815e-05,
                 "c13 + p --> n14 <nacr_reaclib__>": 0.00012943869407433355,
                 "n13 --> c13 <wc12_reaclib_weak_>": 2.5475016632596765e-07,
                 "n13 + p --> o14 <lg06_reaclib__>": 4.8517620910445875e-06,
                 "n14 + p --> o15 <im05_reaclib__>": 9.8137074572314962e-07,
                 "n15 + p --> he4 + c12 <nacr_reaclib__>": 0.087518552257659241,
                 "o14 --> n14 <wc12_reaclib_weak_>": 2.0036691481625654e-06,
                 "o15 --> n15 <wc12_reaclib_weak_>": 1.0822012944765837e-06}

        rv = rc.evaluate_rates(1.e4, 1.e8, c)

        for r in rv:
            assert rv[r] == approx(rates[r.get_rate_id()])

    def test_overview(self, rc):

        ostr = """
p
  consumed by:
     C12 + p ⟶ N13 + 𝛾
     C13 + p ⟶ N14 + 𝛾
     N13 + p ⟶ O14 + 𝛾
     N14 + p ⟶ O15 + 𝛾
     N15 + p ⟶ He4 + C12
  produced by:

he4
  consumed by:
  produced by:
     N15 + p ⟶ He4 + C12

c12
  consumed by:
     C12 + p ⟶ N13 + 𝛾
  produced by:
     N15 + p ⟶ He4 + C12

c13
  consumed by:
     C13 + p ⟶ N14 + 𝛾
  produced by:
     N13 ⟶ C13 + e⁺ + 𝜈

n13
  consumed by:
     N13 ⟶ C13 + e⁺ + 𝜈
     N13 + p ⟶ O14 + 𝛾
  produced by:
     C12 + p ⟶ N13 + 𝛾

n14
  consumed by:
     N14 + p ⟶ O15 + 𝛾
  produced by:
     C13 + p ⟶ N14 + 𝛾
     O14 ⟶ N14 + e⁺ + 𝜈

n15
  consumed by:
     N15 + p ⟶ He4 + C12
  produced by:
     O15 ⟶ N15 + e⁺ + 𝜈

o14
  consumed by:
     O14 ⟶ N14 + e⁺ + 𝜈
  produced by:
     N13 + p ⟶ O14 + 𝛾

o15
  consumed by:
     O15 ⟶ N15 + e⁺ + 𝜈
  produced by:
     N14 + p ⟶ O15 + 𝛾
"""
        assert rc.network_overview().replace(" ", "").strip() == ostr.replace(" ", "").strip()
