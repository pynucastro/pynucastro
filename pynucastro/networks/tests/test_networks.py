# unit tests for networks
import pytest
from pytest import approx

from pynucastro import networks
from pynucastro.nucdata import Nucleus


class TestRateCollection:
    @pytest.fixture(scope="class")
    def rc(self, reaclib_library):
        rate_names = ["c12(p,g)n13",
                      "c13(p,g)n14",
                      "n13(,)c13",
                      "n13(p,g)o14",
                      "n14(p,g)o15",
                      "n15(p,a)c12",
                      "o14(,)n14",
                      "o15(,)n15"]
        rates = reaclib_library.get_rate_by_name(rate_names)
        return networks.RateCollection(rates=rates)

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
