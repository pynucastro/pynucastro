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

        rates = {"C12 + p --> N13 <reaclib_ls09>": 4.3825344233265815e-05,
                 "C13 + p --> N14 <reaclib_nacr>": 0.00012943869407433355,
                 "N13 --> C13 <reaclib_wc12_weak_beta_decay>": 2.5475016632596765e-07,
                 "N13 + p --> O14 <reaclib_lg06>": 4.8517620910445875e-06,
                 "N14 + p --> O15 <reaclib_im05>": 9.8137074572314962e-07,
                 "N15 + p --> He4 + C12 <reaclib_nacr>": 0.087518552257659241,
                 "O14 --> N14 <reaclib_wc12_weak_electron_capture>": 2.0036691481625654e-06,
                 "O15 --> N15 <reaclib_wc12_weak_electron_capture>": 1.0822012944765837e-06}

        rv = rc.evaluate_rates(1.e4, 1.e8, c)

        for r in rv:
            assert rv[r] == approx(rates[r.id])

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

He4
  consumed by:
  produced by:
     N15 + p ⟶ He4 + C12

C12
  consumed by:
     C12 + p ⟶ N13 + 𝛾
  produced by:
     N15 + p ⟶ He4 + C12

C13
  consumed by:
     C13 + p ⟶ N14 + 𝛾
  produced by:
     N13 ⟶ C13 + e⁺ + 𝜈

N13
  consumed by:
     N13 ⟶ C13 + e⁺ + 𝜈
     N13 + p ⟶ O14 + 𝛾
  produced by:
     C12 + p ⟶ N13 + 𝛾

N14
  consumed by:
     N14 + p ⟶ O15 + 𝛾
  produced by:
     C13 + p ⟶ N14 + 𝛾
     O14 ⟶ N14 + e⁺ + 𝜈

N15
  consumed by:
     N15 + p ⟶ He4 + C12
  produced by:
     O15 ⟶ N15 + e⁺ + 𝜈

O14
  consumed by:
     O14 ⟶ N14 + e⁺ + 𝜈
  produced by:
     N13 + p ⟶ O14 + 𝛾

O15
  consumed by:
     O15 ⟶ N15 + e⁺ + 𝜈
  produced by:
     N14 + p ⟶ O15 + 𝛾
"""
        assert rc.network_overview().replace(" ", "").strip() == ostr.replace(" ", "").strip()
