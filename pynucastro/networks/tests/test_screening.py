import pytest
from pytest import approx

from pynucastro import networks
from pynucastro.screening import (chugunov_2007, chugunov_2009,
                                  get_screening_map)


class TestScreening:
    @pytest.fixture(scope="class")
    def rc(self, reaclib_library):
        rate_names = ["c12(a,g)o16",
                      "c12(c12,a)ne20",
                      "c12(c12,n)mg23",
                      "c12(c12,p)na23",
                      "he4(aa,g)c12"]
        rates = reaclib_library.get_rate_by_name(rate_names)

        return networks.RateCollection(rates=rates)

    def test_screening_map(self, rc):

        screening_map = get_screening_map(rc.get_rates())

        assert len(screening_map) == 4
        assert len(screening_map[0].rates) == 1
        assert len(screening_map[1].rates) == 3
        assert len(screening_map[2].rates) == 1
        assert len(screening_map[3].rates) == 1
        # two triple-alpha screening steps
        assert screening_map[2].rates[0] == screening_map[3].rates[0] == rc.rates[4]

    def test_screening_chugunov_2007(self, rc):
        c = networks.Composition(rc.unique_nuclei)
        c.set_solar_like()

        rates = {"C12 + He4 --> O16 <nac2_reaclib__>": 5.794539791829924,
                 "C12 + C12 --> He4 + Ne20 <cf88_reaclib__>": 103.21274049093526,
                 "C12 + C12 --> n + Mg23 <cf88_reaclib__derived_from_inverse>": 103.21274049093526,
                 "C12 + C12 --> p + Na23 <cf88_reaclib__>": 103.21274049093526,
                 "3 He4 --> C12 <fy05_reaclib__>": 6.502599619793744}

        factors = {}
        for r in rc.get_rates():
            factors[r] = r.evaluate_screening(1.e6, 1.e8, c, screen_func=chugunov_2007)

        for r, factor in factors.items():
            assert factor == approx(rates[r.id])

    def test_screening_chugunov_2009(self, rc):
        c = networks.Composition(rc.unique_nuclei)
        c.set_solar_like()

        rates = {"C12 + He4 --> O16 <nac2_reaclib__>": 4.405674333522246,
                 "C12 + C12 --> He4 + Ne20 <cf88_reaclib__>": 89.6640543016441,
                 "C12 + C12 --> n + Mg23 <cf88_reaclib__derived_from_inverse>": 89.6640543016441,
                 "C12 + C12 --> p + Na23 <cf88_reaclib__>": 89.6640543016441,
                 "3 He4 --> C12 <fy05_reaclib__>": 4.380701422122169}

        factors = {}
        for r in rc.get_rates():
            factors[r] = r.evaluate_screening(1.e6, 1.e8, c, screen_func=chugunov_2009)

        for r, factor in factors.items():
            assert factor == approx(rates[r.id])
