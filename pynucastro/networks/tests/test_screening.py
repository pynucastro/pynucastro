import pytest
from pytest import approx

from pynucastro import networks
from pynucastro.screening import chugunov_2007, chugunov_2009


class TestScreening:
    @pytest.fixture(scope="class")
    def rc(self):
        files = ["c12-ag-o16-nac2",
                 "c12-c12a-ne20-cf88",
                 "c12-c12n-mg23-cf88",
                 "c12-c12p-na23-cf88",
                 "he4-aag-c12-fy05"]

        return networks.RateCollection(files)

    def test_screening_map(self, rc):

        screening_map = rc.get_screening_map()

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

        rates = {"c12 + he4 --> o16 <nac2_reaclib__>": 5.794539791829924,
                 "c12 + c12 --> he4 + ne20 <cf88_reaclib__>": 103.21274049093526,
                 "c12 + c12 --> n + mg23 <cf88_reaclib__reverse>": 103.21274049093526,
                 "c12 + c12 --> p + na23 <cf88_reaclib__>": 103.21274049093526,
                 "he4 + he4 + he4 --> c12 <fy05_reaclib__>": 6.502599619793744}

        factors = rc.evaluate_screening(1.e6, 1.e8, c, screen_func=chugunov_2007)

        for r, factor in factors.items():
            assert factor == approx(rates[r.get_rate_id()])

    def test_screening_chugunov_2009(self, rc):
        c = networks.Composition(rc.unique_nuclei)
        c.set_solar_like()

        rates = {"c12 + he4 --> o16 <nac2_reaclib__>": 4.405674333522246,
                 "c12 + c12 --> he4 + ne20 <cf88_reaclib__>": 89.6640543016441,
                 "c12 + c12 --> n + mg23 <cf88_reaclib__reverse>": 89.6640543016441,
                 "c12 + c12 --> p + na23 <cf88_reaclib__>": 89.6640543016441,
                 "he4 + he4 + he4 --> c12 <fy05_reaclib__>": 4.380701422122169}

        factors = rc.evaluate_screening(1.e6, 1.e8, c, screen_func=chugunov_2009)

        for r, factor in factors.items():
            assert factor == approx(rates[r.get_rate_id()])
