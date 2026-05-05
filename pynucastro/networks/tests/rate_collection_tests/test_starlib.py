# Test that a StarLib rate is stored separately from a
# TemperatureTabularRate

import pytest

from pynucastro.networks import RateCollection
from pynucastro.rates.alternate_rates import IliadisO16pgF17


class TestStarLibNet:
    @pytest.fixture(scope="class")
    def rc(self, reaclib_library, starlib_library):
        nuc = ["p", "n15", "o16"]
        lib = reaclib_library.linking_nuclei(nuc, with_reverse=False)

        r = IliadisO16pgF17()
        lib.add_rate(r)

        r2 = starlib_library.get_rate_by_name("c12(a,g)o16")
        lib.add_rate(r2)

        return RateCollection(libraries=[lib])

    def test_num_rates(self, rc):
        assert len(rc.all_rates) == 3

    def test_num_starlib_rates(self, rc):
        assert len(rc.starlib_rates) == 1

    def test_num_temperature_tabular_rates(self, rc):
        assert len(rc.temperature_tabular_rates) == 1

    def test_num_reaclib_rates(self, rc):
        assert len(rc.reaclib_rates) == 1
