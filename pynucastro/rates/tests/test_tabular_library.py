import pytest
from pytest import approx

import pynucastro as pyna


class TestTabularLibrary:
    @pytest.fixture(scope="class")
    def tl_default(self):
        return pyna.TabularLibrary()

    def test_number_of_rates(self, tl_default):
        assert tl_default.num_rates == 453

        suzuki_rates = [r for r in tl_default.get_rates() if r.rfile.startswith("suzuki")]
        assert len(suzuki_rates) == 61

        langanke_rates = [r for r in tl_default.get_rates() if r.rfile.startswith("langanke")]
        assert len(langanke_rates) == 212

        ffn_rates = [r for r in tl_default.get_rates() if r.rfile.startswith("ffn")]
        assert len(ffn_rates) == 180

        # make sure the sum of all the different sources equals the
        # total number of rates
        assert len(suzuki_rates) + len(langanke_rates) + len(ffn_rates) == len(tl_default.get_rates())

    def test_ordering(self):

        tl_new = pyna.TabularLibrary(ordering=["suzuki", "langanke", "ffn"])

        assert tl_new.num_rates == 453

        suzuki_rates = [r for r in tl_new.get_rates() if r.rfile.startswith("suzuki")]
        assert len(suzuki_rates) == 15

        langanke_rates = [r for r in tl_new.get_rates() if r.rfile.startswith("langanke")]
        assert len(langanke_rates) == 62

        ffn_rates = [r for r in tl_new.get_rates() if r.rfile.startswith("ffn")]
        assert len(ffn_rates) == 376

        # make sure the sum of all the different sources equals the
        # total number of rates
        assert len(suzuki_rates) + len(langanke_rates) + len(ffn_rates) == len(tl_new.get_rates())

