import pytest

import pynucastro as pyna


class TestTabularLibrary:
    @pytest.fixture(scope="class")
    def tl_default(self):
        return pyna.TabularLibrary()

    def test_number_of_rates(self, tl_default):
        assert tl_default.num_rates == 458

        suzuki_rates = [r for r in tl_default.get_rates() if r.rfile.startswith("suzuki")]
        assert len(suzuki_rates) == 61

        langanke_rates = [r for r in tl_default.get_rates() if r.rfile.startswith("langanke")]
        assert len(langanke_rates) == 212

        ffn_rates = [r for r in tl_default.get_rates() if r.rfile.startswith("ffn")]
        assert len(ffn_rates) == 88

        oda_rates = [r for r in tl_default.get_rates() if r.rfile.startswith("oda")]
        assert len(oda_rates) == 97

        # make sure the sum of all the different sources equals the
        # total number of rates
        assert len(suzuki_rates) + len(langanke_rates) + len(ffn_rates) + len(oda_rates) == len(tl_default.get_rates())

    def test_ordering(self):

        tl_new = pyna.TabularLibrary(ordering=["suzuki", "oda", "langanke", "ffn"])

        assert tl_new.num_rates == 458

        suzuki_rates = [r for r in tl_new.get_rates() if r.rfile.startswith("suzuki")]
        assert len(suzuki_rates) == 0

        langanke_rates = [r for r in tl_new.get_rates() if r.rfile.startswith("langanke")]
        assert len(langanke_rates) == 62

        ffn_rates = [r for r in tl_new.get_rates() if r.rfile.startswith("ffn")]
        assert len(ffn_rates) == 376

        oda_rates = [r for r in tl_new.get_rates() if r.rfile.startswith("oda")]
        assert len(oda_rates) == 20

        # make sure the sum of all the different sources equals the
        # total number of rates
        assert len(suzuki_rates) + len(langanke_rates) + len(ffn_rates) + len(oda_rates) == len(tl_new.get_rates())

    def test_sources(self):

        ffn_lib = pyna.FFNLibrary()
        assert ffn_lib.get_rates()[0].source["Label"] == "ffn"

        suzuki_lib = pyna.SuzukiLibrary()
        assert suzuki_lib.get_rates()[0].source["Label"] == "suzuki"

        langanke_lib = pyna.LangankeLibrary()
        assert langanke_lib.get_rates()[0].source["Label"] == "langanke"

        oda_lib = pyna.OdaLibrary()
        assert oda_lib.get_rates()[0].source["Label"] == "oda"
