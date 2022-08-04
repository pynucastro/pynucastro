from pynucastro import networks
import pytest


class TestScreening:
    @pytest.fixture(scope="class")
    def rc(self):
        files = ["c12-ag-o16-nac2",
                 "c12-c12a-ne20-cf88",
                 "c12-c12n-mg23-cf88",
                 "c12-c12p-na23-cf88"]

        return networks.RateCollection(files)

    def test_screening(self, rc):

        screening_map = rc.get_screening_map()

        assert len(screening_map[0].rates) == 1
        assert len(screening_map[1].rates) == 3
