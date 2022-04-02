import pynucastro.networks as networks
import pynucastro.rates as rates

from pytest import approx


class TestScreening:
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """
        pass

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        pass

    def setup_method(self):
        """ this is run before each test """

        files = ["c12-ag-o16-nac2",
                 "c12-c12a-ne20-cf88",
                 "c12-c12n-mg23-cf88",
                 "c12-c12p-na23-cf88"]

        self.rc = networks.RateCollection(files)

    def teardown_method(self):
        """ this is run after each test """
        self.rc = None

    def test_screening(self):

        screening_map = self.rc.get_screening_map()

        assert len(screening_map[0].rates) == 1
        assert len(screening_map[1].rates) == 3
