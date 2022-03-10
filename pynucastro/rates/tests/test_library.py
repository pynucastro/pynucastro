# unit tests for rates
import math

import pynucastro as pyna


class TestLibrary:
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

        trates = []
        trates.append(pyna.Rate("c12-pg-n13-ls09"))
        trates.append(pyna.Rate("c13-pg-n14-nacr"))
        trates.append(pyna.Rate("n13--c13-wc12"))
        trates.append(pyna.Rate("n13-pg-o14-lg06"))
        trates.append(pyna.Rate("n14-pg-o15-im05"))
        trates.append(pyna.Rate("n15-pa-c12-nacr"))
        trates.append(pyna.Rate("o14--n14-wc12"))
        trates.append(pyna.Rate("o15--n15-wc12"))

        self.library = pyna.Library(rates=trates)

    def teardown_method(self):
        """ this is run after each test """
        self.library = None

    def test_heaviest(self):
        assert self.library.heaviest() == pyna.Nucleus("n15")

    def test_lightest(self):
        assert self.library.lightest() == pyna.Nucleus("p")

    def test_get_num_rates(self):
        assert self.library.get_num_rates() == 8

    def test_get_rate(self):
        assert self.library.get_rate("c12 + p --> n13 <ls09_reaclib__>") == pyna.Rate("c12-pg-n13-ls09")
