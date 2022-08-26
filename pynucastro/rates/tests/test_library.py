# unit tests for rates

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
        trates.append(pyna.load_rate("c12-pg-n13-ls09"))
        trates.append(pyna.load_rate("c13-pg-n14-nacr"))
        trates.append(pyna.load_rate("n13--c13-wc12"))
        trates.append(pyna.load_rate("n13-pg-o14-lg06"))
        trates.append(pyna.load_rate("n14-pg-o15-im05"))
        trates.append(pyna.load_rate("n15-pa-c12-nacr"))
        trates.append(pyna.load_rate("o14--n14-wc12"))
        trates.append(pyna.load_rate("o15--n15-wc12"))

        self.library = pyna.Library(rates=trates)

        self.removed_rates = []
        self.removed_rates.append(trates.pop())
        self.removed_rates.append(trates.pop())

        self.smaller_lib = pyna.Library(rates=trates)

    def teardown_method(self):
        """ this is run after each test """
        self.library = None
        self.smaller_lib = None
        self.removed_rates = None

    def test_heaviest(self):
        assert self.library.heaviest() == pyna.Nucleus("n15")

    def test_lightest(self):
        assert self.library.lightest() == pyna.Nucleus("p")

    def test_get_num_rates(self):
        assert self.library.get_num_rates() == 8

    def test_get_rate(self):
        assert self.library.get_rate("c12 + p --> n13 <ls09_reaclib__>") == pyna.load_rate("c12-pg-n13-ls09")

    def test_diff(self):
        diff_lib = self.library - self.smaller_lib
        assert sorted(diff_lib.get_rates()) == sorted(self.removed_rates)

    def test_linking_nuclei(self):
        new_lib = self.library.linking_nuclei(["p", "c12", "n13", "c13"])

        assert sorted(new_lib.get_rates()) == sorted([pyna.load_rate("c12-pg-n13-ls09"),
                                                      pyna.load_rate("n13--c13-wc12")])

    def test_forward_backward(self):
        assert self.library.backward() is None
