# unit tests for rates

import pytest

import pynucastro as pyna


class TestLibrary:
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """

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

    def test_num_rates(self):
        assert self.library.num_rates == 8

    def test_get_rate(self):
        # get by rate id
        assert self.library.get_rate("c12 + p --> n13 <reaclib_ls09>") == pyna.load_rate("c12-pg-n13-ls09")

        # get by fname
        assert self.library.get_rate("p_N14_to_O15_reaclib") == pyna.load_rate("n14-pg-o15-im05")

        # get by fname
        assert self.library.get_rate("p_n15_to_he4_c12_reaclib") == pyna.load_rate("n15-pa-c12-nacr")

        # get by fname without label, i.e. use base name
        assert self.library.get_rate("p_N14_to_O15") == pyna.load_rate("n14-pg-o15-im05")

        # get by fname without label, lowercase
        assert self.library.get_rate("p_n15_to_he4_c12") == pyna.load_rate("n15-pa-c12-nacr")

    def test_get_rate_failure(self):
        # missing rate id
        with pytest.raises(LookupError):
            self.library.get_rate("N15 + p --> O16 <reaclib_li10>")

        # missing fname
        with pytest.raises(LookupError):
            self.library.get_rate("F18_to_He4_N14")

        # invalid rate
        with pytest.raises(LookupError):
            self.library.get_rate("this is not a rate")

    def test_get_rate_by_nuclei(self):
        assert self.library.get_rate_by_nuclei(
            [pyna.Nucleus("p"), pyna.Nucleus("c13")], [pyna.Nucleus("n14")]
        ) == pyna.load_rate("c13-pg-n14-nacr")

        assert self.library.get_rate_by_nuclei(
            [pyna.Nucleus("p")], [pyna.Nucleus("n14")]
        ) is None

        dup_rates = [pyna.load_rate("f17--o17-wc12"), pyna.load_rate("suzuki-17f-17o_electroncapture.dat")]
        dup_lib = self.library + pyna.Library(rates=dup_rates)

        assert dup_lib.get_rate_by_nuclei(
            [pyna.Nucleus("f17")], [pyna.Nucleus("o17")]
        ) == dup_rates

    def test_diff(self):
        diff_lib = self.library - self.smaller_lib
        assert sorted(diff_lib.get_rates()) == sorted(self.removed_rates)

    def test_linking_nuclei(self):
        new_lib = self.library.linking_nuclei(["p", "c12", "n13", "c13"])

        assert sorted(new_lib.get_rates()) == sorted([pyna.load_rate("c12-pg-n13-ls09"),
                                                      pyna.load_rate("n13--c13-wc12")])

    def test_forward_backward(self):
        assert self.library.backward() is None
