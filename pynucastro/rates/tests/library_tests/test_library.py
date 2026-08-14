# unit tests for rates

import pytest

import pynucastro as pyna


class TestLibrary:

    @pytest.fixture(scope="class")
    @classmethod
    def regular_library(cls, reaclib_library):

        trates = []
        trates.append(reaclib_library.get_rate_by_name("c12(p,g)n13"))
        trates.append(reaclib_library.get_rate_by_name("c13(p,g)n14"))
        trates.append(reaclib_library.get_rate_by_name("n13(,)c13"))
        trates.append(reaclib_library.get_rate_by_name("n13(p,g)o14"))
        trates.append(reaclib_library.get_rate_by_name("n14(p,g)o15"))
        trates.append(reaclib_library.get_rate_by_name("n15(p,a)c12"))
        trates.append(reaclib_library.get_rate_by_name("o14(,)n14"))
        trates.append(reaclib_library.get_rate_by_name("o15(,)n15"))

        return pyna.Library(rates=trates)

    @pytest.fixture(scope="class")
    @classmethod
    def smaller_library(cls, regular_library):

        rates = regular_library.get_rates()
        rates.pop()
        rates.pop()
        return pyna.Library(rates=rates)

    @pytest.fixture(scope="class")
    @classmethod
    def removed_rates(cls, regular_library):
        trates = regular_library.get_rates()
        rates = []
        rates.append(trates.pop())
        rates.append(trates.pop())
        return pyna.Library(rates=rates)

    def test_heaviest(self, regular_library):
        assert regular_library.heaviest() == pyna.Nucleus("n15")

    def test_lightest(self, regular_library):
        assert regular_library.lightest() == pyna.Nucleus("p")

    def test_num_rates(self, regular_library):
        assert regular_library.num_rates == 8

    def test_get_rate(self, regular_library, reaclib_library):
        # get by rate id
        assert (regular_library.get_rate("c12 + p --> n13 <reaclib_ls09>") ==
                reaclib_library.get_rate_by_name("c12(p,g)n13"))

        # get by fname
        assert (regular_library.get_rate("p_N14_to_O15_reaclib") ==
                reaclib_library.get_rate_by_name("n14(p,g)o15"))

        # get by fname
        assert (regular_library.get_rate("p_n15_to_he4_c12_reaclib") ==
                reaclib_library.get_rate_by_name("n15(p,a)c12"))

        # get by fname without label, i.e. use base name
        assert (regular_library.get_rate("p_N14_to_O15") ==
                reaclib_library.get_rate_by_name("n14(p,g)o15"))

        # get by fname without label, lowercase
        assert (regular_library.get_rate("p_n15_to_he4_c12") ==
                reaclib_library.get_rate_by_name("n15(p,a)c12"))

    def test_get_rate_failure(self, regular_library):
        # missing rate id
        with pytest.raises(LookupError):
            regular_library.get_rate("N15 + p --> O16 <reaclib_li10>")

        # missing fname
        with pytest.raises(LookupError):
            regular_library.get_rate("F18_to_He4_N14")

        # invalid rate
        with pytest.raises(LookupError):
            regular_library.get_rate("this is not a rate")

    def test_get_rate_by_nuclei(self, regular_library, reaclib_library, suzuki_library):
        assert (regular_library.get_rate_by_nuclei([pyna.Nucleus("p"), pyna.Nucleus("c13")],
                                                   [pyna.Nucleus("n14")]) ==
                reaclib_library.get_rate_by_name("c13(p,g)n14"))

        assert regular_library.get_rate_by_nuclei([pyna.Nucleus("p")],
                                                  [pyna.Nucleus("n14")]) is None

        dup_rates = [reaclib_library.get_rate_by_name("f17(,)o17"),
                     suzuki_library.get_rate_by_name("f17(,)o17")]
        dup_lib = regular_library + pyna.Library(rates=dup_rates)

        assert (dup_lib.get_rate_by_nuclei([pyna.Nucleus("f17")],
                                           [pyna.Nucleus("o17")]) ==
                dup_rates)

    def test_diff(self, regular_library, smaller_library, removed_rates):
        diff_lib = regular_library - smaller_library
        assert sorted(diff_lib.get_rates()) == sorted(removed_rates.get_rates())

    def test_linking_nuclei(self, regular_library, reaclib_library):
        new_lib = regular_library.linking_nuclei(["p", "c12", "n13", "c13"])

        assert (sorted(new_lib.get_rates()) ==
                sorted([reaclib_library.get_rate_by_name("c12(p,g)n13"),
                        reaclib_library.get_rate_by_name("n13(,)c13")]))

    def test_forward_backward(self, regular_library):
        assert regular_library.backward() is None
