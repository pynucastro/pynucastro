# test our custom copying of a Rate object

import copy

import pytest

from pynucastro.nucdata import Nucleus


class TestRateCopy:
    @pytest.fixture(scope="class")
    def rl_rate(self, reaclib_library):
        return reaclib_library.get_rate_by_name("c12(a,g)o16")

    @pytest.fixture(scope="class")
    def tl_rate(self, tabular_library):
        return tabular_library.get_rate_by_name("ni56(,)co56")

    def test_copy_rl(self, rl_rate):

        rl_copy = copy.copy(rl_rate)

        # first check that the copy preserves things

        assert isinstance(rl_copy, type(rl_rate))
        assert rl_copy == rl_rate
        assert rl_copy.reactants == rl_rate.reactants
        assert rl_copy.id == rl_rate.id

        # now change the copy check that they are different

        rl_copy.stoichiometry = {Nucleus("he4"): 2, Nucleus("o16"): 1.25}
        rl_copy._set_print_representation()  # pylint: disable=protected-access
        assert rl_copy.stoichiometry != rl_rate.stoichiometry
        assert rl_copy.id != rl_rate.id

    def test_copy_tl(self, tl_rate):

        tl_copy = copy.copy(tl_rate)

        # first check that the copy preserves things

        assert isinstance(tl_copy, type(tl_rate))
        assert tl_copy == tl_rate
        assert tl_copy.reactants == tl_rate.reactants
        assert tl_copy.id == tl_rate.id

        # now change the copy and check that they are different

        tl_copy.reactants = [Nucleus("co56")]
        tl_copy.products = [Nucleus("fe56")]
        tl_copy._set_print_representation()  # pylint: disable=protected-access
        assert tl_copy != tl_rate
