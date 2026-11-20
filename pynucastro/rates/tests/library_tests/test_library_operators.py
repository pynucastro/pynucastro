# unit operators (+, -, in) for Library

import pytest

from pynucastro.rates import Library, Rate


class TestLibrary:

    @pytest.fixture(scope="class")
    @classmethod
    def lib(cls):

        a = Rate(reactants=["c12", "he4"], products=["o16"], label="one")
        b = Rate(reactants=["c12", "he4"], products=["o16"], label="two")
        c = Rate(reactants=["c12", "c12"], products=["mg24"], label="cc")
        return Library(rates=[a, b, c])

    def test_add(self, lib):
        assert lib.num_rates == 3

        d = Rate(reactants=["c12", "o16"], products=["si28"], label="co")
        new_lib = lib + Library(rates=[d])

        assert d in new_lib
        assert new_lib.num_rates == 4

        assert d not in lib

    def test_sub(self, lib):

        a = Rate(reactants=["c12", "he4"], products=["o16"], label="one")
        new_lib = lib - Library(rates=[a])

        assert new_lib.num_rates == 2
        assert a not in new_lib
