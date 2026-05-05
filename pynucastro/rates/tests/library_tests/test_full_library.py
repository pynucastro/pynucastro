import pytest

import pynucastro as pyna


class TestFullLibrary:
    @pytest.fixture(scope="class")
    def full_lib(self):
        return pyna.full_library()

    def test_ni56_rates(self, full_lib):
        rf = pyna.RateFilter(reactants=["ni56"], exact=False)
        new_lib = full_lib.filter(rf)
        dupes = new_lib.find_duplicate_links()

        assert len(dupes) == 1

        # there are 3 electron capture rates
        assert len(dupes[0]) == 3
