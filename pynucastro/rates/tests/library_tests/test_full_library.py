import pytest

import pynucastro as pyna
from pynucastro.rates import ReacLibRate, StarLibRate


class TestFullLibrary:
    @pytest.fixture(scope="class")
    @classmethod
    def full_lib(cls):
        return pyna.full_library()

    def test_ni56_rates(self, full_lib):
        rf = pyna.RateFilter(reactants=["ni56"], exact=False)
        new_lib = full_lib.filter(rf)
        dupes = new_lib.find_duplicate_links()

        assert len(dupes) == 13

        # most of the dupes are strong-mediated rates for which there
        # are just ReacLib and StarLib versions
        for dd in dupes:
            if len(dd) == 2:
                if isinstance(dd[0], ReacLibRate):
                    assert isinstance(dd[1], StarLibRate)
                elif isinstance(dd[0], StarLibRate):
                    assert isinstance(dd[1], ReacLibRate)
                else:
                    # shouldn't get here
                    raise AssertionError("invalid duple")

        # there are 4 electron capture rates
        assert len(dupes[0]) == 4
