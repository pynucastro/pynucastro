# unit tests for rates
import pytest

import pynucastro as pyna


class TestDuplicates:

    @pytest.fixture(scope="class")
    def ecsn_lib(self, reaclib_library, tabular_library):
        all_nuclei = ["p", "he4",
                      "ne20", "o20", "f20",
                      "mg24", "al27", "o16",
                      "si28", "s32", "p31"]
        ecsn_rl_lib = reaclib_library.linking_nuclei(all_nuclei, with_reverse=True)

        ecsn_tabular_lib = tabular_library.linking_nuclei(["f20", "o20", "ne20"])

        return ecsn_rl_lib + ecsn_tabular_lib

    def test_rate_collection(self, ecsn_lib):

        with pytest.raises(pyna.networks.RateDuplicationError):
            pyna.RateCollection(libraries=[ecsn_lib])

    def test_library_find_duplicate_links(self, ecsn_lib):

        dupes = ecsn_lib.find_duplicate_links()

        assert len(dupes) == 2

        d1 = dupes[0][0]
        assert d1.reactants == [pyna.Nucleus("o20")]
        assert d1.products == [pyna.Nucleus("f20")]

        d2 = dupes[1][0]
        assert d2.reactants == [pyna.Nucleus("f20")]
        assert d2.products == [pyna.Nucleus("ne20")]
