# unit tests for rates
import pytest

import pynucastro as pyna


class TestDuplicates:

    @pytest.fixture(scope="class")
    def ecsn_lib(self, reaclib_library):
        all_nuclei = ["p", "he4",
                      "ne20", "o20", "f20",
                      "mg24", "al27", "o16",
                      "si28", "s32", "p31"]
        ecsn_rl_lib = reaclib_library.linking_nuclei(all_nuclei, with_reverse=True)

        tabular_lib = pyna.TabularLibrary()
        ecsn_tabular_lib = tabular_lib.linking_nuclei(["f20", "o20", "ne20"])

        return ecsn_rl_lib + ecsn_tabular_lib

    def test_validate(self, ecsn_lib):

        with pytest.raises(pyna.networks.RateDuplicationError):
            pyna.RateCollection(libraries=[ecsn_lib])
