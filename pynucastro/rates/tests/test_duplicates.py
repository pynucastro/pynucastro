import pytest

import pynucastro as pyna


class TestDuplicates:
    @pytest.fixture(scope="class")
    def nuclei(self):
        nuclei = (["p", "he4"] +
                  pyna.get_nuclei_in_range("fe", A_range=(52, 56)) +
                  pyna.get_nuclei_in_range("co", A_range=(54, 58)) +
                  pyna.get_nuclei_in_range("ni", A_range=(56, 60)))
        return nuclei

    @pytest.fixture(scope="class")
    def all_lib(self, reaclib_library, tabular_library):
        return reaclib_library + tabular_library

    def test_find_duplicate_links(self, nuclei, all_lib):
        lib = all_lib.linking_nuclei(nuclei)
        assert len(lib.find_duplicate_links()) == 5

    def test_eliminate_duplicates1(self, nuclei, all_lib):
        lib = all_lib.linking_nuclei(nuclei)
        assert lib.num_rates == 59

        lib.eliminate_duplicates(rate_type_preference="tabular")

        assert lib.num_rates == 54

        num_tabular = len([r for r in lib.get_rates()
                           if isinstance(r, pyna.rates.TabularRate)])

        assert num_tabular == 12

    def test_eliminate_duplicates2(self, nuclei, all_lib):
        lib = all_lib.linking_nuclei(nuclei)
        assert lib.num_rates == 59

        lib.eliminate_duplicates(rate_type_preference="reaclib")

        assert lib.num_rates == 54

        num_tabular = len([r for r in lib.get_rates()
                           if isinstance(r, pyna.rates.TabularRate)])

        assert num_tabular == 7

