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
    def all_lib(self):
        # we need to be careful about mutating the reaclib_library or
        # tabular_library fixtures, so we'l create our own private
        # copies here
        rl = pyna.ReacLibLibrary()
        tl = pyna.TabularLibrary()
        return rl + tl

    @pytest.fixture(scope="class")
    def pp_reaclib(self, reaclib_library):
        return reaclib_library.linking_nuclei(["p", "d"])

    @pytest.fixture(scope="class")
    def pp_starlib(self, starlib_library):
        return starlib_library.linking_nuclei(["p", "d"])

    def check_pp_dupes_reaclib(self, pp_reaclib):
        assert len(pp_reaclib.find_duplicate_links()) == 0

    def check_pp_dupes_starlib(self, pp_starlib):
        assert len(pp_starlib.find_duplicate_links()) == 0

    def check_pp_dupes_bothlibs(self, pp_reaclib, pp_starlib):
        newlib = pp_reaclib + pp_starlib
        assert len(newlib.find_duplicate_links()) == 1

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
