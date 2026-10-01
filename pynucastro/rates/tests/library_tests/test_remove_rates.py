# focused tests on the interface for remove_rates

# unit tests for rates

import pytest

from pynucastro.rates import Library


class TestLibrary:

    @pytest.fixture(scope="class")
    @classmethod
    def rates(cls, reaclib_library):
        rates = []
        rates.append(reaclib_library.get_rate_by_name("a(aa,g)c12"))
        rates.append(reaclib_library.get_rate_by_name("c12(a,g)o16"))
        return rates

    def test_remove_by_id(self, rates):
        lib = Library(rates=rates)

        lib.remove_rate(rates[0].id)
        assert len(lib.get_rates()) == 1

    def test_remove_by_rid(self, rates):
        lib = Library(rates=rates)

        # this will fail because rid is not unique it is simply meant
        # to be human readable for comments.
        with pytest.raises(LookupError):
            lib.remove_rate(rates[0].rid)

        assert len(lib.get_rates()) == 2

    def test_remove_by_shortname(self, rates):
        lib = Library(rates=rates)

        lib.remove_rate("a(aa,g)c12")

        assert len(lib.get_rates()) == 1
