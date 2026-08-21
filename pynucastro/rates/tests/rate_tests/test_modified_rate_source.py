# test the `rate_source=` keyword argument to `ModifiedRate`

import pytest

import pynucastro as pyna


class TestModifiedRateSource:

    @pytest.fixture(scope="class")
    @classmethod
    def pprates(cls, reaclib_library):
        return reaclib_library.get_rate_by_name("p(p,)d")

    def test_without_source(self, pprates):

        rpp, rpep = pprates

        rppp_he3 = pyna.ModifiedRate(rpp,
                             new_products=[pyna.Nucleus("he3")],
                             stoichiometry={pyna.Nucleus("p"): 3})

        rpepp_he3 = pyna.ModifiedRate(rpep,
                              new_products=[pyna.Nucleus("he3")],
                              stoichiometry={pyna.Nucleus("p"): 3})

        # these two rates have the same id, so an error will be raised
        with pytest.raises(ValueError,
                           match="supplied a Rate object already in the Library: 3 p --> He3 <modified>"):
            _ = pyna.Library(rates=[rppp_he3, rpepp_he3])
        
    def test_with_source(self, pprates):

        rpp, rpep = pprates

        rppp_he3 = pyna.ModifiedRate(rpp,
                             new_products=[pyna.Nucleus("he3")],
                             stoichiometry={pyna.Nucleus("p"): 3},
                             rate_source=rpp.src)

        assert rppp_he3.id == "3 p --> He3 <modified_bet+>"

        rpepp_he3 = pyna.ModifiedRate(rpep,
                              new_products=[pyna.Nucleus("he3")],
                              stoichiometry={pyna.Nucleus("p"): 3},
                              rate_source=rpep.src)

        assert rpepp_he3.id == "3 p --> He3 <modified_ec>"

        lib = pyna.Library(rates=[rppp_he3, rpepp_he3])

        assert len(lib.get_rates()) == 2
