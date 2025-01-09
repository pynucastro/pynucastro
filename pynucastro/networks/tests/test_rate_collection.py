# unit tests for a rate collection
import pytest

import pynucastro as pyna


class TestRateCollection:
    @pytest.fixture(scope="class")
    def rc(self, reaclib_library):
        mylib = reaclib_library.linking_nuclei(["he4", "c12", "o16"])
        return pyna.RateCollection(libraries=[mylib])

    def test_get_ratesg(self, rc):
        rr = rc.get_rates()

        assert len(rr) == 4
        assert rr[0].fname == "O16__He4_C12"
        assert rr[1].fname == "C12__He4_He4_He4"
        assert rr[2].fname == "He4_C12__O16"
        assert rr[3].fname == "He4_He4_He4__C12"

    def test_get_rate(self, rc):
        r = rc.get_rate("He4_He4_He4__C12")
        assert r.fname == "He4_He4_He4__C12"

    def test_find_reverse(self, rc):
        rr = rc.find_reverse(rc.get_rate("He4_C12__O16"))
        assert rr.fname == "O16__He4_C12"

    def test_evaluate_energy_gen(self, rc):
        # define a composition
        comp = pyna.Composition(rc.unique_nuclei)

        comp.set_all(0.3)
        comp.normalize()

        rho = 1e5
        T = 1e8
        assert rc.evaluate_energy_generation(rho, T, comp) == 32.247985011082925

    def test_add_rates(self, rc, reaclib_library):
        new_rate_names = ["o16(a,g)ne20",
                          "ne20(a,g)mg24",
                          "mg24(a,g)si28"]

        new_rates = reaclib_library.get_rate_by_name(new_rate_names)

        # test adding only a single rate
        rc.add_rates(new_rates[0])

        # check to see if the rate is added
        assert new_rates[0] in rc.get_rates()

        # now add them all -- and make sure we don't duplicate a rate
        rc.add_rates(new_rates)

        for r in new_rates:
            assert r in rc.get_rates()

        assert rc.get_rates().count(new_rates[0]) == 1

        # note: this modifies the fixture rc so tests that follow
        # will see the additional rates

    def test_remove_rates(self, rc):
        rate = rc.get_rates()[-1]
        rc.remove_rates(rate)

        assert rate not in rc.get_rates()
        assert len(rc.get_rates()) == 6


class TestUnimportantRates:
    @pytest.fixture(scope="class")
    def rc(self, reaclib_library):
        rate_names = ["c12(p,g)n13",
                      "c13(p,g)n14",
                      "n13(,)c13",
                      "n13(p,g)o14",
                      "n14(p,g)o15",
                      "n15(p,a)c12",
                      "o14(,)n14",
                      "o15(,)n15"]
        rates = reaclib_library.get_rate_by_name(rate_names)
        return pyna.RateCollection(rates=rates)

    @pytest.fixture(scope="class")
    def comp(self, rc):
        # define a composition
        comp = pyna.Composition(rc.unique_nuclei)
        comp.set_all(1.0)
        comp.normalize()
        return comp

    def test_temp_1e8(self, rc, comp):
        rho = 1e5
        T = 1e8

        expected = {rc.rates[i] for i in [2, 3, 4, 6, 7]}
        unimportant = rc.find_unimportant_rates([(rho, T, comp)], cutoff_ratio=1e-4)
        assert rc.rates[0] not in unimportant
        assert unimportant.keys() == expected

    def test_temp_1e10(self, rc, comp):
        rho = 1e5
        T = 1e10

        expected = {rc.rates[i] for i in [0, 2, 4, 6, 7]}
        unimportant = rc.find_unimportant_rates([(rho, T, comp)], cutoff_ratio=1e-4)
        assert rc.rates[3] not in unimportant
        assert unimportant.keys() == expected

    def test_temp_both(self, rc, comp):
        states = [(1.e5, 1.e8, comp), (1.e5, 1.e10, comp)]

        expected = {rc.rates[i] for i in [2, 4, 6, 7]}
        unimportant = rc.find_unimportant_rates(states, cutoff_ratio=1e-4)
        # C12(p,g)N13 (rate 0) is important at T=1e8, but not T=1e10.
        # N13(p,g)O14 (rate 3) is important at T=1e10, but not T=1e8.
        # If we include both temperatures, then both rates should be considered important.
        assert rc.rates[0] not in unimportant
        assert rc.rates[3] not in unimportant
        assert unimportant.keys() == expected
