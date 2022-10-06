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
        assert rr[0].fname == "o16__he4_c12"
        assert rr[1].fname == "c12__he4_he4_he4"
        assert rr[2].fname == "he4_c12__o16"
        assert rr[3].fname == "he4_he4_he4__c12"

    def test_get_rate(self, rc):
        r = rc.get_rate("he4_he4_he4__c12")
        assert r.fname == "he4_he4_he4__c12"

    def test_find_reverse(self, rc):
        rr = rc.find_reverse(rc.get_rate("he4_c12__o16"))
        assert rr.fname == "o16__he4_c12"

    def test_evaluate_energy_gen(self, rc):
        # define a composition
        comp = pyna.Composition(rc.unique_nuclei)

        comp.set_all(0.3)
        comp.normalize()

        rho = 1e5
        T = 1e8
        assert rc.evaluate_energy_generation(rho, T, comp) == 32.24796008826701


class TestUnimportantRates:
    @pytest.fixture(scope="class")
    def rc(self):
        files = ["c12-pg-n13-ls09",
                 "c13-pg-n14-nacr",
                 "n13--c13-wc12",
                 "n13-pg-o14-lg06",
                 "n14-pg-o15-im05",
                 "n15-pa-c12-nacr",
                 "o14--n14-wc12",
                 "o15--n15-wc12"]
        return pyna.RateCollection(files)

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
