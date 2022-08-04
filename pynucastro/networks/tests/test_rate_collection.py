# unit tests for a rate collection
import pynucastro as pyna
import pytest


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
