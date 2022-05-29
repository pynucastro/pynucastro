# unit tests for a rate collection
import pynucastro as pyna


class TestRateCollection:
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """
        pass

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        pass

    def setup_method(self):
        """ this is run before each test """
        lib = pyna.ReacLibLibrary()
        mylib = lib.linking_nuclei(["he4", "c12", "o16"])
        self.rc = pyna.RateCollection(libraries=[mylib])

    def teardown_method(self):
        """ this is run after each test """
        self.rc = 0

    def test_get_ratesg(self):
        rr = self.rc.get_rates()

        assert len(rr) == 4
        assert rr[0].fname == "o16__he4_c12"
        assert rr[1].fname == "c12__he4_he4_he4"
        assert rr[2].fname == "he4_c12__o16"
        assert rr[3].fname == "he4_he4_he4__c12"

    def test_get_rate(self):
        r = self.rc.get_rate("he4_he4_he4__c12")
        assert r.fname == "he4_he4_he4__c12"

    def test_find_reverse(self):
        rr = self.rc.find_reverse(self.rc.get_rate("he4_c12__o16"))
        assert rr.fname == "o16__he4_c12"
