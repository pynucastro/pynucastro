# unit tests for rates


class TestGetRatesByName:

    def test_reaclib_rates(self, reaclib_library):

        r1 = "c12(a,g)o16"
        r2 = "c12(c12,a)ne20"
        r3 = "a(aa,g)c12"
        r4 = "f20(e,nu)o20"   # not found in ReacLib
        r5 = "f20(,e)ne20"

        assert reaclib_library.get_rate_by_name(r1).fname == "he4_c12__o16"
        assert reaclib_library.get_rate_by_name(r2).fname == "c12_c12__he4_ne20"
        assert reaclib_library.get_rate_by_name(r3).fname == "he4_he4_he4__c12"
        assert reaclib_library.get_rate_by_name(r4) is None
        assert reaclib_library.get_rate_by_name(r5).fname == "f20__ne20__weak__wc12"

    def test_tabular_rates(self, tabular_library):

        r4 = "f20(e,nu)o20"   # note found in ReacLib
        r5 = "f20(,e)ne20"

        assert tabular_library.get_rate_by_name(r4).fname == "f20__o20"
        assert tabular_library.get_rate_by_name(r5).fname == "f20__ne20"
