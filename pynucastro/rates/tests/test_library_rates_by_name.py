# unit tests for rates


class TestGetRatesByName:

    def test_reaclib_rates(self, reaclib_library):

        r1 = "c12(a,g)o16"
        r2 = "c12(c12,a)ne20"
        r3 = "a(aa,g)c12"
        r4 = "f20(e,nu)o20"   # not found in ReacLib
        r5 = "f20(,e)ne20"
        r6 = "he3(he3,pp)he4"

        assert reaclib_library.get_rate_by_name(r1).fname == "He4_C12_to_O16"
        assert reaclib_library.get_rate_by_name(r2).fname == "C12_C12_to_He4_Ne20"
        assert reaclib_library.get_rate_by_name(r3).fname == "He4_He4_He4_to_C12"
        assert reaclib_library.get_rate_by_name(r4) is None
        assert reaclib_library.get_rate_by_name(r5).fname == "F20_to_Ne20_weak_wc12"
        assert reaclib_library.get_rate_by_name(r6).fname == "He3_He3_to_p_p_He4"

    def test_tabular_rates(self, tabular_library):

        r4 = "f20(e,nu)o20"
        r5 = "f20(,e)ne20"
        r6 = "ti45(e,nu)sc45"
        r7 = "ti45(,e)v45"

        assert tabular_library.get_rate_by_name(r4).fname == "F20_to_O20"
        assert tabular_library.get_rate_by_name(r5).fname == "F20_to_Ne20"
        assert tabular_library.get_rate_by_name(r6).fname == "Ti45_to_Sc45"
        assert tabular_library.get_rate_by_name(r7).fname == "Ti45_to_V45"
