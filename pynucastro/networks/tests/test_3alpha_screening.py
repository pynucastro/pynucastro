"""
This tests 3-alpha screenning when we have both
a(aa,g)c12 and a(aa,p)b11
"""

import io

import pytest

import pynucastro as pyna
from pynucastro.screening import get_screening_pair_set


class Test3alphaScreening:
    @pytest.fixture(scope="class")
    def net(self):
        return pyna.network_helper(["n", "p", "he4", "be9", "b11", "c12"])

    def test_ion_screen(self, net):
        r1 = net.get_rate_by_name("a(aa,g)c12")
        r2 = net.get_rate_by_name("a(aa,p)b11")
        r3 = net.get_rate_by_nuclei(["n","p","he4","he4"],
                                    ["p","be9"])

        ans1 = pyna.Nucleus.cast_list(["he4", "he4", "he4"])
        ans2 = pyna.Nucleus.cast_list(["p", "he4", "he4"])

        assert ans1 == r1.ion_screen == r2.ion_screen
        assert ans2 == r3.ion_screen

    def test_screening_pairs(self, net):
        r1 = net.get_rate_by_name("a(aa,g)c12")
        r2 = net.get_rate_by_name("a(aa,p)b11")
        r3 = net.get_rate_by_nuclei(["n","p","he4","he4"],
                                    ["p","be9"])

        ans1 = [(pyna.Nucleus("He4"), pyna.Nucleus("He4")),
                (pyna.Nucleus("He4"), pyna.Nucleus("Be8"))]
        ans2 = [(pyna.Nucleus("He4"), pyna.Nucleus("He4")),
                (pyna.Nucleus("p"), pyna.Nucleus("Be8"))]

        assert ans1 == r1.screening_pairs == r2.screening_pairs
        assert ans2 == r3.screening_pairs

    def test_screening_pair_set(self, net):
        ans = set([(pyna.Nucleus("He4"), pyna.Nucleus("He4")),
                   (pyna.Nucleus("He4"), pyna.Nucleus("Be8")),
                   (pyna.Nucleus("He4"), pyna.Nucleus("Be9")),
                   (pyna.Nucleus("p"), pyna.Nucleus("Be8")),
                   (pyna.Nucleus("p"), pyna.Nucleus("Be9")),
                   (pyna.Nucleus("p"), pyna.Nucleus("B11"))])

        assert ans == get_screening_pair_set(net.get_rates())

    def test_python_network(self, net):
        assert isinstance(net, pyna.PythonNetwork)

        scr_str = net.screening_string()
        good = \
"""log_scor_He4_He4 = 0.0
log_scor_He4_Be9 = 0.0
log_scor_p_Be9 = 0.0
log_scor_p_B11 = 0.0
log_scor_He4_Be8 = 0.0
log_scor_p_Be8 = 0.0

if screen_func is not None:
    plasma_state = PlasmaState(T, rho, Y, Z)

    scn_fac = ScreenFactors(2, 4, 2, 4)
    log_scor_He4_He4 = screen_func(plasma_state, scn_fac)
    scn_fac = ScreenFactors(2, 4, 4, 9)
    log_scor_He4_Be9 = screen_func(plasma_state, scn_fac)
    scn_fac = ScreenFactors(1, 1, 4, 9)
    log_scor_p_Be9 = screen_func(plasma_state, scn_fac)
    scn_fac = ScreenFactors(1, 1, 5, 11)
    log_scor_p_B11 = screen_func(plasma_state, scn_fac)
    scn_fac = ScreenFactors(2, 4, 4, 8)
    log_scor_He4_Be8 = screen_func(plasma_state, scn_fac)
    scn_fac = ScreenFactors(1, 1, 4, 8)
    log_scor_p_Be8 = screen_func(plasma_state, scn_fac)
"""

        assert scr_str == good

    def test_amrex_astro_cxx_network(self, net):
        net2 = net.export_as(pyna.AmrexAstroCxxNetwork)

        output = io.StringIO()
        net2._compute_screening_factors(0, output)  # pylint: disable=protected-access
        res = output.getvalue()
        output.close()

        good = \
"""{
    constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
    static_assert(scn_fac.z1 == 2.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
    rate_eval.log_screen(k_He4_He4) = log_scor;
    if constexpr (do_T_derivatives) {
        rate_eval.dlog_screen_dT(k_He4_He4) = dlog_scor_dT;
    }
}

{
    constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 9.0_rt);
    static_assert(scn_fac.z1 == 2.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
    rate_eval.log_screen(k_He4_Be9) = log_scor;
    if constexpr (do_T_derivatives) {
        rate_eval.dlog_screen_dT(k_He4_Be9) = dlog_scor_dT;
    }
}

{
    constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 4.0_rt, 9.0_rt);
    static_assert(scn_fac.z1 == 1.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
    rate_eval.log_screen(k_p_Be9) = log_scor;
    if constexpr (do_T_derivatives) {
        rate_eval.dlog_screen_dT(k_p_Be9) = dlog_scor_dT;
    }
}

{
    constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 5.0_rt, 11.0_rt);
    static_assert(scn_fac.z1 == 1.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
    rate_eval.log_screen(k_p_B11) = log_scor;
    if constexpr (do_T_derivatives) {
        rate_eval.dlog_screen_dT(k_p_B11) = dlog_scor_dT;
    }
}

{
    constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
    static_assert(scn_fac.z1 == 2.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
    rate_eval.log_screen(k_He4_Be8) = log_scor;
    if constexpr (do_T_derivatives) {
        rate_eval.dlog_screen_dT(k_He4_Be8) = dlog_scor_dT;
    }
}

{
    constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 4.0_rt, 8.0_rt);
    static_assert(scn_fac.z1 == 1.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
    rate_eval.log_screen(k_p_Be8) = log_scor;
    if constexpr (do_T_derivatives) {
        rate_eval.dlog_screen_dT(k_p_Be8) = dlog_scor_dT;
    }
}

"""

        assert res == good

    def test_simple_cxx_network(self, net):
        net3 = net.export_as(pyna.SimpleCxxNetwork)

        output = io.StringIO()
        net3._compute_screening_factors(0, output)  # pylint: disable=protected-access
        res = output.getvalue()
        output.close()

        good = \
"""{
    auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_He4_He4) = log_scor;
}

{
    auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 9.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_He4_Be9) = log_scor;
}

{
    auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 4.0_rt, 9.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_p_Be9) = log_scor;
}

{
    auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 5.0_rt, 11.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_p_B11) = log_scor;
}

{
    auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_He4_Be8) = log_scor;
}

{
    auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 4.0_rt, 8.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_p_Be8) = log_scor;
}

"""

        assert res == good
