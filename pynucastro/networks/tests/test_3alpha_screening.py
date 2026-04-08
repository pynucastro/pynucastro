"""
This tests 3-alpha screenning when we have both
a(aa,g)c12 and a(aa,p)b11
"""

import io

import pytest

import pynucastro as pyna


class Test3alphaScreening:
    @pytest.fixture(scope="class")
    def net(self):
        return pyna.network_helper(["p", "he4", "b11", "c12"])

    def test_python_network(self, net):
        assert isinstance(net, pyna.PythonNetwork)

        scr_str = net.screening_string()
        good = \
"""log_scor_p_B11 = 0.0
log_scor_He4_He4 = 0.0
log_scor_He4_Be8 = 0.0

if screen_func is not None:
    plasma_state = PlasmaState(T, rho, Y, Z)

    scn_fac = ScreenFactors(1, 1, 5, 11)
    log_scor_p_B11 = screen_func(plasma_state, scn_fac)
    scn_fac = ScreenFactors(2, 4, 2, 4)
    log_scor_He4_He4 = screen_func(plasma_state, scn_fac)
    scn_fac = ScreenFactors(2, 4, 4, 8)
    log_scor_He4_Be8 = screen_func(plasma_state, scn_fac)
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
    constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 5.0_rt, 11.0_rt);
    static_assert(scn_fac.z1 == 1.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
    rate_eval.log_screen(k_p_B11) = log_scor;
    if constexpr (do_T_derivatives) {
        rate_eval.dlog_screen_dT(k_p_B11) = dlog_scor_dT;
    }
}

{
    constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
    static_assert(scn_fac.z1 == 2.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
    rate_eval.log_screen(k_He4_He4) = log_scor;
    if constexpr (do_T_derivatives) {
        rate_eval.dlog_screen_dT(k_He4_He4) = dlog_scor_dT;
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
    auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 5.0_rt, 11.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_p_B11) = log_scor;
}

{
    auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_He4_He4) = log_scor;
}

{
    auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
    actual_log_screen(pstate, scn_fac, log_scor);
    rate_eval.log_screen(k_He4_Be8) = log_scor;
}

"""

        assert res == good
