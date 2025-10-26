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
"""plasma_state = PlasmaState(T, rho, Y, Z)

scn_fac = ScreenFactors(1, 1, 5, 11)
scor = screen_func(plasma_state, scn_fac)
rate_eval.p_B11__C12 *= scor
rate_eval.p_B11__He4_He4_He4 *= scor

scn_fac = ScreenFactors(2, 4, 2, 4)
scor = screen_func(plasma_state, scn_fac)
scn_fac2 = ScreenFactors(2, 4, 4, 8)
scor2 = screen_func(plasma_state, scn_fac2)
rate_eval.He4_He4_He4__C12 *= scor * scor2
rate_eval.He4_He4_He4__p_B11__derived *= scor * scor2
"""

        assert scr_str == good

    def test_amrex_astro_cxx_network(self, net):
        net2 = net.export_as(pyna.AmrexAstroCxxNetwork)

        output = io.StringIO()
        net2._compute_screening_factors(0, output)  # pylint: disable=protected-access
        res = output.getvalue()
        output.close()

        good = \
"""
{
    constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 5.0_rt, 11.0_rt);
    static_assert(scn_fac.z1 == 1.0_rt);
    actual_screen(pstate, scn_fac, scor, dscor_dt);
}

ratraw = rate_eval.screened_rates(k_p_B11_to_C12);
rate_eval.screened_rates(k_p_B11_to_C12) *= scor;
if constexpr (std::is_same_v<T, rate_derivs_t>) {
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_B11_to_C12);
    rate_eval.dscreened_rates_dT(k_p_B11_to_C12) = ratraw * dscor_dt + dratraw_dT * scor;
}

ratraw = rate_eval.screened_rates(k_p_B11_to_He4_He4_He4);
rate_eval.screened_rates(k_p_B11_to_He4_He4_He4) *= scor;
if constexpr (std::is_same_v<T, rate_derivs_t>) {
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_B11_to_He4_He4_He4);
    rate_eval.dscreened_rates_dT(k_p_B11_to_He4_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
}


{
    constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
    static_assert(scn_fac.z1 == 2.0_rt);
    actual_screen(pstate, scn_fac, scor, dscor_dt);
}


{
    constexpr auto scn_fac2 = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
    static_assert(scn_fac2.z1 == 2.0_rt);
    actual_screen(pstate, scn_fac2, scor2, dscor2_dt);
}

ratraw = rate_eval.screened_rates(k_He4_He4_He4_to_C12);
rate_eval.screened_rates(k_He4_He4_He4_to_C12) *= scor * scor2;
if constexpr (std::is_same_v<T, rate_derivs_t>) {
    dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12);
    rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
}

ratraw = rate_eval.screened_rates(k_He4_He4_He4_to_p_B11_derived);
rate_eval.screened_rates(k_He4_He4_He4_to_p_B11_derived) *= scor * scor2;
if constexpr (std::is_same_v<T, rate_derivs_t>) {
    dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_p_B11_derived);
    rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_p_B11_derived) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
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
"""
{
    auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 5.0_rt, 11.0_rt);
    actual_screen(pstate, scn_fac, scor);
}

rate_eval.screened_rates(k_p_B11_to_C12) *= scor;
rate_eval.screened_rates(k_p_B11_to_He4_He4_He4) *= scor;


{
    auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
    actual_screen(pstate, scn_fac, scor);
}


{
    auto scn_fac2 = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
    actual_screen(pstate, scn_fac2, scor2);
}

rate_eval.screened_rates(k_He4_He4_He4_to_C12) *= scor * scor2;
rate_eval.screened_rates(k_He4_He4_He4_to_p_B11_derived) *= scor * scor2;

"""

        assert res == good
