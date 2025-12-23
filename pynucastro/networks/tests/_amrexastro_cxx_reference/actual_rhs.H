#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <burn_type.H>
#include <jacobian_utilities.H>
#ifdef SCREENING
#include <screen.H>
#endif
#include <microphysics_autodiff.H>
#ifdef NEUTRINOS
#include <neutrino.H>
#endif
#include <reaclib_rates.H>
#include <table_rates.H>

using namespace amrex::literals;
using namespace ArrayUtil;

using namespace Species;
using namespace Rates;

using namespace rate_tables;


template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate(T const& dydt, amrex::Real& enuc)
{

    // Computes the instantaneous energy generation rate (from the nuclei)

    // This is basically e = m c**2

    enuc = 0.0_rt;

    amrex::constexpr_for<1, NumSpec+1>([&] (auto n)
    {
        enuc += dydt(n) * network::mion<n>();
    });

    enuc *= C::enuc_conv2;
}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state,
                    const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                    T& rate_eval) {

    // here Y is consistent with state.xn

    [[maybe_unused]] amrex::Real rhoy = state.rho * state.y_e;

    // Calculate Reaclib rates

    using number_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    number_t temp = state.T;
    if constexpr (do_T_derivatives) {
        // seed the dual number for temperature before calculating anything with it
        autodiff::seed(temp);
    }

#ifdef SCREENING
    plasma_state_t<number_t> pstate{};
    fill_plasma_state(pstate, temp, state.rho, Y);
#endif

    const tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate Derived Rates. This should go last but before approx rates.
    fill_derived_rates<do_T_derivatives, T>(tfactors, rate_eval);


#ifdef SCREENING
    // Evaluate screening factors

    amrex::Real ratraw, dratraw_dT;
    amrex::Real scor, dscor_dt;
    [[maybe_unused]] amrex::Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_n_Mg23_reaclib);
    rate_eval.screened_rates(k_C12_C12_to_n_Mg23_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib);
    rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_C12_to_O16_reaclib);
    rate_eval.screened_rates(k_He4_C12_to_O16_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_O16_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

#endif

    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0_rt;

    amrex::Real log_temp = std::log10(state.T);
    amrex::Real log_rhoy = std::log10(rhoy);

    tabular_evaluate(j_Na23_Ne23_meta, j_Na23_Ne23_rhoy, j_Na23_Ne23_temp, j_Na23_Ne23_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Na23_to_Ne23_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na23_to_Ne23_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Na23) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne23_Na23_meta, j_Ne23_Na23_rhoy, j_Ne23_Na23_temp, j_Ne23_Na23_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne23_to_Na23_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne23_to_Na23_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Ne23) * (edot_nu + edot_gamma);


}


// backwards-compatible wrapper
template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state,
                    T& rate_eval) {

    // create molar fractions
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    evaluate_rates<do_T_derivatives>(state, Y, rate_eval);
}


#ifdef NSE_NET
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_ydot_weak(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             amrex::Real& enuc_weak,
             [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {
    ///
    /// Calculate Ydots contribute only from weak reactions.
    /// This is used to calculate dyedt and energy generation from
    /// weak reactions for self-consistent NSE
    ///


    // initialize ydot_nuc to 0

    for (int i = 1; i <= neqs; ++i) {
        ydot_nuc(i) = 0.0_rt;
    }

    rate_t rate_eval;

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;
    [[maybe_unused]] amrex::Real rhoy = state.rho * state.y_e;

    rate_eval.enuc_weak = 0.0_rt;

    // Calculate tabular rates and get ydot_weak

    amrex::Real log_temp = std::log10(state.T);
    amrex::Real log_rhoy = std::log10(rhoy);

    tabular_evaluate(j_Na23_Ne23_meta, j_Na23_Ne23_rhoy, j_Na23_Ne23_temp, j_Na23_Ne23_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Na23_to_Ne23_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Na23) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne23_Na23_meta, j_Ne23_Na23_rhoy, j_Ne23_Na23_temp, j_Ne23_Na23_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne23_to_Na23_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Ne23) * (edot_nu + edot_gamma);

    auto screened_rates = rate_eval.screened_rates;

    ydot_nuc(N) =
        -screened_rates(k_n_to_p_reaclib)*Y(N);

    ydot_nuc(H1) =
        screened_rates(k_n_to_p_reaclib)*Y(N);

    ydot_nuc(He4) = 0.0_rt;

    ydot_nuc(C12) = 0.0_rt;

    ydot_nuc(O16) = 0.0_rt;

    ydot_nuc(Ne20) = 0.0_rt;

    ydot_nuc(Ne23) =
        (-screened_rates(k_Ne23_to_Na23_weaktab)*Y(Ne23) + screened_rates(k_Na23_to_Ne23_weaktab)*Y(Na23));

    ydot_nuc(Na23) =
        (screened_rates(k_Ne23_to_Na23_weaktab)*Y(Ne23) + -screened_rates(k_Na23_to_Ne23_weaktab)*Y(Na23));

    ydot_nuc(Mg23) = 0.0_rt;

    enuc_weak = rate_eval.enuc_weak;
}
#endif


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        -screened_rates(k_n_to_p_reaclib)*Y(N) +
        0.5*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho;

    ydot_nuc(H1) =
        0.5*screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        screened_rates(k_n_to_p_reaclib)*Y(N);

    ydot_nuc(He4) =
        0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho;

    ydot_nuc(C12) =
        -screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        -screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_C12_C12_to_n_Mg23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho;

    ydot_nuc(O16) =
        screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho;

    ydot_nuc(Ne20) =
        0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho;

    ydot_nuc(Ne23) =
        (-screened_rates(k_Ne23_to_Na23_weaktab)*Y(Ne23) + screened_rates(k_Na23_to_Ne23_weaktab)*Y(Na23));

    ydot_nuc(Na23) =
        0.5*screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (screened_rates(k_Ne23_to_Na23_weaktab)*Y(Ne23) + -screened_rates(k_Na23_to_Ne23_weaktab)*Y(Na23));

    ydot_nuc(Mg23) =
        0.5*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho;

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (const burn_t& state, amrex::Array1D<amrex::Real, 1, neqs>& ydot)
{

    // Set molar abundances
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;

    evaluate_rates<do_T_derivatives, rate_t>(state, Y, rate_eval);

    rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

    // ion binding energy contributions

    amrex::Real enuc;
    ener_gener_rate(ydot, enuc);

    // include any weak rate neutrino losses
    enuc += rate_eval.enuc_weak;

    // Get the thermal neutrino losses

    amrex::Real sneut{};
#ifdef NEUTRINOS
    constexpr int do_derivatives{0};
    amrex::Real dsneutdt{}, dsneutdd{}, dsnuda{}, dsnudz{};

    neutrino_cooling<do_derivatives>(state.T, state.rho, state.abar, state.zbar,
                                     sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);
#endif

    // Append the energy equation (this is erg/g/s)

    ydot(net_ienuc) = enuc - sneut;

}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jac_nuc(const burn_t& state,
             MatrixType& jac,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates)
{

    amrex::Real scratch;

    scratch = -screened_rates(k_n_to_p_reaclib);
    jac.set(N, N, scratch);

    jac.set(N, H1, 0.0);

    jac.set(N, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho;
    jac.set(N, C12, scratch);

    jac.set(N, O16, 0.0);

    jac.set(N, Ne20, 0.0);

    jac.set(N, Ne23, 0.0);

    jac.set(N, Na23, 0.0);

    jac.set(N, Mg23, 0.0);

    scratch = screened_rates(k_n_to_p_reaclib);
    jac.set(H1, N, scratch);

    jac.set(H1, H1, 0.0);

    jac.set(H1, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho;
    jac.set(H1, C12, scratch);

    jac.set(H1, O16, 0.0);

    jac.set(H1, Ne20, 0.0);

    jac.set(H1, Ne23, 0.0);

    jac.set(H1, Na23, 0.0);

    jac.set(H1, Mg23, 0.0);

    jac.set(He4, N, 0.0);

    jac.set(He4, H1, 0.0);

    scratch = -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    jac.set(He4, O16, 0.0);

    jac.set(He4, Ne20, 0.0);

    jac.set(He4, Ne23, 0.0);

    jac.set(He4, Na23, 0.0);

    jac.set(He4, Mg23, 0.0);

    jac.set(C12, N, 0.0);

    jac.set(C12, H1, 0.0);

    scratch = -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(C12, C12, scratch);

    jac.set(C12, O16, 0.0);

    jac.set(C12, Ne20, 0.0);

    jac.set(C12, Ne23, 0.0);

    jac.set(C12, Na23, 0.0);

    jac.set(C12, Mg23, 0.0);

    jac.set(O16, N, 0.0);

    jac.set(O16, H1, 0.0);

    scratch = screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    jac.set(O16, O16, 0.0);

    jac.set(O16, Ne20, 0.0);

    jac.set(O16, Ne23, 0.0);

    jac.set(O16, Na23, 0.0);

    jac.set(O16, Mg23, 0.0);

    jac.set(Ne20, N, 0.0);

    jac.set(Ne20, H1, 0.0);

    jac.set(Ne20, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    jac.set(Ne20, O16, 0.0);

    jac.set(Ne20, Ne20, 0.0);

    jac.set(Ne20, Ne23, 0.0);

    jac.set(Ne20, Na23, 0.0);

    jac.set(Ne20, Mg23, 0.0);

    jac.set(Ne23, N, 0.0);

    jac.set(Ne23, H1, 0.0);

    jac.set(Ne23, He4, 0.0);

    jac.set(Ne23, C12, 0.0);

    jac.set(Ne23, O16, 0.0);

    jac.set(Ne23, Ne20, 0.0);

    scratch = -screened_rates(k_Ne23_to_Na23_weaktab);
    jac.set(Ne23, Ne23, scratch);

    scratch = screened_rates(k_Na23_to_Ne23_weaktab);
    jac.set(Ne23, Na23, scratch);

    jac.set(Ne23, Mg23, 0.0);

    jac.set(Na23, N, 0.0);

    jac.set(Na23, H1, 0.0);

    jac.set(Na23, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    jac.set(Na23, O16, 0.0);

    jac.set(Na23, Ne20, 0.0);

    scratch = screened_rates(k_Ne23_to_Na23_weaktab);
    jac.set(Na23, Ne23, scratch);

    scratch = -screened_rates(k_Na23_to_Ne23_weaktab);
    jac.set(Na23, Na23, scratch);

    jac.set(Na23, Mg23, 0.0);

    jac.set(Mg23, N, 0.0);

    jac.set(Mg23, H1, 0.0);

    jac.set(Mg23, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho;
    jac.set(Mg23, C12, scratch);

    jac.set(Mg23, O16, 0.0);

    jac.set(Mg23, Ne20, 0.0);

    jac.set(Mg23, Ne23, 0.0);

    jac.set(Mg23, Na23, 0.0);

    jac.set(Mg23, Mg23, 0.0);


}



template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac(const burn_t& state, MatrixType& jac)
{

    // Set molar abundances
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }


    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives = 1;

    evaluate_rates<do_T_derivatives, rate_derivs_t>(state, Y, rate_eval);

    // Species Jacobian elements with respect to other species

    jac_nuc(state, jac, Y, rate_eval.screened_rates);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> amrex::Real { return jac.get(i, j); };
        amrex::Real tmp{};
        ener_gener_rate(jac_slice_2, tmp);
        jac(net_ienuc,j) = tmp;
    }

    // Account for the thermal neutrino losses

    amrex::Real dsneutdt{};
#ifdef NEUTRINOS
    amrex::Real sneut, dsneutdd, dsnuda, dsnudz;
    constexpr int do_derivatives{1};
    neutrino_cooling<do_derivatives>(state.T, state.rho, state.abar, state.zbar,
                                     sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

    for (int j = 1; j <= NumSpec; ++j) {
       amrex::Real b1 = (-state.abar * state.abar * dsnuda + (zion[j-1] - state.zbar) * state.abar * dsnudz);
       jac.add(net_ienuc, j, -b1);
    }
#endif

    // Evaluate the Jacobian elements with respect to energy by
    // calling the RHS using d(rate) / dT and then transform them
    // to our energy integration variable.

    amrex::Array1D<amrex::Real, 1, neqs>  yderivs;

    rhs_nuc(state, yderivs, Y, rate_eval.dscreened_rates_dT);

    for (int k = 1; k <= NumSpec; k++) {
        jac.set(k, net_ienuc, temperature_to_energy_jacobian(state, yderivs(k)));
    }


    // finally, d(de/dt)/de

    amrex::Real jac_e_T;
    ener_gener_rate(yderivs, jac_e_T);
    jac_e_T -= dsneutdt;
    jac.set(net_ienuc, net_ienuc, temperature_to_energy_jacobian(state, jac_e_T));

}


AMREX_INLINE
void actual_rhs_init () {

    init_tabular();

}


#endif
