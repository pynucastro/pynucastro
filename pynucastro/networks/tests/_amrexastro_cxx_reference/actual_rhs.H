#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <burn_type.H>
#include <jacobian_utilities.H>
#include <screen.H>
#include <sneut5.H>
#include <reaclib_rates.H>
#include <table_rates.H>
#include <extern_parameters.H>

using namespace amrex;
using namespace ArrayUtil;

using namespace Species;
using namespace Rates;

using namespace rate_tables;


template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate(T const& dydt, Real& enuc)
{

    // Computes the instantaneous energy generation rate (from the nuclei)

    // This is basically e = m c**2

    enuc = 0.0_rt;

    for (int n = 1; n <= NumSpec; ++n) {
        enuc += dydt(n) * network::mion(n);
    }

    enuc *= C::Legacy::enuc_conv2;
}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state, T& rate_eval) {


    // create molar fractions

    Array1D<Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    [[maybe_unused]] Real rhoy = state.rho * state.y_e;

    // Calculate Reaclib rates

    plasma_state_t pstate;
    fill_plasma_state(pstate, state.T, state.rho, Y);

    tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);



    // Evaluate screening factors

    Real ratraw, dratraw_dT;
    Real scor, dscor_dt;
    Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_c12_c12__he4_ne20);
    rate_eval.screened_rates(k_c12_c12__he4_ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_c12__he4_ne20);
        rate_eval.dscreened_rates_dT(k_c12_c12__he4_ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_c12_c12__n_mg23);
    rate_eval.screened_rates(k_c12_c12__n_mg23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_c12__n_mg23);
        rate_eval.dscreened_rates_dT(k_c12_c12__n_mg23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_c12_c12__p_na23);
    rate_eval.screened_rates(k_c12_c12__p_na23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_c12__p_na23);
        rate_eval.dscreened_rates_dT(k_c12_c12__p_na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_c12__o16);
    rate_eval.screened_rates(k_he4_c12__o16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_c12__o16);
        rate_eval.dscreened_rates_dT(k_he4_c12__o16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] Real rate, drate_dt, edot_nu, edot_gamma;

    tabular_evaluate(j_na23_ne23_meta, j_na23_ne23_rhoy, j_na23_ne23_temp, j_na23_ne23_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_na23__ne23) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_na23__ne23) = drate_dt;
    }
    rate_eval.add_energy_rate(k_na23__ne23) = edot_nu + edot_gamma;

    tabular_evaluate(j_ne23_na23_meta, j_ne23_na23_rhoy, j_ne23_na23_temp, j_ne23_na23_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_ne23__na23) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ne23__na23) = drate_dt;
    }
    rate_eval.add_energy_rate(k_ne23__na23) = edot_nu + edot_gamma;


}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             Array1D<Real, 1, neqs>& ydot_nuc,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        -screened_rates(k_n__p__weak__wc12)*Y(N) +
        0.5*screened_rates(k_c12_c12__n_mg23)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(H1) =
        0.5*screened_rates(k_c12_c12__p_na23)*std::pow(Y(C12), 2)*state.rho +
        screened_rates(k_n__p__weak__wc12)*Y(N);

    ydot_nuc(He4) =
        0.5*screened_rates(k_c12_c12__he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho;

    ydot_nuc(C12) =
        -screened_rates(k_c12_c12__he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_c12_c12__p_na23)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_c12_c12__n_mg23)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(O16) =
        screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho;

    ydot_nuc(Ne20) =
        0.5*screened_rates(k_c12_c12__he4_ne20)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(Ne23) =
        (-screened_rates(k_ne23__na23)*Y(Ne23) + screened_rates(k_na23__ne23)*Y(Na23));

    ydot_nuc(Na23) =
        0.5*screened_rates(k_c12_c12__p_na23)*std::pow(Y(C12), 2)*state.rho +
        (screened_rates(k_ne23__na23)*Y(Ne23) + -screened_rates(k_na23__ne23)*Y(Na23));

    ydot_nuc(Mg23) =
        0.5*screened_rates(k_c12_c12__n_mg23)*std::pow(Y(C12), 2)*state.rho;

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    for (int i = 1; i <= neqs; ++i) {
        ydot(i) = 0.0_rt;
    }


    // Set molar abundances
    Array1D<Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;
    evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);

    rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

    // ion binding energy contributions

    Real enuc;
    ener_gener_rate(ydot, enuc);

    // include reaction neutrino losses (non-thermal) and gamma heating
    enuc += C::Legacy::n_A * Y(Na23) * rate_eval.add_energy_rate(k_na23__ne23);
    enuc += C::Legacy::n_A * Y(Ne23) * rate_eval.add_energy_rate(k_ne23__na23);

    // Get the thermal neutrino losses

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;

    sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    // Append the energy equation (this is erg/g/s)

    ydot(net_ienuc) = enuc - sneut;

}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jac_nuc(const burn_t& state,
             MatrixType& jac,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates)
{

    Real scratch;

    scratch = -screened_rates(k_n__p__weak__wc12);
    jac.set(N, N, scratch);

    scratch = 1.0*screened_rates(k_c12_c12__n_mg23)*Y(C12)*state.rho;
    jac.set(N, C12, scratch);

    scratch = screened_rates(k_n__p__weak__wc12);
    jac.set(H1, N, scratch);

    scratch = 1.0*screened_rates(k_c12_c12__p_na23)*Y(C12)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(C12)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12__he4_ne20)*Y(C12)*state.rho - screened_rates(k_he4_c12__o16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(C12)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_c12_c12__he4_ne20)*Y(C12)*state.rho - 2.0*screened_rates(k_c12_c12__n_mg23)*Y(C12)*state.rho - 2.0*screened_rates(k_c12_c12__p_na23)*Y(C12)*state.rho - screened_rates(k_he4_c12__o16)*Y(He4)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = screened_rates(k_he4_c12__o16)*Y(C12)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_he4_c12__o16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = 1.0*screened_rates(k_c12_c12__he4_ne20)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = -screened_rates(k_ne23__na23);
    jac.set(Ne23, Ne23, scratch);

    scratch = screened_rates(k_na23__ne23);
    jac.set(Ne23, Na23, scratch);

    scratch = 1.0*screened_rates(k_c12_c12__p_na23)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    scratch = screened_rates(k_ne23__na23);
    jac.set(Na23, Ne23, scratch);

    scratch = -screened_rates(k_na23__ne23);
    jac.set(Na23, Na23, scratch);

    scratch = 1.0*screened_rates(k_c12_c12__n_mg23)*Y(C12)*state.rho;
    jac.set(Mg23, C12, scratch);


}



template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac(const burn_t& state, MatrixType& jac)
{

    // Set molar abundances
    Array1D<Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }


    jac.zero();

    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives = 1;
    evaluate_rates<do_T_derivatives, rate_derivs_t>(state, rate_eval);

    // Species Jacobian elements with respect to other species

    jac_nuc(state, jac, Y, rate_eval.screened_rates);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc,j));
    }

    // Account for the thermal neutrino losses

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    for (int j = 1; j <= NumSpec; ++j) {
       Real b1 = (-state.abar * state.abar * snuda + (zion[j-1] - state.zbar) * state.abar * snudz);
       jac.add(net_ienuc, j, -b1);
    }


    // Evaluate the Jacobian elements with respect to energy by
    // calling the RHS using d(rate) / dT and then transform them
    // to our energy integration variable.

    Array1D<Real, 1, neqs>  yderivs;

    rhs_nuc(state, yderivs, Y, rate_eval.dscreened_rates_dT);

    for (int k = 1; k <= NumSpec; k++) {
        jac.set(k, net_ienuc, temperature_to_energy_jacobian(state, yderivs(k)));
    }


    // finally, d(de/dt)/de

    Real jac_e_T;
    ener_gener_rate(yderivs, jac_e_T);
    jac_e_T -= dsneutdt;
    jac.set(net_ienuc, net_ienuc, temperature_to_energy_jacobian(state, jac_e_T));

}


AMREX_INLINE
void actual_rhs_init () {

    init_tabular();

}


#endif
