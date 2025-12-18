#ifndef actual_rhs_H
#define actual_rhs_H

#include <amrex_bridge.H>

#include <actual_network.H>
#include <burn_type.H>
#include <screen.H>
#include <reaclib_rates.H>

using namespace Species;
using namespace Rates;


template<class T>
inline
void ener_gener_rate(T const& dydt, Real& enuc)
{

    // Computes the instantaneous energy generation rate (from the nuclei)

    // This is basically e = m c**2

    enuc = 0.0_rt;

    for (int n = 1; n <= NumSpec; ++n) {
        enuc += dydt(n) * network::mion(n);
    }

    enuc *= C::enuc_conv2;
}


template <int do_T_derivatives, typename T>
inline
void evaluate_rates(const burn_t& state, T& rate_eval) {

    // create molar fractions

    Array1D<Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    // Calculate Reaclib rates

    plasma_state_t pstate{};
    fill_plasma_state(pstate, state.T, state.rho, Y);

    tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);

#ifndef DISABLE_SCREENING

    // Evaluating screening factors

    Real scor;
    [[maybe_unused]] Real scor2;


    {
        auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);
        actual_screen(pstate, scn_fac, scor);
    }

    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib) *= scor;
    rate_eval.screened_rates(k_C12_C12_to_n_Mg23_reaclib) *= scor;
    rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib) *= scor;


    {
        auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);
        actual_screen(pstate, scn_fac, scor);
    }

    rate_eval.screened_rates(k_He4_C12_to_O16_reaclib) *= scor;


    {
        auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);
        actual_screen(pstate, scn_fac, scor);
    }

    rate_eval.screened_rates(k_He4_Ne20_to_Mg24_reaclib) *= scor;

#endif

    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

}

inline
void rhs_nuc(const burn_t& state,
             Array1D<Real, 1, NumSpec>& ydot_nuc,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        -screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(Mg23)*Y(N)*state.rho +
        -screened_rates(k_n_to_p_reaclib)*Y(N) +
        0.5*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(H1) =
        0.5*screened_rates(k_C12_C12_to_p_Na23_reaclib)*std::pow(Y(C12), 2)*state.rho +
        screened_rates(k_n_to_p_reaclib)*Y(N);

    ydot_nuc(He4) =
        0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho;

    ydot_nuc(C12) =
        -screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_C12_C12_to_p_Na23_reaclib)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_C12_C12_to_n_Mg23_reaclib)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(O16) =
        screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho;

    ydot_nuc(Ne20) =
        0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho;

    ydot_nuc(Na23) =
        0.5*screened_rates(k_C12_C12_to_p_Na23_reaclib)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(Mg23) =
        -screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(Mg23)*Y(N)*state.rho +
        0.5*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(Mg24) =
        screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho +
        screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(Mg23)*Y(N)*state.rho;

}


inline
void actual_rhs (burn_t& state, Array1D<Real, 1, NumSpec>& ydot)
{

    compute_ye(state);

    for (int i = 1; i <= NumSpec; ++i) {
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

}


template<class MatrixType>
inline
void jac_nuc(const burn_t& state,
             MatrixType& jac,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates)
{

    Real scratch;

    scratch = -screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(Mg23)*state.rho - screened_rates(k_n_to_p_reaclib);
    jac.set(N, N, scratch);

    jac.set(N, H1, 0.0);

    jac.set(N, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho;
    jac.set(N, C12, scratch);

    jac.set(N, O16, 0.0);

    jac.set(N, Ne20, 0.0);

    jac.set(N, Na23, 0.0);

    scratch = -screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(N)*state.rho;
    jac.set(N, Mg23, scratch);

    jac.set(N, Mg24, 0.0);

    scratch = screened_rates(k_n_to_p_reaclib);
    jac.set(H1, N, scratch);

    jac.set(H1, H1, 0.0);

    jac.set(H1, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho;
    jac.set(H1, C12, scratch);

    jac.set(H1, O16, 0.0);

    jac.set(H1, Ne20, 0.0);

    jac.set(H1, Na23, 0.0);

    jac.set(H1, Mg23, 0.0);

    jac.set(H1, Mg24, 0.0);

    jac.set(He4, N, 0.0);

    jac.set(He4, H1, 0.0);

    scratch = -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    jac.set(He4, O16, 0.0);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Ne20, scratch);

    jac.set(He4, Na23, 0.0);

    jac.set(He4, Mg23, 0.0);

    jac.set(He4, Mg24, 0.0);

    jac.set(C12, N, 0.0);

    jac.set(C12, H1, 0.0);

    scratch = -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(C12, C12, scratch);

    jac.set(C12, O16, 0.0);

    jac.set(C12, Ne20, 0.0);

    jac.set(C12, Na23, 0.0);

    jac.set(C12, Mg23, 0.0);

    jac.set(C12, Mg24, 0.0);

    jac.set(O16, N, 0.0);

    jac.set(O16, H1, 0.0);

    scratch = screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    jac.set(O16, O16, 0.0);

    jac.set(O16, Ne20, 0.0);

    jac.set(O16, Na23, 0.0);

    jac.set(O16, Mg23, 0.0);

    jac.set(O16, Mg24, 0.0);

    jac.set(Ne20, N, 0.0);

    jac.set(Ne20, H1, 0.0);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    jac.set(Ne20, O16, 0.0);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho;
    jac.set(Ne20, Ne20, scratch);

    jac.set(Ne20, Na23, 0.0);

    jac.set(Ne20, Mg23, 0.0);

    jac.set(Ne20, Mg24, 0.0);

    jac.set(Na23, N, 0.0);

    jac.set(Na23, H1, 0.0);

    jac.set(Na23, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    jac.set(Na23, O16, 0.0);

    jac.set(Na23, Ne20, 0.0);

    jac.set(Na23, Na23, 0.0);

    jac.set(Na23, Mg23, 0.0);

    jac.set(Na23, Mg24, 0.0);

    scratch = -screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(Mg23)*state.rho;
    jac.set(Mg23, N, scratch);

    jac.set(Mg23, H1, 0.0);

    jac.set(Mg23, He4, 0.0);

    scratch = 1.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho;
    jac.set(Mg23, C12, scratch);

    jac.set(Mg23, O16, 0.0);

    jac.set(Mg23, Ne20, 0.0);

    jac.set(Mg23, Na23, 0.0);

    scratch = -screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(N)*state.rho;
    jac.set(Mg23, Mg23, scratch);

    jac.set(Mg23, Mg24, 0.0);

    scratch = screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(Mg23)*state.rho;
    jac.set(Mg24, N, scratch);

    jac.set(Mg24, H1, 0.0);

    scratch = screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    jac.set(Mg24, C12, 0.0);

    jac.set(Mg24, O16, 0.0);

    scratch = screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    jac.set(Mg24, Na23, 0.0);

    scratch = screened_rates(k_n_Mg23_to_Mg24_reaclib)*Y(N)*state.rho;
    jac.set(Mg24, Mg23, scratch);

    jac.set(Mg24, Mg24, 0.0);


}



template<class MatrixType>
inline
void actual_jac(burn_t& state, MatrixType& jac)
{

    compute_ye(state);

    // Set molar abundances
    Array1D<Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }


    jac.zero();

    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives = 0;
    evaluate_rates<do_T_derivatives, rate_derivs_t>(state, rate_eval);

    // Species Jacobian elements with respect to other species

    jac_nuc(state, jac, Y, rate_eval.screened_rates);

}

#endif
