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

    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20) *= scor;
    rate_eval.screened_rates(k_C12_C12_to_n_Mg23) *= scor;
    rate_eval.screened_rates(k_C12_C12_to_p_Na23) *= scor;


    {
        auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);
        actual_screen(pstate, scn_fac, scor);
    }

    rate_eval.screened_rates(k_He4_C12_to_O16) *= scor;

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
        -screened_rates(k_n_to_p_weak_wc12)*Y(N) +
        0.5*screened_rates(k_C12_C12_to_n_Mg23)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(H1) =
        0.5*screened_rates(k_C12_C12_to_p_Na23)*std::pow(Y(C12), 2)*state.rho +
        screened_rates(k_n_to_p_weak_wc12)*Y(N);

    ydot_nuc(He4) =
        0.5*screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho;

    ydot_nuc(C12) =
        -screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_C12_C12_to_p_Na23)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_C12_C12_to_n_Mg23)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(O16) =
        screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho;

    ydot_nuc(Ne20) =
        0.5*screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(Na23) =
        0.5*screened_rates(k_C12_C12_to_p_Na23)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(Mg23) =
        0.5*screened_rates(k_C12_C12_to_n_Mg23)*std::pow(Y(C12), 2)*state.rho;

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

    scratch = -screened_rates(k_n_to_p_weak_wc12);
    jac.set(N, N, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_n_Mg23)*Y(C12)*state.rho;
    jac.set(N, C12, scratch);

    scratch = screened_rates(k_n_to_p_weak_wc12);
    jac.set(H1, N, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_n_Mg23)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_n_Mg23)*Y(C12)*state.rho;
    jac.set(Mg23, C12, scratch);


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
