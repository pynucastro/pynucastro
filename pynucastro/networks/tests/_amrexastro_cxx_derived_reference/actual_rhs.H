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

    for (int n = 1; n <= NumSpec; ++n) {
        enuc += dydt(n) * network::mion(n);
    }

    enuc *= C::enuc_conv2;
}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state, T& rate_eval) {


    // create molar fractions

    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

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

    tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);


#ifdef SCREENING
    // Evaluate screening factors

    amrex::Real ratraw, dratraw_dT;
    amrex::Real scor, dscor_dt;
    [[maybe_unused]] amrex::Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_Ni56);
    rate_eval.screened_rates(k_He4_Fe52_to_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_p_Co55);
    rate_eval.screened_rates(k_He4_Fe52_to_p_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Co55_to_Ni56);
    rate_eval.screened_rates(k_p_Co55_to_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56);
        rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived);
    rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived);
        rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }

#endif

    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0_rt;


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


    ydot_nuc(H1) = 0.0_rt;

    ydot_nuc(He4) = 0.0_rt;

    ydot_nuc(Fe52) = 0.0_rt;

    ydot_nuc(Co55) = 0.0_rt;

    ydot_nuc(Ni56) = 0.0_rt;

    enuc_weak = rate_eval.enuc_weak;
}
#endif


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(H1) =
        (-screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_Ni56_to_p_Co55_derived)*Y(Ni56)) +
        (screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(H1)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        (-screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(H1)*state.rho);

    ydot_nuc(Fe52) =
        (-screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        (-screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(H1)*state.rho);

    ydot_nuc(Co55) =
        (-screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_Ni56_to_p_Co55_derived)*Y(Ni56)) +
        (screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(H1)*state.rho);

    ydot_nuc(Ni56) =
        (screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        (screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*Y(H1)*state.rho + -screened_rates(k_Ni56_to_p_Co55_derived)*Y(Ni56));

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, amrex::Array1D<amrex::Real, 1, neqs>& ydot)
{
    for (int i = 1; i <= neqs; ++i) {
        ydot(i) = 0.0_rt;
    }


    // Set molar abundances
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;

    evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);

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
    amrex::Real pair{}, phot{}, plas{}, brem{};

    neutrino_cooling<do_derivatives>(state.T, state.rho, state.abar, state.zbar,
                                     sneut, dsneutdt, dsneutdd, dsnuda, dsnudz,
                                     pair, phot, plas, brem);
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

    scratch = -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho - screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = screened_rates(k_He4_Fe52_to_p_Co55)*Y(He4)*state.rho;
    jac.set(H1, Fe52, scratch);

    scratch = -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(H1)*state.rho - screened_rates(k_p_Co55_to_Ni56)*Y(H1)*state.rho;
    jac.set(H1, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_p_Co55_derived);
    jac.set(H1, Ni56, scratch);

    scratch = screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = -screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55)*Y(He4)*state.rho;
    jac.set(He4, Fe52, scratch);

    scratch = screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(H1)*state.rho;
    jac.set(He4, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_He4_Fe52_derived);
    jac.set(He4, Ni56, scratch);

    scratch = screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho;
    jac.set(Fe52, H1, scratch);

    scratch = -screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    scratch = -screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55)*Y(He4)*state.rho;
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(H1)*state.rho;
    jac.set(Fe52, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_He4_Fe52_derived);
    jac.set(Fe52, Ni56, scratch);

    scratch = -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho - screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*state.rho;
    jac.set(Co55, H1, scratch);

    scratch = screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*state.rho;
    jac.set(Co55, He4, scratch);

    scratch = screened_rates(k_He4_Fe52_to_p_Co55)*Y(He4)*state.rho;
    jac.set(Co55, Fe52, scratch);

    scratch = -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(H1)*state.rho - screened_rates(k_p_Co55_to_Ni56)*Y(H1)*state.rho;
    jac.set(Co55, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_p_Co55_derived);
    jac.set(Co55, Ni56, scratch);

    scratch = screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*state.rho;
    jac.set(Ni56, H1, scratch);

    scratch = screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho;
    jac.set(Ni56, He4, scratch);

    scratch = screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    scratch = screened_rates(k_p_Co55_to_Ni56)*Y(H1)*state.rho;
    jac.set(Ni56, Co55, scratch);

    scratch = -screened_rates(k_Ni56_to_He4_Fe52_derived) - screened_rates(k_Ni56_to_p_Co55_derived);
    jac.set(Ni56, Ni56, scratch);


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


    jac.zero();

    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives = 1;

    evaluate_rates<do_T_derivatives, rate_derivs_t>(state, rate_eval);

    // Species Jacobian elements with respect to other species

    jac_nuc(state, jac, Y, rate_eval.screened_rates);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> amrex::Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc,j));
    }

    // Account for the thermal neutrino losses

    amrex::Real dsneutdt{};
#ifdef NEUTRINOS
    amrex::Real sneut, dsneutdd, dsnuda, dsnudz;
    amrex::Real pair{}, phot{}, plas{}, brem{};
    constexpr int do_derivatives{1};
    neutrino_cooling<do_derivatives>(state.T, state.rho, state.abar, state.zbar,
                                     sneut, dsneutdt, dsneutdd, dsnuda, dsnudz,
                                     pair, phot, plas, brem);

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
