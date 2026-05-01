#ifndef APPROXIMATE_RATES_H
#define APPROXIMATE_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

using namespace Rates;
using namespace Species;

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Mg24_He4_to_Si28_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Al27_to_Si28_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Al27_to_He4_Mg24_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Mg24_to_Si28_reaclib);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Mg24_to_p_Al27_reaclib);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28_reaclib);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27_reaclib);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Si28_to_Mg24_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Al27_to_Si28_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Al27_to_He4_Mg24_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ga = rate_eval.screened_rates(k_Si28_to_He4_Mg24_reaclib);
    amrex::Real r_gp = rate_eval.screened_rates(k_Si28_to_p_Al27_reaclib);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24_reaclib);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27_reaclib);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Si28_He4_to_S32_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_P31_to_S32_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_P31_to_He4_Si28_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Si28_to_S32_reaclib);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Si28_to_p_P31_reaclib);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32_reaclib);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31_reaclib);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_P31_to_S32_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_to_Si28_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_P31_to_S32_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_P31_to_He4_Si28_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ga = rate_eval.screened_rates(k_S32_to_He4_Si28_reaclib);
    amrex::Real r_gp = rate_eval.screened_rates(k_S32_to_p_P31_reaclib);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28_reaclib);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_S32_to_p_P31_reaclib);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_P31_to_S32_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}



template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_approx_rates([[maybe_unused]] const tf_t& tfactors,
                  [[maybe_unused]] const amrex::Real rho,
                  [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                  [[maybe_unused]] T& rate_eval)
{

    [[maybe_unused]] amrex::Real rate{};
    [[maybe_unused]] amrex::Real drate_dT{};

    rate_Mg24_He4_to_Si28_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Mg24_He4_to_Si28_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_He4_to_Si28_approx) = drate_dT;
    }

    rate_Si28_to_Mg24_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Si28_to_Mg24_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_to_Mg24_He4_approx) = drate_dT;
    }

    rate_Si28_He4_to_S32_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Si28_He4_to_S32_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_He4_to_S32_approx) = drate_dT;
    }

    rate_S32_to_Si28_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_S32_to_Si28_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_to_Si28_He4_approx) = drate_dT;
    }


}

#endif
