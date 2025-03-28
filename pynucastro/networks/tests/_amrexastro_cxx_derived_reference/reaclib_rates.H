#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <partition_functions.H>

using namespace Rates;
using namespace Species;

struct rate_t {
    amrex::Array1D<amrex::Real, 1, NumRates>  screened_rates;
    amrex::Real enuc_weak;
};

struct rate_derivs_t {
    amrex::Array1D<amrex::Real, 1, NumRates>  screened_rates;
    amrex::Array1D<amrex::Real, 1, NumRates>  dscreened_rates_dT;
    amrex::Real enuc_weak;
};


template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe52_to_Ni56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + He4 --> Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  66.6417 + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                         + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                  + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co55_to_Ni56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + p --> Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.3736 + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                         + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                  + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe52_to_p_Co55(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + He4 --> p + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.2207 + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                         + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                  + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_to_He4_Fe52_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 --> He4 + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  91.62258922626623 + -92.80109932767355 * tfactors.T9i + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                         + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  92.80109932767355 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                  + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors, pf_cache, Ni56_pf, dNi56_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors, pf_cache, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = He4_pf * Fe52_pf;
    amrex::Real z_p = Ni56_pf;

    amrex::Real dz_r_dT = Fe52_pf * dHe4_pf_dT + He4_pf * dFe52_pf_dT;
    amrex::Real dz_p_dT = dNi56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_to_p_Co55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 --> p + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  63.131770606802746 + -83.14741674753726 * tfactors.T9i + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                         + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  83.14741674753726 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                  + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors, pf_cache, Ni56_pf, dNi56_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors, pf_cache, Co55_pf, dCo55_pf_dT);

    amrex::Real z_r = p_pf * Co55_pf;
    amrex::Real z_p = Ni56_pf;

    amrex::Real dz_r_dT = Co55_pf * dp_pf_dT + p_pf * dCo55_pf_dT;
    amrex::Real dz_p_dT = dNi56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co55_to_He4_Fe52_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 + p --> He4 + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.443418619463486 + -9.653647766581932 * tfactors.T9i + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                         + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  9.653647766581932 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                  + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors, pf_cache, Co55_pf, dCo55_pf_dT);

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors, pf_cache, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = He4_pf * Fe52_pf;
    amrex::Real z_p = p_pf * Co55_pf;

    amrex::Real dz_r_dT = Fe52_pf * dHe4_pf_dT + He4_pf * dFe52_pf_dT;
    amrex::Real dz_p_dT = Co55_pf * dp_pf_dT + p_pf * dCo55_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}



template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::Real rate;
    amrex::Real drate_dT;

    part_fun::pf_cache_t pf_cache{};

    rate_He4_Fe52_to_Ni56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe52_to_Ni56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56) = drate_dT;

    }
    rate_p_Co55_to_Ni56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co55_to_Ni56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56) = drate_dT;

    }
    rate_He4_Fe52_to_p_Co55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe52_to_p_Co55) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55) = drate_dT;

    }
    rate_Ni56_to_He4_Fe52_derived<do_T_derivatives>(tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ni56_to_He4_Fe52_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_derived) = drate_dT;

    }
    rate_Ni56_to_p_Co55_derived<do_T_derivatives>(tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ni56_to_p_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_derived) = drate_dT;

    }
    rate_p_Co55_to_He4_Fe52_derived<do_T_derivatives>(tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived) = drate_dT;

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


}

#endif
