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
    Array1D<Real, 1, NumRates>  screened_rates;
    Real enuc_weak;
};

struct rate_derivs_t {
    Array1D<Real, 1, NumRates>  screened_rates;
    Array1D<Real, 1, NumRates>  dscreened_rates_dT;
    Real enuc_weak;
};


template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_fe52_to_ni56(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // fe52 + he4 --> ni56

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

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
void rate_p_co55_to_ni56(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // co55 + p --> ni56

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

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
void rate_he4_fe52_to_p_co55(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // fe52 + he4 --> p + co55

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

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
void rate_ni56_to_he4_fe52_derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ni56 --> he4 + fe52

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  91.62258922810439 + -92.801099329237 * tfactors.T9i + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                         + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  92.801099329237 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                  + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real ni56_pf, dni56_pf_dT;
    // interpolating ni56 partition function
    get_partition_function(Ni56, tfactors, ni56_pf, dni56_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real fe52_pf, dfe52_pf_dT;
    // interpolating fe52 partition function
    get_partition_function(Fe52, tfactors, fe52_pf, dfe52_pf_dT);

    Real z_r = he4_pf * fe52_pf;
    Real z_p = ni56_pf;

    Real dz_r_dT = fe52_pf * dhe4_pf_dT + he4_pf * dfe52_pf_dT;
    Real dz_p_dT = dni56_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ni56_to_p_co55_derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ni56 --> p + co55

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  63.131770608640906 + -83.14741674893808 * tfactors.T9i + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                         + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  83.14741674893808 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                  + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real ni56_pf, dni56_pf_dT;
    // interpolating ni56 partition function
    get_partition_function(Ni56, tfactors, ni56_pf, dni56_pf_dT);

    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real co55_pf, dco55_pf_dT;
    // interpolating co55 partition function
    get_partition_function(Co55, tfactors, co55_pf, dco55_pf_dT);

    Real z_r = p_pf * co55_pf;
    Real z_p = ni56_pf;

    Real dz_r_dT = co55_pf * dp_pf_dT + p_pf * dco55_pf_dT;
    Real dz_p_dT = dni56_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_co55_to_he4_fe52_derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // co55 + p --> he4 + fe52

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.443418619463486 + -9.65364776674457 * tfactors.T9i + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                         + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  9.65364776674457 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                  + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real co55_pf, dco55_pf_dT;
    // interpolating co55 partition function
    get_partition_function(Co55, tfactors, co55_pf, dco55_pf_dT);

    Real fe52_pf, dfe52_pf_dT;
    // interpolating fe52 partition function
    get_partition_function(Fe52, tfactors, fe52_pf, dfe52_pf_dT);

    Real z_r = he4_pf * fe52_pf;
    Real z_p = p_pf * co55_pf;

    Real dz_r_dT = fe52_pf * dhe4_pf_dT + he4_pf * dfe52_pf_dT;
    Real dz_p_dT = co55_pf * dp_pf_dT + p_pf * dco55_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}



template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, T& rate_eval)
{

    Real rate;
    Real drate_dT;

    rate_he4_fe52_to_ni56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_fe52_to_ni56) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_fe52_to_ni56) = drate_dT;

    }
    rate_p_co55_to_ni56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_co55_to_ni56) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_co55_to_ni56) = drate_dT;

    }
    rate_he4_fe52_to_p_co55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_fe52_to_p_co55) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_fe52_to_p_co55) = drate_dT;

    }
    rate_ni56_to_he4_fe52_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ni56_to_he4_fe52_derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ni56_to_he4_fe52_derived) = drate_dT;

    }
    rate_ni56_to_p_co55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ni56_to_p_co55_derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ni56_to_p_co55_derived) = drate_dT;

    }
    rate_p_co55_to_he4_fe52_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_co55_to_he4_fe52_derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_co55_to_he4_fe52_derived) = drate_dT;

    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_approx_rates([[maybe_unused]] const tf_t& tfactors, [[maybe_unused]] T& rate_eval)
{

    [[maybe_unused]] Real rate{};
    [[maybe_unused]] Real drate_dT{};


}

#endif
