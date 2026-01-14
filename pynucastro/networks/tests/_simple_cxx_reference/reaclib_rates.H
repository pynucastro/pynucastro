#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <amrex_bridge.H>

#include <tfactors.H>
#include <actual_network.H>
#include <partition_function.H>

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
inline
void rate_C12_C12_to_He4_Ne20_reaclib(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // C12 + C12 --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  61.2863 + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                         + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                  + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
inline
void rate_C12_C12_to_n_Mg23_reaclib(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // C12 + C12 --> n + Mg23

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  -12.8056 + -30.1498 * tfactors.T9i + 11.4826 * tfactors.T913
                         + 1.82849 * tfactors.T9 + -0.34844 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  30.1498 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.4826 * tfactors.T923i
                                  + 1.82849 + (5.0/3.0) * -0.34844 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
inline
void rate_C12_C12_to_p_Na23_reaclib(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // C12 + C12 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  60.9649 + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                         + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                                  + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
inline
void rate_He4_C12_to_O16_reaclib(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // C12 + He4 --> O16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nac2 
    ln_set_rate =  254.634 + -1.84097 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 137.303 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.84097 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                  + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 137.303 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // nac2 
    ln_set_rate =  69.6526 + -1.39254 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                         + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 70.3554 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.39254 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                  + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 70.3554 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
inline
void rate_n_to_p_reaclib(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n --> p

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -6.78161;
    amrex::ignore_unused(tfactors);

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}


template <typename T>
inline
void rate_He4_Ne20_to_C12_C12_derived(const tf_t& tfactors, Real& rate, Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 + He4 --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  61.474151468919175 + -53.57824966896331 * tfactors.T9i + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                         + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  53.57824966896331 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                  + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <typename T>
inline
void rate_p_Na23_to_C12_C12_derived(const tf_t& tfactors, Real& rate, Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na23 + p --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  60.92541574740554 + -26.004360836807493 * tfactors.T9i + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                         + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  26.004360836807493 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                                  + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <typename T>
inline
void rate_O16_to_He4_C12_derived(const tf_t& tfactors, Real& rate, Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 --> He4 + C12

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nac2 
    ln_set_rate =  279.29694929711803 + -84.95157686791681 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 138.803 * tfactors.lnT9;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.95157686791681 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                  + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 138.803 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // nac2 
    ln_set_rate =  94.31554929711804 + -84.50314686791681 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                         + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 71.8554 * tfactors.lnT9;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.50314686791681 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                  + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 71.8554 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}



template <int do_T_derivatives, typename T>
inline
void
fill_reaclib_rates(const tf_t& tfactors, T& rate_eval)
{

    Real rate;
    Real drate_dT;

    rate_C12_C12_to_He4_Ne20_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib) = drate_dT;

    }
    rate_C12_C12_to_n_Mg23_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_n_Mg23_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23_reaclib) = drate_dT;

    }
    rate_C12_C12_to_p_Na23_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib) = drate_dT;

    }
    rate_He4_C12_to_O16_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_C12_to_O16_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16_reaclib) = drate_dT;

    }
    rate_n_to_p_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_to_p_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_to_p_reaclib) = drate_dT;

    }

}

template <int do_T_derivatives, typename T>
inline
void
fill_derived_rates(const tf_t& tfactors, T& rate_eval)
{

    Real rate;
    Real drate_dT;

    part_fun::pf_cache_t pf_cache{};

    rate_He4_Ne20_to_C12_C12_derived<T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived) = drate_dT;

    }
    rate_p_Na23_to_C12_C12_derived<T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_Na23_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12_derived) = drate_dT;

    }
    rate_O16_to_He4_C12_derived<T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;

    }

}

template <int do_T_derivatives, typename T>
inline
void
fill_approx_rates([[maybe_unused]] const tf_t& tfactors,
                  [[maybe_unused]] const Real rho,
                  [[maybe_unused]] const Array1D<Real, 1, NumSpec>& Y,
                  [[maybe_unused]] T& rate_eval)
{

    Real rate{};
    Real drate_dT{};


}

#endif
