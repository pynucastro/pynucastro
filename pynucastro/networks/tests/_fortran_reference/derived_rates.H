#ifndef DERIVED_RATES_H
#define DERIVED_RATES_H

#include <amrex_bridge.H>

#include <tfactors.H>
#include <actual_network.H>
#include <partition_function.H>

using namespace Rates;

template <int do_T_derivatives, typename T>
inline
void rate_He4_Ne20_to_C12_C12_derived(const tf_t& tfactors, Real& rate, Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 + He4 --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Assume no partition function effects
    Real net_log_pf{0.0};
    [[maybe_unused]] Real net_dlog_pf_dT9{0.0};

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  61.474151468919175 + -53.57824966896331 * tfactors.T9i + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                         + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  53.57824966896331 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                  + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
inline
void rate_p_Na23_to_C12_C12_derived(const tf_t& tfactors, Real& rate, Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na23 + p --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Assume no partition function effects
    Real net_log_pf{0.0};
    [[maybe_unused]] Real net_dlog_pf_dT9{0.0};

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  60.92541574740554 + -26.004360836807496 * tfactors.T9i + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                         + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  26.004360836807496 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                                  + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
inline
void rate_O16_to_He4_C12_derived(const tf_t& tfactors, Real& rate, Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 --> He4 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Assume no partition function effects
    Real net_log_pf{0.0};
    [[maybe_unused]] Real net_dlog_pf_dT9{0.0};

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nac2 
    ln_set_rate =  279.29694929711803 + -84.95157686791683 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 138.803 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.95157686791683 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                  + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 138.803 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // nac2 
    ln_set_rate =  94.31554929711803 + -84.50314686791683 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                         + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 71.8554 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.50314686791683 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                  + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 71.8554 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
inline
void rate_Mg24_to_He4_Ne20_derived(const tf_t& tfactors, Real& rate, Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    // Assume no partition function effects
    Real net_log_pf{0.0};
    [[maybe_unused]] Real net_dlog_pf_dT9{0.0};

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  16.022538839949377 + -120.89510873006793 * tfactors.T9i + 16.9229 * tfactors.T913
                         + -2.57325 * tfactors.T9 + 0.208997 * tfactors.T953;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  120.89510873006793 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 16.9229 * tfactors.T923i
                                  + -2.57325 + (5.0/3.0) * 0.208997 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // il10r
    ln_set_rate =  26.803878839949377 + -117.33446873006793 * tfactors.T9i;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  117.33446873006793 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // il10r
    ln_set_rate =  -13.884691160050625 + -110.62025873006793 * tfactors.T9i;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  110.62025873006793 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // il10n
    ln_set_rate =  49.32660883994937 + -108.11420873006793 * tfactors.T9i + -46.2525 * tfactors.T913i + 5.58901 * tfactors.T913
                         + 7.61843 * tfactors.T9 + -3.683 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  108.11420873006793 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.2525 * tfactors.T943i + (1.0/3.0) * 5.58901 * tfactors.T923i
                                  + 7.61843 + (5.0/3.0) * -3.683 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
inline
void rate_Mg24_to_n_Mg23_derived(const tf_t& tfactors, Real& rate, Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> n + Mg23

    rate = 0.0;
    drate_dT = 0.0;

    // Assume no partition function effects
    Real net_log_pf{0.0};
    [[maybe_unused]] Real net_dlog_pf_dT9{0.0};

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  32.05683058133977 + -191.83691284889673 * tfactors.T9i + 2.66964 * tfactors.T913
                         + -0.448904 * tfactors.T9 + 0.0326505 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  191.83691284889673 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.66964 * tfactors.T923i
                                  + -0.448904 + (5.0/3.0) * 0.0326505 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}


template <int do_T_derivatives, typename T>
inline
void
fill_derived_rates(const tf_t& tfactors, T& rate_eval)
{
    amrex::ignore_unused(tfactors);
    amrex::ignore_unused(rate_eval);

    [[maybe_unused]] Real rate;
    [[maybe_unused]] Real drate_dT;

    part_fun::pf_cache_t pf_cache{};

    rate_He4_Ne20_to_C12_C12_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived) = drate_dT;

    }
    rate_p_Na23_to_C12_C12_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_Na23_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12_derived) = drate_dT;

    }
    rate_O16_to_He4_C12_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;

    }
    rate_Mg24_to_He4_Ne20_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;

    }
    rate_Mg24_to_n_Mg23_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_n_Mg23_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_n_Mg23_derived) = drate_dT;

    }

}
#endif
