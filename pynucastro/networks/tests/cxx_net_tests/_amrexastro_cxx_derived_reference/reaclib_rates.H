#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>


namespace reaclib_rates {

    using namespace Rates;
    using namespace Species;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_He4_to_Ni56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe52 + He4 --> Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  66.6417 + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                             + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                      + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe52_He4_to_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_He4_to_Ni56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_p_to_Ni56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co55 + p --> Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  37.3736 + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                             + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                      + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co55_p_to_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_p_to_Ni56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_He4_to_p_Co55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe52 + He4 --> p + Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.2207 + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                             + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                      + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe52_He4_to_p_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_He4_to_p_Co55_reaclib) = drate_dT;
        }
    }


}
#endif
