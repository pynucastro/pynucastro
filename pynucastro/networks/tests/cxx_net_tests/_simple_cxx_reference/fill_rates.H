#ifndef FILL_RATES_H
#define FILL_RATES_H

#include <amrex_bridge.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

namespace rate_util {

    using namespace Rates;
    using namespace Species;

    template <typename T>
    inline
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] const Real temp,
               [[maybe_unused]] const Real rho,
               [[maybe_unused]] const Real rhoy,
               [[maybe_unused]] const Array1D<Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] Real rate{};
        [[maybe_unused]] Real drate_dT{};

        [[maybe_unused]] const Real log_temp = std::log10(temp);
        [[maybe_unused]] const Real log_rhoy = std::log10(rhoy);

        rate_eval.enuc_weak = 0.0_rt;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            for (int n = 1; n <= NumSpec; ++n) {
                rate_eval.dweak_ydot_dYe(n) = 0.0_rt;
                rate_eval.denuc_weak_dY(n) = 0.0_rt;
            }
            rate_eval.denuc_weak_dT = 0.0_rt;
            rate_eval.denuc_weak_dYe = 0.0_rt;

        }

        // rates are filled in a topologically-sorted order

        part_fun::pf_cache_t pf_cache{};

        {
            // C12_C12_to_He4_Ne20_reaclib
            Real log_scor {0.0_rt};
            Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            reaclib_rates::rate_C12_C12_to_He4_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C12_C12_to_n_Mg23_reaclib
            Real log_scor {0.0_rt};
            Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            reaclib_rates::rate_C12_C12_to_n_Mg23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C12_C12_to_p_Na23_reaclib
            Real log_scor {0.0_rt};
            Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            reaclib_rates::rate_C12_C12_to_p_Na23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C12_He4_to_O16_reaclib
            Real log_scor {0.0_rt};
            Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C12);
            }
#endif
            reaclib_rates::rate_C12_He4_to_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // n_to_p_beta_neg_reaclib
            reaclib_rates::rate_n_to_p_beta_neg_reaclib(rate_eval, tfactors);
        }

        {
            // Ne20_He4_to_C12_C12_derived
            Real log_scor {0.0_rt};
            Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            derived_rates::rate_Ne20_He4_to_C12_C12_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Na23_p_to_C12_C12_derived
            Real log_scor {0.0_rt};
            Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            derived_rates::rate_Na23_p_to_C12_C12_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // O16_to_He4_C12_derived
            derived_rates::rate_O16_to_He4_C12_derived(rate_eval, tfactors, pf_cache);
        }


    }

}

#endif
