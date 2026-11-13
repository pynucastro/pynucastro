#ifndef FILL_RATES_H
#define FILL_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

namespace rate_util {

    using namespace Rates;
    using namespace Species;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] const amrex::Real temp,
               [[maybe_unused]] const amrex::Real rho,
               [[maybe_unused]] const amrex::Real rhoy,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate{};
        [[maybe_unused]] amrex::Real drate_dT{};

        [[maybe_unused]] const amrex::Real log_temp = std::log10(temp);
        [[maybe_unused]] const amrex::Real log_rhoy = std::log10(rhoy);

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

        {
            // Mg24_He4_to_p_Al27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            reaclib_rates::rate_Mg24_He4_to_p_Al27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_He4_to_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            reaclib_rates::rate_Mg24_He4_to_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_p_to_He4_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            reaclib_rates::rate_Al27_p_to_He4_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_p_to_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            reaclib_rates::rate_Al27_p_to_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si28_to_He4_Mg24_reaclib
            reaclib_rates::rate_Si28_to_He4_Mg24_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_to_p_Al27_reaclib
            reaclib_rates::rate_Si28_to_p_Al27_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_He4_to_p_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            reaclib_rates::rate_Si28_He4_to_p_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si28_He4_to_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            reaclib_rates::rate_Si28_He4_to_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_p_to_He4_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            reaclib_rates::rate_P31_p_to_He4_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_p_to_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            reaclib_rates::rate_P31_p_to_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S32_to_He4_Si28_reaclib
            reaclib_rates::rate_S32_to_He4_Si28_reaclib(rate_eval, tfactors);
        }

        {
            // S32_to_p_P31_reaclib
            reaclib_rates::rate_S32_to_p_P31_reaclib(rate_eval, tfactors);
        }

        {
            // Mg24_He4_to_Si28_approx
            approximate_rates::rate_Mg24_He4_to_Si28_approx(rate_eval);
        }

        {
            // Si28_to_Mg24_He4_approx
            approximate_rates::rate_Si28_to_Mg24_He4_approx(rate_eval);
        }

        {
            // Si28_He4_to_S32_approx
            approximate_rates::rate_Si28_He4_to_S32_approx(rate_eval);
        }

        {
            // S32_to_Si28_He4_approx
            approximate_rates::rate_S32_to_Si28_He4_approx(rate_eval);
        }


    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_weak_rates([[maybe_unused]] const tf_t& tfactors,
                    [[maybe_unused]] const amrex::Real temp,
                    [[maybe_unused]] const amrex::Real rho,
                    [[maybe_unused]] const amrex::Real rhoy,
                    [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                    [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate{};
        [[maybe_unused]] amrex::Real drate_dT{};

        [[maybe_unused]] const amrex::Real log_temp = std::log10(temp);
        [[maybe_unused]] const amrex::Real log_rhoy = std::log10(rhoy);

        // we are not doing derivatives
        static_assert(std::is_same_v<T, rate_t>);

        rate_eval.enuc_weak = 0.0_rt;

        // rates are filled in a topologically-sorted order


    }

}

#endif
