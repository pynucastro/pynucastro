#ifndef TABULAR_WEAK_RATES_H
#define TABULAR_WEAK_RATES_H

#include <cmath>

#include <AMReX.H>

#include <actual_network.H>
#include <fundamental_constants.H>
#include <network_properties.H>
#include <rate_type.H>
#include <table_rates.H>

namespace tabular_weak_rates {

    using namespace Rates;
    using namespace Species;

    using namespace rate_tables;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na23_to_Ne23_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Na23 --> Ne23

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Na23_Ne23_meta, j_Na23_Ne23_rhoy, j_Na23_Ne23_temp, j_Na23_Ne23_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Na23_to_Ne23_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_to_Ne23_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Na23 and Ne23
            // this is the derivative of ∂Y(Na23)/∂t = -Y(Na23) λ and
            //                           ∂Y(Ne23)/∂t = +Y(Na23) λ
            rate_eval.dweak_ydot_dYe(Na23) -= rho * Y(Na23) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ne23) += rho * Y(Na23) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} on T and Y        rate_eval.denuc_weak_dY(Na23) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Na23) * table_values.denu_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Na23) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne23_to_Na23_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ne23 --> Na23

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ne23_Na23_meta, j_Ne23_Na23_rhoy, j_Ne23_Na23_temp, j_Ne23_Na23_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ne23_to_Na23_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne23_to_Na23_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ne23 and Na23
            // this is the derivative of ∂Y(Ne23)/∂t = -Y(Ne23) λ and
            //                           ∂Y(Na23)/∂t = +Y(Ne23) λ
            rate_eval.dweak_ydot_dYe(Ne23) -= rho * Y(Ne23) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Na23) += rho * Y(Ne23) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} on T and Y        rate_eval.denuc_weak_dY(Ne23) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ne23) * table_values.denu_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ne23) * (table_values.enu + table_values.gamma);
    }


    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const amrex::Real temp,
               [[maybe_unused]] const amrex::Real rho,
               [[maybe_unused]] const amrex::Real rhoy,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        // Calculate tabular weak rates

        const amrex::Real log_temp = std::log10(temp);
        const amrex::Real log_rhoy = std::log10(rhoy);

        rate_eval.enuc_weak = 0.0_rt;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            for (int n = 1; n <= NumSpec; ++n) {
                rate_eval.dweak_ydot_dYe(n) = 0.0_rt;
                rate_eval.denuc_weak_dY(n) = 0.0_rt;
            }
            rate_eval.denuc_weak_dT = 0.0_rt;

        }

        {
            // Na23_to_Ne23_electron_capture_weaktab
            rate_Na23_to_Ne23_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ne23_to_Na23_beta_neg_weaktab
            rate_Ne23_to_Na23_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

    }

}
#endif
