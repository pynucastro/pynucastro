#ifndef DERIVED_RATES_H
#define DERIVED_RATES_H

#include <tfactors.H>
#include <actual_network.H>
#include <partition_functions.H>
#include <temperature_table_rates.H>
#include <interp_tools.H>

using namespace Rates;

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_F17_to_O16_p_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F17 --> O16 + p

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 9.674279321716613e-05;
    amrex::Real Q_kBT = -0.6002689639985703 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_O16_p_to_F17_iliadis);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_O16_p_to_F17_iliadis);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F17_pf, dF17_pf_dT;
    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors, pf_cache, F17_pf, dF17_pf_dT);

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real z_r = O16_pf * p_pf;
    amrex::Real z_p = F17_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = p_pf * dO16_pf_dT + O16_pf * dp_pf_dT;
        amrex::Real dz_p_dT = dF17_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_derived_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::ignore_unused(tfactors);
    amrex::ignore_unused(rate_eval);

    [[maybe_unused]] amrex::Real rate;
    [[maybe_unused]] amrex::Real drate_dT;

    part_fun::pf_cache_t pf_cache{};

    pf_cache.index_temp_array_1 = interp_net::find_index(tfactors.T9, part_fun::temp_array_1);

    rate_F17_to_O16_p_derived<T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_F17_to_O16_p_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F17_to_O16_p_derived) = drate_dT;

    }

}
#endif
