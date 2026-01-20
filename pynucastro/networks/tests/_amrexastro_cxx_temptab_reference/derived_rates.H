#ifndef DERIVED_RATES_H
#define DERIVED_RATES_H

#include <tfactors.H>
#include <actual_network.H>
#include <partition_functions.H>
#include <temperature_table_rates.H>
#include <interp_tools.H>

using namespace Rates;

template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_F17_to_O16_p_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F17 --> O16 + p

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F17_log_pf, dF17_log_pf_dT9;

    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors.T9, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = O16_log_pf + p_log_pf - F17_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dO16_log_pf_dT9 + dp_log_pf_dT9 - dF17_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::cubic_interp_uneven<do_T_derivatives>(
                                               tfactors.lnT9,
                                               O16_p_to_F17_iliadis_data::log_t9,
                                               O16_p_to_F17_iliadis_data::log_rate);

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -0.6002689639985703 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 21.841444037879352 + Q_kBT + net_log_pf;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

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

    rate_F17_to_O16_p_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_F17_to_O16_p_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F17_to_O16_p_derived) = drate_dT;

    }

}
#endif
