#ifndef TEMPERATURE_TABLE_RATES_H
#define TEMPERATURE_TABLE_RATES_H

#include <AMReX_Array.H>

#include <tfactors.H>
#include <actual_network.H>


using namespace amrex::literals;

using namespace Rates;


namespace temp_tabular {

    // this can be replaced by std::numbers::ln10_v when we go to C++ 20
    constexpr amrex::Real ln10_inv = 1.0_rt / 2.302585092994046_rt;


    template <int do_T_derivatives, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    fill_rates(const tf_t& tfactors, T& rate_eval)
    {

        amrex::Real rate;
        amrex::Real drate_dT;


    }


}
#endif
