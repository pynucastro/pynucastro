#ifndef BRANCHED_RATES_H
#define BRANCHED_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>
#include <temperature_table_rates.H>
#include <reaclib_rates.H>

namespace branched_rates {

    using namespace Rates;

    using namespace temp_tabular;


    template <int do_T_derivatives, typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate;
        [[maybe_unused]] amrex::Real drate_dT;


    }
}
#endif
