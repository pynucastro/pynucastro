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
void
fill_derived_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::Real rate;
    amrex::Real drate_dT;


}
#endif
