#ifndef MODIFIED_RATES_H
#define MODIFIED_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

using namespace Rates;
using namespace Species;


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_modified_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::Real rate;
    amrex::Real drate_dT;


}

#endif
