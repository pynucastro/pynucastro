#ifndef BETA_LIMITED_RATES_H
#define BETA_LIMITED_RATES_H

#include <amrex_bridge.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

namespace beta_limited_rates {

    using namespace Rates;


    template <int do_T_derivatives, typename T>
    inline
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] const Real rho,
               [[maybe_unused]] const Array1D<Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] Real rate;
        [[maybe_unused]] Real drate_dT;


    }
}
#endif
