#ifndef MODIFIED_RATES_H
#define MODIFIED_RATES_H

#include <amrex_bridge.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>
#include <reaclib_rates.H>
#include <beta_limited_rates.H>

namespace modified_rates {

    using namespace reaclib_rates;
    using namespace Rates;
    using namespace beta_limited_rates;


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
