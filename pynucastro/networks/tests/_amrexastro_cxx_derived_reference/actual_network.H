#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates
{

    enum NetworkRates
    {
        k_he4_fe52_to_ni56 = 1,
        k_p_co55_to_ni56 = 2,
        k_he4_fe52_to_p_co55 = 3,
        k_ni56_to_he4_fe52_derived = 4,
        k_ni56_to_p_co55_derived = 5,
        k_p_co55_to_he4_fe52_derived = 6,
        NumRates = k_p_co55_to_he4_fe52_derived
    };

    // number of reaclib rates

    const int NrateReaclib = 6;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "he4_fe52_to_ni56",  // 1,
        "p_co55_to_ni56",  // 2,
        "he4_fe52_to_p_co55",  // 3,
        "ni56_to_he4_fe52_derived",  // 4,
        "ni56_to_p_co55_derived",  // 5,
        "p_co55_to_he4_fe52_derived"  // 6,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int h1_index = 0;
    constexpr int n_index = -1;
    constexpr int he4_index = 1;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
