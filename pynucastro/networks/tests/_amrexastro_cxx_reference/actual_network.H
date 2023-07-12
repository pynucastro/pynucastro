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
        k_c12_c12_to_he4_ne20 = 1,
        k_c12_c12_to_n_mg23 = 2,
        k_c12_c12_to_p_na23 = 3,
        k_he4_c12_to_o16 = 4,
        k_n_to_p_weak_wc12 = 5,
        k_na23_to_ne23 = 6,
        k_ne23_to_na23 = 7,
        NumRates = k_ne23_to_na23
    };

    // number of reaclib rates

    const int NrateReaclib = 5;

    // number of tabular rates

    const int NrateTabular = 2;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "c12_c12_to_he4_ne20",  // 1,
        "c12_c12_to_n_mg23",  // 2,
        "c12_c12_to_p_na23",  // 3,
        "he4_c12_to_o16",  // 4,
        "n_to_p_weak_wc12",  // 5,
        "na23_to_ne23",  // 6,
        "ne23_to_na23"  // 7,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int h1_index = 1;
    constexpr int n_index = 0;
    constexpr int he4_index = 2;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
    extern AMREX_GPU_MANAGED bool initialized;
}
#endif

#endif
