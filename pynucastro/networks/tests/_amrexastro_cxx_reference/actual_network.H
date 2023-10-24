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
        k_C12_C12_to_He4_Ne20 = 1,
        k_C12_C12_to_n_Mg23 = 2,
        k_C12_C12_to_p_Na23 = 3,
        k_He4_C12_to_O16 = 4,
        k_n_to_p_weak_wc12 = 5,
        k_Na23_to_Ne23 = 6,
        k_Ne23_to_Na23 = 7,
        NumRates = k_Ne23_to_Na23
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
        "C12_C12_to_He4_Ne20",  // 1,
        "C12_C12_to_n_Mg23",  // 2,
        "C12_C12_to_p_Na23",  // 3,
        "He4_C12_to_O16",  // 4,
        "n_to_p_weak_wc12",  // 5,
        "Na23_to_Ne23",  // 6,
        "Ne23_to_Na23"  // 7,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = 0;
    constexpr int He4_index = 2;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
