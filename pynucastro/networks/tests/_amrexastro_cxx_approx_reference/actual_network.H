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
        k_he4_mg24_to_si28_removed = 1,
        k_he4_mg24_to_p_al27_removed = 2,
        k_p_al27_to_si28_removed = 3,
        k_si28_to_he4_mg24_removed = 4,
        k_si28_to_p_al27_removed = 5,
        k_p_al27_to_he4_mg24_removed = 6,
        k_he4_si28_to_s32_removed = 7,
        k_he4_si28_to_p_p31_removed = 8,
        k_p_p31_to_s32_removed = 9,
        k_s32_to_he4_si28_removed = 10,
        k_s32_to_p_p31_removed = 11,
        k_p_p31_to_he4_si28_removed = 12,
        k_mg24_he4_to_si28_approx = 13,
        k_si28_to_mg24_he4_approx = 14,
        k_si28_he4_to_s32_approx = 15,
        k_s32_to_si28_he4_approx = 16,
        NumRates = k_s32_to_si28_he4_approx
    };

    // number of reaclib rates

    const int NrateReaclib = 12;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "he4_mg24_to_si28_removed",  // 1,
        "he4_mg24_to_p_al27_removed",  // 2,
        "p_al27_to_si28_removed",  // 3,
        "si28_to_he4_mg24_removed",  // 4,
        "si28_to_p_al27_removed",  // 5,
        "p_al27_to_he4_mg24_removed",  // 6,
        "he4_si28_to_s32_removed",  // 7,
        "he4_si28_to_p_p31_removed",  // 8,
        "p_p31_to_s32_removed",  // 9,
        "s32_to_he4_si28_removed",  // 10,
        "s32_to_p_p31_removed",  // 11,
        "p_p31_to_he4_si28_removed",  // 12,
        "mg24_he4_to_si28_approx",  // 13,
        "si28_to_mg24_he4_approx",  // 14,
        "si28_he4_to_s32_approx",  // 15,
        "s32_to_si28_he4_approx"  // 16,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    extern AMREX_GPU_MANAGED int h1_index;
    extern AMREX_GPU_MANAGED int n_index;
    extern AMREX_GPU_MANAGED int he4_index;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
    extern AMREX_GPU_MANAGED bool initialized;
}
#endif

#endif
