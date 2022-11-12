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
        k_he4_mg24__si28__removed = 1,
        k_he4_mg24__p_al27__removed = 2,
        k_p_al27__si28__removed = 3,
        k_si28__he4_mg24__removed = 4,
        k_si28__p_al27__removed = 5,
        k_p_al27__he4_mg24__removed = 6,
        k_he4_si28__s32__removed = 7,
        k_he4_si28__p_p31__removed = 8,
        k_p_p31__s32__removed = 9,
        k_s32__he4_si28__removed = 10,
        k_s32__p_p31__removed = 11,
        k_p_p31__he4_si28__removed = 12,
        k_mg24_he4__si28__approx = 13,
        k_si28__mg24_he4__approx = 14,
        k_si28_he4__s32__approx = 15,
        k_s32__si28_he4__approx = 16,
        NumRates = k_s32__si28_he4__approx
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
        "he4_mg24__si28__removed",  // 1,
        "he4_mg24__p_al27__removed",  // 2,
        "p_al27__si28__removed",  // 3,
        "si28__he4_mg24__removed",  // 4,
        "si28__p_al27__removed",  // 5,
        "p_al27__he4_mg24__removed",  // 6,
        "he4_si28__s32__removed",  // 7,
        "he4_si28__p_p31__removed",  // 8,
        "p_p31__s32__removed",  // 9,
        "s32__he4_si28__removed",  // 10,
        "s32__p_p31__removed",  // 11,
        "p_p31__he4_si28__removed",  // 12,
        "mg24_he4__si28__approx",  // 13,
        "si28__mg24_he4__approx",  // 14,
        "si28_he4__s32__approx",  // 15,
        "s32__si28_he4__approx"  // 16,
    };

}

#endif
