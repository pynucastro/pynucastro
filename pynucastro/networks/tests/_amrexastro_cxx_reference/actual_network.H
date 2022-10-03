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
        k_c12_c12__he4_ne20 = 1,
        k_c12_c12__n_mg23 = 2,
        k_c12_c12__p_na23 = 3,
        k_he4_c12__o16 = 4,
        k_n__p__weak__wc12 = 5,
        k_na23__ne23 = 6,
        k_ne23__na23 = 7,
        NumRates = k_ne23__na23
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
        "c12_c12__he4_ne20",  // 1,
        "c12_c12__n_mg23",  // 2,
        "c12_c12__p_na23",  // 3,
        "he4_c12__o16",  // 4,
        "n__p__weak__wc12",  // 5,
        "na23__ne23",  // 6,
        "ne23__na23"  // 7,
    };

}

#endif
