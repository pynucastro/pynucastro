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
    const int NumReaclibSets = 6;

    // number of tabular rates

    const int NrateTabular = 2;

}

namespace reaclib_rates
{

    // Temperature coefficient arrays (numbers correspond to reaction
    // numbers in net_info)

    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, 7, 1, Rates::NumReaclibSets> ctemp_rate;

    // Index into ctemp_rate, dimension 2, where each rate's
    // coefficients start

    extern AMREX_GPU_MANAGED amrex::Array1D<int, 1, Rates::NrateReaclib> rate_start_idx;

    // Reaction multiplicities-1 (how many rates contribute - 1)

    extern AMREX_GPU_MANAGED amrex::Array1D<int, 1, Rates::NrateReaclib> rate_extra_mult;

}

#endif
