#ifndef actual_network_H
#define actual_network_H

#include <amrex_bridge.H>

#include <fundamental_constants.H>
#include <network_properties.H>

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{
    extern Array1D<Real, 1, NumSpec> bion;
    extern Array1D<Real, 1, NumSpec> mion;
}

namespace Rates
{

    enum NetworkRates
    {
        k_C12_C12_to_He4_Ne20_reaclib = 1,
        k_C12_C12_to_n_Mg23_reaclib = 2,
        k_C12_C12_to_p_Na23_reaclib = 3,
        k_He4_C12_to_O16_reaclib = 4,
        k_n_to_p_reaclib = 5,
        NumRates = k_n_to_p_reaclib
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "C12_C12_to_He4_Ne20_reaclib",  // 1,
        "C12_C12_to_n_Mg23_reaclib",  // 2,
        "C12_C12_to_p_Na23_reaclib",  // 3,
        "He4_C12_to_O16_reaclib",  // 4,
        "n_to_p_reaclib"  // 5,
    };

}

#endif
