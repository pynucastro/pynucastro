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
        k_c12_c12_to_he4_ne20 = 1,
        k_c12_c12_to_n_mg23 = 2,
        k_c12_c12_to_p_na23 = 3,
        k_he4_c12_to_o16 = 4,
        k_n_to_p_weak_wc12 = 5,
        NumRates = k_n_to_p_weak_wc12
    };

    // number of reaclib rates

    const int NrateReaclib = 5;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "c12_c12_to_he4_ne20",  // 1,
        "c12_c12_to_n_mg23",  // 2,
        "c12_c12_to_p_na23",  // 3,
        "he4_c12_to_o16",  // 4,
        "n_to_p_weak_wc12"  // 5,
    };

}

#endif
