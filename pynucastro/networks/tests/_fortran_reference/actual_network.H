#ifndef actual_network_H
#define actual_network_H

#include <amrex_bridge.H>

#include <fundamental_constants.H>
#include <network_properties.H>

inline void actual_network_init() {}

const std::string network_name = "pynucastro-cxx";

namespace network
{

    // nuclei masses (in grams)

    inline Array1D<Real, 1, NumSpec> mion{
        1.674927498034172e-24_rt,  // n
        1.6735328377636005e-24_rt,  // p
        6.646479071584587e-24_rt,  // He4
        1.99264687992e-23_rt,  // C12
        2.6560180592333686e-23_rt,  // O16
        3.3198227947612416e-23_rt,  // Ne20
        3.817541002484691e-23_rt,  // Na23
        3.8182640828719474e-23_rt,  // Mg23
        3.9828098739467446e-23_rt,  // Mg24
    };
}

namespace Rates
{

    enum NetworkRates
    {
        k_C12_C12_to_He4_Ne20_reaclib = 1,
        k_C12_C12_to_n_Mg23_reaclib = 2,
        k_C12_C12_to_p_Na23_reaclib = 3,
        k_He4_C12_to_O16_reaclib = 4,
        k_He4_Ne20_to_Mg24_reaclib = 5,
        k_n_Mg23_to_Mg24_reaclib = 6,
        k_n_to_p_reaclib = 7,
        k_He4_Ne20_to_C12_C12_derived = 8,
        k_p_Na23_to_C12_C12_derived = 9,
        k_O16_to_He4_C12_derived = 10,
        k_Mg24_to_He4_Ne20_derived = 11,
        k_Mg24_to_n_Mg23_derived = 12,
        NumRates = k_Mg24_to_n_Mg23_derived
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
        "He4_Ne20_to_Mg24_reaclib",  // 5,
        "n_Mg23_to_Mg24_reaclib",  // 6,
        "n_to_p_reaclib",  // 7,
        "He4_Ne20_to_C12_C12_derived",  // 8,
        "p_Na23_to_C12_C12_derived",  // 9,
        "O16_to_He4_C12_derived",  // 10,
        "Mg24_to_He4_Ne20_derived",  // 11,
        "Mg24_to_n_Mg23_derived"  // 12,
    };

}

#endif
