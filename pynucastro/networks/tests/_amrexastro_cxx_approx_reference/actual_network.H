#ifndef actual_network_H
#define actual_network_H

#include <string_view>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;

void actual_network_init();

constexpr std::string_view network_name = "pynucastro-cxx";

namespace network
{

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element

        if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.53684539599638_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real mion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the mass of the element

        if constexpr (spec == He4) {
            return 6.646479071584587e-24_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
        }
        else if constexpr (spec == Si28) {
            return 4.6456779473820677e-23_rt;
        }
        else if constexpr (spec == S32) {
            return 5.309087322384128e-23_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real bion (int spec) {
        using namespace Species;

        amrex::Real b = 0.0_rt;

        // Set the binding energy of the element
        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                b = bion<n>();
            }
        });

        return b;
    }

}

namespace Rates
{

    enum NetworkRates :
    std::uint8_t
    {
        k_He4_Mg24_to_Si28_removed = 1,
        k_He4_Mg24_to_p_Al27_removed = 2,
        k_p_Al27_to_Si28_removed = 3,
        k_Si28_to_He4_Mg24_removed = 4,
        k_Si28_to_p_Al27_removed = 5,
        k_p_Al27_to_He4_Mg24_removed = 6,
        k_He4_Si28_to_S32_removed = 7,
        k_He4_Si28_to_p_P31_removed = 8,
        k_p_P31_to_S32_removed = 9,
        k_S32_to_He4_Si28_removed = 10,
        k_S32_to_p_P31_removed = 11,
        k_p_P31_to_He4_Si28_removed = 12,
        k_Mg24_He4_to_Si28_approx = 13,
        k_Si28_to_Mg24_He4_approx = 14,
        k_Si28_He4_to_S32_approx = 15,
        k_S32_to_Si28_He4_approx = 16,
        NumRates = k_S32_to_Si28_He4_approx
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "He4_Mg24_to_Si28_removed",  // 1,
        "He4_Mg24_to_p_Al27_removed",  // 2,
        "p_Al27_to_Si28_removed",  // 3,
        "Si28_to_He4_Mg24_removed",  // 4,
        "Si28_to_p_Al27_removed",  // 5,
        "p_Al27_to_He4_Mg24_removed",  // 6,
        "He4_Si28_to_S32_removed",  // 7,
        "He4_Si28_to_p_P31_removed",  // 8,
        "p_P31_to_S32_removed",  // 9,
        "S32_to_He4_Si28_removed",  // 10,
        "S32_to_p_P31_removed",  // 11,
        "p_P31_to_He4_Si28_removed",  // 12,
        "Mg24_He4_to_Si28_approx",  // 13,
        "Si28_to_Mg24_He4_approx",  // 14,
        "Si28_He4_to_S32_approx",  // 15,
        "S32_to_Si28_He4_approx"  // 16,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = -1;
    constexpr int N_index = -1;
    constexpr int He4_index = 0;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices;
}
#endif

#endif
