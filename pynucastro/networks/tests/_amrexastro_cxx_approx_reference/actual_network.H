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

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real mion (int spec) {
        using namespace Species;

        amrex::Real m = 0.0_rt;

        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                m = mion<n>();
            }
        });

        return m;
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
    // A RatePair data table, 2D array of shape (NumNSERatePairs, 10).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Use enum names following NetworkSpecies and NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Columns 7-8: 1-based index of non (n, H1, He4) nuclei.
    //  - Column  9:   1-based index of the forward rate
    //  - Column  10:  1-based index of the reverse rate

    using namespace Species;
    using namespace Rates;

    constexpr int H1_index = -1;
    constexpr int N_index = -1;
    constexpr int He4_index = 1;
    constexpr int NumNSERatePairs = 2;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        Mg24, He4, -1, Si28, -1, -1, Mg24, Si28, k_Mg24_He4_to_Si28_approx, k_Si28_to_Mg24_He4_approx,
        Si28, He4, -1, S32, -1, -1, Si28, S32, k_Si28_He4_to_S32_approx, k_S32_to_Si28_He4_approx
    };
}
#endif

#endif
