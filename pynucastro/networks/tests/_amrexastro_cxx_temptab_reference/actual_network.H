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

        if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == F17) {
            return 128.21958437599824_rt;
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

        if constexpr (spec == H1) {
            return 1.6735328377636005e-24_rt;
        }
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == F17) {
            return 2.8232643353472436e-23_rt;
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
        k_O16_p_to_F17_iliadis = 1,
        k_F17_to_O16_p_derived = 2,
        NumRates = k_F17_to_O16_p_derived
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "O16_p_to_F17_iliadis",  // 1,
        "F17_to_O16_p_derived"  // 2,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = -1;
    constexpr int He4_index = -1;

    // A RatePair data table, 2D array of shape (NumNSERatePairs, 8).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Nuclei indices follow NetworkSpecies and rate indices follow NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Column 7:    1-based index of the forward rate
    //  - Column 8:    1-based index of the reverse rate

    constexpr int NumNSERatePairs = 1;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 8, amrex::Order::C> rate_pair_data {
        -1, 1, 2, -1, -1, 3, 1, 2  // fr: O16_p_to_F17_iliadis, rr: F17_to_O16_p_derived
    };
}
#endif

#endif
